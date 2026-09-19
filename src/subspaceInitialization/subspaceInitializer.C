/*---------------------------------------------------------------------------*\
                  Class subspaceInitializer Implementation
\*---------------------------------------------------------------------------*/

#include "subspaceInitializer.H"
#include "EigenMatrix.H"
#include "Switch.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subspaceInitializer::subspaceInitializer
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const word& fieldName,
    const dictionary& controlDict,
    const label maxLenHistory,
    const label maxNumProbes,
    const UList<scalar>& decayRates
)
:
    matrix_(matrix),
    interfaceBouCoeffs_(interfaceBouCoeffs),
    interfaces_(interfaces),
    meshPtr_(isA<const fvMesh>(matrix.mesh())),
    fieldName_(fieldName),
    maxLenHistory_
    (
        maxLenHistory >= 0
      ? maxLenHistory
      : controlDict.getOrDefault<label>("lenHistory", 0)
    ),
    maxNumProbes_
    (
        maxNumProbes >= 0
      ? maxNumProbes
      : controlDict.getOrDefault<label>("numProbes", 4)
    ),
    decayRates_(),
    decayRate_(controlDict.getOrDefault<scalar>("decayRate", 0)),
    truncTol_(controlDict.getOrDefault<scalar>("truncTol", 1e-14)),
    projection_
    (
        subspaceInitializer::readProjection(controlDict.getOrDefault<word>("projection", "galerkin"))
    ),
    randomSeed_(0),
    persistState_(Switch(controlDict.getOrDefault<word>("persistState", "no"))),
    window_(nullptr),
    ewmaSketches_(),
    pushed_(false)
{
    if (Foam::isNull(decayRates))
    {
        checkExclusive
        (
            controlDict.getOrDefault<label>("lenHistory", 0),
            decayRate_
        );

        if (decayRate_ > 0)
        {
            decayRates_ = List<scalar>(1, decayRate_);
        }
    }
    else
    {
        decayRates_ = decayRates;
    }

    if (meshPtr_)
    {
        randomSeed_ =
            meshPtr_->time().controlDict().getOrDefault<label>("randomSeed", 0);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::subspaceInitializer::projection
Foam::subspaceInitializer::readProjection(const word& name)
{
    if (name == "leastSquares")
    {
        return LEAST_SQUARES;
    }
    if (name != "galerkin")
    {
        FatalErrorInFunction
            << "projection must be galerkin or leastSquares, got " << name
            << exit(FatalError);
    }
    return GALERKIN;
}



void Foam::subspaceInitializer::checkExclusive
(
    const label lenHistory,
    const scalar decayRate
)
{
    if (lenHistory > 0 && decayRate > 0)
    {
        FatalErrorInFunction
            << "lenHistory (" << lenHistory << ") and decayRate ("
            << decayRate << ") are alternatives; set exactly one."
            << "  lenHistory keeps a window of past iterates and re-sketches"
            << " it; decayRate keeps no iterates and accumulates the sketch"
            << " itself.  There is no sketch that is both."
            << exit(FatalError);
    }
}


Foam::label
Foam::subspaceInitializer::rateIndex(const scalar decayRate) const
{
    // Rates that render to the same registry key share one sketch.
    const word rateKey = Foam::name(decayRate);

    forAll(decayRates_, i)
    {
        if (Foam::name(decayRates_[i]) == rateKey)
        {
            return i;
        }
    }

    return -1;
}


Foam::label Foam::subspaceInitializer::correct
(
    const direction cmpt,
    const label lenHistory,
    const label numProbes,
    const EWMASketch* sketch,
    solveScalarField& psi,
    solveScalarField& rA
) const
{
    const label nCells = psi.size();
    const bool streaming = (sketch != nullptr);

    // Limit the projection size to the configured capacity and stored samples.
    const label windowDepth =
        streaming ? 0 : min(min(lenHistory, maxLenHistory_), window_->depth());

    const label numDirections =
        streaming
      ? min(min(numProbes, maxNumProbes_), min(sketch->numProbes(), sketch->depth()))
      : min(min(numProbes, maxNumProbes_), windowDepth);

    if (numDirections <= 0)
    {
        return 0;
    }

    if (projection_ == GALERKIN && !matrix_.symmetric())
    {
        // The Galerkin Gram matrix is stored and solved as symmetric.
        FatalErrorInFunction
            << "projection galerkin requires a symmetric matrix;"
            << " use leastSquares" << exit(FatalError);
    }

    const label comm = matrix_.mesh().comm();

    // Build the correction basis from past iterates shifted by the entry psi.
    List<solveScalarField> basis(numDirections);
    List<solveScalar*> basisData(numDirections);
    for (label j = 0; j < numDirections; ++j)
    {
        // Window contributions accumulate; EWMA assigns each element.
        if (streaming)
        {
            basis[j].resize_nocopy(nCells);
        }
        else
        {
            basis[j].resize_fill(nCells, Zero);
        }
        basisData[j] = basis[j].begin();
    }

    if (streaming)
    {
        // Shift the stored EWMA sketch: Omega = Y - psi0*w^T.
        const solveScalar* __restrict__ initialPsi = psi.begin();

        for (label j = 0; j < numDirections; ++j)
        {
            const solveScalar* __restrict__ sketchColumn = sketch->column(j).begin();
            const solveScalar sketchWeight = sketch->weight(j);
            solveScalar* __restrict__ basisColumn = basisData[j];

            for (label cell = 0; cell < nCells; ++cell)
            {
                basisColumn[cell] = sketchColumn[cell] - sketchWeight*initialPsi[cell];
            }
        }
    }
    else
    {
        const RectangularMatrix<scalar>& Z = window_->Z();
        const solveScalar* __restrict__ initialPsi = psi.begin();

        for (label k = 0; k < windowDepth; ++k)
        {
            const solveScalar* __restrict__ historyColumn = window_->column(k).begin();

            for (label cell = 0; cell < nCells; ++cell)
            {
                // Difference before sketching to reduce cancellation.
                const solveScalar difference = historyColumn[cell] - initialPsi[cell];
                for (label j = 0; j < numDirections; ++j)
                {
                    basisData[j][cell] += solveScalar(Z(k, j))*difference;
                }
            }
        }
    }

    // Store the entry iterate after reading the old state, before correcting psi.
    pushState(psi);

    // matrixBasis = A*basis. Sequential Amul calls share interface MPI buffers.
    List<solveScalarField> matrixBasis(numDirections);
    for (label j = 0; j < numDirections; ++j)
    {
        matrixBasis[j].resize_nocopy(nCells);
        matrix_.Amul(matrixBasis[j], basis[j], interfaceBouCoeffs_, interfaces_, cmpt);
    }

    // Galerkin projects against the basis; least squares against A*basis.
    const List<solveScalarField>& testBasis = (projection_ == GALERKIN) ? basis : matrixBasis;

    const label packedSize = numDirections*(numDirections + 1)/2 + numDirections + 1;
    solveScalarField reductionBuffer(packedSize, Zero);

    // Pack local entries of the Gram matrix's upper triangle, row by row.
    label packedIndex = 0;
    for (label i = 0; i < numDirections; ++i)
    {
        for (label j = i; j < numDirections; ++j)
        {
            reductionBuffer[packedIndex++] = sumProd(testBasis[i], matrixBasis[j]);
        }
    }
    // Append the reduced RHS: testBasis^T rA.
    for (label i = 0; i < numDirections; ++i)
    {
        reductionBuffer[packedIndex++] = sumProd(testBasis[i], rA);
    }
    // The final entry is trace(A), used to identify the Galerkin matrix sign.
    reductionBuffer[packedIndex++] = solveScalar(sum(matrix_.diag()));

    // Sum all three parts across ranks in one collective.
    Foam::reduce
    (
        reductionBuffer.data(), int(reductionBuffer.size()), sumOp<solveScalar>(),
        UPstream::msgType(),
        comm
    );

    SquareMatrix<solveScalar> gram(numDirections, Zero);
    solveScalarField rhs(numDirections, Zero);

    // Unpack the global Gram matrix, mirroring its symmetric lower triangle.
    packedIndex = 0;
    for (label i = 0; i < numDirections; ++i)
    {
        for (label j = i; j < numDirections; ++j)
        {
            gram(i, j) = reductionBuffer[packedIndex];
            gram(j, i) = reductionBuffer[packedIndex];
            ++packedIndex;
        }
    }
    // Unpack the global RHS and read trace(A) from the remaining entry.
    for (label i = 0; i < numDirections; ++i)
    {
        rhs[i] = reductionBuffer[packedIndex++];
    }
    const solveScalar definitenessSign =
        (projection_ == GALERKIN) ? ((reductionBuffer[packedIndex] < 0) ? -1 : 1) : 1;

    // Threshold eigenvalue magnitudes to support either sign of definite matrix.
    EigenMatrix<solveScalar> eigenSystem(gram, true);
    const DiagonalMatrix<solveScalar>& eigenvalues = eigenSystem.EValsRe();
    const SquareMatrix<solveScalar>& eigenvectors = eigenSystem.EVecs();     // columns are vectors

    solveScalar maxEigenvalueMag = 0;
    for (label k = 0; k < numDirections; ++k)
    {
        maxEigenvalueMag = max(maxEigenvalueMag, mag(eigenvalues[k]));
    }
    if (maxEigenvalueMag <= 0)
    {
        return 0;
    }

    label retainedRank = 0;
    solveScalarField coefficients(numDirections, Zero);
    // Solve in the retained eigenmodes: sum_k v_k*(v_k^T rhs)/lambda_k.
    for (label k = 0; k < numDirections; ++k)
    {
        if (mag(eigenvalues[k]) <= truncTol_*maxEigenvalueMag || definitenessSign*eigenvalues[k] < 0)
        {
            // Discard small modes and modes with the wrong definiteness sign.
            continue;
        }

        solveScalar projectedRhs = 0;
        for (label i = 0; i < numDirections; ++i)
        {
            projectedRhs += eigenvectors(i, k)*rhs[i];
        }

        const solveScalar modeWeight = projectedRhs/eigenvalues[k];
        for (label i = 0; i < numDirections; ++i)
        {
            coefficients[i] += modeWeight*eigenvectors(i, k);
        }
        ++retainedRank;
    }

    if (retainedRank == 0)
    {
        return 0;
    }

    // Reject a correction with negative predicted improvement from roundoff.
    if (definitenessSign*sumProd(coefficients, rhs) < 0)
    {
        WarningInFunction
            << "subspace correction would raise the energy norm; skipping"
            << endl;
        return 0;
    }

    // Update in place: the solver retains pointers into psi and rA.
    solveScalar* __restrict__ psiData = psi.begin();
    solveScalar* __restrict__ residualData = rA.begin();

    for (label j = 0; j < numDirections; ++j)
    {
        const solveScalar coefficient = coefficients[j];
        if (coefficient == 0)
        {
            continue;
        }

        const solveScalar* __restrict__ basisColumn = basis[j].begin();
        const solveScalar* __restrict__ matrixBasisColumn = matrixBasis[j].begin();

        for (label cell = 0; cell < nCells; ++cell)
        {
            psiData[cell] += coefficient*basisColumn[cell];
            residualData[cell] -= coefficient*matrixBasisColumn[cell];
        }
    }

    return retainedRank;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::subspaceInitializer::pushState(const solveScalarField& psi) const
{
    // Advance state at most once per solve; require an fvMesh registry.
    if (pushed_ || !meshPtr_)
    {
        return;
    }
    pushed_ = true;

    const label timeIndex = meshPtr_->time().timeIndex();

    if (window_)
    {
        window_->update(psi, timeIndex);
    }

    // Keep every configured rate current, including unselected arms.
    forAll(ewmaSketches_, i)
    {
        if (ewmaSketches_[i])
        {
            ewmaSketches_[i]->update(psi, timeIndex);
        }
    }
}


Foam::label Foam::subspaceInitializer::update
(
    const solveScalarField& psi,
    const direction cmpt
) const
{
    window_ = nullptr;
    ewmaSketches_.clear();
    pushed_ = false;

    if (!meshPtr_ || (maxLenHistory_ <= 0 && decayRates_.empty()))
    {
        return -1;
    }

    // Check for time discontinuities before using stored state.
    const label timeIndex = meshPtr_->time().timeIndex();
    const label comm = matrix_.mesh().comm();

    // Find or create the mesh-owned state; constructors read checkpoints.
    DynamicList<subspaceStateBase*> states(1 + decayRates_.size());

    contextWindow* window = nullptr;

    if (maxLenHistory_ > 0)
    {
        // p and pFinal share state for the same field component.
        const word windowKey = contextWindow::makeKey(fieldName_, cmpt);

        // Derive the same per-object seed on every rank.
        const unsigned windowSeedOffset = string::hasher()(windowKey) & 0x00ffffffu;

        window = &contextWindow::get
        (
            *meshPtr_,
            windowKey,
            maxLenHistory_,
            psi.size(),
            randomSeed_ + label(windowSeedOffset),
            persistState_
        );

        window->checkRestart(timeIndex);
        states.append(window);
    }

    if (!decayRates_.empty())
    {
        ewmaSketches_.setSize(decayRates_.size(), nullptr);

        forAll(decayRates_, i)
        {
            const word sketchKey =
                EWMASketch::makeKey(fieldName_, cmpt, decayRates_[i]);
            const unsigned sketchSeedOffset = string::hasher()(sketchKey) & 0x00ffffffu;

            EWMASketch& sketch = EWMASketch::get
            (
                *meshPtr_,
                sketchKey,
                maxNumProbes_,
                decayRates_[i],
                psi.size(),
                randomSeed_ + label(sketchSeedOffset),
                persistState_
            );

            sketch.checkRestart(timeIndex);
            ewmaSketches_[i] = &sketch;
            states.append(&sketch);
        }
    }

    // Validate restored state once across ranks. Matching depths and timestamps
    // keep projection collectives and restart resets synchronized.
    bool allAgreed = true;
    forAll(states, k)
    {
        allAgreed = allAgreed && states[k]->agreed();
    }

    if (!allAgreed)
    {
        // Pack [depth, -depth, timestamp, -timestamp] for one min reduction.
        // Negative depths mark cell-count mismatches; -min(-x) gives max(x).
        labelList stateBounds(4*states.size());
        forAll(states, k)
        {
            const bool sizeMatches = (states[k]->nCells() == psi.size());
            stateBounds[4*k]     = sizeMatches ?  states[k]->depth() : -1;
            stateBounds[4*k + 1] = sizeMatches ? -states[k]->depth() : -labelMax;
            stateBounds[4*k + 2] =  states[k]->stamp();
            stateBounds[4*k + 3] = -states[k]->stamp();
        }

        Foam::reduce
        (
            stateBounds.data(), int(stateBounds.size()), minOp<label>(),
            UPstream::msgType(),
            comm
        );

        forAll(states, k)
        {
            const bool sizeMatches = (stateBounds[4*k] >= 0);

            if
            (
                !sizeMatches
             || stateBounds[4*k] != -stateBounds[4*k + 1]
             || stateBounds[4*k + 2] != -stateBounds[4*k + 3]
            )
            {
                // Clear inconsistent histories on every rank.
                states[k]->clear();
            }

            if (sizeMatches)
            {
                states[k]->setAgreed();
            }
            else
            {
                // Exclude objects whose stored cell count does not match psi.
                if (states[k] == window)
                {
                    window = nullptr;
                }
                forAll(ewmaSketches_, i)
                {
                    if (ewmaSketches_[i] == states[k])
                    {
                        ewmaSketches_[i] = nullptr;
                    }
                }
            }
        }
    }

    // Borrow the registry-owned window for this solve.
    window_ = window;

    label depth = -1;
    if (window_)
    {
        depth = window_->depth();
    }
    forAll(ewmaSketches_, i)
    {
        if (ewmaSketches_[i])
        {
            depth = max(depth, ewmaSketches_[i]->depth());
        }
    }

    return depth;
}


Foam::label Foam::subspaceInitializer::initialize
(
    const direction cmpt,
    const label lenHistory,
    const label numProbes,
    const scalar decayRate,
    const solveScalar normFactor,
    solveScalarField& psi,
    solveScalarField& rA,
    solverPerformance& solverPerf
) const
{
    // Each initialization request selects either a window or a stream.
    checkExclusive(lenHistory, decayRate);

    // Select the state validated by update() for this solve.
    const EWMASketch* sketch = nullptr;

    bool canCorrect = (numProbes > 0);

    if (canCorrect && decayRate > 0)
    {
        const label i = rateIndex(decayRate);

        if (i < 0)
        {
            // Undeclared rates are errors; declared but empty sketches are valid.
            FatalErrorInFunction
                << "no sketch for decayRate " << decayRate
                << "; the rates this initializer maintains are " << decayRates_
                << ".  Every rate an arm can name must be passed to the"
                << " constructor, or its sketch is never updated."
                << exit(FatalError);
        }

        sketch = (i < ewmaSketches_.size()) ? ewmaSketches_[i] : nullptr;
        canCorrect = (sketch && sketch->depth() > 0);
    }
    else if (canCorrect)
    {
        canCorrect = (window_ && window_->depth() > 0);
    }

    // Advance state even when no correction is possible.
    // If correct() already stored the entry psi, pushState() does nothing.
    const label rank =
        canCorrect
      ? correct(cmpt, lenHistory, numProbes, sketch, psi, rA)
      : 0;

    pushState(psi);

    if (rank <= 0)
    {
        return -1;
    }

    // Preserve initialResidual so the relative-convergence target is unchanged.
    solverPerf.finalResidual() =
        gSumMag(rA, matrix_.mesh().comm())/normFactor;

    return rank;
}

// ************************************************************************* //
