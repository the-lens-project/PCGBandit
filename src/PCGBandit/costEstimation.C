/*---------------------------------------------------------------------------*\
                     Class costEstimation Implementation
\*---------------------------------------------------------------------------*/

#include "costEstimation.H"
#include "configurationSpace.H"
#include "GAMGAgglomeration.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    // --- Deterministic cost estimation constants (calibrated on 4-core 2x-resolved pitzDaily)
    const scalar ICTC_SETUP_WEIGHT = 10.0;
    const scalar GAMG_ASSEMBLY_WEIGHT = 6.0;
    const scalar GAMG_AGGLOMERATION_WEIGHT = 45.0;
    const scalar ICTC_SMOOTHER_SETUP_WEIGHT = 18.0;
    const scalar ICTC_SMOOTHER_APPLY_WEIGHT = 2.5;
    const scalar ICTC_SMOOTHER_BASELINE = 10.0;
    const scalar COARSE_CG_WEIGHT = 0.08;
    const scalar COMM_EVENT_FLOPS = 1.0e4;
    const scalar DIRECT_LU_WEIGHT = 0.3;
    const scalar SUBSPACE_WEIGHT = 0.24;
    const scalar SUBSPACE_GRAM_WEIGHT = 0.72;

    //- Allreduce latency grows as log2 of the rank count on any tree
    //  implementation; a nearest-neighbour halo exchange does not.
    static inline scalar log2Procs(const label nProcs)
    {
        return Foam::log(scalar(nProcs))/Foam::log(scalar(2));
    }

    static inline bool isGAMG(const word& p)
    {
        return p == "GAMG" || p == "FGAMG";
    }

    static scalar smootherApplyCost
    (
        const word& smoother,
        const label nnzL,
        const label nCells,
        const label nSweeps,
        const scalar fillFactor = 1.0,
        const scalar nnzBaseline = 0.0
    )
    {
        if (nSweeps <= 0) {
            return 0.0;
        }

        scalar c = 0.0;
        const scalar matvec = scalar(2 * nnzL + nCells);
        if (smoother.find("ICTC") != string::npos) {
            c = matvec + (ICTC_SMOOTHER_APPLY_WEIGHT * fillFactor + nnzBaseline) * scalar(nnzL) + scalar(2 * nCells);
            if (smoother.find("GaussSeidel") != string::npos) {
                c += matvec;
            }
        } else if (smoother == "symGaussSeidel") {
            c = scalar(4 * nnzL + 2 * nCells);
        } else {
            if (smoother == "GaussSeidel"
                || smoother == "DICGaussSeidel"
                || smoother.find("SOR_") == 0
                || smoother.find("DICSOR_") == 0) {
                c += scalar(2 * nnzL + nCells);
            }
            if (smoother == "DIC"
                || smoother == "DICGaussSeidel"
                || smoother.find("DICSOR_") == 0) {
                c += scalar(4 * nnzL + nCells);
            }
            if (c == 0.0) {
                c = matvec;
            }
        }
        return c * scalar(nSweeps);
    }

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::costEstimation::communicationCostEstimate
(
    const dictionary& armDict,
    const word& preconditioner
) const
{
    const label comm = matrix_.mesh().comm();

    if (Pstream::nProcs(comm) <= 1 || !isGAMG(preconditioner)) {
        return 0.0;
    }

    const GAMGAgglomeration& agg = GAMGAgglomeration::New(matrix_, armDict);
    const label L = agg.size();

    // --- Each level: restrict/prolong + smoother halo exchange + scale reduction
    scalar perVcycle = scalar(3 * L);

    // --- DIC-PCG iterations (n) if directSolveCoarsest is off
    if (L > 0 && !Switch(armDict.getOrDefault<word>("directSolveCoarsest", "no"))) {
        const scalar ncG = returnReduce(scalar(max(agg.nCells(L - 1), label(1))), sumOp<scalar>(), UPstream::msgType(), comm);
        perVcycle += COARSE_CG_WEIGHT * ncG;
    }

    return COMM_EVENT_FLOPS * perVcycle * scalar(armDict.getOrDefault<label>("nVcycles", 2));
}

Foam::scalar Foam::costEstimation::perIterationCostEstimate
(
    const dictionary& armDict,
    const word& preconditioner
) const
{
    const label nCells = matrix_.diag().size();
    const label nnzL = matrix_.lower().size();

    // --- ICTC scattered-triangular-solve baseline (density-independent, 3D only).
    const scalar ictcBaseline = (nGeometricD_ >= 3) ? ICTC_SMOOTHER_BASELINE: 0.0;

    // --- One CG step has a matvec (2 * nnzL + nCells) and five vector operations
    const scalar cgStep = scalar(2 * nnzL + 6 * nCells);

    // --- And, in parallel, one halo exchange (the Amul) and three allreduces
    const label comm = matrix_.mesh().comm();
    const label nProcs = Pstream::nProcs(comm);
    const scalar cgComm = (nProcs > 1) ? COMM_EVENT_FLOPS*(1.0 + 3.0*log2Procs(nProcs)) : 0.0;

    if (preconditioner == "ICTC") {
        return returnReduce(cgStep + scalar(2 * (debug::controlDict().get<label>("ICTC_NNZ") + nCells)), maxOp<scalar>(), UPstream::msgType(), comm) + cgComm;
    }
    if (preconditioner == "DIC") {
        return returnReduce(cgStep + scalar(4 * nnzL + nCells), maxOp<scalar>(), UPstream::msgType(), comm) + cgComm;
    }

    const label nPreSweeps  = armDict.getOrDefault<label>("nPreSweeps", 0);
    const label nPostSweeps = armDict.getOrDefault<label>("nPostSweeps", 2);
    const label maxPreSweeps = 4;
    const label maxPostSweeps = 4;
    const label preSweepsLevelMultiplier = 1;
    const label postSweepsLevelMultiplier = 1;
    const label nVcycles = armDict.getOrDefault<label>("nVcycles", 2);
    const label nFinestSweeps = armDict.getOrDefault<label>("nFinestSweeps", 2);
    const word smoother = armDict.get<word>("smoother");
    const bool interpolateCorrection = Switch(armDict.getOrDefault<word>("interpolateCorrection", "no"));
    const bool directSolveCoarsest = Switch(armDict.getOrDefault<word>("directSolveCoarsest", "no"));

    const bool scaleCorrection = matrix_.symmetric();
    const GAMGAgglomeration& agg = GAMGAgglomeration::New(matrix_, armDict);
    const label L = agg.size();

    // --- Compute ratio of ICTC smoother fills to agglomeration matrix fills
    scalar fillFactor = 1.0;
    if (smoother.contains("ICTC")) {
        scalar nnzAgg = scalar(nnzL);
        // Built-in GAMG still constructs an unused final coarse smoother.
        const label nFactorLevels = preconditioner == "FGAMG" ? L - 1 : L;
        for (label i = 0; i < nFactorLevels; i++) {
            nnzAgg += scalar(agg.nFaces(i));
        }
        fillFactor = scalar(debug::controlDict().get<label>("ICTC_SMOOTHER_NNZ")) / nnzAgg;
    }

    scalar perVcycle = 0.0;

    // --- Finest level: smoothing + prolongation + corrections
    perVcycle += smootherApplyCost(smoother, nnzL, nCells, nFinestSweeps, fillFactor, ictcBaseline);
    perVcycle += scalar(2 * nCells);
    if (interpolateCorrection) {
        perVcycle += scalar(2 * nnzL + 3 * nCells);
    }
    if (scaleCorrection) {
        perVcycle += scalar(2 * nnzL + 4 * nCells);
    }

    // --- Per vcycle restriction + smoothing + residual computation + prolongation
    for (label i = 0; i < L; i++) {

        // Transfers span all coarse levels; smoothing excludes the final matrix.
        const label nc  = agg.nCells(i);
        perVcycle += scalar(4 * nc);
        if (i == L - 1) {
            continue;
        }
        const label nf  = agg.nFaces(i);
        const label nPre  = (nPreSweeps  > 0)
            ? min(nPreSweeps  + preSweepsLevelMultiplier  * i, maxPreSweeps)  : 0;
        const label nPost = (nPostSweeps > 0)
            ? min(nPostSweeps + postSweepsLevelMultiplier * i, maxPostSweeps) : 0;

        perVcycle += smootherApplyCost(smoother, nf, nc, nPre + nPost, fillFactor, ictcBaseline);
        if (nPre > 0) {
            perVcycle += scalar(2 * nf + nc);
        }
        if (interpolateCorrection) {
            perVcycle += scalar(2 * nf + 3 * nc);
        }
        if (scaleCorrection) {
            perVcycle += scalar(2 * nf + 4 * nc);
        }
    }

    // --- Coarsest-level direct (n^2) or Poisson-like DIC-PCG (sqrt(n)) solve
    if (L > 0) {
        const scalar nc = max(scalar(agg.nCells(L - 1)), scalar(1));
        const label  nf = agg.nFaces(L - 1);
        if (directSolveCoarsest) {
            const scalar ncG = returnReduce(scalar(agg.nCells(L - 1)), sumOp<scalar>(), UPstream::msgType(), comm);
            perVcycle += 2.0 * ncG * ncG;
        } else {
            perVcycle += sqrt(nc) * scalar(6 * nf + 3 * label(nc));
        }
    }

    return returnReduce(cgStep + perVcycle * scalar(nVcycles), maxOp<scalar>(), UPstream::msgType(), comm)
         + communicationCostEstimate(armDict, preconditioner) + cgComm;
}


Foam::scalar Foam::costEstimation::setupCostEstimate
(
    const dictionary& armDict,
    const word& preconditioner
) const
{
    const scalar pICE = perIterationCostEstimate(armDict, preconditioner);
    const label comm = matrix_.mesh().comm();

    if (preconditioner == "ICTC") {
        return ICTC_SETUP_WEIGHT * pICE;
    }
    if (!isGAMG(preconditioner)) {
        return pICE;
    }

    // --- GAMG: rebuild agglomeration (unless cached) and smoothers every solve
    const GAMGAgglomeration& agg = GAMGAgglomeration::New(matrix_, armDict);
    scalar hierarchy = scalar(2 * matrix_.lower().size() + matrix_.diag().size());
    for (label i = 0; i < agg.size(); i++) {
        hierarchy += scalar(2 * agg.nFaces(i) + agg.nCells(i));
    }
    const bool cached = armDict.getOrDefault<label>("cacheAgglomeration", 1);
    scalar weight = GAMG_ASSEMBLY_WEIGHT + (cached ? 0.0 : GAMG_AGGLOMERATION_WEIGHT);
    scalar setup = weight * hierarchy;

    // --- FGAMG: precomputes smoother factors once per solve
    const word smoother = armDict.getOrDefault<word>("smoother", "");
    const word coarsestSmoother = armDict.getOrDefault<word>("coarsestSmoother", smoother);
    if (smoother.find("ICTC") != string::npos
     || coarsestSmoother.find("ICTC") != string::npos) {
        const label smootherNNZ = debug::controlDict().get<label>("ICTC_SMOOTHER_NNZ");
        setup += ICTC_SMOOTHER_SETUP_WEIGHT * scalar(smootherNNZ);
    }

    // --- Dense LU factorisation (n^3) if directSolveCoarsest is on
    if (agg.size() > 0 && Switch(armDict.getOrDefault<word>("directSolveCoarsest", "no"))) {
        const scalar ncG = returnReduce(scalar(agg.nCells(agg.size() - 1)), sumOp<scalar>(), UPstream::msgType(), comm);
        setup += DIRECT_LU_WEIGHT * ncG * ncG * ncG;
    }

    return returnReduce(setup, maxOp<scalar>(), UPstream::msgType(), comm);
}



Foam::scalar Foam::costEstimation::subspaceCostEstimate
(
    const dictionary& armDict
) const
{

    const label lenHistory = armDict.getOrDefault<label>("lenHistory", 0);
    const label numProbes = armDict.getOrDefault<label>("numProbes", 4);
    const scalar decayRate = armDict.getOrDefault<scalar>("decayRate", 0);

    if (subspaceOff(armDict)) {
        return 0.0;                 // the off arm must cost exactly zero
    }

    const scalar n = scalar(matrix_.diag().size());
    const label nnzL = matrix_.lower().size();

    // Forming Omega is the only term the two sketches disagree on. The context
    // window's push and the EWMA stream's rank-1 update are NOT charged, as
    // they run every solve regardless of arm drawn.
    const scalar omegaCost =
        (decayRate > 0)
      ? 3.0*n*scalar(numProbes)
      : n*scalar(lenHistory)*(1.0 + 2.0*scalar(numProbes));

    // leastSquares dots W against itself where galerkin dots it against Omega
    const scalar fixedCost = galerkin_ ? 3.0*n : 2.0*n;

    scalar cost =
        SUBSPACE_WEIGHT*
        (
            omegaCost                                                   // Omega
          + scalar(numProbes)*scalar(4*nnzL + label(n))                 // W = A*Omega
          + 4.0*n*scalar(numProbes)                                     // psi += Omega y ; rA -= W y
          + fixedCost                                                   // gSumMag refresh
          + 10.0*scalar(numProbes)*scalar(numProbes)*scalar(numProbes)  // eigensolve
        )
      + SUBSPACE_GRAM_WEIGHT*n*scalar(numProbes)*scalar(numProbes + 3); // G and t = L^T r

    const label comm = matrix_.mesh().comm();
    const label nProcs = Pstream::nProcs(comm);
    if (nProcs > 1) {
        // Everything above is a per-rank count reduced with maxOp, so only
        // collectives are left. The numProbes halo exchanges inside Amul are
        // nearest-neighbour, while the fused reduce and the finalResidual
        // gSumMag are allreduces, latency log2(nProcs).
        cost += COMM_EVENT_FLOPS*(scalar(numProbes) + 2.0*log2Procs(nProcs));
    }

    return returnReduce(cost, maxOp<scalar>(), UPstream::msgType(), comm);
}


Foam::scalar Foam::costEstimation::totalCostEstimate
(
    const dictionary& armDict,
    const label nIterations,
    const label maxIter,
    const label backstop
) const
{
    scalar cost = subspaceCostEstimate(armDict);

    if (nIterations == 0) {
        return cost;
    }

    const word preconditioner = armDict.get<word>("preconditioner");
    const scalar pICE = perIterationCostEstimate(armDict, preconditioner);

    cost += setupCostEstimate(armDict, preconditioner);

    label backstopIter = maxIter;
    if (backstop == -1) {
        backstopIter = label(scalar(backstopIter) * perIterationCostEstimate(armDict, "DIC") / pICE);
    }
    if (nIterations > backstopIter) {
        cost += pICE * scalar(backstopIter)
                + perIterationCostEstimate(armDict, "DIC") * scalar(nIterations - backstopIter + label(preconditioner != "DIC"));
    } else {
        cost += pICE * scalar(nIterations);
    }

    return returnReduce
    (
        cost, maxOp<scalar>(), UPstream::msgType(), matrix_.mesh().comm()
    );
}

// ************************************************************************* //
