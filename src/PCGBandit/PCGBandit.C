/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

//
#include "PCGBandit.H"
#include "PrecisionAdaptor.H"

#include "clockValue.H"
#include "fvMesh.H"
#include "GAMGAgglomeration.H"
#include "Pstream.H"
#include "Random.H"

//#define PCGB_DEBUG
//#define DUMP_ABSOL

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    clockValue PCGTime = clockValue();

    defineTypeNameAndDebug(PCGBandit, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PCGBandit>
        addPCGBanditSymMatrixConstructorToTable_;

    Random rndGen;

    dictionary preconditionerDict;
    dictionary subDict;
    dictionary learningDicts;

    // --- Configuration space specification
    List<word> wordGAMGParams = {"smoother", "agglomerator", "directSolveCoarsest"};
    List<word> labelGAMGParams = {"nCellsInCoarsestLevel", "mergeLevels", "nPreSweeps", "nPostSweeps", "nFinestSweeps", "nVcycles"};

    boolField wordGAMGTune(wordGAMGParams.size(), false);
    boolField labelGAMGTune(labelGAMGParams.size(), false);

    labelField wordGAMGSizes(wordGAMGParams.size());
    labelField labelGAMGSizes(labelGAMGParams.size());

    #ifdef PCGB_DEBUG
    dictionary initialDiags;
    dictionary initialLowers;
    dictionary initialTargets;
    dictionary previousDiags;
    dictionary previousLowers;
    dictionary previousTargets;
    #endif

    #ifdef DUMP_ABSOL
    #include "initializeDumping.H"
    #endif

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PCGBandit::PCGBandit
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    ),
    numDroptols_(0),
    dGAMG_(0),
    cacheAgglomeration_(false)
{

    // --- Contextual information specification
    word preconditioner = solverControls.get<word>("preconditioner");
    const fvMesh& mesh = dynamicCast<const fvMesh>(matrix.mesh());
    if (preconditioner == "separate") {
        banditName_ = mesh.name() + "." + fieldName;
    } else if (preconditioner == "joint") {
        banditName_ = "joint";
    } else {
        banditName_ = preconditioner;
    }
    if (relTol_ == 0.0 and Switch(solverControls.getOrDefault<word>("residualContext", "no"))) {
        banditName_ += "Final";
    }

    // --- Learning algorithm specification
    lossEstimator_ = solverControls.getOrDefault<word>("lossEstimator", "IW");
    deterministic_ = Switch(solverControls.getOrDefault<word>("deterministic", "no"));
    seed_ = mesh.time().controlDict().getOrDefault<label>("randomSeed", 0);
    backstop_ = solverControls.getOrDefault<label>("backstop", -1);
    static_ = label(solverControls.getOrDefault<label>("static", -1));
    randomUniform_ = Switch(solverControls.getOrDefault<word>("randomUniform", "no"));
    banditAlgorithm_ = solverControls.getOrDefault<word>("banditAlgorithm", "TsallisINF");

    dictionary& learningDict = learningDicts.subDictOrAdd(banditName_);
    dictionary& GAMGOptions = learningDict.subDictOrAdd("GAMGOptions");

    if (!learningDict.found("dGAMG")) {
        // --- Configuration space initialization
        if (Pstream::myProcNo() == 0)
        {
            rndGen.reset(seed_);
        }

        // --- ICTC specification
        maxLogDroptol_ = solverControls.getOrDefault<scalar>("maxLogDroptol", -0.5);
        minLogDroptol_ = solverControls.getOrDefault<scalar>("minLogDroptol", -4.0);
        numDroptols_ = solverControls.getOrDefault<label>("numDroptols", 0);

    
        // Label Parameters
        auto addLabelListOrDefault = [&](const word& fvName, const word& keyName, const List<label>& defList)
        {
        List<label> vals = solverControls.getOrDefault<List<label>>(fvName, defList);
            if (vals.empty())
                vals = defList;
            GAMGOptions.set(keyName, vals);
        };

        addLabelListOrDefault("nPreSweepsOptions", "nPreSweeps", {0, 2});
        addLabelListOrDefault("nPostSweepsOptions", "nPostSweeps", {1, 2});
        addLabelListOrDefault("nFinestSweepsOptions", "nFinestSweeps", {2});
        addLabelListOrDefault("nCellsInCoarsestLevelOptions", "nCellsInCoarsestLevel", {10, 100, 1000});
        addLabelListOrDefault("mergeLevelsOptions", "mergeLevels", {1, 2});
        addLabelListOrDefault("nVcyclesOptions", "nVcycles", {1, 2});

        // Word Parameters (string options)
        auto addWordListOrDefault = [&](const word& fvName, const word& keyName, const List<word>& defList)
        {
            List<word> vals = solverControls.getOrDefault<List<word>>(fvName, defList);
            if (vals.empty())
               vals = defList;
            GAMGOptions.set(keyName, vals);
        };

        addWordListOrDefault("smootherOptions", "smoother", {"GaussSeidel", "DIC", "DICGaussSeidel", "symGaussSeidel"});
        addWordListOrDefault("agglomeratorOptions", "agglomerator", {"faceAreaPair", "algebraicPair"});
        addWordListOrDefault("directSolveCoarsestOptions", "directSolveCoarsest", {"no", "yes"});

        // --- GAMG specification
        dGAMG_ = label(Switch(solverControls.getOrDefault<word>("GAMGTune", "no")));

        // labelGAMG
        List<List<label>> labelGAMGlist(labelGAMGParams.size());
        for (label j = 0; j < labelGAMGParams.size(); ++j) {
            const word optKey = labelGAMGParams[j] + "Options";
            List<label> vals;
            if (solverControls.found(optKey))
            {
                vals = solverControls.getOrDefault<List<label>>(optKey, List<label>());
            }
            labelGAMGlist[j] = vals;
            labelGAMGSizes[j] = GAMGOptions.get<List<label>>(labelGAMGParams[j]).size();
            if (Switch(solverControls.getOrDefault<word>(labelGAMGParams[j] + "Tune", "no")))
            {
                labelGAMGTune[j] = true;
                dGAMG_ = max(dGAMG_, 1) * labelGAMGSizes[j];
            }
        }

        // wordGAMG
        List<List<word>> wordGAMGlist(wordGAMGParams.size());
        for (label j = 0; j < wordGAMGParams.size(); ++j) {
            const word optKey = wordGAMGParams[j] + "Options";
            List<word> vals;
            if (solverControls.found(optKey))
            {
                vals = solverControls.getOrDefault<List<word>>(optKey, List<word>());
            }
            wordGAMGlist[j] = vals;
            wordGAMGSizes[j] = GAMGOptions.get<List<word>>(wordGAMGParams[j]).size();
            if (Switch(solverControls.getOrDefault<word>(wordGAMGParams[j] + "Tune", "no")))
            {
                wordGAMGTune[j] = true;
                dGAMG_ = max(dGAMG_, 1) * wordGAMGSizes[j];
            }
        }

        learningDict.set("labelGAMG", labelGAMGlist);
        learningDict.set("wordGAMG", wordGAMGlist);
        learningDict.set("dGAMG", dGAMG_);
    } else {
        // *** Load the previously computed dimension back in ***
        dGAMG_ = learningDict.get<label>("dGAMG");
    }

    {
        // --- No cache agglomeration if GAMG is tuned
        cacheAgglomeration_ = (static_ > -1 || !(wordGAMGTune[1] || labelGAMGTune[0] || labelGAMGTune[1]));
        if (cacheAgglomeration_)
        {
            cacheAgglomeration_ = Switch(solverControls.getOrDefault<word>("cacheAgglomeration", "yes"));
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PCGBandit::queryLearner
(
    const scalar initialResidual
) const
{
    label i = static_;
    dictionary& learningDict = learningDicts.subDictOrAdd(banditName_);
    dictionary& GAMGOptions = learningDict.subDict("GAMGOptions");

    if (i == -1) {

        if (Pstream::myProcNo() == 0) {

            label d = numDroptols_ + 1 + dGAMG_;
            #ifdef PCGB_DEBUG
            Info<< "PCGBandit INFO:" << endl;
            #endif

            if (randomUniform_) {
                i = floor(scalar(d) * rndGen.sample01<scalar>());
            } else if (banditAlgorithm_ == "ThompsonSampling") {
                #include "ThompsonSampling.H"
            } else {
                #include "TsallisINF.H"
            }
        }

        Pstream::scatter(i);

        #ifdef PCGB_DEBUG
        Info<< "\tselected: ";
    } else {

        Info<< "Static INFO: ";
        #endif
    }

    if (i == numDroptols_) {
        subDict.set("preconditioner", "DIC");
        #ifdef PCGB_DEBUG
        Info<< "preconditioner=DIC" << endl;
        #endif
    } else if (i > numDroptols_) {
        i -= numDroptols_ + 1;
        subDict.set("preconditioner", "GAMG");
        #ifdef PCGB_DEBUG
        Info<< "preconditioner=GAMG";
        #endif
        if (not wordGAMGTune[0]) {
            subDict.set("smoother", "DICGaussSeidel");
            #ifdef PCGB_DEBUG
            Info<< ", smoother=DICGaussSeidel";
            #endif
        }
        subDict.set("cacheAgglomeration", cacheAgglomeration_);
        #ifdef PCGB_DEBUG
        Info<< ", cacheAgglomeration=" << cacheAgglomeration_;
        #endif
        for (label j = labelGAMGParams.size()-1; j >= 0; j--) {
            if (labelGAMGTune[j]) {
                word param = labelGAMGParams[j];
                label size = labelGAMGSizes[j];
                label option = GAMGOptions.get<List<label>>(param)[i % size];
                i /= size;
                subDict.set(param, option);
                #ifdef PCGB_DEBUG
                Info<< ", " << param << "=" << option;
                #endif
            }
        }
        for (label j = wordGAMGParams.size()-1; j >= 0; j--) {
            if (wordGAMGTune[j]) {
                word param = wordGAMGParams[j];
                label size = wordGAMGSizes[j];
                word option = GAMGOptions.get<List<word>>(param)[i % size];
                i /= size;
                subDict.set(param, option);
                #ifdef PCGB_DEBUG
                Info<< ", " << param << "=" << option;
                #endif
            }
        }
        #ifdef PCGB_DEBUG
        Info<< endl;
        #endif
    } else {
        subDict.set("preconditioner", "ICTC");
        scalar droptol;
        if (i == 0) {
            droptol = pow(10.0, minLogDroptol_);
        } else {
            droptol = pow(10.0, minLogDroptol_+(maxLogDroptol_-minLogDroptol_)*scalar(i)/scalar(numDroptols_-1));
        }
        subDict.set("droptol", droptol);
        #ifdef PCGB_DEBUG
        Info<< "preconditioner=ICTC, droptol=" << droptol << endl;
        #endif
    }
    preconditionerDict.set("preconditioner", subDict);
}


Foam::scalar Foam::PCGBandit::perIterationCostEstimate
(
    const word preconditioner
) const
{

    label nCells = matrix_.diag().size();
    label nnz = matrix_.lower().size();
    label cost = 2 * nnz + 6 * nCells;
    if (preconditioner == "ICTC") {
        nnz = Foam::debug::controlDict().get<label>("ICTC_NNZ");
        return returnReduce(scalar(cost + 2 * (nnz + nCells)), maxOp<scalar>());
    } else {
        if (preconditioner == "DIC") {
            return returnReduce(scalar(cost + 4 * nnz + nCells), maxOp<scalar>());
        }
    }

    label nPreSweeps = subDict.getOrDefault<label>("nPreSweeps", 0);
    label nPostSweeps = subDict.getOrDefault<label>("nPostSweeps", 2);
    label maxPreSweeps = 4;
    label maxPostSweeps = 4;
    label preSweepsLevelMultiplier = 1;
    label postSweepsLevelMultiplier = 1;
    label nVcycles = subDict.getOrDefault<label>("nVcycles", 2);
    label nSweeps;
    word smoother = subDict.get<word>("smoother");
    const GAMGAgglomeration *agglomeration = &GAMGAgglomeration::New(matrix_, subDict);

    cost += 2 * nnz + nCells;
    for (label i = 0; i <= agglomeration->size(); i++) {

        if (i > 0) {
            nCells = agglomeration->nCells(i-1);
            nnz = agglomeration->nFaces(i-1);
            cost += (2 * nnz + nCells) * nVcycles;
            nSweeps = 0;
            if (nPreSweeps > 0) {
                nSweeps += min(nPreSweeps+preSweepsLevelMultiplier*(i-1), maxPreSweeps);
            }
            if (nPostSweeps > 0) {
                nSweeps += min(nPostSweeps+postSweepsLevelMultiplier*(i-1), maxPostSweeps);
            }
        } else {
            nSweeps = subDict.getOrDefault<label>("nFinestSweeps", 2);
        }

        nSweeps *= nVcycles;
        if (smoother == "symGaussSeidel") {
            cost += (4 * nnz + 2 * nCells) * nSweeps;
        } else {
            if (smoother == "GaussSeidel" || smoother == "DICGaussSeidel") {
                cost += (2 * nnz + nCells) * nSweeps;
            }
            if (smoother == "DIC" || smoother == "DICGaussSeidel") {
                cost += (4 * nnz + nCells) * nSweeps;
            }
        }
    }

    return returnReduce(scalar(cost), maxOp<scalar>());
}

Foam::scalar Foam::PCGBandit::totalCostEstimate
(
    const label nIterations
) const
{

    word preconditioner = subDict.get<word>("preconditioner");
    scalar pICE = perIterationCostEstimate(preconditioner);
    scalar cost;

    if (preconditioner == "GAMG") {
        cost = scalar(2 * matrix_.lower().size() + matrix_.diag().size());
    } else if (preconditioner == "ICTC") {
        cost = 10.0 * pICE;
    } else {
        cost = pICE;
    }

    label backstopIter = maxIter_;
    if (backstop_ == -1) {
        backstopIter = label(scalar(backstopIter) * perIterationCostEstimate("DIC") / pICE);
    }
    if (nIterations > backstopIter) {
        cost += pICE * scalar(backstopIter) 
                + perIterationCostEstimate("DIC") * scalar(nIterations - backstopIter + label(preconditioner != "DIC"));
    } else {
        cost += pICE * scalar(nIterations);
    }

    return returnReduce(cost, maxOp<scalar>());
}

Foam::solverPerformance Foam::PCGBandit::scalarSolve
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt
) const
{

    #ifdef DUMP_ABSOL
    #include "startDump.H"
    #endif

    #ifdef PCGB_DEBUG
    solveScalarField currentDiag = solveScalarField(matrix_.diag());
    solveScalarField currentLower = solveScalarField(matrix_.lower());
    solveScalarField currentTarget = solveScalarField(source);
    Info<< "\tsq Frob norm of current diag: " << gSum(sqr(currentDiag)) << endl;
    Info<< "\tsq Frob norm of current lower: " << gSum(sqr(currentLower)) << endl;
    Info<< "\tsq Eucl norm of current target: " << gSum(sqr(currentTarget)) << endl;
    if (initialDiags.found(banditName_)) {
        Info<< "\tsq Frob dist from initial diag: " << gSum(sqr(currentDiag - initialDiags.get<solveScalarField>(banditName_))) << endl;
        Info<< "\tsq Frob dist from initial lower: " << gSum(sqr(currentLower - initialLowers.get<solveScalarField>(banditName_))) << endl;
        Info<< "\tsq Eucl norm from initial target: " << gSum(sqr(currentTarget - initialTargets.get<solveScalarField>(banditName_))) << endl;
        Info<< "\tsq Frob dist from previous diag: " << gSum(sqr(currentDiag - previousDiags.get<solveScalarField>(banditName_))) << endl;
        Info<< "\tsq Frob dist from previous lower: " << gSum(sqr(currentLower - previousLowers.get<solveScalarField>(banditName_))) << endl;
        Info<< "\tsq Eucl norm from previous target: " << gSum(sqr(currentTarget - previousTargets.get<solveScalarField>(banditName_))) << endl;
    } else {
        initialDiags.add<solveScalarField>(banditName_, currentDiag);
        initialLowers.add<solveScalarField>(banditName_, currentLower);
        initialTargets.add<solveScalarField>(banditName_, currentTarget);
    }
    previousDiags.set<solveScalarField>(banditName_, currentDiag);
    previousLowers.set<solveScalarField>(banditName_, currentLower);
    previousTargets.set<solveScalarField>(banditName_, currentTarget);
    #endif

    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );
    clockValue preconstructTime;
    clockValue iterationTime;
    clockValue learningTime;
    clockValue solverTime = clockValue::now();

    label maxIter = maxIter_;
    label backstopIter = maxIter_;
    dictionary backstopDict;
    autoPtr<lduMatrix::preconditioner> preconPtr;

    label nCells = psi.size();
    solveScalar* __restrict__ psiPtr = psi.begin();

    solveScalarField pA(nCells);
    solveScalar* __restrict__ pAPtr = pA.begin();

    solveScalarField wA(nCells);
    solveScalar* __restrict__ wAPtr = wA.begin();

    solveScalar wArA = solverPerf.great_;
    solveScalar wArAold = wArA;

    // --- Calculate A.psi
    matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    solveScalarField rA(source - wA);
    solveScalar* __restrict__ rAPtr = rA.begin();

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        true
    );

    // --- Calculate normalisation factor
    solveScalar normFactor = this->normFactor(psi, source, wA, pA);

    if ((log_ >= 2) || (lduMatrix::debug >= 2))
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())
       /normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    for (label backstop = 0; backstop <= label(backstop_ != 0); backstop++) {

        // --- Check convergence, solve if not converged
        if
        (
            minIter_ > 0
         || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
        )
        {

            // --- Select and construct the preconditioner
            if (backstop) {

                // --- Revert to backstopping preconditioner
                preconstructTime += preconstructTime.now();
                if (subDict.get<word>("preconditioner") != "DIC") {
                    preconPtr = lduMatrix::preconditioner::New(*this, backstopDict);
                }
                preconstructTime -= clockValue::now();
                iterationTime += clockValue::now();
            } else {
                // --- Get preconditioner from learning algorithm
                learningTime = learningTime.now();
                queryLearner(solverPerf.initialResidual());

                learningTime -= clockValue::now();
                preconstructTime = preconstructTime.now();

                preconPtr = lduMatrix::preconditioner::New(*this, preconditionerDict);

                // --- Default backstop iteration computed via a cost estimate ratio
                if (backstop_ == -1) {
                    backstopIter = label(scalar(maxIter)
                                         * perIterationCostEstimate("DIC")
                                         / perIterationCostEstimate(subDict.get<word>("preconditioner")));
                    maxIter = backstopIter;
                }
                preconstructTime -= clockValue::now();
                iterationTime = iterationTime.now();
            }

            // --- Solver iteration
            do
            {

                // --- Store previous wArA
                wArAold = wArA;

                // --- Precondition residual
                preconPtr->precondition(wA, rA, cmpt);

                // --- Update search directions:
                wArA = gSumProd(wA, rA, matrix().mesh().comm());

                if (solverPerf.nIterations() == 0 || solverPerf.nIterations() == backstopIter)
                {
                    for (label cell=0; cell<nCells; cell++)
                    {
                        pAPtr[cell] = wAPtr[cell];
                    }
                }
                else
                {
                    solveScalar beta = wArA/wArAold;

                    for (label cell=0; cell<nCells; cell++)
                    {
                        pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                    }
                }

                // --- Update preconditioned residual
                matrix_.Amul(wA, pA, interfaceBouCoeffs_, interfaces_, cmpt);

                solveScalar wApA = gSumProd(wA, pA, matrix().mesh().comm());

                // --- Test for singularity
                if (solverPerf.checkSingularity(mag(wApA)/normFactor)) break;

                // --- Update solution and residual:

                solveScalar alpha = wArA/wApA;

                for (label cell=0; cell<nCells; cell++) 
                {
                    psiPtr[cell] += alpha*pAPtr[cell];
                    rAPtr[cell] -= alpha*wAPtr[cell];
                }

                solverPerf.finalResidual() =
                    gSumMag(rA, matrix().mesh().comm())
                   /normFactor;

            } while
            (
                (
                  ++solverPerf.nIterations() < maxIter
                && !solverPerf.checkConvergence(tolerance_, relTol_, log_)
                )
             || solverPerf.nIterations() < minIter_
            );
            iterationTime -= clockValue::now();
        }

        // --- Exit if converged or if already tried backstopping
        if (backstop == 1 || solverPerf.checkConvergence(tolerance_, relTol_, log_)) 
        { 
            matrix().setResidualField
            (
                ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
                fieldName_,
                false
            );
            break;
        } 

        Info << "PCG backstopping at iteration " << backstopIter << endl;
        if (backstop_ == -1) {
            maxIter += maxIter_;
        } else {
            maxIter += backstop_;
        }
        if (subDict.get<word>("preconditioner") == "DIC") {
            backstopIter = maxIter;
        } else {
            backstopDict.set("preconditioner", "DIC");
        }
    }

    solverTime -= clockValue::now();
    learningTime += clockValue::now();
    scalar costEstimate = 0.0;
    if (solverPerf.nIterations() > 0) {

        // --- Compute solver cost
        if (deterministic_) {
            costEstimate = 1e-9 * totalCostEstimate(solverPerf.nIterations());
        } else {
            costEstimate = -solverTime;
        }

        // --- Pass cost to learning algorithm
        if (static_ == -1 && !randomUniform_ && Pstream::myProcNo() == 0) {
            dictionary& learningDict = learningDicts.subDict(banditName_);
            learningDict.set<scalar>("loss", costEstimate);
        }
    }
    learningTime -= clockValue::now();
    PCGTime -= solverTime;

    Info<< "INFO: banditName=" << banditName_;
    Info<< ", fieldName=" << fieldName_;
    Info<< ", relativeTolerance=" << relTol_;
    Info<< ", tolerance=" << tolerance_;
    Info<< ", initialResidual=" << solverPerf.initialResidual();
    Info<< ", finalResidual=" << solverPerf.finalResidual();
    Info<< ", nIterations=" << solverPerf.nIterations();
    Info<< ", preconstructTime=" << -preconstructTime;
    Info<< ", iterationTime=" << -iterationTime;
    Info<< ", learningTime=" << -learningTime;
    Info<< ", solverTime=" << -solverTime;
    Info<< ", PCGTime=" << PCGTime;
    if (deterministic_ && solverPerf.nIterations() > 0) {
        Info<< ", costEstimate=" << costEstimate;
    }
    Info<< endl;

    #ifdef DUMP_ABSOL
    #include "finishDump.H"
    #endif

    return solverPerf;
}

Foam::solverPerformance Foam::PCGBandit::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    return scalarSolve
    (
        tpsi.ref(),
        ConstPrecisionAdaptor<solveScalar, scalar>(source)(),
        cmpt
    );
}

// ************************************************************************* //
