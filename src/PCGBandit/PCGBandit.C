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

#include "HashPtrTable.H"

//#define PCGB_DEBUG
//#define DUMP_ABSOL

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    clockValue PCGTime = clockValue();
    scalar PCGCost = 0.0;

    defineTypeNameAndDebug(PCGBandit, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PCGBandit>
        addPCGBanditSymMatrixConstructorToTable_;
    word GAMG_or_FGAMG = (lduMatrix::preconditioner::symMatrixConstructorTablePtr_->sortedToc()).found("FGAMG") ? "FGAMG" : "GAMG";

    Random rndGen;

    dictionary preconditionerDict;
    dictionary subDict;
    HashTable<List<dictionary>> preconditionerDictsMap;
    dictionary learningDicts;

    // --- GAMG configuration space specification and defaults
    const List<Tuple2<word, List<word>>> GAMGDefaultLists = {
        Tuple2<word, List<word>>("smoother",                {"GaussSeidel", "DIC", "DICGaussSeidel", "symGaussSeidel"}),
        Tuple2<word, List<word>>("agglomerator",            {"faceAreaPair", "algebraicPair"}),
        Tuple2<word, List<word>>("directSolveCoarsest",     {"no", "yes"}),
        Tuple2<word, List<word>>("nCellsInCoarsestLevel",   {"10", "100", "1000"}),
        Tuple2<word, List<word>>("mergeLevels",             {"1", "2"}),
        Tuple2<word, List<word>>("nPreSweeps",              {"0", "2"}),
        Tuple2<word, List<word>>("nPostSweeps",             {"1", "2"}),
        Tuple2<word, List<word>>("nFinestSweeps",           {"2"}),
        Tuple2<word, List<word>>("nVcycles",                {"1", "2"})
    };
    const List<word> ICTCSuffixes = {"m4", "m3p5", "m3", "m2p5", "m2", "m1p5", "m1", "m0p5"};
    const HashSet<word> noCacheAgglomeration = {"agglomerator", "nCellsInCoarsestLevel", "mergeLevels"};

    #ifdef PCGB_DEBUG
    #include "initializeTracking.H"
    #endif

    #ifdef DUMP_ABSOL
    #include "initializeDumping.H"
    #endif

    HashPtrTable<DecomposedLaplacian> nonSerializableObjects_;

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
    )
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
    lossEstimator_ = solverControls.getOrDefault<word>("lossEstimator", "RV");
    deterministic_ = Switch(solverControls.getOrDefault<word>("deterministic", "no"));
    seed_ = mesh.time().controlDict().getOrDefault<label>("randomSeed", 0);
    backstop_ = solverControls.getOrDefault<label>("backstop", -1);
    static_ = label(solverControls.getOrDefault<label>("static", -1));
    randomUniform_ = Switch(solverControls.getOrDefault<word>("randomUniform", "no"));
    banditAlgorithm_ = solverControls.getOrDefault<word>("banditAlgorithm", "TsallisINF");

    if (!preconditionerDictsMap.found(banditName_))
    {
        if (Pstream::myProcNo() == 0) {
            rndGen.reset(seed_);
        }

        // --- Read GAMG tune flags and option lists
        label minDroptolIdx = ICTCSuffixes.find(solverControls.getOrDefault<word>("minSmootherLogDroptol", "m4"));
        if (minDroptolIdx == -1) {
            FatalErrorInFunction << "minSmootherLogDroptol must be one of " << ICTCSuffixes << exit(FatalError);
        }
        label maxDroptolIdx = ICTCSuffixes.find(solverControls.getOrDefault<word>("maxSmootherLogDroptol", "m0p5"));
        if (maxDroptolIdx == -1) {
            FatalErrorInFunction << "maxSmootherLogDroptol must be one of " << ICTCSuffixes << exit(FatalError);
        }
        if (minDroptolIdx > maxDroptolIdx) {
            FatalErrorInFunction << "minSmootherLogDroptol cannot be greater than maxSmootherLogDroptol" << exit(FatalError);
        }
        bool coarsestSmootherTune = Switch(solverControls.getOrDefault<word>("coarsestSmootherTune", "no"));
        if (coarsestSmootherTune && GAMG_or_FGAMG == "GAMG") {
            FatalErrorInFunction << "coarsestSmootherTune requires FGAMG" << exit(FatalError);
        }
        label inc = max(ICTCSuffixes.size() / solverControls.getOrDefault<label>("numSmootherDroptols", ICTCSuffixes.size()), 1);

        bool cacheAgglomeration = true;
        label dGAMG = 0;
        List<List<word>> GAMGOptions(GAMGDefaultLists.size());
        for (label j = 0; j < GAMGDefaultLists.size(); ++j) {
            const word param = GAMGDefaultLists[j].first();
            const List<word>& defaults = GAMGDefaultLists[j].second();
            const word tuneKey = param + "Tune";
            if (solverControls.found(tuneKey)) {
                ITstream& is = solverControls.lookup(tuneKey);
                token tok(is);
                if (tok.isPunctuation(token::BEGIN_LIST)) {
                    const List<token>& toks = is;
                    DynamicList<word> opts;
                    for (const token& t : toks) {
                        if (!t.isPunctuation()) {
                            OStringStream os;
                            os << t;
                            word config = os.str();
                            if (param == "smoother" && (config == "ICTC" || config == "ICTCGaussSeidel")) { // add ICTC smoothers
                                if (GAMG_or_FGAMG == "GAMG") {
                                    WarningInFunction<< "Set " << config << " smoother but FGAMG unavailable; using GAMG (may be slow)" << endl;
                                }
                                for (label finestIdx = minDroptolIdx; finestIdx <= maxDroptolIdx; finestIdx += inc) {
                                    word finest = config + "_" + ICTCSuffixes[finestIdx];
                                    if (coarsestSmootherTune) { // pair up smoother and coarsestSmoother
                                        for (label coarsestIdx = minDroptolIdx; coarsestIdx <= maxDroptolIdx; coarsestIdx += inc) {
                                            opts.append(finest + "; coarsestSmoother " + config + "_" + ICTCSuffixes[coarsestIdx]);
                                        }
                                    } else {
                                        opts.append(finest);
                                    }
                                }
                            } else {
                                opts.append(config);
                            }
                        }
                    }

                    GAMGOptions[j] = opts;
                } else if (Switch::found(tok.wordToken())) {
                    if (Switch(tok.wordToken())) {
                        GAMGOptions[j] = defaults;
                    }
                } else {
                    FatalErrorInFunction << tuneKey << " must be a bool or list" << exit(FatalError);
                }
            }
            dGAMG = max(dGAMG, 1) * max(GAMGOptions[j].size(), min(dGAMG, 1));
            if (noCacheAgglomeration.found(param) && GAMGOptions[j].size() > 1) {
                cacheAgglomeration = false; // turn off cacheAgglomeration if tuning any param that affects agglomeration
            }
        }

        if (cacheAgglomeration || static_ > -1) {
            cacheAgglomeration = Switch(solverControls.getOrDefault<word>("cacheAgglomeration", "yes"));
        }

        // --- Build preconditioner dictionary list
        const scalar maxLogDroptol_ = solverControls.getOrDefault<scalar>("maxLogDroptol", -0.5);
        const scalar minLogDroptol_ = solverControls.getOrDefault<scalar>("minLogDroptol", -4.0);
        const label numDroptols_ = solverControls.getOrDefault<label>("numDroptols", 0);
        bool DICTune = Switch(solverControls.getOrDefault<word>("DICTune", "yes"));
        label d = numDroptols_ + label(DICTune) + dGAMG;
        List<dictionary> preconditionerDicts(d);
        for (label i = 0; i < numDroptols_; ++i) {
            preconditionerDicts[i].set("preconditioner", "ICTC");
            scalar droptol;
            if (i == 0) {
                droptol = pow(10.0, minLogDroptol_); // handles case of numDroptols_ = 1
            } else {
                droptol = pow(10.0, minLogDroptol_ + (maxLogDroptol_ - minLogDroptol_) * scalar(i) / scalar(numDroptols_ - 1));
            }
            preconditionerDicts[i].set("droptol", droptol);
        }

        if (DICTune){
            preconditionerDicts[numDroptols_].set("preconditioner", "DIC");
        }

        for (label i = 0; i < dGAMG; ++i) {
            dictionary& pd = preconditionerDicts[d - dGAMG + i];
            pd.set("preconditioner", GAMG_or_FGAMG);
            if (GAMGOptions[0].size() == 0) {
                pd.set("smoother", "DICGaussSeidel");
            }
            pd.set("cacheAgglomeration", cacheAgglomeration);

            label remaining = i;
            for (label j = GAMGDefaultLists.size()-1; j >= 0; j--) {
                if (GAMGOptions[j].size() > 0) {
                    label size = GAMGOptions[j].size();
                    word config = GAMGOptions[j][remaining % size];
                    pd.set(GAMGDefaultLists[j].first(), config);
                    label idx = config.find("; coarsestSmoother ");
                    if (idx != -1) {
                        pd.set("coarsestSmoother", word(config.substr(idx + 19)));
                    }
                    remaining /= size;
                }
            }
 
            OStringStream oss;
            pd.write(oss, false);
            pd = dictionary(IStringStream(oss.str())());
        }

        preconditionerDictsMap.set(banditName_, preconditionerDicts);
        #ifdef PCGB_DEBUG
        Info<< "Preconditioner configurations for " << banditName_ << " : " << preconditionerDictsMap[banditName_] << endl;
        #endif

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
    const List<dictionary>& preconditionerDicts = preconditionerDictsMap[banditName_];

    if (i == -1) {

        if (Pstream::myProcNo() == 0) {
            label d = preconditionerDicts.size();
            if (d == 1) {
                i = 0;
            } else if (randomUniform_) {
                i = floor(scalar(d) * rndGen.sample01<scalar>());
            } else if (banditAlgorithm_ == "ThompsonSampling") {
                #include "ThompsonSampling.H"
            } else if (banditAlgorithm_ == "simTsallisINF") {
                #include "simTsallisINF.H"
            } else if (banditAlgorithm_ == "SpeKL") {
                #include "SpeKL.H"
            } else {
                #include "TsallisINF.H"
            }
        }

        Pstream::broadcast(i);

        #ifdef PCGB_DEBUG
        Info<< banditAlgorithm_ << " Selection: ";
    } else {

        Info<< "Static Preconditioner: ";
        #endif
    }


    subDict = preconditionerDicts[i];
    preconditionerDict.set("preconditioner", subDict);

    #ifdef PCGB_DEBUG
    Info<< subDict << endl;
    #endif
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
    #include "printTracking.H"
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
        PCGCost += costEstimate;
        Info<< ", PCGCost=" << PCGCost;
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
