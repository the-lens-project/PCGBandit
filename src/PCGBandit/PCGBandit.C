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
    const List<word> ICTCSuffixes =
        {"m5", "m4p5", "m4", "m3p5", "m3", "m2p5", "m2", "m1p5", "m1", "m0p5"};
    const List<word> SORSuffixes =
        {"p0p1", "p0p2", "p0p3", "p0p4", "p0p5", "p0p6", "p0p7", "p0p8", "p0p9",
         "p1p0", "p1p1", "p1p2", "p1p3", "p1p4", "p1p5", "p1p6", "p1p7", "p1p8", "p1p9"};
    const HashSet<word> noCacheAgglomeration = {"agglomerator", "nCellsInCoarsestLevel", "mergeLevels"};

    // A strided window [minIdx, maxIdx] step inc into a suffix list; lets the 
    // ICTC droptol axis and the SOR omega axis be tuned through one mechanism.
    struct SmootherRange
    {
        label minIdx;
        label maxIdx;
        label inc;
    };

    // Build a SmootherRange over `suffixes` from solverControls: minKey/maxKey name
    // the endpoint suffixes (defaulting to minDefault/maxDefault) and numKey sets
    // how many entries to select evenly across that window (0 or omitted = do not
    // tune this axis, i.e. select no ICTC/SOR smoothers).
    static SmootherRange readSmootherRange
    (
        const dictionary& solverControls,
        const List<word>& suffixes,
        const word& minKey,
        const word& maxKey,
        const word& numKey,
        const word& minDefault,
        const word& maxDefault,
        const label& numDefault
    )
    {
        SmootherRange range;
        range.minIdx = suffixes.find(solverControls.getOrDefault<word>(minKey, minDefault));
        range.maxIdx = suffixes.find(solverControls.getOrDefault<word>(maxKey, maxDefault));
        if (range.minIdx == -1 || range.maxIdx == -1) {
            FatalErrorInFunction
                << minKey << " and " << maxKey << " must each be one of "
                << suffixes << exit(FatalError);
        }
        if (range.minIdx > range.maxIdx) {
            FatalErrorInFunction
                << minKey << " cannot come before " << maxKey << " in "
                << suffixes << exit(FatalError);
        }
        const label num = solverControls.getOrDefault<label>(numKey, numDefault);
        range.inc = max((range.maxIdx - range.minIdx) / max(num - 1, label(1)), label(1));
        return range;
    }

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

    #ifdef DUMP_ABSOL
    #include "Absol/initializeDumping.H"
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
    nGeometricD_ = mesh.nGeometricD();          // 2 (planar) or 3
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
    backstop_ = solverControls.getOrDefault<label>("backstop", -1);
    static_ = label(solverControls.getOrDefault<label>("static", -1));
    randomUniform_ = Switch(solverControls.getOrDefault<word>("randomUniform", "no"));
    banditAlgorithm_ = solverControls.getOrDefault<word>("banditAlgorithm", "TsallisINF");

    if (!preconditionerDictsMap.found(banditName_))
    {
        if (Pstream::myProcNo() == 0) {
            rndGen.reset(mesh.time().controlDict().getOrDefault<label>("randomSeed", 0));
        }

        // --- Read GAMG tune flags and option lists
        const SmootherRange ictcRange = readSmootherRange
        (
            solverControls, ICTCSuffixes,
            "minSmootherLogDroptol", "maxSmootherLogDroptol", "numSmootherDroptols",
            "m4", "m0p5", 4
        );
        const SmootherRange sorRange = readSmootherRange
        (
            solverControls, SORSuffixes,
            "minSmootherOmega", "maxSmootherOmega", "numSmootherOmegas",
            "p0p8", "p1p2", 5
        );

        bool cacheAgglomeration = true;
        label dGAMG = 0;
        label nCellsMin = -1;
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
                    // Deduplicate options so e.g. an explicit GaussSeidel and
                    // SOR's omega=1.0 (also GaussSeidel) collapse to one arm.
                    auto appendUnique = [&opts](const word& option) {
                        if (!opts.found(option)) opts.append(option);
                    };
                    for (const token& t : toks) {
                        if (!t.isPunctuation()) {
                            OStringStream os;
                            os << t;
                            word config = os.str();
                            if (param == "smoother" && (config == "ICTC" || config == "ICTCGaussSeidel")) { // add ICTC smoothers
                                if (GAMG_or_FGAMG == "GAMG") {
                                    WarningInFunction<< "Set " << config << " smoother but FGAMG unavailable; using GAMG (may be slow)" << endl;
                                }
                                for (label idx = ictcRange.minIdx; idx <= ictcRange.maxIdx; idx += ictcRange.inc) {
                                    appendUnique(config + "_" + ICTCSuffixes[idx]);
                                }
                            } else if (param == "smoother" && (config == "SOR" || config == "DICSOR")) { // add SOR smoothers
                                for (label idx = sorRange.minIdx; idx <= sorRange.maxIdx; idx += sorRange.inc) {
                                    const word& suffix = SORSuffixes[idx];
                                    if (suffix == "p1p0") { // omega = 1 is the built-in (DIC)GaussSeidel smoother
                                        appendUnique(config == "SOR" ? word("GaussSeidel") : word("DICGaussSeidel"));
                                    } else {
                                        appendUnique(config + "_" + suffix);
                                    }
                                }
                            } else if (param == "nCellsInCoarsestLevel") {
                                if (nCellsMin == -1) {
                                    nCellsMin = returnReduce(matrix.diag().size(), minOp<label>());
                                }
                                label resolved = readLabel(config);
                                label clamped = max(label(1), min(resolved, nCellsMin));
                                if (clamped != resolved) {
                                    WarningInFunction
                                        << "nCellsInCoarsestLevel value " << resolved
                                        << " clamped to " << clamped
                                        << " (per-processor min cells: " << nCellsMin << ")" << endl;
                                }
                                appendUnique(Foam::name(clamped));
                            } else {
                                appendUnique(config);
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

namespace Foam 
{

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
            if (smoother == "GaussSeidel" || smoother == "DICGaussSeidel") {
                c += scalar(2 * nnzL + nCells);
            }
            if (smoother == "DIC" || smoother == "DICGaussSeidel") {
                c += scalar(4 * nnzL + nCells);
            }
            if (c == 0.0) {
                c = matvec;
            }
        }
        return c * scalar(nSweeps);
    }

}

Foam::scalar Foam::PCGBandit::communicationCostEstimate
(
  const word preconditioner
) const
{
    if (Pstream::nProcs(matrix_.mesh().comm()) <= 1 || !isGAMG(preconditioner)) {
        return 0.0;
    }

    const GAMGAgglomeration& agg = GAMGAgglomeration::New(matrix_, subDict);
    const label L = agg.size();

    // --- Each level: restrict/prolong + smoother halo exchange + scale reduction
    scalar perVcycle = scalar(3 * L);

    // --- DIC-PCG iterations (n) if directSolveCoarsest is off
    if (L > 0 && !Switch(subDict.getOrDefault<word>("directSolveCoarsest", "no"))) {
        const scalar ncG = returnReduce(scalar(max(agg.nCells(L - 1), label(1))), sumOp<scalar>());
        perVcycle += COARSE_CG_WEIGHT * ncG;
    }

    return COMM_EVENT_FLOPS * perVcycle * scalar(subDict.getOrDefault<label>("nVcycles", 2));
}

Foam::scalar Foam::PCGBandit::perIterationCostEstimate
(
    const word preconditioner
) const
{
    const label nCells = matrix_.diag().size();
    const label nnzL = matrix_.lower().size();

    // ICTC scattered-triangular-solve baseline (density-independent, 3D only).
    const scalar ictcBaseline = (nGeometricD_ >= 3) ? ICTC_SMOOTHER_BASELINE: 0.0;

    // --- One CG step has a matvec (2 * nnzL + nCells) and five vector operations
    const scalar cgStep = scalar(2 * nnzL + 6 * nCells);

    if (preconditioner == "ICTC") {
        return returnReduce(cgStep + scalar(2 * (debug::controlDict().get<label>("ICTC_NNZ") + nCells)), maxOp<scalar>());
    }
    if (preconditioner == "DIC") {
        return returnReduce(cgStep + scalar(4 * nnzL + nCells), maxOp<scalar>());
    }

    const label nPreSweeps  = subDict.getOrDefault<label>("nPreSweeps", 0);
    const label nPostSweeps = subDict.getOrDefault<label>("nPostSweeps", 2);
    const label maxPreSweeps = 4;
    const label maxPostSweeps = 4;
    const label preSweepsLevelMultiplier = 1;
    const label postSweepsLevelMultiplier = 1;
    const label nVcycles = subDict.getOrDefault<label>("nVcycles", 2);
    const label nFinestSweeps = subDict.getOrDefault<label>("nFinestSweeps", 2);
    const word smoother = subDict.get<word>("smoother");
    const bool interpolateCorrection = Switch(subDict.getOrDefault<word>("interpolateCorrection", "no"));
    const bool directSolveCoarsest = Switch(subDict.getOrDefault<word>("directSolveCoarsest", "no"));

    const bool scaleCorrection = matrix_.symmetric();
    const GAMGAgglomeration& agg = GAMGAgglomeration::New(matrix_, subDict);
    const label L = agg.size();

    // --- Compute ratio of ICTC smoother fills to agglomeration matrix fills
    scalar fillFactor = 1.0;
    if (smoother.find("ICTC") != std::string::npos) {
        scalar nnzAgg = scalar(nnzL);
        for (label i = 0; i < L; i++) {
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
        const label nc  = agg.nCells(i);
        const label nf  = agg.nFaces(i);
        const label nPre  = (nPreSweeps  > 0)
            ? min(nPreSweeps  + preSweepsLevelMultiplier  * i, maxPreSweeps)  : 0;
        const label nPost = (nPostSweeps > 0)
            ? min(nPostSweeps + postSweepsLevelMultiplier * i, maxPostSweeps) : 0;

        perVcycle += scalar(2 * nc);
        perVcycle += smootherApplyCost(smoother, nf, nc, nPre + nPost, fillFactor, ictcBaseline);
        if (nPre > 0) {
            perVcycle += scalar(2 * nf + nc);
        }
        perVcycle += scalar(2 * nc);
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
            const scalar ncG = returnReduce(scalar(agg.nCells(L - 1)), sumOp<scalar>());
            perVcycle += 2.0 * ncG * ncG;
        } else {
            perVcycle += sqrt(nc) * scalar(6 * nf + 3 * label(nc));
        }
    }

    return returnReduce(cgStep + perVcycle * scalar(nVcycles), maxOp<scalar>()) + communicationCostEstimate(preconditioner);
}


Foam::scalar Foam::PCGBandit::setupCostEstimate
(
    const word preconditioner
) const
{
    const scalar pICE = perIterationCostEstimate(preconditioner);

    if (preconditioner == "ICTC") {
        return ICTC_SETUP_WEIGHT * pICE;
    }
    if (!isGAMG(preconditioner)) {
        return pICE;
    }

    // --- GAMG: rebuild agglomeration (unless cached) and smoothers every solve
    const GAMGAgglomeration& agg = GAMGAgglomeration::New(matrix_, subDict);
    scalar hierarchy = scalar(2 * matrix_.lower().size() + matrix_.diag().size());
    for (label i = 0; i < agg.size(); i++) {
        hierarchy += scalar(2 * agg.nFaces(i) + agg.nCells(i));
    }
    const bool cached = subDict.getOrDefault<label>("cacheAgglomeration", 1);
    scalar weight = GAMG_ASSEMBLY_WEIGHT + (cached ? 0.0 : GAMG_AGGLOMERATION_WEIGHT);
    scalar setup = weight * hierarchy;

    // --- FGAMG: precomputes smoother factors once per solve
    const word smoother = subDict.getOrDefault<word>("smoother", "");
    const word coarsestSmoother = subDict.getOrDefault<word>("coarsestSmoother", smoother);
    if (smoother.find("ICTC") != string::npos
     || coarsestSmoother.find("ICTC") != string::npos) {
        const label smootherNNZ = debug::controlDict().get<label>("ICTC_SMOOTHER_NNZ");
        setup += ICTC_SMOOTHER_SETUP_WEIGHT * scalar(smootherNNZ);
    }

    // --- Dense LU factorisation (n^3) if directSolveCoarsest is on
    if (agg.size() > 0 && Switch(subDict.getOrDefault<word>("directSolveCoarsest", "no"))) {
        const scalar ncG = returnReduce(scalar(agg.nCells(agg.size() - 1)), sumOp<scalar>());
        setup += DIRECT_LU_WEIGHT * ncG * ncG * ncG;
    }

    return returnReduce(setup, maxOp<scalar>());
}



Foam::scalar Foam::PCGBandit::totalCostEstimate
(
    const label nIterations
) const
{
    const word preconditioner = subDict.get<word>("preconditioner");
    const scalar pICE = perIterationCostEstimate(preconditioner);
    scalar cost = setupCostEstimate(preconditioner);

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
    #include "Absol/startDump.H"
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

                // --- Resets the ICTCSmoother NNZ accumulator for deterministic mode
                debug::controlDict().set<label>("ICTC_SMOOTHER_NNZ", 0);

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
    #include "Absol/finishDump.H"
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
