/*---------------------------------------------------------------------------*\
                  Class configurationSpace Implementation
\*---------------------------------------------------------------------------*/

#include "configurationSpace.H"
#include "Switch.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    word GAMG_or_FGAMG = (lduMatrix::preconditioner::symMatrixConstructorTablePtr_->sortedToc()).found("FGAMG") ? "FGAMG" : "GAMG";

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
    struct smootherRange
    {
        label minIdx;
        label maxIdx;
        label inc;
    };

    // Build a smootherRange over `suffixes` from solverControls: minKey/maxKey name
    // the endpoint suffixes (defaulting to minDefault/maxDefault) and numKey sets
    // how many to select evenly across that window (0 = do not tune this axis)
    static smootherRange readSmootherRange
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
        smootherRange range;
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

    // Read an axis; an unset non-positive default gives an empty list.
    template<class T>
    static List<T> readSubspaceAxis
    (
        const dictionary& solverControls,
        const word& param,
        const List<T>& defaults
    )
    {
        List<T> values;
        T defaultValue = solverControls.getOrDefault<T>(param, defaults[0]);
        if (defaultValue > 0) {
            values = List<T>(1, defaultValue);
        }

        const word tuneKey = param + "Tune";
        if (solverControls.found(tuneKey)) {

            ITstream& is = solverControls.lookup(tuneKey);
            token tok(is);

            if (tok.isPunctuation(token::BEGIN_LIST)) {
                is.putBack(tok);
                values = List<T>(is);
                for (const T value : values) {
                    if (!(value >= 0)) {
                        FatalErrorInFunction
                            << tuneKey << " requires non-negative "
                            << pTraits<T>::typeName << " values, found "
                            << value << exit(FatalError);
                    }
                }
            } else if (tok.isWord() && Switch::found(tok.wordToken())) {
                if (Switch(tok.wordToken())) {
                    values = defaults;
                }
            } else {
                FatalErrorInFunction
                    << tuneKey << " must be yes, no, or a list, found "
                    << tok << exit(FatalError);
            }
        }

        return values;
    }


    // One point of the subspace configuration space. lenHistory and
    // decayRate are alternatives, so at most one can be non-zero
    struct subspaceConfig
    {
        label numProbes;
        label lenHistory;
        scalar decayRate;

        bool operator==(const subspaceConfig& o) const
        {
            return numProbes == o.numProbes
                && lenHistory == o.lenHistory
                && decayRate == o.decayRate;
        }
    };


    // The subspace configs, pruned to be unique by clamping numProbes to be at
    // most lenHistory and setting numProbes to zero if the decayRate is zero
    static List<subspaceConfig> readSubspaceConfigs
    (
        const dictionary& solverControls,
        DynamicList<scalar>& decayRates
    )
    {
        DynamicList<subspaceConfig> cfgs;
        List<label> nums = readSubspaceAxis(solverControls, "numProbes", numProbesDefault);

        List<label> lens = readSubspaceAxis(solverControls, "lenHistory", lenHistoryDefault);
        for (const label lenHistory: lens) {
            for (const label numProbes : nums) {
                const subspaceConfig cfg{min(numProbes, lenHistory), numProbes > 0 ? lenHistory : 0, 0.0};
                if (!cfgs.found(cfg)) {
                    cfgs.append(cfg);
                }
            }
        }

        List<scalar> rates = readSubspaceAxis(solverControls, "decayRate", decayRateDefault);
        const bool hasProbes = !nums.empty() && max(nums) > 0;
        for (const scalar decayRate: rates) {
            // Collect each input rate before expanding probe and preconditioner arms.
            if (decayRate > 0 && hasProbes) {
                decayRates.append(decayRate);
            }
            for (const label numProbes : nums) {
                const subspaceConfig cfg{decayRate > 0.0 ? numProbes : 0, 0, numProbes > 0 ? decayRate : 0.0};
                if (!cfgs.found(cfg)) {
                    cfgs.append(cfg);
                }
            }
        }

        return List<subspaceConfig>(cfgs);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::configurationSpace::appendICTC(const dictionary& solverControls)
{
    const scalar maxLogDroptol_ = solverControls.getOrDefault<scalar>("maxLogDroptol", -0.5);
    const scalar minLogDroptol_ = solverControls.getOrDefault<scalar>("minLogDroptol", -4.0);
    const label numDroptols_ = solverControls.getOrDefault<label>("numDroptols", 0);

    for (label i = 0; i < numDroptols_; ++i) {
        dictionary pd;
        pd.set("preconditioner", "ICTC");
        scalar droptol;
        if (i == 0) {
            droptol = pow(10.0, minLogDroptol_); // handles case of numDroptols_ = 1
        } else {
            droptol = pow(10.0, minLogDroptol_ + (maxLogDroptol_ - minLogDroptol_) * scalar(i) / scalar(numDroptols_ - 1));
        }
        pd.set("droptol", droptol);
        dicts_.append(pd);
    }
}


void Foam::configurationSpace::appendDIC(const dictionary& solverControls)
{
    if (Switch(solverControls.getOrDefault<word>("DICTune", "yes"))) {
        dictionary pd;
        pd.set("preconditioner", "DIC");
        dicts_.append(pd);
    }
}


void Foam::configurationSpace::appendGAMG
(
    const dictionary& solverControls,
    const lduMatrix& matrix,
    const label staticArm
)
{
    // --- Read GAMG tune flags and option lists
    const smootherRange ictcRange = readSmootherRange
    (
        solverControls, ICTCSuffixes,
        "minSmootherLogDroptol", "maxSmootherLogDroptol", "numSmootherDroptols",
        "m4", "m0p5", 4
    );
    const smootherRange sorRange = readSmootherRange
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
                // Collapse repeated options to one arm.
                auto appendUnique = [&opts](const word& option) {
                    if (!opts.found(option)) opts.append(option);
                };
                for (const token& t : toks) {
                    if (!t.isPunctuation()) {
                        OStringStream os;
                        os << t;
                        word config = os.str();
                        if (param == "smoother" && GAMG_or_FGAMG == "GAMG" && config.startsWith("ICTC")) {
                            if (Switch(solverControls.getOrDefault<word>("deterministic", "no"))
                                || solverControls.getOrDefault<label>("backstop", -1) == -1) {
                                FatalErrorInFunction
                                    << "cost estimation not implemented for ICTC smoothers for GAMG; "
                                    << "load libFGAMG" << exit(FatalError);
                            }
                            WarningInFunction
                                << "Set " << config  << " smoother but FGAMG unavailable; "
                                << "using GAMG (may be slow)" << endl;
                        }
                        if (param == "smoother" && (config == "ICTC" || config == "ICTCGaussSeidel")) { // add ICTC smoothers
                            for (label idx = ictcRange.minIdx; idx <= ictcRange.maxIdx; idx += ictcRange.inc) {
                                appendUnique(config + "_" + ICTCSuffixes[idx]);
                            }
                        } else if (param == "smoother" && (config == "SOR" || config == "DICSOR")) { // add SOR smoothers
                            for (label idx = sorRange.minIdx; idx <= sorRange.maxIdx; idx += sorRange.inc) {
                                appendUnique(config + "_" + SORSuffixes[idx]);
                            }
                        } else if (param == "nCellsInCoarsestLevel") {
                            if (nCellsMin == -1) {
                                nCellsMin = returnReduce
                                (
                                    matrix.diag().size(), minOp<label>(),
                                    UPstream::msgType(), matrix.mesh().comm()
                                );
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

    if (cacheAgglomeration || staticArm > -1) {
        cacheAgglomeration = Switch(solverControls.getOrDefault<word>("cacheAgglomeration", "yes"));
    }
    for (label i = 0; i < dGAMG; ++i) {
        dictionary pd;
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
        dicts_.append(dictionary(IStringStream(oss.str())()));
    }
}


void Foam::configurationSpace::expandInitializations
(
    const dictionary& solverControls
)
{
    // --- Cross-multiply the subspace configs onto every arm.
    const List<subspaceConfig> subspaceConfigs = readSubspaceConfigs(solverControls, decayRates_);
    if (subspaceConfigs.size() > 0) {
        const label dP = dicts_.size();
        List<dictionary> expanded(dP*subspaceConfigs.size());
        forAll(subspaceConfigs, k) {
            for (label i = 0; i < dP; ++i) {
                dictionary pd = dicts_[i];
                pd.set("lenHistory", subspaceConfigs[k].lenHistory);
                pd.set("numProbes", subspaceConfigs[k].numProbes);
                pd.set("decayRate", subspaceConfigs[k].decayRate);
                expanded[k*dP + i] = pd;
            }
        }
        dicts_.transfer(expanded);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::configurationSpace::configurationSpace
(
    const dictionary& solverControls,
    const lduMatrix& matrix,
    const label staticArm
)
{
    // Order is the contract; see the header.
    appendICTC(solverControls);
    appendDIC(solverControls);
    appendGAMG(solverControls, matrix, staticArm);
    expandInitializations(solverControls);

    // Ceilings over the assembled arms, not over the raw tune lists: the
    // window and the test matrix are sized from these and shared by every arm,
    // and pruned pairs must not inflate them.  Taking them from the space also
    // means every solver dictionary that shares a bandit shares one ceiling,
    // so p and pFinal cannot end up with arms one of them is unable to run.
    for (const dictionary& pd : dicts_) {
        maxNumProbes_ = max(maxNumProbes_, pd.getOrDefault<label>("numProbes", 4));
        maxLenHistory_ = max(maxLenHistory_, pd.getOrDefault<label>("lenHistory", 0));
    }
}

// ************************************************************************* //
