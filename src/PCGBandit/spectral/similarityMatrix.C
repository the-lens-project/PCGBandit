/*---------------------------------------------------------------------------*\
        Functions to build matrices of preconditioner similarities.
\*---------------------------------------------------------------------------*/

#include "similarityMatrix.H"
#include "configurationSpace.H"

namespace Foam
{

// Parse the droptol scalar from a word-encoded smoother name.
// e.g. "ICTCGaussSeidel_m3p5" -> 10^-3.5, "DICGaussSeidel" -> 1.0
static scalar smootherToDroptol(const word& smootherName) {
    if (smootherName == "DICSOR_p1p0") return 1.0;

    label underscoreIdx = smootherName.rfind('_');
    if (underscoreIdx == -1) return 1.0;

    word suffix = smootherName.substr(underscoreIdx + 1);
    word magnitudeStr = suffix.substr(1);

    label decimalIdx = magnitudeStr.find('p');
    if (decimalIdx == -1) {
        return pow(10.0, -1.0 * std::stod(magnitudeStr));
    } else {
        scalar wholePart   = std::stod(magnitudeStr.substr(0, decimalIdx));
        scalar fracPart    = std::stod(magnitudeStr.substr(decimalIdx + 1));
        scalar fracDivisor = pow(10.0, scalar(magnitudeStr.substr(decimalIdx + 1).size()));
        return pow(10.0, -1.0 * (wholePart + fracPart / fracDivisor));
    }
}

static bool smootherIsICTCLike(const word& smootherName) {
    return smootherName.startsWith("ICTC_")
        || smootherName == "ICTC"
        || smootherName == "DIC";
}

static bool smootherIsICTCGaussSeidelLike(const word& smootherName) {
    return smootherName.startsWith("ICTCGaussSeidel_")
        || smootherName == "ICTCGaussSeidel"
        || smootherName == "DICGaussSeidel"
        || smootherName == "DICSOR_p1p0";
}

// Parse the SOR relaxation factor omega from a word-encoded smoother name.
// e.g. "SOR_p0p5" -> 0.5, "DICSOR_p1p3" -> 1.3, "DICGaussSeidel" -> 1.0 
static scalar smootherToOmega(const word& smootherName) {

    label underscoreIdx = smootherName.rfind('_');
    if (underscoreIdx == -1) return 1.0;

    word suffix = smootherName.substr(underscoreIdx + 1);
    label decimalIdx = suffix.find('p', 1);

    scalar wholePart   = std::stod(suffix.substr(1, decimalIdx - 1));
    word fracStr       = suffix.substr(decimalIdx + 1);
    scalar fracDivisor = pow(10.0, scalar(fracStr.size()));
    return wholePart + std::stod(fracStr) / fracDivisor;
}

static bool smootherIsSORLike(const word& smootherName) {
    return smootherName.startsWith("SOR_") 
        || smootherName == "GaussSeidel";
}

static bool smootherIsDICSORLike(const word& smootherName) {
    return smootherName.startsWith("DICSOR_") 
        || smootherName == "DICGaussSeidel";
}


static scalar omegaSimilarity(const scalar omega_i, const scalar omega_j)
{
    return 1.0 / (1.0 + mag(omega_i - omega_j));
}

static scalar ICTCSimilarity(const scalar droptol_i, const scalar droptol_j)
{
    return 1.0 / (1.0 + mag(log10(droptol_i / droptol_j)));
}

static scalar similarityIC
(
    const dictionary& preconDict_i,
    const dictionary& preconDict_j
)
{
    word type_i = preconDict_i.get<word>("preconditioner");
    word type_j = preconDict_j.get<word>("preconditioner");

    // --- DIC has no droptol, treat as 1.0
    scalar droptol_i = (type_i == "DIC") ? 1.0 : preconDict_i.get<scalar>("droptol");
    scalar droptol_j = (type_j == "DIC") ? 1.0 : preconDict_j.get<scalar>("droptol");

    scalar logRatioDroptol = mag(log10(droptol_i / droptol_j));
    return 1.0 / (1.0 + logRatioDroptol);
}

// Similarity between two GAMG configurations. Score is the average of three
// components: smoother similarity, mergeLevels match, and nCellsInCoarsest proximity.
// Each component contributes equally, giving a score in [0, 1].
static scalar similarityGAMG
(
    const dictionary& preconDict_i,
    const dictionary& preconDict_j
)
{
    scalar score = 0.0;

    word smoother_i = preconDict_i.get<word>("smoother");
    word smoother_j = preconDict_j.get<word>("smoother");

    // --- Smoother similarity
    if (
        (smootherIsICTCLike(smoother_i) && smootherIsICTCLike(smoother_j))
     || (smootherIsICTCGaussSeidelLike(smoother_i) && smootherIsICTCGaussSeidelLike(smoother_j))
    )
    {
        score += ICTCSimilarity(smootherToDroptol(preconDict_i.get<word>("smoother")),
                                smootherToDroptol(preconDict_j.get<word>("smoother")));
    }
    else if (
        (smootherIsSORLike(smoother_i) && smootherIsSORLike(smoother_j))
     || (smootherIsDICSORLike(smoother_i) && smootherIsDICSORLike(smoother_j))
    )
    {
        score += omegaSimilarity(smootherToOmega(preconDict_i.get<word>("smoother")), 
                                 smootherToOmega(preconDict_j.get<word>("smoother")));
    }
    else if (smoother_i == smoother_j)
    {
        score += 1.0;
    }

    // --- mergeLevels similarity: binary match
    label mergeLevels_i = preconDict_i.getOrDefault<label>("mergeLevels", 1);
    label mergeLevels_j = preconDict_j.getOrDefault<label>("mergeLevels", 1);
    if (mergeLevels_i == mergeLevels_j)
        score += 1.0;

    // --- nCellsInCoarsest similarity: log-ratio based
    scalar nCells_i = preconDict_i.getOrDefault<label>("nCellsInCoarsestLevel", 10);
    scalar nCells_j = preconDict_j.getOrDefault<label>("nCellsInCoarsestLevel", 10);
    scalar logRatioNCells = mag(log10(nCells_i / nCells_j));
    score += 1.0 / (1.0 + logRatioNCells);

    return score / 3.0;
}

// Similarity along the three subspace-initialization axes, in [0, 1].
static scalar subspaceSimilarity
(
    const dictionary& preconDict_i,
    const dictionary& preconDict_j
)
{
    label lenHistory_i = preconDict_i.getOrDefault<label>("lenHistory", 0);
    label lenHistory_j = preconDict_j.getOrDefault<label>("lenHistory", 0);
    label numProbes_i  = preconDict_i.getOrDefault<label>("numProbes", 4);
    label numProbes_j  = preconDict_j.getOrDefault<label>("numProbes", 4);
    scalar decay_i     = preconDict_i.getOrDefault<scalar>("decayRate", 0);
    scalar decay_j     = preconDict_j.getOrDefault<scalar>("decayRate", 0);

    const bool off_i = subspaceOff(preconDict_i);
    const bool off_j = subspaceOff(preconDict_j);

    if (off_i || off_j) return (off_i && off_j) ? 1.0 : 0.0;

    if ((decay_i > 0) != (decay_j > 0)) return 0.0;

    scalar logRatioNumProbes  = mag(log2(scalar(numProbes_i)/scalar(numProbes_j)));

    const scalar depth_i =
        (decay_i > 0) ? 1.0/max(SMALL, 1.0 - decay_i) : scalar(lenHistory_i);
    const scalar depth_j =
        (decay_j > 0) ? 1.0/max(SMALL, 1.0 - decay_j) : scalar(lenHistory_j);

    scalar logRatioDepth = mag(log2(depth_i/depth_j));

    return 0.5/(1.0 + logRatioDepth) + 0.5/(1.0 + logRatioNumProbes);
}


// Build the full similarity matrix over all preconditioner configurations.
// S[i][j] is the similarity between arm i and arm j, in [0, 1].
// Cross-type similarity (IC vs GAMG) is always 0.
SquareMatrix<scalar> similarityMatrix(
    const List<dictionary>& preconditionerDicts
) {

    label numConfigs = preconditionerDicts.size();
    SquareMatrix<scalar> S(numConfigs, 0.0);

    for (label i = 0; i < numConfigs; ++i){

        S(i, i) = 1.0;
        
        for (label j = i + 1; j < numConfigs; j++) {
            
            const dictionary& dict_i = preconditionerDicts[i];
            const dictionary& dict_j = preconditionerDicts[j];
            scalar similarity = 0.0;

            word type_i = dict_i.get<word>("preconditioner");
            word type_j = dict_j.get<word>("preconditioner");

            bool iGAMG = (type_i == "GAMG" || type_i == "FGAMG");
            bool jGAMG = (type_j == "GAMG" || type_j == "FGAMG");
            bool iIC   = (type_i == "ICTC" || type_i == "DIC");
            bool jIC   = (type_j == "ICTC" || type_j == "DIC");

            if (iGAMG && jGAMG) {
                similarity = similarityGAMG(dict_i, dict_j);
            } else if (iIC && jIC) {
                similarity = similarityIC(dict_i, dict_j);
            }

            similarity *= subspaceSimilarity(dict_i, dict_j);

            S(i, j) = similarity;
            S(j, i) = similarity;
        }
    }

    return S;

}

template<class T>
static dictionary rankDict(DynamicList<T>& values) {
    Foam::sort(values);
    dictionary dict;
    label rank = 0;
    forAll(values, k) {
        word key = Foam::name(values[k]);
        if (!dict.found(key)) {
            dict.add(key, rank++);
        }
    }
    return dict;
}

// Build the full similarity matrix over all preconditioner configurations.
// S[i][j] is the similarity between arm i and arm j, in {0, 1}.
// Cross-type similarity (IC vs GAMG) is always 0.
SquareMatrix<scalar> pathMatrix(
    const List<dictionary>& preconditionerDicts
) {

    label numConfigs = preconditionerDicts.size();
    DynamicList<scalar> droptolList;
    DynamicList<label> nCellsList;
    DynamicList<scalar> smootherList;
    DynamicList<scalar> omegaList;
    DynamicList<label> lenHistoryList;
    DynamicList<label> numProbesList;
    DynamicList<scalar> decayRateList;

    for (label i = 0; i < numConfigs; ++i) {
        const dictionary& dict = preconditionerDicts[i];
        lenHistoryList.append(dict.getOrDefault<label>("lenHistory", 0));
        numProbesList.append(dict.getOrDefault<label>("numProbes", 4));
        decayRateList.append(dict.getOrDefault<scalar>("decayRate", 0));
        word type = dict.get<word>("preconditioner");
        if (type == "ICTC") {
            droptolList.append(dict.get<scalar>("droptol"));
        } else if (type == "DIC") {
            droptolList.append(1.0);
        } else if (type == "GAMG" || type == "FGAMG") {
            nCellsList.append(dict.getOrDefault<label>("nCellsInCoarsestLevel", 10));
            word smoother = dict.get<word>("smoother");
            if (smootherIsICTCLike(smoother) || smootherIsICTCGaussSeidelLike(smoother)) {
                smootherList.append(smootherToDroptol(smoother));
            }
            if (smootherIsSORLike(smoother) || smootherIsDICSORLike(smoother)) {
                omegaList.append(smootherToOmega(smoother));
            }
        }
    }

    dictionary droptolRanks = rankDict(droptolList);
    dictionary nCellsRanks = rankDict(nCellsList);
    dictionary smootherRanks = rankDict(smootherList);
    dictionary omegaRanks = rankDict(omegaList);
    dictionary lenHistoryRanks = rankDict(lenHistoryList);
    dictionary numProbesRanks = rankDict(numProbesList);
    dictionary decayRateRanks = rankDict(decayRateList);

    SquareMatrix<scalar> S(numConfigs, 0.0);

    for (label i = 0; i < numConfigs; ++i) {

        S(i, i) = 1.0;

        for (label j = i + 1; j < numConfigs; j++) {

            const dictionary& dict_i = preconditionerDicts[i];
            const dictionary& dict_j = preconditionerDicts[j];
            scalar adjacent = 0.0;
            label diff = 0;

            auto axisDiff =
                [&](const dictionary& ranks, const word& key, auto defaultValue)
                {
                    typedef decltype(defaultValue) axisType;
                    return mag
                    (
                        ranks.get<label>(name(dict_i.getOrDefault<axisType>(key, defaultValue)))
                      - ranks.get<label>(name(dict_j.getOrDefault<axisType>(key, defaultValue)))
                    );
                };

            diff += axisDiff(lenHistoryRanks, "lenHistory", label(0));
            diff += axisDiff(numProbesRanks,  "numProbes",  label(4));
            diff += axisDiff(decayRateRanks,  "decayRate",  scalar(0));

            word type_i = dict_i.get<word>("preconditioner");
            word type_j = dict_j.get<word>("preconditioner");

            bool iGAMG = (type_i == "GAMG" || type_i == "FGAMG");
            bool jGAMG = (type_j == "GAMG" || type_j == "FGAMG");
            bool iIC   = (type_i == "ICTC" || type_i == "DIC");
            bool jIC   = (type_j == "ICTC" || type_j == "DIC");

            if (iGAMG && jGAMG) {
                label iNC = nCellsRanks.get<label>(name(dict_i.getOrDefault<label>("nCellsInCoarsestLevel", 10)));
                label jNC = nCellsRanks.get<label>(name(dict_j.getOrDefault<label>("nCellsInCoarsestLevel", 10)));
                diff += mag(iNC - jNC);
                word smoother_i = dict_i.get<word>("smoother");
                word smoother_j = dict_j.get<word>("smoother");
                if (
                    (smootherIsICTCLike(smoother_i) && smootherIsICTCLike(smoother_j))
                 || (smootherIsICTCGaussSeidelLike(smoother_i) && smootherIsICTCGaussSeidelLike(smoother_j))
                ) {
                    diff += mag(
                        smootherRanks.get<label>(name(smootherToDroptol(smoother_i))) -
                        smootherRanks.get<label>(name(smootherToDroptol(smoother_j)))
                    );
                } else if (
                    (smootherIsSORLike(smoother_i) && smootherIsSORLike(smoother_j))
                 || (smootherIsDICSORLike(smoother_i) && smootherIsDICSORLike(smoother_j))
                ) {
                    diff += mag(
                        omegaRanks.get<label>(name(smootherToOmega(smoother_i))) -
                        omegaRanks.get<label>(name(smootherToOmega(smoother_j)))
                    );
                } else if (smoother_i != smoother_j) {
                    diff++;
                }
                if (dict_i.getOrDefault<label>("mergeLevels", 1) != dict_j.getOrDefault<label>("mergeLevels", 1)) {
                    diff++;
                }
                if (diff < 2) {
                    adjacent = 1.0;
                }
            } else if (iIC && jIC) {
                label iD = droptolRanks.get<label>(name(dict_i.getOrDefault<scalar>("droptol", 1.0)));
                label jD = droptolRanks.get<label>(name(dict_j.getOrDefault<scalar>("droptol", 1.0)));
                diff += mag(iD - jD);
                if (diff < 2) {
                    adjacent = 1.0;
                }
            }

            S(i, j) = adjacent;
            S(j, i) = adjacent;

        }
    }

    return S;

}


void elementwisePower(
    SquareMatrix<scalar>& S, 
    const scalar power
) {

    if (power != 1.0) {
        for (label i = 0; i < S.n(); i++) {
            for (label j = 0; j < S.n(); j++) {
                S(i, j) = pow(S(i, j), power);
            }
        }
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
