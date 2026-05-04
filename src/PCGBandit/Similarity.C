/*---------------------------------------------------------------------------*\
        Functions to build matrices of preconditioner similarities.
\*---------------------------------------------------------------------------*/

#include "Similarity.H"

namespace Foam
{

// Parse the droptol scalar from a word-encoded smoother name.
// e.g. "ICTCGaussSeidel_m3p5" -> 10^-3.5, "DICGaussSeidel" -> 1.0 (no suffix)
static scalar smootherToDroptol(const word& smootherName) {
    
    label underscoreIdx = smootherName.rfind('_');
    if (underscoreIdx == -1) return 1.0;

    word suffix = smootherName.substr(underscoreIdx + 1); // e.g. "m3p5"
    word magnitudeStr = suffix.substr(1);                 // remove leading 'm', e.g. "3p5"

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

// Returns true if the smoother is ICTC-like (i.e. has a droptol parameter).
static bool smootherIsICTCLike(const word& smootherName) {
    return smootherName.startsWith("ICTC_") || smootherName == "DIC";
}

// Extract the effective (droptol, coarsestDroptol) pair from a GAMG dict.
// For DIC, both are 1.0. If no coarsestSmoother is specified, defaults
// to the same droptol as the main smoother.
static Pair<scalar> effectiveICTCParams (
    const dictionary& preconDict
) {
    word smoother = preconDict.get<word>("smoother");

    if (smoother == "DIC")
        return Pair<scalar>(1.0, 1.0);

    scalar smootherDroptol = smootherToDroptol(smoother);
    scalar coarsestDroptol = smootherDroptol;

    if (preconDict.found("coarsestSmoother"))
        coarsestDroptol = smootherToDroptol(preconDict.get<word>("coarsestSmoother"));

    return Pair<scalar>(smootherDroptol, coarsestDroptol);
}

// Similarity between two ICTC-like smoothers based on log-ratio of droptols.
// Each component contributes 0.5, giving a score in [0, 1].
static scalar ICTCSimilarity
(
    const scalar droptol_i,
    const scalar coarsestDroptol_i,
    const scalar droptol_j,
    const scalar coarsestDroptol_j
)
{
    scalar logRatioDroptol  = mag(log10(droptol_i / droptol_j));
    scalar logRatioCoarsest = mag(log10(coarsestDroptol_i / coarsestDroptol_j));
    return 0.5 / (1.0 + logRatioDroptol) + 0.5 / (1.0 + logRatioCoarsest);
}

// Similarity between two standalone IC preconditioners based on log-ratio of droptols.
// Score in [0, 1], equal to 1 when droptols are identical.
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
    if (smootherIsICTCLike(smoother_i) && smootherIsICTCLike(smoother_j))
    {
        Pair<scalar> ICTCParams_i = effectiveICTCParams(preconDict_i);
        Pair<scalar> ICTCParams_j = effectiveICTCParams(preconDict_j);
        score += ICTCSimilarity
        (
            ICTCParams_i.first(), ICTCParams_i.second(),
            ICTCParams_j.first(), ICTCParams_j.second()
        );
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
    DynamicList<scalar> coarsestSmootherList;
    
    for (label i = 0; i < numConfigs; ++i) {
        const dictionary& dict = preconditionerDicts[i];
        word type = dict.get<word>("preconditioner");
        if (type == "ICTC") {
            droptolList.append(dict.get<scalar>("droptol"));
        } else if (type == "DIC") {
            droptolList.append(1.0);
        } else if (type == "GAMG" || type == "FGAMG") {
            nCellsList.append(dict.getOrDefault<label>("nCellsInCoarsestLevel", 10));
            word smoother = dict.get<word>("smoother");
            if (smootherIsICTCLike(smoother)) {
                Pair<scalar> pair = effectiveICTCParams(dict);
                smootherList.append(pair.first());
                coarsestSmootherList.append(pair.second());
            }
        }
    }

    dictionary droptolRanks = rankDict(droptolList);
    dictionary nCellsRanks = rankDict(nCellsList);
    dictionary smootherRanks = rankDict(smootherList);
    dictionary coarsestSmootherRanks = rankDict(coarsestSmootherList);

    SquareMatrix<scalar> S(numConfigs, 0.0);

    for (label i = 0; i < numConfigs; ++i) {

        S(i, i) = 1.0;

        for (label j = i + 1; j < numConfigs; j++) {

            const dictionary& dict_i = preconditionerDicts[i];
            const dictionary& dict_j = preconditionerDicts[j];
            scalar adjacent = 0.0;
            label diff = 0;

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
                if (smootherIsICTCLike(smoother_i) && smootherIsICTCLike(smoother_j)) {
                    Pair<scalar> pair_i = effectiveICTCParams(dict_i);
                    Pair<scalar> pair_j = effectiveICTCParams(dict_j);
                    diff += mag(
                        smootherRanks.get<label>(name(pair_i.first())) - 
                        smootherRanks.get<label>(name(pair_j.first()))
                    );
                    diff += mag(
                        coarsestSmootherRanks.get<label>(name(pair_i.second())) - 
                        coarsestSmootherRanks.get<label>(name(pair_j.second()))
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
