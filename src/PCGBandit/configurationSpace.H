/*---------------------------------------------------------------------------*\
Class
    Foam::configurationSpace

Description
    The arm list for one bandit: every preconditioner configuration the learner
    may draw, as a dictionary the solver can hand straight to
    lduMatrix::preconditioner::New.

SourceFiles
    configurationSpace.C
\*---------------------------------------------------------------------------*/

#ifndef configurationSpace_H
#define configurationSpace_H

#include "lduMatrix.H"
#include "dictionary.H"
#include "DynamicList.H"

namespace Foam
{

// --- subspaceInitialization configuration space defaults
inline const List<label> numProbesDefault({4, 8, 12, 16});
inline const List<label> lenHistoryDefault({0, 4, 8, 12, 16});
inline const List<scalar> decayRateDefault({0.0, 0.25, 0.5, 0.75});

//  Whether the arm's subspace initialization is on or off
inline bool subspaceOff(const dictionary& armDict)
{
    return armDict.getOrDefault<label>("numProbes", numProbesDefault[0]) <= 0
        || (armDict.getOrDefault<label>("lenHistory", lenHistoryDefault[0]) <= 0
            && armDict.getOrDefault<scalar>("decayRate", decayRateDefault[0]) <= 0);
}


/*---------------------------------------------------------------------------*\
                     Class configurationSpace Declaration
\*---------------------------------------------------------------------------*/

class configurationSpace
{
    List<dictionary> dicts_;

    //- Largest subspace request any arm can make, so one window serves all
    label maxLenHistory_ = 0;
    label maxNumProbes_ = 0;

    //- Positive rates before arm expansion; input rates must identify distinct sketches
    DynamicList<scalar> decayRates_;

    //- One ICTC arm per drop tolerance, log-spaced
    void appendICTC(const dictionary& solverControls);

    //- The single unparameterised DIC arm
    void appendDIC(const dictionary& solverControls);

    //- The cross-product over whichever GAMG <param>Tune axes are set
    void appendGAMG
    (
        const dictionary& solverControls,
        const lduMatrix& matrix,
        const label staticArm
    );

    //- Multiply the subspace initialization axes onto every arm above
    void expandInitializations(const dictionary& solverControls);

public:

    configurationSpace() = default;

    configurationSpace
    (
        const dictionary& solverControls,
        const lduMatrix& matrix,
        const label staticArm
    );

    const List<dictionary>& dicts() const noexcept { return dicts_; }

    label maxLenHistory() const noexcept { return maxLenHistory_; }
    label maxNumProbes() const noexcept { return maxNumProbes_; }
    const UList<scalar>& decayRates() const noexcept { return decayRates_; }
};

} // End namespace Foam

#endif

// ************************************************************************* //
