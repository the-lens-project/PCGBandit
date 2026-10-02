/*---------------------------------------------------------------------------*\
Class
    Foam::subspaceInitializer

Description
    Correct an initial iterate using a subspace supplied by contextWindow
    or EWMASketch. Galerkin minimizes the error's energy norm for a symmetric
    definite matrix; least squares minimizes the residual's Euclidean norm.
    Here e = psiExact - psi and r = b - A*psi = A*e.

SourceFiles
    subspaceInitializer.C
\*---------------------------------------------------------------------------*/

#ifndef subspaceInitializer_H
#define subspaceInitializer_H

#include "lduMatrix.H"
#include "contextWindow.H"
#include "EWMASketch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class subspaceInitializer Declaration
\*---------------------------------------------------------------------------*/

class subspaceInitializer
{
    //- Objective of the reduced solve
    enum projection
    {
        GALERKIN,       // Minimize error energy norm; symmetric definite A
        LEAST_SQUARES   // Minimize residual Euclidean norm
    };

    const lduMatrix& matrix_;
    const FieldField<Field, scalar>& interfaceBouCoeffs_;
    const lduInterfaceFieldPtrsList& interfaces_;

    //- Null disables initialization for matrices without an fvMesh registry
    const fvMesh* meshPtr_;

    const word fieldName_;

    //- State capacities, also used by the default initialize() overload
    label maxLenHistory_;
    label maxNumProbes_;

    //- Caller-supplied positive rates; each must identify a different sketch
    List<scalar> decayRates_;

    //- Dictionary rate for the default initialize() overload
    scalar decayRate_;

    scalar truncTol_;
    projection projection_;
    label randomSeed_;

    //- Write and resume binary state at simulation write times
    bool persistState_;

    //- Mesh-owned window selected by update()
    mutable contextWindow* window_;

    //- Mesh-owned sketches in decayRates_ order, selected by update()
    mutable List<EWMASketch*> ewmaSketches_;

    //- Prevent duplicate state updates between calls to update()
    mutable bool pushed_;

    //- Parse the projection keyword
    static projection readProjection(const word& name);

    //- Reject initialization requests that name both a window and a stream.
    static void checkExclusive(const label lenHistory, const scalar decayRate);

    //- Find the rate by registry-key rendering; return -1 if undeclared
    label rateIndex(const scalar decayRate) const;

    //- Store the entry iterate once per solve, before correcting psi
    void pushState(const solveScalarField& psi) const;

    //- Apply the projection using a sketch, or the window if sketch is null.
    //  Requires state validated by update(); returns the retained rank or 0.
    label correct
    (
        const direction cmpt,
        const label lenHistory,
        const label numProbes,
        const EWMASketch* sketch,
        solveScalarField& psi,
        solveScalarField& rA
    ) const;

    void operator=(const subspaceInitializer&) = delete;

public:

    subspaceInitializer
    (
        const lduMatrix& matrix,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const word& fieldName,
        const dictionary& controlDict,
        const label maxLenHistory = -1,
        const label maxNumProbes = -1,
        const UList<scalar>& decayRates = UList<scalar>::null()
    );

    //- Find and validate state before each solve without adding psi.
    //  Returns the maximum available depth, or -1 when inactive.
    label update(const solveScalarField& psi, const direction cmpt) const;

    //- Correct psi and rA in place and refresh finalResidual after update().
    void initialize
    (
        const direction cmpt,
        const label lenHistory,
        const label numProbes,
        const scalar decayRate,
        const solveScalar normFactor,
        solveScalarField& psi,
        solveScalarField& rA,
        solverPerformance& solverPerf
    ) const;

    //- Initialize using the configured capacities and dictionary decay rate
    void initialize
    (
        const direction cmpt,
        const solveScalar normFactor,
        solveScalarField& psi,
        solveScalarField& rA,
        solverPerformance& solverPerf
    ) const
    {
        initialize
        (
            cmpt, maxLenHistory_, maxNumProbes_, decayRate_,
            normFactor, psi, rA, solverPerf
        );
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
