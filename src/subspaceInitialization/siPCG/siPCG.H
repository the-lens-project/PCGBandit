/*---------------------------------------------------------------------------*\
Class
    Foam::siPCG

Group
    grpLduMatrixSolvers

Description
    OpenFOAM PCG with subspace initialization before the first iteration
    (arXiv:2309.02156, arXiv:2511.18071). Reads numProbes, lenHistory,
    decayRate, and persistState from its fvSolution entry.

SourceFiles
    siPCG.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_siPCG_H
#define Foam_siPCG_H

#include "lduMatrix.H"
#include "subspaceInitializer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class siPCG Declaration
\*---------------------------------------------------------------------------*/

class siPCG
:
    public lduMatrix::solver
{
    // Private Member Data

        //- Cached preconditioner
        mutable autoPtr<lduMatrix::preconditioner> preconPtr_;


    // Private Member Functions

        //- No copy construct
        siPCG(const siPCG&) = delete;

        //- No copy assignment
        void operator=(const siPCG&) = delete;

        //- Uses field-keyed mesh state shared by p and pFinal across solves.
        //  Disabled unless lenHistory or decayRate is set.
        subspaceInitializer subspace_;


public:

    //- Runtime type information
    TypeName("siPCG");


    // Constructors

        //- Construct from matrix components and solver controls
        siPCG
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~siPCG() = default;


    // Member Functions

        //- Solve the matrix with this solver
        virtual solverPerformance scalarSolve
        (
            solveScalarField& psi,
            const solveScalarField& source,
            const direction cmpt=0
        ) const;

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
