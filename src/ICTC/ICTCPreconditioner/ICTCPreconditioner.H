/*---------------------------------------------------------------------------*\

Class
    Foam::ICTCPreconditioner

Group
    grpLduMatrixPreconditioners

Description
    A preconditioner for symmetric matrices that uses column-based 
    incomplete Cholesky factorization by threshold. Based in part on 
    Edmond Chow's ICT package, ported from Fortran and to work with 
    OpenFOAM's sparse matrix format.

SourceFiles
    ICTCPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef ICTCPreconditioner_H
#define ICTCPreconditioner_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ICTCPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class ICTCPreconditioner
:
    public lduMatrix::preconditioner
{
    // Private data

        //- The reciprocal preconditioned diagonal
        solveScalarField diagL_;
        solveScalarField lowerL_;
        labelField rowAddrL_;
        labelField colPtrL_;


public:

    //- Runtime type information
    TypeName("ICTC");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        ICTCPreconditioner
        (
            const lduMatrix::solver&,
            const dictionary& 
        );


    //- Destructor
    virtual ~ICTCPreconditioner() = default;


    // Member Functions

        //- Compute factorization
        static label calcL(solveScalarField&, solveScalarField&, labelField&, labelField&, const lduMatrix&, const scalar, const scalar);

        //- Return wA the preconditioned form of residual rA
        virtual void precondition
        (
            solveScalarField& wA,
            const solveScalarField& rA,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
