/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#ifndef SORSmootherBase_H
#define SORSmootherBase_H

#include "GaussSeidelSmoother.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SORSmootherBase Declaration
\*---------------------------------------------------------------------------*/

class SORSmootherBase
:
    public GaussSeidelSmoother
{
    // Private data

        //- Successive over-relaxation factor
        const scalar omega_;

        //- omega/diag, precomputed once here at construction rather than
        //- re-derived (via a division and an extra multiply-add) on every
        //- cell of every sweep of every application. The matrix and omega
        //- are both fixed for the lifetime of the smoother, so this and the
        //- (1-omega) blend factor are the only omega-dependent quantities
        //- the sweep needs; mirrors the reciprocal-diagonal caching already
        //- used by ICTCSmootherBase/DICPreconditioner.
        solveScalarField rDiagOmega_;


public:

    // Constructors

        //- Construct from matrix components and a fixed relaxation factor
        SORSmootherBase
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls,
            const scalar omega
        );


    // Member Functions

        //- Smooth for the given number of sweeps
        static void smooth
        (
            solveScalarField& psi,
            const lduMatrix& matrix,
            const solveScalarField& source,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const direction cmpt,
            const label nSweeps,
            const scalar omega,
            const solveScalarField& rDiagOmega
        );

        //- Smooth the solution for a given number of sweeps
        virtual void smooth
        (
            solveScalarField& psi,
            const scalarField& source,
            const direction cmpt,
            const label nSweeps
        ) const;

        //- Smooth the solution for a given number of sweeps
        virtual void scalarSmooth
        (
            solveScalarField& psi,
            const solveScalarField& source,
            const direction cmpt,
            const label nSweeps
        ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
