/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#ifndef ICTCSmootherBase_H
#define ICTCSmootherBase_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ICTCSmootherBase Declaration
\*---------------------------------------------------------------------------*/

class ICTCSmootherBase
:
    public lduMatrix::smoother
{
     // Private data
  
        //- The reciprocal preconditioned diagonal
        solveScalarField diagL_;
        solveScalarField lowerL_;
        labelField rowAddrL_;
        labelField colPtrL_;
	solveScalarField rD_; // for debugging purposes
  
public:
  
     // Constructors
  
         //- Construct from matrix components
         ICTCSmootherBase
         (
             const word& fieldName,
             const lduMatrix& matrix,
             const FieldField<Field, scalar>& interfaceBouCoeffs,
             const FieldField<Field, scalar>& interfaceIntCoeffs,
             const lduInterfaceFieldPtrsList& interfaces,
	     const scalar droptol
         );

    // Member Functions
  
        //- Smooth the solution for a given number of sweeps
	void smooth
        (
             solveScalarField& psi,
             const scalarField& source,
             const direction cmpt,
             const label nSweeps
        ) const;

	void scalarSmooth 
	(
  	     solveScalarField& psi,
             const solveScalarField& source,
             const direction cmpt,
             const label nSweeps
         ) const;

};
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
} // End namespace Foam
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
#endif
  
// ************************************************************************* //
