/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "ICTCGaussSeidelSmoother.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER(ClassName, VariableName)		 \
namespace Foam									 \
{ 										 \
    defineTypeNameAndDebug(ClassName, 0); 					 \
 										 \
    lduMatrix::smoother::addsymMatrixConstructorToTable<ClassName> 		 \
        add##ClassName##SymMatrixConstructorToTable_; 				 \
} 										 \
             									 \
 										 \
 										 \
Foam::ClassName::ClassName 							 \
( 										 \
    const word& fieldName, 							 \
    const lduMatrix& matrix, 							 \
    const FieldField<Field, scalar>& interfaceBouCoeffs, 			 \
    const FieldField<Field, scalar>& interfaceIntCoeffs, 			 \
    const lduInterfaceFieldPtrsList& interfaces, 				 \
    const dictionary& solverControls 						 \
) 										 \
: 										 \
    lduMatrix::smoother 							 \
    ( 										 \
        fieldName, 								 \
        matrix, 								 \
        interfaceBouCoeffs, 							 \
        interfaceIntCoeffs, 							 \
        interfaces 								 \
    ), 										 \
    VariableName								 \
    ( 										 \
    	fieldName, 								 \
    	matrix, 								 \
    	interfaceBouCoeffs, 							 \
    	interfaceIntCoeffs,  							 \
    	interfaces, 								 \
    	solverControls 								 \
    ), 										 \
    gsSmoother_ 								 \
    ( 										 \
        fieldName, 								 \
        matrix, 								 \
        interfaceBouCoeffs, 							 \
        interfaceIntCoeffs, 							 \
        interfaces, 								 \
        solverControls 								 \
    ) 										 \
{} 										 \
 										 \
void Foam::ClassName::smooth 							 \
( 										 \
    solveScalarField& psi, 							 \
    const scalarField& source, 							 \
    const direction cmpt,							 \
    const label nSweeps								 \
) const										 \
{										 \
    VariableName.smooth(psi, source, cmpt, nSweeps);				 \
    gsSmoother_.smooth(psi, source, cmpt, nSweeps);				 \
}										 \
										 \
										 \
void Foam::ClassName::scalarSmooth						 \
(										 \
    solveScalarField& psi,							 \
    const solveScalarField& source,						 \
    const direction cmpt,							 \
    const label nSweeps								 \
) const										 \
{ 										 \
    VariableName.scalarSmooth(psi, source, cmpt, nSweeps);			 \
    gsSmoother_.scalarSmooth(psi, source, cmpt, nSweeps);			 \
}										 \

// Instantiate parameterized ICTCGaussSeidel Variants
MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m4, ictc_m4Smoother_);
MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m3p5, ictc_m3p5Smoother_);
MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m3, ictc_m3Smoother_);
MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m2p5, ictc_m2p5Smoother_);
MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m2, ictc_m2Smoother_);
MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m1p5, ictc_m1p5Smoother_);
MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m1, ictc_m1Smoother_);
MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m0p5, ictc_m0p5Smoother_);

// Instantiate default ICTCGaussSeidel Variant
MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel, ictcSmoother_);

#undef MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER

// ************************************************************************* //
