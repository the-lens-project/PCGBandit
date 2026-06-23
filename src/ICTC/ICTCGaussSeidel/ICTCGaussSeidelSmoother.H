/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#ifndef ICTCGaussSeidelSmoother_H
#define ICTCGaussSeidelSmoother_H

#include "../ICTCSmoother/ICTCSmoother.H"
#include "GaussSeidelSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DECLARE_ICTC_GAUSS_SEIDEL_SMOOTHER(ClassName, Type, VariableName)       \
namespace Foam								        \
{									        \
class ClassName 							        \
: 										\
    public lduMatrix::smoother 							\
{ 										\
    Type VariableName;								\
    GaussSeidelSmoother gsSmoother_; 						\
 										\
public: 									\
 										\
        TypeName(#ClassName); 							\
 										\
	ClassName 								\
        ( 									\
            const word& fieldName,					 	\
            const lduMatrix& matrix, 						\
            const FieldField<Field, scalar>& interfaceBouCoeffs, 		\
            const FieldField<Field, scalar>& interfaceIntCoeffs,		\
            const lduInterfaceFieldPtrsList& interfaces, 			\
            const dictionary& solverControls 					\
        ); 									\
 										\
 										\
        virtual void smooth 							\
        ( 									\
            solveScalarField& psi, 						\
            const scalarField& source, 						\
            const direction cmpt, 						\
            const label nSweeps 						\
        ) const; 								\
				 						\
        virtual void scalarSmooth 						\
        ( 									\
            solveScalarField& psi, 						\
            const solveScalarField& source, 					\
            const direction cmpt, 						\
            const label nSweeps	 						\
        ) const; 								\
}; 										\
}										\

DECLARE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m4, ICTC_m4, ictc_m4Smoother_);
DECLARE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m3p5, ICTC_m3p5, ictc_m3p5Smoother_);
DECLARE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m3, ICTC_m3, ictc_m3Smoother_);
DECLARE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m2p5, ICTC_m2p5, ictc_m2p5Smoother_);
DECLARE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m2, ICTC_m2, ictc_m2Smoother_);
DECLARE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m1p5, ICTC_m1p5, ictc_m1p5Smoother_);
DECLARE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m1, ICTC_m1, ictc_m1Smoother_);
DECLARE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel_m0p5, ICTC_m0p5, ictc_m0p5Smoother_);

DECLARE_ICTC_GAUSS_SEIDEL_SMOOTHER(ICTCGaussSeidel, ICTC, ictcSmoother_);


#undef MAKE_ICTC_GAUSS_SEIDEL_SMOOTHER

#endif
// ************************************************************************* //
