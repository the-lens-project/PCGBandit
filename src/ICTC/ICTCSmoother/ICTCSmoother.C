/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "ICTCSmootherBase.H"
#include "ICTCSmoother.H"
#include "PrecisionAdaptor.H"
#include <algorithm>

#define MAKE_PARAMATERIZED_ICTC_SMOOTHER(ClassName, dropTol)                   \
namespace Foam                                                                 \
{                                                                              \
    defineTypeNameAndDebug(ClassName, 0);                                      \
                                                                               \
    lduMatrix::smoother::                                                      \
    addsymMatrixConstructorToTable<ClassName>                                  \
        add##ClassName##SymMatrixConstructorToTable_;                          \
}                                                                              \
                                                                               \
Foam::ClassName::ClassName                                                     \
(                                                                              \
    const word& fieldName,                                                     \
    const lduMatrix& matrix,                                                   \
    const FieldField<Field, scalar>& interfaceBouCoeffs,                       \
    const FieldField<Field, scalar>& interfaceIntCoeffs,                       \
    const lduInterfaceFieldPtrsList& interfaces,                               \
    const dictionary& 	                                                       \
)                                                                              \
:                                                                              \
    ICTCSmootherBase                                                           \
    (                                                                          \
        fieldName,                                                             \
        matrix,                                                                \
        interfaceBouCoeffs,                                                    \
        interfaceIntCoeffs,                                                    \
        interfaces,                                                            \
        dropTol                                                                \
    )                                                                          \
{}                                                                             \


// Instantiates ICTCSmoother Variants

MAKE_PARAMATERIZED_ICTC_SMOOTHER(ICTC_m4,    1e-4);
MAKE_PARAMATERIZED_ICTC_SMOOTHER(ICTC_m3p5,  pow(10.0, -3.5));
MAKE_PARAMATERIZED_ICTC_SMOOTHER(ICTC_m3,    1e-3);
MAKE_PARAMATERIZED_ICTC_SMOOTHER(ICTC_m2p5,  pow(10.0, -2.5));
MAKE_PARAMATERIZED_ICTC_SMOOTHER(ICTC_m2,    1e-2);
MAKE_PARAMATERIZED_ICTC_SMOOTHER(ICTC_m1p5,  pow(10.0, -1.5));
MAKE_PARAMATERIZED_ICTC_SMOOTHER(ICTC_m1,    1e-1);
MAKE_PARAMATERIZED_ICTC_SMOOTHER(ICTC_m0p5,  pow(10.0, -0.5));

// Default Smoother Instantiation
MAKE_PARAMATERIZED_ICTC_SMOOTHER(ICTC,  1);

#undef MAKE_PARAMATERIZED_ICTC_SMOOTHER

// ************************************************************************* //

