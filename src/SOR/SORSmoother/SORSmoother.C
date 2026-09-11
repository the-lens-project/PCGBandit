/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "SORSmoother.H"

#define MAKE_PARAMETERIZED_SOR_SMOOTHER(ClassName, Omega)                    \
namespace Foam                                                               \
{                                                                            \
    defineTypeNameAndDebug(ClassName, 0);                                    \
                                                                             \
    lduMatrix::smoother::addsymMatrixConstructorToTable<ClassName>           \
        add##ClassName##SymMatrixConstructorToTable_;                        \
                                                                             \
    lduMatrix::smoother::addasymMatrixConstructorToTable<ClassName>          \
        add##ClassName##AsymMatrixConstructorToTable_;                       \
}                                                                            \
                                                                             \
Foam::ClassName::ClassName                                                   \
(                                                                            \
    const word& fieldName,                                                   \
    const lduMatrix& matrix,                                                 \
    const FieldField<Field, scalar>& interfaceBouCoeffs,                     \
    const FieldField<Field, scalar>& interfaceIntCoeffs,                     \
    const lduInterfaceFieldPtrsList& interfaces,                             \
    const dictionary& solverControls                                        \
)                                                                            \
:                                                                            \
    SORSmootherBase                                                          \
    (                                                                        \
        fieldName,                                                           \
        matrix,                                                              \
        interfaceBouCoeffs,                                                  \
        interfaceIntCoeffs,                                                  \
        interfaces,                                                         \
        solverControls,                                                      \
        Omega                                                                \
    )                                                                        \
{}

MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p0p1, 0.1);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p0p2, 0.2);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p0p3, 0.3);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p0p4, 0.4);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p0p5, 0.5);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p0p6, 0.6);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p0p7, 0.7);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p0p8, 0.8);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p0p9, 0.9);

MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p1p1, 1.1);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p1p2, 1.2);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p1p3, 1.3);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p1p4, 1.4);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p1p5, 1.5);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p1p6, 1.6);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p1p7, 1.7);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p1p8, 1.8);
MAKE_PARAMETERIZED_SOR_SMOOTHER(SOR_p1p9, 1.9);

#undef MAKE_PARAMETERIZED_SOR_SMOOTHER

// ************************************************************************* //
