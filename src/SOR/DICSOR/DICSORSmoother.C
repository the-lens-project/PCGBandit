/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "DICSORSmoother.H"

/*---------------------------------------------------------------------------*\
                   Macro to Define DICSORSmoother Variants
\*---------------------------------------------------------------------------*/

#define MAKE_DIC_SOR_SMOOTHER(ClassName, VariableName)                         \
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
    const dictionary& solverControls                                           \
)                                                                              \
:                                                                              \
    lduMatrix::smoother                                                        \
    (                                                                          \
        fieldName,                                                             \
        matrix,                                                                \
        interfaceBouCoeffs,                                                    \
        interfaceIntCoeffs,                                                    \
        interfaces                                                             \
    ),                                                                         \
    dicSmoother_                                                               \
    (                                                                          \
        fieldName,                                                             \
        matrix,                                                                \
        interfaceBouCoeffs,                                                    \
        interfaceIntCoeffs,                                                    \
        interfaces,                                                            \
        solverControls                                                         \
    ),                                                                         \
    VariableName                                                               \
    (                                                                          \
        fieldName,                                                             \
        matrix,                                                                \
        interfaceBouCoeffs,                                                    \
        interfaceIntCoeffs,                                                    \
        interfaces,                                                            \
        solverControls                                                         \
    )                                                                          \
{}                                                                             \
                                                                               \
void Foam::ClassName::smooth                                                   \
(                                                                              \
    solveScalarField& psi,                                                     \
    const scalarField& source,                                                 \
    const direction cmpt,                                                      \
    const label nSweeps                                                        \
) const                                                                        \
{                                                                              \
    dicSmoother_.smooth(psi, source, cmpt, nSweeps);                           \
    VariableName.smooth(psi, source, cmpt, nSweeps);                           \
}                                                                              \
                                                                               \
void Foam::ClassName::scalarSmooth                                             \
(                                                                              \
    solveScalarField& psi,                                                     \
    const solveScalarField& source,                                            \
    const direction cmpt,                                                      \
    const label nSweeps                                                        \
) const                                                                        \
{                                                                              \
    dicSmoother_.scalarSmooth(psi, source, cmpt, nSweeps);                     \
    VariableName.scalarSmooth(psi, source, cmpt, nSweeps);                     \
}

MAKE_DIC_SOR_SMOOTHER(DICSOR_p0p1, sor_p0p1Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p0p2, sor_p0p2Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p0p3, sor_p0p3Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p0p4, sor_p0p4Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p0p5, sor_p0p5Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p0p6, sor_p0p6Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p0p7, sor_p0p7Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p0p8, sor_p0p8Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p0p9, sor_p0p9Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p1p0, sor_p1p0Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p1p1, sor_p1p1Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p1p2, sor_p1p2Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p1p3, sor_p1p3Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p1p4, sor_p1p4Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p1p5, sor_p1p5Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p1p6, sor_p1p6Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p1p7, sor_p1p7Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p1p8, sor_p1p8Smoother_)
MAKE_DIC_SOR_SMOOTHER(DICSOR_p1p9, sor_p1p9Smoother_)

#undef MAKE_DIC_SOR_SMOOTHER

// ************************************************************************* //
