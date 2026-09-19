/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#ifndef DICSORSmoother_H
#define DICSORSmoother_H

#include "../SORSmoother/SORSmoother.H"
#include "DICSmoother.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Macro to Declare DICSORSmoother Variants
\*---------------------------------------------------------------------------*/

#define DECLARE_DIC_SOR_SMOOTHER(ClassName, SORType, VariableName)             \
class ClassName                                                               \
:                                                                             \
    public lduMatrix::smoother                                                 \
{                                                                             \
    DICSmoother dicSmoother_;                                                  \
    SORType VariableName;                                                      \
                                                                              \
public:                                                                       \
                                                                              \
    TypeName(#ClassName);                                                      \
                                                                              \
    ClassName                                                                 \
    (                                                                         \
        const word& fieldName,                                                 \
        const lduMatrix& matrix,                                               \
        const FieldField<Field, scalar>& interfaceBouCoeffs,                   \
        const FieldField<Field, scalar>& interfaceIntCoeffs,                   \
        const lduInterfaceFieldPtrsList& interfaces,                           \
        const dictionary& solverControls                                      \
    );                                                                        \
                                                                              \
    virtual void smooth                                                       \
    (                                                                         \
        solveScalarField& psi,                                                 \
        const scalarField& source,                                             \
        const direction cmpt,                                                  \
        const label nSweeps                                                    \
    ) const;                                                                  \
                                                                              \
    virtual void scalarSmooth                                                 \
    (                                                                         \
        solveScalarField& psi,                                                 \
        const solveScalarField& source,                                        \
        const direction cmpt,                                                  \
        const label nSweeps                                                    \
    ) const;                                                                  \
};

DECLARE_DIC_SOR_SMOOTHER(DICSOR_p0p1, SOR_p0p1, sor_p0p1Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p0p2, SOR_p0p2, sor_p0p2Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p0p3, SOR_p0p3, sor_p0p3Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p0p4, SOR_p0p4, sor_p0p4Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p0p5, SOR_p0p5, sor_p0p5Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p0p6, SOR_p0p6, sor_p0p6Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p0p7, SOR_p0p7, sor_p0p7Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p0p8, SOR_p0p8, sor_p0p8Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p0p9, SOR_p0p9, sor_p0p9Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p1p0, SOR_p1p0, sor_p1p0Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p1p1, SOR_p1p1, sor_p1p1Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p1p2, SOR_p1p2, sor_p1p2Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p1p3, SOR_p1p3, sor_p1p3Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p1p4, SOR_p1p4, sor_p1p4Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p1p5, SOR_p1p5, sor_p1p5Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p1p6, SOR_p1p6, sor_p1p6Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p1p7, SOR_p1p7, sor_p1p7Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p1p8, SOR_p1p8, sor_p1p8Smoother_)
DECLARE_DIC_SOR_SMOOTHER(DICSOR_p1p9, SOR_p1p9, sor_p1p9Smoother_)

#undef DECLARE_DIC_SOR_SMOOTHER

} // End namespace Foam

#endif

// ************************************************************************* //
