/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#ifndef SORSmoother_H
#define SORSmoother_H

#include "SORSmootherBase.H"

namespace Foam
{

#define DECLARE_SOR_SMOOTHER(ClassName)                                      \
class ClassName                                                              \
:                                                                            \
    public SORSmootherBase                                                   \
{                                                                            \
public:                                                                      \
                                                                             \
    TypeName(#ClassName);                                                    \
                                                                             \
    ClassName                                                                \
    (                                                                        \
        const word& fieldName,                                               \
        const lduMatrix& matrix,                                             \
        const FieldField<Field, scalar>& interfaceBouCoeffs,                 \
        const FieldField<Field, scalar>& interfaceIntCoeffs,                 \
        const lduInterfaceFieldPtrsList& interfaces,                         \
        const dictionary& solverControls                                    \
    );                                                                       \
};

DECLARE_SOR_SMOOTHER(SOR_p0p1)
DECLARE_SOR_SMOOTHER(SOR_p0p2)
DECLARE_SOR_SMOOTHER(SOR_p0p3)
DECLARE_SOR_SMOOTHER(SOR_p0p4)
DECLARE_SOR_SMOOTHER(SOR_p0p5)
DECLARE_SOR_SMOOTHER(SOR_p0p6)
DECLARE_SOR_SMOOTHER(SOR_p0p7)
DECLARE_SOR_SMOOTHER(SOR_p0p8)
DECLARE_SOR_SMOOTHER(SOR_p0p9)
DECLARE_SOR_SMOOTHER(SOR_p1p0)
DECLARE_SOR_SMOOTHER(SOR_p1p1)
DECLARE_SOR_SMOOTHER(SOR_p1p2)
DECLARE_SOR_SMOOTHER(SOR_p1p3)
DECLARE_SOR_SMOOTHER(SOR_p1p4)
DECLARE_SOR_SMOOTHER(SOR_p1p5)
DECLARE_SOR_SMOOTHER(SOR_p1p6)
DECLARE_SOR_SMOOTHER(SOR_p1p7)
DECLARE_SOR_SMOOTHER(SOR_p1p8)
DECLARE_SOR_SMOOTHER(SOR_p1p9)

#undef DECLARE_SOR_SMOOTHER

} // End namespace Foam

#endif

// ************************************************************************* //
