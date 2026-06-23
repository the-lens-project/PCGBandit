/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#ifndef ICTCSmoother_H
#define ICTCSmoother_H

#include "lduMatrix.H"
#include "ICTCSmootherBase.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Macro to Declare ICTCSmoother Variants
\*---------------------------------------------------------------------------*/

#define DECLARE_ICTC_SMOOTHER(ClassName)                                       \
class ClassName : public ICTCSmootherBase                                      \
{                                                                              \
public:                                                                        \
                                                                               \
    TypeName(#ClassName);                                                      \
                                                                               \
    ClassName                                                                  \
    (                                                                          \
        const word& fieldName,                                                 \
        const lduMatrix& matrix,                                               \
        const FieldField<Field, scalar>& interfaceBouCoeffs,                   \
        const FieldField<Field, scalar>& interfaceIntCoeffs,                   \
        const lduInterfaceFieldPtrsList& interfaces,                           \
        const dictionary&                                                      \
    );                                                                         \
};

// Declaring ICTCSmoother variants

DECLARE_ICTC_SMOOTHER(ICTC_m4)
DECLARE_ICTC_SMOOTHER(ICTC_m3p5)
DECLARE_ICTC_SMOOTHER(ICTC_m3)
DECLARE_ICTC_SMOOTHER(ICTC_m2p5)
DECLARE_ICTC_SMOOTHER(ICTC_m2)
DECLARE_ICTC_SMOOTHER(ICTC_m1p5)
DECLARE_ICTC_SMOOTHER(ICTC_m1)
DECLARE_ICTC_SMOOTHER(ICTC_m0p5)

// Declaring default ICTCSmoother
DECLARE_ICTC_SMOOTHER(ICTC)

#undef DECLARE_ICTC_SMOOTHER

} // End namespace Foam

#endif
// ************************************************************************* //

