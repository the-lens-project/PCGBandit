/*---------------------------------------------------------------------------*\
    Preconditioned conjugate gradient with a callable matrix operator.
\*---------------------------------------------------------------------------*/

#ifndef conjugateGradient_H
#define conjugateGradient_H

#include "scalarField.H"
#include "error.H"

namespace Foam
{

template<class MatrixOperator, class Preconditioner>
bool conjugateGradient
(
    const MatrixOperator& A,
    const scalarField& b,
    const Preconditioner& precondition,
    scalarField& x,
    const scalar tolerance = 1e-6,
    const label maxCGIter = 100
)
{

    const scalar toleranceSqr = tolerance * tolerance;
    x.setSize(b.size());
    x = 0;
    scalarField q(b.size());
    A(x, q);
    scalarField r(b - q);
    scalar residualSqr = sumSqr(r);
    if (residualSqr <= toleranceSqr) return true;

    scalarField s(b.size());
    scalar rs = 0;
    scalarField d(b.size());
    for (label i = 0; i < maxCGIter; ++i)
    {
        precondition(r, s);
        const scalar denom = rs;
        rs = sumProd(r, s);
        if (i == 0) d = s;
        else {
            d *= rs / denom;
            d += s;
        }
        A(d, q);
        const scalar alpha = rs / sumProd(d, q);
        x += alpha * d;
        if ((i+1) % 25 == 0) {
            A(x, q);
            r = b - q;
        } else {
            r -= alpha * q;
        }
        residualSqr = sumSqr(r);
        if (residualSqr <= toleranceSqr) return true;
    }

    WarningInFunction
        << "CG did not converge"
        << "; residual = " << sqrt(residualSqr)
        << "; tolerance = " << tolerance << endl;
    return false;
}

} // End namespace Foam

#endif

// ************************************************************************* //
