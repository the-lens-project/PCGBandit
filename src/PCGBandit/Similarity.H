/*---------------------------------------------------------------------------*\
        Functions to build matrices of preconditioner similarities.
\*---------------------------------------------------------------------------*/

#ifndef Similarity_H
#define Similarity_H
#include "SquareMatrix.H"

namespace Foam
{

SquareMatrix<scalar> similarityMatrix(
    const List<dictionary>& preconditionerDicts
);

SquareMatrix<scalar> pathMatrix(
    const List<dictionary>& preconditionerDicts
);

void elementwisePower(
    SquareMatrix<scalar>& S, 
    const scalar power
);

} // End namespace Foam

#endif

// ************************************************************************* //