/*---------------------------------------------------------------------------*\
Function
    Foam::preconditionerSimilarityMatrix

Description
    Build the similarity matrix over all preconditioner configurations.
\*---------------------------------------------------------------------------*/

#ifndef similarityMatrix_H
#define similarityMatrix_H

#include "SquareMatrix.H"

namespace Foam
{

SquareMatrix<scalar> preconditionerSimilarityMatrix
(
    const List<dictionary>& preconditionerDicts
);

} // End namespace Foam

#endif

// ************************************************************************* //