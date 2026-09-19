/*---------------------------------------------------------------------------*\
                Functions for entropic mirror descent updates.
\*---------------------------------------------------------------------------*/

#ifndef MirrorDescent_H
#define MirrorDescent_H
#include "scalarField.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


scalarField softmax(const scalarField& logits);

scalar tsallisINF(
    const scalarField& loss, 
    const scalar eta, 
    scalarField& probs, 
    scalar x
);

scalar tsallisINF(
    const scalarField& loss, 
    const scalar eta, 
    scalarField& probs, 
    const scalar x,
    const scalar alpha
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
// ************************************************************************* //