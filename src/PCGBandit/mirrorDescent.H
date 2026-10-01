/*---------------------------------------------------------------------------*\
                Functions for entropic mirror descent updates.
\*---------------------------------------------------------------------------*/

#ifndef MirrorDescent_H
#define MirrorDescent_H
#include "scalarField.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


//- Return exp(logits - max(logits)), normalized to sum to one.
scalarField softmax(const scalarField& logits);

//- Fill probs=(2/(eta*(loss-x)))^2; warm-start and return the Newton offset x.
scalar tsallisINF(
    const scalarField& loss, 
    const scalar eta, 
    scalarField& probs, 
    scalar x
);

//- Fill probs=(eta*(loss-x))^(1/(alpha-1)); alpha=0.5 delegates above.
//  As above, x is a Newton warm start and the return value.
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