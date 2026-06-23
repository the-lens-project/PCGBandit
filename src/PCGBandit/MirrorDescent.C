/*---------------------------------------------------------------------------*\
                Functions for entropic mirror descent updates.
\*---------------------------------------------------------------------------*/

#include "MirrorDescent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

scalarField softmax(const scalarField& logits) {

    scalarField probs = exp(logits - max(logits));
    probs /= sum(probs);
    return probs;

}

const label maxNewtonIter = 100;
const scalar tolerance = 1e-8;

scalar tsallisINF(
    const scalarField& loss, 
    const scalar eta, 
    scalarField& probs, 
    scalar x
) {

    scalar minLoss = min(loss);
    scalar update;
    label i;
    for (i = 0; i < maxNewtonIter; i++) {

        if (x >= minLoss) {
            x = minLoss - 1.0;
        }
        
        probs = 2.0 / (eta * (loss - x));
        probs *= probs;
        update = (sum(probs) - 1.0) / (eta * sum(pow(probs, 1.5)));
        x -= update;
        
        if (mag(update) < tolerance) {
            break;
        }

    }

    if (i == maxNewtonIter) {
        Info<< "Newton solver did not converge: " << update << endl;
    }
    return x;

}

scalar tsallisINF(
    const scalarField& loss, 
    const scalar eta, 
    scalarField& probs, 
    scalar x, 
    const scalar alpha
) {

    if (alpha == 0.5) {
        return tsallisINF(loss, eta, probs, x);
    }

    scalar minLoss = min(loss);
    scalar update;
    label i;
    for (i = 0; i < maxNewtonIter; i++) {
        
        if (x >= minLoss) {
            x = minLoss - 1.0;
        }

        probs = pow(eta * (loss - x), 1.0 / (alpha - 1.0));
        update = (1.0 - alpha) * (sum(probs) - 1.0) / (eta * sum(pow(probs, 2.0 - alpha)));
        x -= update;
        
        if (mag(update) < tolerance) {
            break;
        }

    }

    if (i == maxNewtonIter) {
        Info<< "Newton solver did not converge: " << update << endl;
    }
    return x;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
