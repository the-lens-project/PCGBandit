/*---------------------------------------------------------------------------*\
                Functions for entropic mirror descent updates.
\*---------------------------------------------------------------------------*/

#include "mirrorDescent.H"

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
        WarningInFunction << "Newton solver did not converge:"  << update << endl;
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

    if (alpha == 1.0) {
        probs = softmax(-eta * loss);
        return x;
    }

    // - The general alpha implementation uses the Tsallis entropy formula from Ito (2024),
    //   which differs from the formula in Zimmert & Seldin (2021) by a 1 / alpha factor.
    //   We thus divide the step-size by alpha when passing to the alpha = 0.5 method above,
    //   which uses the Zimmert & Seldin (2021) formula.
    if (alpha == 0.5) {
        return tsallisINF(loss, eta / alpha, probs, x);
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
        WarningInFunction << "Newton solver did not converge:"  << update << endl;
    }
    return x;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
