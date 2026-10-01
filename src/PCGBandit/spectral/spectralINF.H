const scalar mu = controlDict_.getOrDefault<scalar>("mu", 0.001);
const scalar eta = controlDict_.getOrDefault<scalar>("eta", 0.1);
// Defaults to setting alpha automatically using Lap.dTr(mu)
const scalar alpha = controlDict_.getOrDefault<scalar>("alpha", 0.0);
const scalar gamma = controlDict_.getOrDefault<scalar>("gamma", 0.0);
const bool useBonus = Switch(controlDict_.getOrDefault<word>("useBonus", "no"));

scalar t = learningDict.getOrAdd<scalar>("t", 0.0);

if (t == 0.0) {

    const decomposedLaplacian& Lap = learningDict.emplace<decomposedLaplacian>
    (
        "laplacian",
        similarityMatrix(configurationSpaces[banditName_].dicts(), controlDict_)
    );
    scalarField& probs = learningDict.emplace<scalarField>("probs", d);
    t += 1.0;
    learningDict.set<scalar>("t", t);
    if (alpha == 0.0) {
        learningDict.emplace<scalar>("alpha", 1.0 - 0.5 / max(1.0, log(scalar(d) / Lap.dTr(mu))));
    } else {
        learningDict.emplace<scalar>("alpha", alpha);
    }
    if (gamma > 0.0) {
        probs = Lap.DOptimalDesign(mu);
        learningDict.emplace<scalarField>("pi", probs);
    } else {
        probs = 1.0 / d;
    }
    learningDict.emplace<scalarField>("F", d, 0.0);
    if (useBonus) {
        learningDict.emplace<scalarField>("B", d, 0.0);
    }

} else {

    const decomposedLaplacian& Lap = learningDict.get<decomposedLaplacian>("laplacian");
    scalarField& probs = learningDict.get<scalarField>("probs");
    scalar scale = learningDict.getOrAdd<scalar>("scale", 1.0);
    scalar loss = learningDict.getOrAdd<scalar>("loss", 0.0);
    scalarField& F = learningDict.get<scalarField>("F");
    
    // --- Process feedback
    if (loss > 0.0) {

        if (useBonus) {
            Pair<scalarField> hb = Lap.getHatAndBonus(probs, mu, learningDict.get<label>("i"));
            F += hb[0] * loss;
            learningDict.get<scalarField>("B") += sqrt(hb[1] / t);
        } else {
            F += loss * Lap.getHat(probs, mu, learningDict.get<label>("i"));
        }
            
        scale = (scale * (t - 1.0) + loss) / t;
        t += 1.0;
        learningDict.set<scalar>("scale", scale);
        learningDict.set<scalar>("t", t);
        
    }

    // --- Update probabilities
    scalarField lhat = F / scale;
    if (useBonus) {
        lhat -= learningDict.get<scalarField>("B");
    }
    scalar x = learningDict.getOrAdd<scalar>("x", -1.0);
    learningDict.set<scalar>("x", tsallisINF(lhat, eta, probs, x, learningDict.get<scalar>("alpha")));
    if (gamma > 0.0) {
        probs *= (1.0 - gamma);
        probs += gamma * learningDict.get<scalarField>("pi");
    }

}

const scalarField& probs = learningDict.get<scalarField>("probs");

#ifdef PCGB_DEBUG
Info<< "SpectralINF Probabilities: " << probs << endl;
#endif

// --- Naive random sampling
scalar r = rndGen.sample01<scalar>() * sum(probs);
scalar cumulative = 0.0;
for (i = 0; i < d; i++) {
    cumulative += probs[i];
    if (r <= cumulative) {
        break;
    }
}
learningDict.set<label>("i", i);
