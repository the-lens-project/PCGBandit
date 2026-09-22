const scalar mu = controlDict_.getOrDefault<scalar>("mu", 0.001);
const scalar eta = controlDict_.getOrDefault<scalar>("eta", 0.1);
const scalar gamma = controlDict_.getOrDefault<scalar>("gamma", 0.0);
const bool useBonus = Switch(controlDict_.getOrDefault<word>("useBonus", "no"));

scalarField probs(d);
scalar t = learningDict.getOrAdd("t", 0.0);

if (t == 0.0) {

    t += 1.0;
    learningDict.set<scalar>("t", t);
    nonSerializableObjects_.set(
        banditName_,
        new decomposedLaplacian(
            similarityMatrix(configurationSpaces[banditName_].dicts(), controlDict_)
        )
    );
    const decomposedLaplacian& Lap = *nonSerializableObjects_[banditName_];
    scalar dEff = Lap.dEff(mu);
    learningDict.set<scalar>("dEff", dEff);
    if (gamma == 0.0) {
        probs = 1.0 / d;
    } else {
        probs = Lap.DOptimalDesign(mu);
    }
    learningDict.set<scalarField>("Pi", probs);
    learningDict.set<scalarField>("F", scalarField(d, 0.0));
    if (useBonus) {
        learningDict.set<scalarField>("B", scalarField(d, 0.0));
    }

} else {

    scalar scale = learningDict.getOrAdd<scalar>("scale", 1.0);
    scalar loss = learningDict.getOrAdd<scalar>("loss", 0.0);
    scalarField F = learningDict.get<scalarField>("F");
    
    if (loss > 0.0) {

        const decomposedLaplacian& Lap = *nonSerializableObjects_[banditName_];
        if (useBonus) {
            Pair<scalarField> hb = Lap.getHatAndBonus(learningDict.get<scalarField>("probs"), mu, learningDict.get<label>("i"));
            F += hb[0] * loss;
            learningDict.set<scalarField>("B", learningDict.get<scalarField>("B") + sqrt(hb[1] / t));
        } else {
            F += loss * Lap.getHat(learningDict.get<scalarField>("probs"), mu, learningDict.get<label>("i"));
        }
            
        scale = (scale * (t - 1.0) + loss) / t;
        t += 1.0;
        learningDict.set<scalarField>("F", F);
        learningDict.set<scalar>("scale", scale);
        learningDict.set<scalar>("t", t);
        
    }

    // --- Update probabilities
    scalar dEff = learningDict.get<scalar>("dEff");
    scalarField lhat = F / scale;
    if (useBonus) {
        lhat -= learningDict.get<scalarField>("B");
    }
    scalar alpha = 1.0 - 0.5 / max(1.0, log(scalar(d) / dEff));
    scalar x = learningDict.getOrAdd<scalar>("x", -1.0);
    learningDict.set<scalar>("x", tsallisINF(lhat, eta, probs, x, alpha));
    if (gamma > 0.0) {
      probs *= (1.0 - gamma);
      probs += gamma * learningDict.get<scalarField>("Pi");
    }

}

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
learningDict.set<scalarField>("probs", probs);
