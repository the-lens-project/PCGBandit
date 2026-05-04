const scalar mu = controlDict_.getOrDefault<scalar>("mu", 1.0);
const scalar power = controlDict_.getOrDefault<scalar>("power", 1.0);
const bool pathSimilarity = Switch(controlDict_.getOrDefault<word>("pathSimilarity", "yes"));
const bool useTsallis = Switch(controlDict_.getOrDefault<word>("useTsallis", "no"));
const bool useBonus = Switch(controlDict_.getOrDefault<word>("useBonus", "no"));
const scalar gammaMax = controlDict_.getOrDefault<scalar>("gammaMax", 0.0);

scalarField probs(d);
scalar t = learningDict.getOrAdd("t", 0.0);

if (t == 0.0) {

    t += 1.0;
    learningDict.set<scalar>("t", t);
    SquareMatrix<scalar> S;
    if (pathSimilarity) {
        S = pathMatrix(preconditionerDictsMap[banditName_]);
    } else {
        S = similarityMatrix(preconditionerDictsMap[banditName_]);
        elementwisePower(S, power);
    }
    nonSerializableObjects_.set(banditName_, new DecomposedLaplacian(S));
    const DecomposedLaplacian& Lap = *nonSerializableObjects_[banditName_];
    scalar dEff = Lap.dEff(mu);
    learningDict.set<scalar>("dEff", dEff);
    if (gammaMax == 0.0) {
        probs = 1.0 / d;
    } else {
        probs = Lap.DOptimalDesign(mu);
    }
    learningDict.set<scalarField>("Pi", probs);
    learningDict.set<scalarField>("F", scalarField(d, 0.0));
    if (useBonus) {
        learningDict.set<scalarField>("B", scalarField(d, 0.0));
    }

    #ifdef PCGB_DEBUG
    Info<< "dEff: " << dEff << endl;
    #endif

} else {

    scalar scale = learningDict.getOrAdd<scalar>("scale", 1.0);
    scalar loss = learningDict.getOrAdd<scalar>("loss", 0.0);
    scalarField F = learningDict.get<scalarField>("F");
    
    if (loss > 0.0) {

        const DecomposedLaplacian& Lap = *nonSerializableObjects_[banditName_];
        if (useBonus || useTsallis) {
            Pair<scalarField> hb = Lap.getHatAndBonus(learningDict.get<scalarField>("probs"), mu, learningDict.get<label>("i"));
            F += hb[0] * loss;
            if (useBonus) {
                learningDict.set<scalarField>("B", learningDict.get<scalarField>("B") + sqrt(hb[1] / t));
            }
            if (useTsallis) {
                learningDict.set<scalar>("gPi", max(hb[1]));
            }
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
    scalar eta;
    scalar gamma;
    scalarField lhat = F / scale;
    if (useBonus) {
        lhat -= learningDict.get<scalarField>("B");
    }
    if (useTsallis) {
        scalar alpha = 1.0 - 0.5 / max(1.0, log(scalar(d) / dEff));
        scalar oma = 1.0 - alpha;
        eta = 0.5 * sqrt(oma * pow(scalar(d), oma) / (2.0 * alpha * (t + 1.0) * pow(dEff, alpha)));
        scalar gPi = learningDict.get<scalar>("gPi");
        gamma = gPi * pow(8.0 * eta / oma * (1.0 + 1.0 / sqrt(gPi * t)), 1.0 / alpha);
        scalar x = learningDict.getOrAdd<scalar>("x", -1.0);
        learningDict.set<scalar>("x", tsallisINF(lhat, eta, probs, x, alpha));
    } else {
        eta = sqrt(log(scalar(d)) / (2.0 * dEff * t));
        gamma = eta * dEff * (1.0 + 1.0 / sqrt(t));
        probs = softmax(-eta * lhat);
    }
    gamma = min(gammaMax, gamma);
    probs *= (1.0 - gamma);
    probs += gamma * learningDict.get<scalarField>("Pi");

    #ifdef PCGB_DEBUG
    Info<< "t: " << t << "; eta: " << eta << "; gamma: " << gamma << endl;
    Info<< "F: " << F << endl;
    if (useBonus) {
        Info<< "B: " << learningDict.get<scalarField>("B") << endl;
    }
    Info<< "loss: " << loss << endl;
    Info<< "scale: " << scale << endl;
    Info<< "entropy: " << -sum(probs * log(probs + SMALL)) << "; max: " << max(probs) << endl;
    #endif

}

#ifdef PCGB_DEBUG
Info<< "SpeKL Probabilities: " << probs << endl;
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
