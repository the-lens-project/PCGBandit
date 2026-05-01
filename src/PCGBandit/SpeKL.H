scalar mu = 1.0;

scalarField probs(d);
scalar t = learningDict.getOrAdd("t", 0.0);

if (t == 0.0) {

    t += 1.0;
    learningDict.set<scalar>("t", t);
    SquareMatrix<scalar> similarityMatrix = preconditionerSimilarityMatrix(preconditionerDictsMap[banditName_]);
    nonSerializableObjects_.set(banditName_, new DecomposedLaplacian(similarityMatrix));
    const DecomposedLaplacian& Lap = *nonSerializableObjects_[banditName_];
    scalar dEff = Lap.dEff(mu);
    learningDict.set<scalar>("dEff", dEff);
    probs = Lap.DOptimalDesign(mu);
    learningDict.set<scalarField>("Pi", probs);
    learningDict.set<scalarField>("F", scalarField(d, 0.0));
    learningDict.set<scalarField>("B", scalarField(d, 0.0));

} else {

    scalar scale = learningDict.getOrAdd<scalar>("scale", 1.0);
    scalar loss = learningDict.getOrAdd<scalar>("loss", 0.0);
    scalarField F = learningDict.get<scalarField>("F");
    scalarField B = learningDict.get<scalarField>("B");
    
    if (loss > 0.0) {
        const DecomposedLaplacian& Lap = *nonSerializableObjects_[banditName_];        
        Pair<scalarField> hb = Lap.HatAndBonus(learningDict.get<scalarField>("probs"), mu, learningDict.get<label>("i"));
        F += hb[0] * loss;
        B += sqrt(hb[1] / t);
        scale = (scale * (t - 1.0) + loss) / t;
        t += 1.0;
        learningDict.set<scalarField>("F", F);
        learningDict.set<scalarField>("B", B);
        learningDict.set<scalar>("scale", scale);
        learningDict.set<scalar>("t", t);
    }

    // --- Compute softmax and probabilities
    scalar dEff = learningDict.get<scalar>("dEff");
    scalar eta = sqrt(log(scalar(d)) / (2.0 * dEff * t));
    scalar gamma = min(1.0, eta * dEff * (1.0 + 1.0 / sqrt(t)));
    scalarField logits = eta * (B - F / scale);
    scalarField e = exp(logits - max(logits));
    probs = (1.0 - gamma) * (e / sum(e)) + gamma * learningDict.get<scalarField>("Pi");

    Info<< "SpeKL t: " << t << "; eta: " << eta << "; gamma: " << gamma << endl;
    Info<< "F: " << F << endl;
    Info<< "B: " << B << endl;
    Info<< "loss: " << loss << endl;
    Info<< "scale: " << scale << endl;
    Info<< "dEff: " << dEff << endl;
    Info<< "Pi: " << learningDict.get<scalarField>("Pi") << endl;

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
