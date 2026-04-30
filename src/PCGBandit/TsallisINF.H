label maxNewtonIter = 100;
scalar tolerance = 1e-8;

scalarField probs(d);
scalar t = learningDict.getOrAdd("t", 0.0);

if (t == 0.0) {

    t += 1.0;
    learningDict.set<scalar>("t", t);
    
    // use numDroptols instead of ictc loop
    label numICTC = controlDict_.get<label>("numDroptols");
    
    // p = 0 gives uniform
    label priorCenter =
        controlDict_.getOrDefault<label>("priorCenter", numICTC/2);
    
    // priorPower default is uniform (0)
    scalar priorPower =
        controlDict_.getOrDefault<scalar>("priorPower", 0.0);
    
    // validate priorCenter is in valid range
    if (priorCenter < 0 || priorCenter >= numICTC)
    {
        WarningInFunction
            << "priorCenter = " << priorCenter
            << " is outside range [0, " << numICTC-1 << "]"
            << nl << "    using nearest valid arm instead" << endl;
        
        priorCenter = max(label(0), min(priorCenter, numICTC-1));
    }
    
    // using: W_i = (1 / (1 + |i-k|))^p
    // where k = priorCenter, i = arm index
    
    // compute weights for ICTC arms
    for (label k = 0; k < numICTC; k++)
    {
        scalar dist = mag(scalar(k) - scalar(priorCenter));
        scalar base = 1.0 / (1.0 + dist);
        probs[k] = pow(base, priorPower);
    }
    
    // DIC and GAMG get same weight as center to prevent downweighting
    for (label k = numICTC; k < d; k++)
    {
        probs[k] = 1.0;
    }
    
    // normalize so weights sum to 1
    probs /= sum(probs);
    
    Info<< "using distance prior: center=" << priorCenter
        << ", power=" << priorPower << endl;
    Info<< "initial probs: " << probs << endl;
    
    scalarField IW(d, 0.0);
    learningDict.set<scalarField>("IW", IW);
    if (lossEstimator_ == "RV") {
        scalarField shift(d, 0.0);
        learningDict.set<scalarField>("shift", shift);
    }

} else {

    scalarField IW = learningDict.get<scalarField>("IW");
    scalar scale = learningDict.getOrAdd<scalar>("scale", 1.0);
    scalar loss = learningDict.getOrAdd<scalar>("loss", 0.0);
    if (loss > 0.0) {
        IW[learningDict.get<label>("i")] += loss / learningDict.get<scalar>("p");
        learningDict.set<scalarField>("IW", IW);
        scale = (scale * (t - 1.0) + loss) / t;
        learningDict.set<scalar>("scale", scale);
        t += 1.0;
        learningDict.set<scalar>("t", t);
    }

    scalarField lhat = IW / scale;
    scalar eta = 2.0 / sqrt(t);
    if (lossEstimator_ == "RV") {
        lhat += learningDict.get<scalarField>("shift");
        eta *= 2.0;
    }

    scalar x = learningDict.getOrAdd<scalar>("x", -1.0);

    scalar update;
    label j;
    for (j = 0; j < maxNewtonIter; j++) {
        probs = 2.0  / (eta * (lhat - x));
        probs *= probs;
        update = (sum(probs) - 1.0) / (eta * sum(pow(probs, 1.5)));
        x -= update;
        if (mag(update) < tolerance) {
            break;
        }
    }

    if (j == maxNewtonIter) {
        Info<< "Newton solver did not converge: " << update << endl;
    }
    learningDict.set<scalar>("x", x);

}

#ifdef PCGB_DEBUG
Info<< "TsallisINF Probabilities: " << probs << endl;
#endif

scalar r = rndGen.sample01<scalar>() * sum(probs);
scalar cumulative = 0.0;
for (i = 0; i < d; i++) {
    cumulative += probs[i];
    if (r <= cumulative) {
        break;
    }
}

learningDict.set<label>("i", i);
learningDict.set<scalar>("p", probs[i]);

if (lossEstimator_ == "RV" && t >= 16.0) {
    scalarField shift = learningDict.get<scalarField>("shift");
    for (label j = 0; j < d; j++) {
        scalar p = probs[j];
        if (t * p >= 16.0) {
            shift[j] += 0.5;
            if (i == j) {
                shift[j] -= 0.5 / p;
            }
        }
    }
    learningDict.set<scalarField>("shift", shift);
}
