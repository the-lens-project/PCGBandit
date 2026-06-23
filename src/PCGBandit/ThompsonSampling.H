scalarField mu(d);
scalarField sigma(d);
scalar t = learningDict.getOrAdd("t", 0.0);

if (t == 0.0) {

    t += 1.0;
    learningDict.set<scalar>("t", t);
    scalarField L(d, 0.0);
    learningDict.set<scalarField>("L", L);
    scalarField N(d, 0.0);
    learningDict.set<scalarField>("N", N);
    mu = 0.0;
    sigma = 1.0;

} else {

    scalarField L = learningDict.get<scalarField>("L");
    scalarField N = learningDict.get<scalarField>("N");
    scalar loss = learningDict.get<scalar>("loss");
    scalar variance = learningDict.getOrAdd<scalar>("variance", 1.0);
    if (loss > 0.0) {
        scalar scale = learningDict.getOrAdd<scalar>("scale", 0.0);
        variance = (variance * (t - 1.0) + sqr(loss - scale)) / t;
        scale = (scale * (t - 1.0) + loss) / t;
        learningDict.set<scalar>("scale", scale);
        learningDict.set<scalar>("variance", variance);
        t += 1.0;
        learningDict.set<scalar>("t", t);
        i = learningDict.get<label>("i");
        L[i] += loss;
        N[i] += 1.0;
        learningDict.set<scalarField>("L", L);
        learningDict.set<scalarField>("N", N);
    }
    mu = L / (N + 1.0);
    sigma = sqrt(variance / (N + 1.0));

}

#ifdef PCGB_DEBUG
Info<< "ThompsonSampling mu: " << mu << endl;
Info<< "ThompsonSampling sigma: " << sigma << endl;
#endif

i = 0;
scalar thetaMin = mu[0] + sigma[0] * rndGen.GaussNormal<scalar>();
for (label j = 1; j < d; j++) {
    scalar theta = mu[j] + sigma[j] * rndGen.GaussNormal<scalar>();
    if (theta < thetaMin) {
        thetaMin = theta;
        i = j;
    }
}
learningDict.set<label>("i", i);
