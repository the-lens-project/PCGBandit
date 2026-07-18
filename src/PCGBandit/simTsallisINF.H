scalarField probs(d);
scalar t = learningDict.getOrAdd("t", 0.0);

SquareMatrix<scalar> S;
scalarField columnSums;

if (t == 0.0) {

    t += 1.0;
    learningDict.set<scalar>("t", t);
    probs = 1.0 / scalar(d);
    scalarField IW(d, 0.0);
    learningDict.set<scalarField>("IW", IW);
    if (lossEstimator_ == "RV") {
        scalarField shift(d, 0.0);
        learningDict.set<scalarField>("shift", shift);
    }

    // Constructs the matrix
    S = similarityMatrix(preconditionerDictsMap[banditName_]);
    learningDict.set("similarityMatrix", S);

    // Finds the sum of each column of the similarity matrix (to be used for RV estimator)
    if (lossEstimator_ == "RV") {
	    columnSums.setSize(d, 0.0);
    	for (label i = 0; i < d; i++) {
            for (label j = 0; j < d; j++) {
                    columnSums[j] += S(i, j);
            }
        }
        learningDict.set("columnSums", columnSums);
    }

} else {

    S = learningDict.get<SquareMatrix<scalar>>("similarityMatrix");

    // --- Process solver feedback
    scalarField IW = learningDict.get<scalarField>("IW");
    scalar scale = learningDict.getOrAdd<scalar>("scale", 1.0);
    scalar loss = learningDict.getOrAdd<scalar>("loss", 0.0);
    
    if (loss > 0.0) {
	    IW += (loss / learningDict.get<scalar>("p")) * scalarField(S.subColumn(learningDict.get<label>("i")));
        learningDict.set<scalarField>("IW", IW);
        scale = (scale * (t - 1.0) + loss) / t;
        learningDict.set<scalar>("scale", scale);
        t += 1.0;
        learningDict.set<scalar>("t", t);
    }

    // --- Compute loss estimator
    scalarField lhat = IW / scale;
    scalar eta = 2.0 / sqrt(t);
    if (lossEstimator_ == "RV") {
        lhat += learningDict.get<scalarField>("shift");
        eta *= 2.0;
    }

    // --- Tsallis-INF update
    scalar x = learningDict.getOrAdd<scalar>("x", -1.0);
    learningDict.set<scalar>("x", tsallisINF(lhat, eta, probs, x));

}

#ifdef PCGB_DEBUG
Info<< "\tprobabilities: " << probs << endl;
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
learningDict.set<scalar>("p", probs[i]);

// --- Computes shift for the reduced variance case
if (lossEstimator_ == "RV" && t >= 16.0) {
    
    scalarField shift = learningDict.get<scalarField>("shift");
    columnSums = learningDict.get<scalarField>("columnSums");
    scalarField lossFactor = S.subColumn(i);

    for (label j = 0; j < d; j++) {

        scalar p_i = probs[i];
        scalar p_j = probs[j];

        if (t * p_j >= 16.0) {
            shift[j] += 0.5*columnSums[j];
            shift[j] -= (0.5 / p_i)*lossFactor[j];
        }
    }

    learningDict.set<scalarField>("shift", shift);
}
