solveScalarField currentDiag = solveScalarField(matrix_.diag());
solveScalarField currentLower = solveScalarField(matrix_.lower());
solveScalarField currentTarget = solveScalarField(source);
Info<< "\tsq Frob norm of current diag: " << gSum(sqr(currentDiag)) << endl;
Info<< "\tsq Frob norm of current lower: " << gSum(sqr(currentLower)) << endl;
Info<< "\tsq Eucl norm of current target: " << gSum(sqr(currentTarget)) << endl;
if (initialDiags.found(banditName_)) {
    Info<< "\tsq Frob dist from initial diag: " << gSum(sqr(currentDiag - initialDiags.get<solveScalarField>(banditName_))) << endl;
    Info<< "\tsq Frob dist from initial lower: " << gSum(sqr(currentLower - initialLowers.get<solveScalarField>(banditName_))) << endl;
    Info<< "\tsq Eucl norm from initial target: " << gSum(sqr(currentTarget - initialTargets.get<solveScalarField>(banditName_))) << endl;
    Info<< "\tsq Frob dist from previous diag: " << gSum(sqr(currentDiag - previousDiags.get<solveScalarField>(banditName_))) << endl;
    Info<< "\tsq Frob dist from previous lower: " << gSum(sqr(currentLower - previousLowers.get<solveScalarField>(banditName_))) << endl;
    Info<< "\tsq Eucl norm from previous target: " << gSum(sqr(currentTarget - previousTargets.get<solveScalarField>(banditName_))) << endl;
} else {
    initialDiags.add<solveScalarField>(banditName_, currentDiag);
    initialLowers.add<solveScalarField>(banditName_, currentLower);
    initialTargets.add<solveScalarField>(banditName_, currentTarget);
}
previousDiags.set<solveScalarField>(banditName_, currentDiag);
previousLowers.set<solveScalarField>(banditName_, currentLower);
previousTargets.set<solveScalarField>(banditName_, currentTarget);
