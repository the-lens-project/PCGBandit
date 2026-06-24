const auto* meshPtr = isA<fvMesh>(matrix().mesh());
const auto& mesh = *meshPtr;

if (mesh.time().writeTime()) {

    label numDumps = numDumpsDict.getOrAdd<label>(banditName_, 0) + 1;
    numDumpsDict.set<label>(banditName_, numDumps);

    std::stringstream subdir;
    subdir << fieldName_ << "-Absol/" << numDumps << "/";

    dump<solveScalarField, IOField<solveScalar>>
    (
        psi,
        subdir.str() + "init",
        mesh
    );

    check_and_dump<solveScalarField, IOField<solveScalar>>
    (
        banditName_,
        diagDict,
        matrix_.diag(),
        subdir.str() + "diag",
        mesh
    );

}
