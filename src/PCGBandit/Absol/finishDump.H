if (mesh.time().writeTime()) {

    const label numDumps = numDumpsDict.get<label>(banditName_);

    std::stringstream subdir;
    subdir << fieldName_ << "-Absol/" << numDumps << "/";

    // Global cell numbering, so the per-processor pieces can be assembled into
    // a single monolithic global matrix in Absol.py. Collective: every rank 
    // constructs it in lock-step inside this writeTime block.
    const globalIndex globalNumbering(nCells);
    labelList globalCells(nCells);
    forAll(globalCells, c)
    {
        globalCells[c] = globalNumbering.toGlobal(c);
    }

    IOdictionary writeInfo
    (
        IOobject
        (
            subdir.str() + "info",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
        )
    );
    writeInfo.set<scalar>("relativeTolerance", relTol_);
    writeInfo.set<scalar>("tolerance", tolerance_);
    writeInfo.set<scalar>("initialResidual", solverPerf.initialResidual());
    writeInfo.set<scalar>("finalResidual", solverPerf.finalResidual());
    writeInfo.set<scalar>("nIterations", solverPerf.nIterations());
    writeInfo.set<scalar>("solverTime", -solverTime);
    writeInfo.Foam::regIOobject::write();

    dump<solveScalarField, IOField<solveScalar>>
    (
        psi,
        subdir.str() + "sol",
        mesh
    );

    // Internal-face structure (GLOBAL indices) and coefficients, deduplicated 
    // against the most recent dump. The addressing is constant for a static 
    // mesh, so it is dumped once and carried forward; the coefficients are re-
    // dumped whenever they change (e.g. each timestep, but skipped across 
    // correctors that reuse the same matrix).
    check_and_dump<List<label>, IOList<label>>
    (
        banditName_,
        lowerAddrDict,
        globalNumbering.toGlobal(matrix_.lduAddr().lowerAddr()),
        subdir.str() + "lowerAddr",
        mesh
    );
    check_and_dump<List<label>, IOList<label>>
    (
        banditName_,
        upperAddrDict,
        globalNumbering.toGlobal(matrix_.lduAddr().upperAddr()),
        subdir.str() + "upperAddr",
        mesh
    );

    check_and_dump<solveScalarField, IOField<solveScalar>>
    (
        banditName_,
        lowerDict,
        matrix_.lower(),
        subdir.str() + "lower",
        mesh
    );

    // Interface couplings as off-diagonal entries: A[row][col] = -bouCoeff,
    // with row/col the global indices of the local/neighbour cells. Neighbour 
    // global cell numbers are swapped across each interface (canonical 
    // lduPrimitiveMesh pattern). The symmetric partner A[col][row] is 
    // contributed by the neighbour subdomain.
    PtrList<labelList> nbrGlobalCells(interfaces_.size());
    {
        const label startOfRequests = UPstream::nRequests();
        forAll(interfaces_, inti)
        {
            if (interfaces_.set(inti))
            {
                interfaces_[inti].interface().initInternalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    globalCells
                );
            }
        }
        UPstream::waitRequests(startOfRequests);
        forAll(interfaces_, inti)
        {
            if (interfaces_.set(inti))
            {
                nbrGlobalCells.set
                (
                    inti,
                    new labelList
                    (
                        interfaces_[inti].interface().internalFieldTransfer
                        (
                            Pstream::commsTypes::nonBlocking,
                            globalCells
                        )
                    )
                );
            }
        }
    }

    DynamicList<label> interfaceRow;
    DynamicList<label> interfaceCol;
    DynamicList<solveScalar> interfaceLower;
    forAll(interfaces_, inti)
    {
        if (interfaces_.set(inti))
        {
            const labelUList& fc = interfaces_[inti].interface().faceCells();
            const scalarField& bc = interfaceBouCoeffs_[inti];
            const labelList& nbr = nbrGlobalCells[inti];
            forAll(fc, f)
            {
                interfaceRow.append(globalCells[fc[f]]);
                interfaceCol.append(nbr[f]);
                interfaceLower.append(-bc[f]);
            }
        }
    }

    // Deduplicated like the internal structure: interfaceRow/interfaceCol 
    // (addressing) are constant and dumped once; interfaceLower (coefficients)
    // track the matrix. check_and_dump's collective decision keeps the 
    // collated write unanimous.
    check_and_dump<List<label>, IOList<label>>
    (
        banditName_,
        interfaceRowDict,
        labelList(interfaceRow),
        subdir.str() + "interfaceRow",
        mesh
    );
    check_and_dump<List<label>, IOList<label>>
    (
        banditName_,
        interfaceColDict,
        labelList(interfaceCol),
        subdir.str() + "interfaceCol",
        mesh
    );
    check_and_dump<solveScalarField, IOField<solveScalar>>
    (
        banditName_,
        interfaceLowerDict,
        solveScalarField(interfaceLower),
        subdir.str() + "interfaceLower",
        mesh
    );

    // Right-hand side: with the interface coupling represented in the matrix, 
    // the RHS is simply the (local) source vector.
    dump<solveScalarField, IOField<solveScalar>>
    (
        source,
        subdir.str() + "b",
        mesh
    );

}
