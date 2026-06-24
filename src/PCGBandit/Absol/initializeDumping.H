dictionary numDumpsDict;

// Most-recent-dump reference values, keyed by banditName, for deduplication.
HashTable<solveScalarField> diagDict;
HashTable<solveScalarField> lowerDict;
HashTable<solveScalarField> interfaceLowerDict;
HashTable<List<label>> lowerAddrDict;
HashTable<List<label>> upperAddrDict;
HashTable<List<label>> interfaceRowDict;
HashTable<List<label>> interfaceColDict;


template <typename DataType, typename IOType>
void dump(
    const DataType& current,
    const word fname,
    const fvMesh& mesh) {

    IOType writeCurrent(
        IOobject(
            fname,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ
            ),
        current
        );
    writeCurrent.write();

}


template <typename DataType, typename IOType>
void check_and_dump(
    const word key,
    HashTable<DataType>& store,
    const DataType& current,
    const word fname,
    const fvMesh& mesh) {

    // Deduplicate against most recent dump. The stored reference is updated on
    // every write, so the comparison is always against the latest dump.
    bool different = true;
    if (store.found(key)) {
        const DataType& previous = store[key];
        if (previous.size() == current.size()) {
            different = false;
            for (label i = 0; i < current.size(); i++) {
                if (previous[i] != current[i]) {
                    different = true;
                    break;
                }
            }
        }
    }

    // Collective: write on every rank iff the data changed on ANY rank, so the
    // decision is unanimous and a collated (single-file) write stays consistent.
    if (returnReduce(label(different), maxOp<label>()) > 0) {
        store.set(key, current);
        dump<DataType, IOType>(current, fname, mesh);
    }

}
