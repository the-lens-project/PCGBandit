/*---------------------------------------------------------------------------*\
                    Class EWMASketch Implementation
\*---------------------------------------------------------------------------*/

#include "EWMASketch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(EWMASketch, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EWMASketch::EWMASketch
(
    const word& objName,
    const fvMesh& mesh,
    const label numProbes,
    const scalar decayRate,
    const label nCells,
    const label seed,
    const bool persist
)
:
    subspaceState<EWMASketch>(objName, mesh, nCells, persist),
    numProbes_(max(numProbes, label(1))),
    decayRate_(decayRate),
    seed_(seed),
    Y_(numProbes_),
    w_(numProbes_, Zero),
    probeGenerator_(seed)
{
    forAll(Y_, j)
    {
        Y_[j].setSize(nCells_, Zero);
    }

    // Restore only after all members and the probe generator are initialized.
    resume();
}


// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

Foam::word Foam::EWMASketch::makeKey
(
    const word& fieldName,
    const direction cmpt,
    const scalar decayRate
)
{
    return word
    (
        "EWMASketch:" + fieldName
      + ":" + Foam::name(label(cmpt))
      + ":" + Foam::name(decayRate),
        false
    );
}


Foam::EWMASketch& Foam::EWMASketch::get
(
    const fvMesh& mesh,
    const word& objName,
    const label numProbes,
    const scalar decayRate,
    const label nCells,
    const label seed,
    const bool persist
)
{
    // The first caller fixes the sketch size and seed.
    EWMASketch& sketch = const_cast<EWMASketch&>
    (
        MeshObject<fvMesh, TopologicalMeshObject, EWMASketch>::New
        (
            objName, mesh, numProbes, decayRate, nCells, seed, persist
        )
    );

    // Warn once if a later caller requests more probes than the sketch holds.
    if (!sketch.ceilingWarned_ && numProbes > sketch.numProbes_)
    {
        sketch.ceilingWarned_ = true;

        WarningInFunction
            << objName << " was built with numProbes=" << sketch.numProbes_
            << " by the first solver dictionary to use it; a later one asks"
            << " for " << numProbes << " and will be clamped.  Give every"
            << " dictionary solving this field the same numProbes." << endl;
    }

    return sketch;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::EWMASketch::clearState()
{
    forAll(Y_, j)
    {
        Y_[j] = Zero;
        w_[j] = 0;
    }

    // Restart replay assumes draw zero after clear. Random::reset would retain
    // the cached Marsaglia sample, so reconstruct the generator.
    probeGenerator_ = Random(seed_);
}


bool Foam::EWMASketch::update
(
    const solveScalarField& psi,
    const label timeIndex
)
{
    if (psi.size() != nCells_)
    {
        clear();
        return false;
    }

    const solveScalar* __restrict__ psiValues = psi.begin();

    for (label j = 0; j < numProbes_; ++j)
    {
        // Every rank advances the generator in the same order.
        const solveScalar probeWeight = probeGenerator_.GaussNormal<scalar>();

        solveScalar* __restrict__ columnValues = Y_[j].begin();
        for (label celli = 0; celli < nCells_; ++celli)
        {
            columnValues[celli] =
                decayRate_*columnValues[celli] + probeWeight*psiValues[celli];
        }
        w_[j] = decayRate_*w_[j] + probeWeight;
    }

    advance(timeIndex);
    return true;
}

// * * * * * * * * * * * * * * * * *  I/O  * * * * * * * * * * * * * * * * * //

bool Foam::EWMASketch::writeData(Ostream& os) const
{
    os  << nCells_ << token::SPACE
        << numProbes_ << token::SPACE
        << nPushes_ << token::SPACE
        << newestTimeIndex_ << nl;

    os  << w_ << nl;

    forAll(Y_, j)
    {
        os  << Y_[j] << nl;
    }

    return os.good();
}


bool Foam::EWMASketch::readData(Istream& is)
{
    label nCells, numProbes, nUpdates, timeIndex;
    is  >> nCells >> numProbes >> nUpdates >> timeIndex;

    if (nCells != nCells_ || numProbes != numProbes_)
    {
        WarningInFunction
            << name() << " was written with nCells=" << nCells
            << " numProbes=" << numProbes << ", this run has "
            << nCells_ << "/" << numProbes_
            << ".  Starting from an empty sketch." << endl;
        clear();
        return false;
    }

    // Validate temporary buffers before replacing the fixed-size state.
    List<solveScalar> weights;
    is >> weights;

    List<solveScalarField> columns(numProbes_);
    forAll(columns, j)
    {
        is >> columns[j];

        if
        (
            !is.good()
         || columns[j].size() != nCells_
         || weights.size() != numProbes_
        )
        {
            clear();
            return false;
        }
    }

    w_ = weights;
    forAll(Y_, j)
    {
        Y_[j].transfer(columns[j]);
    }

    nPushes_ = max(nUpdates, label(0));
    newestTimeIndex_ = timeIndex;

    // Replay from draw zero, including the cached Marsaglia sample.
    const label nDraws = nPushes_*numProbes_;
    for (label k = 0; k < nDraws; ++k)
    {
        probeGenerator_.GaussNormal<scalar>();
    }

    return true;
}

// ************************************************************************* //
