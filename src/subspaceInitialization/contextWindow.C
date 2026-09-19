/*---------------------------------------------------------------------------*\
                    Class contextWindow Implementation
\*---------------------------------------------------------------------------*/

#include "contextWindow.H"
#include "Random.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(contextWindow, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::contextWindow::contextWindow
(
    const word& objName,
    const fvMesh& mesh,
    const label maxLenHistory,
    const label nCells,
    const label seed,
    const bool persist
)
:
    subspaceState<contextWindow>(objName, mesh, nCells, persist),
    maxLenHistory_(max(maxLenHistory, label(1))),
    solutions_(maxLenHistory_),
    newestIdx_(0),
    Z_(maxLenHistory_, maxLenHistory_, Zero)
{
    forAll(solutions_, i)
    {
        solutions_[i].setSize(nCells_, Zero);
    }

    Random probeGenerator(seed);
    for (label k = 0; k < maxLenHistory_; ++k)
    {
        for (label j = 0; j < maxLenHistory_; ++j)
        {
            Z_(k, j) = probeGenerator.GaussNormal<scalar>();
        }
    }

    // Restore only after all members are initialized.
    resume();
}

// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

Foam::word Foam::contextWindow::makeKey
(
    const word& fieldName,
    const direction cmpt
)
{
    return word
    (
        "contextWindow:" + fieldName + ":" + Foam::name(label(cmpt)),
        false
    );
}


Foam::contextWindow& Foam::contextWindow::get
(
    const fvMesh& mesh,
    const word& objName,
    const label maxLenHistory,
    const label nCells,
    const label seed,
    const bool persist
)
{
    // The first caller fixes the window size and seed.
    contextWindow& window = const_cast<contextWindow&>
    (
        MeshObject<fvMesh, TopologicalMeshObject, contextWindow>::New
        (
            objName, mesh, maxLenHistory, nCells, seed, persist
        )
    );

    // Warn once if a later caller requests more columns than the window holds.
    if (!window.ceilingWarned_ && maxLenHistory > window.maxLenHistory_)
    {
        window.ceilingWarned_ = true;

        WarningInFunction
            << objName << " was built with lenHistory=" << window.maxLenHistory_
            << " by the first solver dictionary to use it; a later one asks"
            << " for " << maxLenHistory
            << " and will be clamped.  Give every solver dictionary solving"
            << " this field the same lenHistory." << endl;
    }

    return window;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::contextWindow::clearState()
{
    newestIdx_ = 0;
}


bool Foam::contextWindow::update
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

    // Move backward so increasing column indices run from newest to oldest.
    // Once full, the next slot contains the oldest column.
    if (nPushes_ > 0)
    {
        newestIdx_ = (newestIdx_ + maxLenHistory_ - 1) % maxLenHistory_;
    }

    solutions_[newestIdx_] = psi;
    advance(timeIndex);
    return true;
}

// * * * * * * * * * * * * * * * * *  I/O  * * * * * * * * * * * * * * * * * //

bool Foam::contextWindow::writeData(Ostream& os) const
{
    const label nColumns = depth();

    os  << nCells_ << token::SPACE
        << maxLenHistory_ << token::SPACE
        << nColumns << token::SPACE
        << newestTimeIndex_ << nl;

    for (label i = 0; i < nColumns; ++i)
    {
        os  << column(i) << nl;
    }

    return os.good();
}


bool Foam::contextWindow::readData(Istream& is)
{
    label nCells, maxLenHistory, nColumns, timeIndex;
    is  >> nCells >> maxLenHistory >> nColumns >> timeIndex;

    if (nCells != nCells_ || maxLenHistory != maxLenHistory_)
    {
        WarningInFunction
            << name() << " was written with nCells=" << nCells
            << " lenHistory=" << maxLenHistory << ", this run has "
            << nCells_ << "/" << maxLenHistory_
            << ".  Starting from an empty window." << endl;
        clear();
        return false;
    }

    // Validate temporary columns before replacing the fixed-size ring buffers.
    nColumns = min(max(nColumns, label(0)), maxLenHistory_);

    List<solveScalarField> columns(nColumns);
    for (label i = 0; i < nColumns; ++i)
    {
        is >> columns[i];

        if (!is.good() || columns[i].size() != nCells_)
        {
            clear();
            return false;
        }
    }

    newestIdx_ = 0;
    for (label i = 0; i < nColumns; ++i)
    {
        solutions_[ringPos(i)].transfer(columns[i]);
    }

    nPushes_ = nColumns;
    newestTimeIndex_ = timeIndex;

    return true;
}

// ************************************************************************* //
