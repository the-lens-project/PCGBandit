/*---------------------------------------------------------------------------*\
Class
    Foam::subspaceState

Description
    Mesh-registered solver state with solve-sequence tracking and optional
    restart persistence. The template gives each derived type its own
    MeshObject registry entry.

    Derived classes provide clearState(), writeData(), and readData(), and
    call resume() after initializing all members that readData() accesses.

SourceFiles
    subspaceState.H
\*---------------------------------------------------------------------------*/

#ifndef subspaceState_H
#define subspaceState_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "Time.H"
#include "primitiveFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class subspaceStateBase Declaration
\*---------------------------------------------------------------------------*/

//- Common interface for checking rank agreement in one collective.
class subspaceStateBase
{
protected:

    const label nCells_;

    //- Pushes since the last clear; a context window's depth may be smaller.
    label nPushes_;

    //- timeIndex of the newest push, -1 when empty
    label newestTimeIndex_;

    //- Rank agreement is checked once; later pushes and clears are collective.
    //  Mesh topology changes destroy the object.
    bool agreed_;

    explicit subspaceStateBase(const label nCells)
    :
        nCells_(nCells),
        nPushes_(0),
        newestTimeIndex_(-1),
        agreed_(false)
    {}

    //- Record a completed update().
    void advance(const label timeIndex) noexcept
    {
        ++nPushes_;
        newestTimeIndex_ = timeIndex;
    }

public:

    virtual ~subspaceStateBase() = default;

    label nCells() const noexcept { return nCells_; }

    //- Newest timeIndex; must agree across ranks for collective restart checks.
    label stamp() const noexcept { return newestTimeIndex_; }

    //- Available samples; must agree across ranks to match Amuls and reductions.
    virtual label depth() const = 0;

    virtual void clear() = 0;

    bool agreed() const noexcept { return agreed_; }
    void setAgreed() noexcept { agreed_ = true; }
};


/*---------------------------------------------------------------------------*\
                       Class subspaceState Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class subspaceState
:
    public MeshObject<fvMesh, TopologicalMeshObject, Type>,
    public subspaceStateBase
{
protected:

    typedef MeshObject<fvMesh, TopologicalMeshObject, Type> MeshObjectType;

    //- Read and write state in time directories.
    const bool persist_;

    //- Warn once if a later solver requests more capacity than was allocated.
    bool ceilingWarned_;


    subspaceState
    (
        const word& objName,
        const fvMesh& mesh,
        const label nCells,
        const bool persist
    )
    :
        MeshObjectType(objName, mesh),
        subspaceStateBase(nCells),
        persist_(persist),
        ceilingWarned_(false)
    {}


    //- Enable persistence and read saved state; return the number of pushes.
    //  Call after initializing the derived members used by readData().
    label resume();


public:

    virtual ~subspaceState() = default;

    //- Clear stored samples and sequence counters.
    virtual void clear()
    {
        static_cast<Type*>(this)->clearState();
        nPushes_ = 0;
        newestTimeIndex_ = -1;
    }

    //- Clear state on backward time steps or gaps. Repeated and consecutive
    //  time indices are valid. Call before reading samples at each solve.
    void checkRestart(const label timeIndex)
    {
        if
        (
            nPushes_ > 0
         && (timeIndex < newestTimeIndex_ || timeIndex > newestTimeIndex_ + 1)
        )
        {
            clear();
        }
    }

    //- Write binary state to avoid rounding at writePrecision on restart.
    //  Preserve the configured compression.
    virtual bool writeObject
    (
        IOstreamOption streamOpt,
        const bool writeOnProc
    ) const
    {
        streamOpt.format(IOstreamOption::BINARY);
        return regIOobject::writeObject(streamOpt, writeOnProc);
    }
};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::label Foam::subspaceState<Type>::resume()
{
    if (!persist_)
    {
        return 0;
    }

    const Time& runTime = this->mesh().time();

    // Use a time instance so writeObject advances it at each write time.
    // The mesh's initial constant/ instance would otherwise be overwritten.
    this->instance() = runTime.timeName();
    this->writeOpt(IOobject::AUTO_WRITE);

    // Find the latest saved state, including the preceding restart time.
    // Exclude constant/, which does not identify a solve time.
    const word instanceName = runTime.findInstance
    (
        fileName::null,
        this->name(),
        IOobject::READ_IF_PRESENT,
        word::null,
        false                       // no constant fallback
    );

    if (instanceName.empty())
    {
        return 0;
    }

    this->instance() = instanceName;
    const IOobject::readOption previousReadOption =
        this->readOpt(IOobject::MUST_READ);

    bool readOk = false;
    if (this->headerOk())
    {
        readOk = this->readData(this->readStream(Type::typeName));
        this->close();
    }

    this->readOpt(previousReadOption);
    this->instance() = runTime.timeName();

    if (!readOk)
    {
        // Discard any partially read state.
        clear();
        return 0;
    }

    // The first solve checks the time index and agrees depth/stamp across ranks.
    return nPushes_;
}

} // End namespace Foam

#endif

// ************************************************************************* //
