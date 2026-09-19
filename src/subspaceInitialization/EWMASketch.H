/*---------------------------------------------------------------------------*\
Class
    Foam::EWMASketch

Description
    Exponentially weighted sketch of previous iterates
    (arXiv:2511.18071, Algorithm 4):

        Y <- r*Y + psi*omega^T,   omega ~ N(0, I) drawn fresh every solve
        w <- r*w + omega
        Omega = Y - psi0*w^T = sum_k r^k (x_k - psi0) omega_k^T

    Stores numProbes vectors and costs O(nCells*numProbes) per update.
    Initialization reads the sketch before adding the entry iterate psi0, so
    all stored iterates are older than psi0.

    Each field component and decay rate has a separate sketch. All sketches
    are updated regardless of the selected arm. Restart reads restore Y and w
    and replay the random draws to restore the probe generator.

SourceFiles
    EWMASketch.C
\*---------------------------------------------------------------------------*/

#ifndef EWMASketch_H
#define EWMASketch_H

#include "subspaceState.H"
#include "Random.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class EWMASketch Declaration
\*---------------------------------------------------------------------------*/

class EWMASketch
:
    public subspaceState<EWMASketch>
{
    const label numProbes_;
    const scalar decayRate_;

    //- Seed used to reset the generator and replay persisted updates
    const label seed_;

    //- The sketch: numProbes_ columns of nCells_
    List<solveScalarField> Y_;

    //- Accumulated probe weights, w_j = sum_k r^k omega_k,j
    List<solveScalar> w_;

    //- All ranks must draw identical probes for the distributed sketch
    Random probeGenerator_;

    EWMASketch(const EWMASketch&) = delete;
    void operator=(const EWMASketch&) = delete;

public:

    TypeName("EWMASketch");

    friend class subspaceState<EWMASketch>;

    EWMASketch
    (
        const word& objName,
        const fvMesh& mesh,
        const label numProbes,
        const scalar decayRate,
        const label nCells,
        const label seed,
        const bool persist
    );

    //- Must not dereference mesh(): fvMesh::~fvMesh clears its registered
    //  objects from the destructor body.
    virtual ~EWMASketch() = default;


    //- Key by field, component and rendered decay rate; equal keys share state
    static word makeKey
    (
        const word& fieldName,
        const direction cmpt,
        const scalar decayRate
    );

    //- Retrieve or create the mesh-registered sketch
    static EWMASketch& get
    (
        const fvMesh& mesh,
        const word& objName,
        const label numProbes,
        const scalar decayRate,
        const label nCells,
        const label seed,
        const bool persist
    );


    label numProbes() const noexcept { return numProbes_; }

    //- Updates since the last clear; rank(Omega) <= min(numProbes, depth)
    virtual label depth() const { return nPushes_; }

    const solveScalarField& column(const label j) const { return Y_[j]; }
    solveScalar weight(const label j) const { return w_[j]; }

    //- Update Y and w after reading the sketch, before correcting psi.
    //  Clears the sketch and returns false if psi has the wrong length.
    bool update(const solveScalarField& psi, const label timeIndex);


    // I/O

        //- Write sizes, update count, time index, weights and sketch columns
        virtual bool writeData(Ostream& os) const;

        //- Validate and restore state, then replay the probe draws
        virtual bool readData(Istream& is);

private:

        //- Reset the sketch and probe generator
        void clearState();
};

} // End namespace Foam

#endif

// ************************************************************************* //
