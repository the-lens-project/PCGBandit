/*---------------------------------------------------------------------------*\
Class
    Foam::contextWindow

Description
    Ring buffer of previous iterates used to form the randomized subspace
    sketch (arXiv:2309.02156).

    For stored columns X, initial iterate psi0, and Gaussian Z:

        Omega = (X - psi0*1^T) Z,     Z: lenHistory x numProbes
        W = A*Omega
        G y = t
        psi += Omega*y,              rA -= W*y

    Galerkin uses G = Omega^T W and t = Omega^T r0.
    Least squares uses G = W^T W and t = W^T r0.

    Initialization reads the window before storing the entry iterate, so all
    columns are older than that iterate and numProbes can reach lenHistory.
    Solves already converged at entry do not update the window.

    One column is stored per solve. Each field component has one window shared
    by its solver dictionaries, including p and pFinal. The mesh registry
    retains the window between solves; subspaceState handles persistence.

SourceFiles
    contextWindow.C
\*---------------------------------------------------------------------------*/

#ifndef contextWindow_H
#define contextWindow_H

#include "subspaceState.H"
#include "RectangularMatrix.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class contextWindow Declaration
\*---------------------------------------------------------------------------*/

class contextWindow
:
    public subspaceState<contextWindow>
{
    const label maxLenHistory_;

    List<solveScalarField> solutions_;
    label newestIdx_;

    //- Gaussian test matrix shared by all window lengths and probe counts.
    //  Drawn identically on every rank, independently of the bandit's RNG.
    RectangularMatrix<scalar> Z_;

    label ringPos(const label i) const { return (newestIdx_ + i) % maxLenHistory_; }

    contextWindow(const contextWindow&) = delete;
    void operator=(const contextWindow&) = delete;

public:

    TypeName("contextWindow");

    friend class subspaceState<contextWindow>;

    contextWindow
    (
        const word& objName,
        const fvMesh& mesh,
        const label maxLenHistory,
        const label nCells,
        const label seed,
        const bool persist
    );

    //- Must not dereference mesh(): fvMesh::~fvMesh clears its registered
    //  objects from the destructor body.
    virtual ~contextWindow() = default;


    //- Registry key shared by solver dictionaries for the same field component
    static word makeKey(const word& fieldName, const direction cmpt);

    //- Retrieve or create the mesh-registered window
    static contextWindow& get
    (
        const fvMesh& mesh,
        const word& objName,
        const label maxLenHistory,
        const label nCells,
        const label seed,
        const bool persist
    );


    //- Number of stored columns; bounds the number of probes
    virtual label depth() const { return min(nPushes_, maxLenHistory_); }

    const RectangularMatrix<scalar>& Z() const noexcept { return Z_; }

    //- The i-th newest iterate, valid for 0 <= i < depth()
    const solveScalarField& column(const label i) const
    {
        return solutions_[ringPos(i)];
    }

    //- Store psi after reading the window, before correcting psi.
    //  Clears the window and returns false if psi has the wrong length.
    bool update(const solveScalarField& psi, const label timeIndex);

    // I/O

        //- Write sizes, depth, time index and columns in newest-first order
        virtual bool writeData(Ostream& os) const;

        //- Validate and restore the stored columns
        virtual bool readData(Istream& is);

private:

        //- Reset the ring position
        void clearState();
};

} // End namespace Foam

#endif

// ************************************************************************* //
