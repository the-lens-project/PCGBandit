/*---------------------------------------------------------------------------*\
Class
    Foam::PCGBandit
Group
    grpLduMatrixSolvers
Description
    Preconditioned conjugate gradient solver for symmetric lduMatrices.
    Selects preconditioner and initial-guess configurations using Tsallis-INF,
    SpectralINF, or Thompson sampling.
SourceFiles
    PCGBandit.C
\*---------------------------------------------------------------------------*/
#ifndef PCGBandit_H
#define PCGBandit_H
#include "lduMatrix.H"
#include "autoPtr.H"
#include "mirrorDescent.H"
#include "similarityMatrix.H"
#include "subspaceInitializer.H"
#include "costEstimation.H"
#include "configurationSpace.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
/*---------------------------------------------------------------------------*\
                             Class PCGBandit Declaration
\*---------------------------------------------------------------------------*/
class PCGBandit
:
    public lduMatrix::solver
{


private:

    // Private Member Functions
        //- No copy construct
        PCGBandit(const PCGBandit&) = delete;
        //- No copy assignment
        void operator=(const PCGBandit&) = delete;

    word banditName_;
    word lossEstimator_;
    bool deterministic_;
    bool randomUniform_;
    word banditAlgorithm_;
    label backstop_;
    label static_;
    autoPtr<subspaceInitializer> subspace_;
    costEstimation costs_;


    void queryLearner(const scalar initialResidual) const;

public:
    //- Runtime type information
    TypeName("PCGBandit");

    // Constructors
        //- Construct from matrix components and solver controls
        PCGBandit
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );

    //- Destructor
    virtual ~PCGBandit() = default;

    // Member Functions
        //- Solve the matrix with this solver
        virtual solverPerformance scalarSolve
        (
            solveScalarField& psi,
            const solveScalarField& source,
            const direction cmpt=0
        ) const;
        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt=0
        ) const;
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
// ************************************************************************* //
