/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "FGAMGPreconditioner.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FGAMGPreconditioner, 0);

    lduMatrix::preconditioner::addsymMatrixConstructorToTable
    <FGAMGPreconditioner> addFGAMGPreconditionerSymMatrixConstructorToTable_;

    lduMatrix::preconditioner::addasymMatrixConstructorToTable
    <FGAMGPreconditioner> addFGAMGPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FGAMGPreconditioner::FGAMGPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary& solverControls
)
:
    FGAMGSolver
    (
        sol.fieldName(),
        sol.matrix(),
        sol.interfaceBouCoeffs(),
        sol.interfaceIntCoeffs(),
        sol.interfaces(),
        solverControls
    ),
    lduMatrix::preconditioner(sol),
    nVcycles_(2)
{
    initVcycle
    (
        coarseCorrFields,
        coarseSources,
        smoothers,
        ApsiScratch,
        finestCorrectionScratch
    );

    readControls();

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FGAMGPreconditioner::readControls()
{
    FGAMGSolver::readControls();
    nVcycles_ = controlDict_.getOrDefault<label>("nVcycles", 2);
}


void Foam::FGAMGPreconditioner::precondition
(
    solveScalarField& wA,
    const solveScalarField& rA_ss,
    const direction cmpt
) const
{

    solveScalarField AwA(wA.size());
    solveScalarField finestCorrection(wA.size());
    solveScalarField finestResidual(rA_ss);

    /*
    // Create coarse grid correction fields
    PtrList<solveScalarField> coarseCorrFields;

    // Create coarse grid sources
    PtrList<solveScalarField> coarseSources;

    // Create the smoothers for all levels
    PtrList<lduMatrix::smoother> smoothers;

    // Scratch fields if processor-agglomerated coarse level meshes
    // are bigger than original. Usually not needed
    solveScalarField ApsiScratch;
    solveScalarField finestCorrectionScratch;

    // Initialise the above data structures
    initVcycle
    (
        coarseCorrFields,
        coarseSources,
        smoothers,
        ApsiScratch,
        finestCorrectionScratch
    );

    */

    // Adapt solveScalarField back to scalarField (as required)
    ConstPrecisionAdaptor<scalar, solveScalar> rA_adaptor(rA_ss);
    const scalarField& rA = rA_adaptor.cref();

    for (label cycle=0; cycle<nVcycles_; cycle++)
    {
        Vcycle
        (
            smoothers,
            wA,
            rA,
            AwA,
            finestCorrection,
            finestResidual,

            (ApsiScratch.size() ? ApsiScratch : AwA),
            (
                finestCorrectionScratch.size()
              ? finestCorrectionScratch
              : finestCorrection
            ),

            coarseCorrFields,
            coarseSources,
            cmpt
        );

        if (cycle < nVcycles_-1)
        {
            // Calculate finest level residual field
            matrix_.Amul(AwA, wA, interfaceBouCoeffs_, interfaces_, cmpt);
            finestResidual = rA_ss;
            finestResidual -= AwA;
        }
    }
}


// ************************************************************************* //
