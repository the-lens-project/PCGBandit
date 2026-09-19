/*---------------------------------------------------------------------------*\
                          Class siPCG Implementation
\*---------------------------------------------------------------------------*/

#include "siPCG.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(siPCG, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<siPCG>
        addsiPCGSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::siPCG::siPCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    ),
    subspace_(matrix, interfaceBouCoeffs, interfaces, fieldName, solverControls)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::siPCG::scalarSolve
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt
) const
{
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    label nCells = psi.size();

    solveScalar* __restrict__ psiPtr = psi.begin();

    solveScalarField pA(nCells);
    solveScalar* __restrict__ pAPtr = pA.begin();

    solveScalarField wA(nCells);
    solveScalar* __restrict__ wAPtr = wA.begin();

    solveScalar wArA = solverPerf.great_;
    solveScalar wArAold = wArA;

    // Attach state before reading it; leave psi unchanged until initialize().
    subspace_.update(psi, cmpt);

    // Calculate A*psi and the initial residual.
    matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    solveScalarField rA(source - wA);
    solveScalar* __restrict__ rAPtr = rA.begin();

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        true
    );

    solveScalar normFactor = this->normFactor(psi, source, wA, pA);

    if ((log_ >= 2) || (lduMatrix::debug >= 2))
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())
       /normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Skip initialization and state updates if already converged. Reducing
    // the energy norm can increase the residual 1-norm used for convergence.
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
    )
    {
        // Preserve the incoming normalization and relative tolerance target.
        subspace_.initialize(cmpt, normFactor, psi, rA, solverPerf);

        // Initialization may converge without constructing a preconditioner.
        if
        (
            minIter_ > 0
         || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
        )
        {
            if (!preconPtr_)
            {
                preconPtr_ = lduMatrix::preconditioner::New
                (
                    *this,
                    controlDict_
                );
            }

            do
            {
                wArAold = wArA;

                preconPtr_->precondition(wA, rA, cmpt);

                // Update the search direction.
                wArA = gSumProd(wA, rA, matrix().mesh().comm());

                if (solverPerf.nIterations() == 0)
                {
                    for (label cell=0; cell<nCells; cell++)
                    {
                        pAPtr[cell] = wAPtr[cell];
                    }
                }
                else
                {
                    const solveScalar beta = wArA/wArAold;

                    for (label cell=0; cell<nCells; cell++)
                    {
                        pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                    }
                }


                // Multiply the search direction by A.
                matrix_.Amul(wA, pA, interfaceBouCoeffs_, interfaces_, cmpt);

                solveScalar wApA = gSumProd(wA, pA, matrix().mesh().comm());

                if (solverPerf.checkSingularity(mag(wApA)/normFactor)) break;


                // Update the solution and residual.

                const solveScalar alpha = wArA/wApA;

                for (label cell=0; cell<nCells; cell++)
                {
                    psiPtr[cell] += alpha*pAPtr[cell];
                    rAPtr[cell] -= alpha*wAPtr[cell];
                }

                solverPerf.finalResidual() =
                    gSumMag(rA, matrix().mesh().comm())
                   /normFactor;

            } while
            (
                (
                  ++solverPerf.nIterations() < maxIter_
                && !solverPerf.checkConvergence(tolerance_, relTol_, log_)
                )
             || solverPerf.nIterations() < minIter_
            );
        }
    }

    if (preconPtr_)
    {
        preconPtr_->setFinished(solverPerf);
    }

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        false
    );

    return solverPerf;
}


Foam::solverPerformance Foam::siPCG::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    return scalarSolve
    (
        tpsi.ref(),
        ConstPrecisionAdaptor<solveScalar, scalar>(source)(),
        cmpt
    );
}


// ************************************************************************* //
