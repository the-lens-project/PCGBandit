/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "SORSmootherBase.H"
#include "PrecisionAdaptor.H"


Foam::SORSmootherBase::SORSmootherBase
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls,
    const scalar omega
)
:
    GaussSeidelSmoother
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    ),
    omega_(omega)
{}


void Foam::SORSmootherBase::smooth
(
    solveScalarField& psi,
    const lduMatrix& matrix,
    const solveScalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt,
    const label nSweeps,
    const scalar omega
)
{
    solveScalar* __restrict__ psiPtr = psi.begin();

    const label nCells = psi.size();

    solveScalarField& bPrime = matrix.work(nCells);
    solveScalar* __restrict__ bPrimePtr = bPrime.begin();

    const scalar* const __restrict__ diagPtr = matrix.diag().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();
    const scalar* const __restrict__ lowerPtr = matrix.lower().begin();

    const label* const __restrict__ upperAddrPtr =
        matrix.lduAddr().upperAddr().begin();
    const label* const __restrict__ ownerStartPtr =
        matrix.lduAddr().ownerStartAddr().begin();

    // This is the OpenFOAM Gauss-Seidel sweep with point-wise relaxation
    // applied before the updated value is used by subsequent rows.
    for (label sweep = 0; sweep < nSweeps; ++sweep)
    {
        bPrime = source;

        const label startRequest = UPstream::nRequests();

        matrix.initMatrixInterfaces
        (
            false,
            interfaceBouCoeffs,
            interfaces,
            psi,
            bPrime,
            cmpt
        );

        matrix.updateMatrixInterfaces
        (
            false,
            interfaceBouCoeffs,
            interfaces,
            psi,
            bPrime,
            cmpt,
            startRequest
        );

        label faceStart;
        label faceEnd = ownerStartPtr[0];

        for (label celli = 0; celli < nCells; ++celli)
        {
            faceStart = faceEnd;
            faceEnd = ownerStartPtr[celli + 1];

            solveScalar psii = bPrimePtr[celli];
            for (label facei = faceStart; facei < faceEnd; ++facei)
            {
                psii -= upperPtr[facei]*psiPtr[upperAddrPtr[facei]];
            }

            psii =
                psiPtr[celli]
              + omega*(psii/diagPtr[celli] - psiPtr[celli]);

            for (label facei = faceStart; facei < faceEnd; ++facei)
            {
                bPrimePtr[upperAddrPtr[facei]] -= lowerPtr[facei]*psii;
            }

            psiPtr[celli] = psii;
        }
    }
}


void Foam::SORSmootherBase::smooth
(
    solveScalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    smooth
    (
        psi,
        matrix_,
        ConstPrecisionAdaptor<solveScalar, scalar>(source),
        interfaceBouCoeffs_,
        interfaces_,
        cmpt,
        nSweeps,
        omega_
    );
}


void Foam::SORSmootherBase::scalarSmooth
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    smooth
    (
        psi,
        matrix_,
        source,
        interfaceBouCoeffs_,
        interfaces_,
        cmpt,
        nSweeps,
        omega_
    );
}

// ************************************************************************* //
