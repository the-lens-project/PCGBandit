/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "FGAMGSolver.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::FGAMGSolver::interpolate
(
    solveScalarField& psi,
    solveScalarField& Apsi,
    const lduMatrix& m,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    solveScalar* __restrict__ psiPtr = psi.begin();

    const label* const __restrict__ uPtr = m.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = m.lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ diagPtr = m.diag().begin();
    const scalar* const __restrict__ upperPtr = m.upper().begin();
    const scalar* const __restrict__ lowerPtr = m.lower().begin();

    Apsi = 0;
    solveScalar* __restrict__ ApsiPtr = Apsi.begin();

    const label startRequest = UPstream::nRequests();

    m.initMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    const label nFaces = m.upper().size();
    for (label face=0; face<nFaces; face++)
    {
        ApsiPtr[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];
        ApsiPtr[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];
    }

    m.updateMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt,
        startRequest
    );

    const label nCells = m.diag().size();
    for (label celli=0; celli<nCells; celli++)
    {
        psiPtr[celli] = -ApsiPtr[celli]/(diagPtr[celli]);
    }
}


void Foam::FGAMGSolver::interpolate
(
    solveScalarField& psi,
    solveScalarField& Apsi,
    const lduMatrix& m,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const labelList& restrictAddressing,
    const solveScalarField& psiC,
    const direction cmpt
) const
{
    interpolate
    (
        psi,
        Apsi,
        m,
        interfaceBouCoeffs,
        interfaces,
        cmpt
    );

    const label nCells = m.diag().size();
    solveScalar* __restrict__ psiPtr = psi.begin();
    const scalar* const __restrict__ diagPtr = m.diag().begin();
    const solveScalar* const __restrict__ psiCPtr = psiC.begin();


    const label nCCells = psiC.size();
    solveScalarField corrC(nCCells, 0);
    solveScalar* __restrict__ corrCPtr = corrC.begin();

    solveScalarField diagC(nCCells, 0);
    solveScalar* __restrict__ diagCPtr = diagC.begin();

    for (label celli=0; celli<nCells; celli++)
    {
        corrCPtr[restrictAddressing[celli]] += diagPtr[celli]*psiPtr[celli];
        diagCPtr[restrictAddressing[celli]] += diagPtr[celli];
    }

    for (label ccelli=0; ccelli<nCCells; ccelli++)
    {
        corrCPtr[ccelli] = psiCPtr[ccelli] - corrCPtr[ccelli]/diagCPtr[ccelli];
    }

    for (label celli=0; celli<nCells; celli++)
    {
        psiPtr[celli] += corrCPtr[restrictAddressing[celli]];
    }
}


// ************************************************************************* //
