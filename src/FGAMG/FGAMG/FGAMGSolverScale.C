/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "FGAMGSolver.H"
#include "FixedList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::FGAMGSolver::scale
(
    solveScalarField& field,
    solveScalarField& Acf,
    const lduMatrix& A,
    const FieldField<Field, scalar>& interfaceLevelBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaceLevel,
    const solveScalarField& source,
    const direction cmpt
) const
{
    A.Amul
    (
        Acf,
        field,
        interfaceLevelBouCoeffs,
        interfaceLevel,
        cmpt
    );


    const label nCells = field.size();
    solveScalar* __restrict__ fieldPtr = field.begin();
    const solveScalar* const __restrict__ sourcePtr = source.begin();
    const solveScalar* const __restrict__ AcfPtr = Acf.begin();


    FixedList<solveScalar, 2> scalingFactor(Zero);

    for (label i=0; i<nCells; i++)
    {
        scalingFactor[0] += fieldPtr[i]*sourcePtr[i];
        scalingFactor[1] += fieldPtr[i]*AcfPtr[i];
    }

    A.mesh().reduce(scalingFactor, sumOp<solveScalar>());

    const solveScalar sf =
    (
        scalingFactor[0]
      / stabilise(scalingFactor[1], pTraits<solveScalar>::vsmall)
    );

    if (debug >= 2)
    {
        Pout<< sf << " ";
    }

    const scalarField& D = A.diag();
    const scalar* const __restrict__ DPtr = D.begin();

    for (label i=0; i<nCells; i++)
    {
        fieldPtr[i] = sf*fieldPtr[i] + (sourcePtr[i] - sf*AcfPtr[i])/DPtr[i];
    }
}


// ************************************************************************* //
