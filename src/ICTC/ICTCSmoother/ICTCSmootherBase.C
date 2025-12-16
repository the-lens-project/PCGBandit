/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "ICTCSmootherBase.H"
#include "../ICTCPreconditioner/ICTCPreconditioner.H"
#include "PrecisionAdaptor.H"
#include <algorithm>


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ICTCSmootherBase::ICTCSmootherBase
(
     const word& fieldName,
     const lduMatrix& matrix,
     const FieldField<Field, scalar>& interfaceBouCoeffs,
     const FieldField<Field, scalar>& interfaceIntCoeffs,
     const lduInterfaceFieldPtrsList& interfaces,
     const scalar droptol
)
:
    lduMatrix::smoother
    (
		fieldName, 
		matrix, 
		interfaceBouCoeffs, 
		interfaceIntCoeffs, 
		interfaces
    ),
    
    diagL_(matrix_.diag().size()),
    lowerL_(200 * matrix_.diag().size()),
    rowAddrL_(200 * matrix_.diag().size()),
    colPtrL_(matrix_.diag().size() + 1)
{
    const scalar pivmin = 1e-16;
    ICTCPreconditioner::calcL(diagL_, lowerL_, rowAddrL_, colPtrL_, matrix_, droptol, pivmin);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ICTCSmootherBase::smooth
(
     solveScalarField& psi,
     const scalarField& source,
     const direction cmpt,
     const label nSweeps
) const
{
    const solveScalar* __restrict__ diagL = diagL_.begin();
    const solveScalar* __restrict__ lowerL = lowerL_.begin();
    const label* __restrict__ rowAddrL = rowAddrL_.begin();
    const label* __restrict__ colPtrL = colPtrL_.begin();

    solveScalarField rA(diagL_.size());

    solveScalar* __restrict__ rAPtr = rA.begin();

    const label n = rA.size();

    bool posdef = (diagL[0] > 0.0);

    for (label sweeps = 0; sweeps < nSweeps; sweeps++) {

		 matrix_.residual
         (
             rA,
             psi,
             source,
             interfaceBouCoeffs_,
             interfaces_,
             cmpt
         );

    	// --- Forward solve with L
    	for (label i = 0; i < n; i++) {

            rAPtr[i] *= diagL[i];
            scalar s = rAPtr[i];
            for (label k = colPtrL[i]; k < colPtrL[i+1]; k++) {
           	rAPtr[rowAddrL[k]] -= lowerL[k] * s;
            }
    	}

    	// --- Backward solve with L^T
    	for (label i = n-1; i >= 0; i--) {
            scalar s = rAPtr[i];
            for (label k = colPtrL[i]; k < colPtrL[i+1]; k++) {
            	s -= lowerL[k] * rAPtr[rowAddrL[k]];
            }

            rAPtr[i] = s * diagL[i];
    	}

	if (posdef) {
		psi += rA; 
	}
	else {
		psi -= rA;
	}


    }

}

void Foam::ICTCSmootherBase::scalarSmooth
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
        ConstPrecisionAdaptor<scalar, solveScalar>(source),
        cmpt,
        nSweeps
    );

}

// ************************************************************************* //
