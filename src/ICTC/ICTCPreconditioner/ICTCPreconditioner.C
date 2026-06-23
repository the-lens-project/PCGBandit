/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "ICTCPreconditioner.H"
#include <algorithm>

//#define ICTC_DEBUG

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ICTCPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<ICTCPreconditioner>
        addICTCPreconditionerSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ICTCPreconditioner::ICTCPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary& controlDict_
)
:
    lduMatrix::preconditioner(sol),
    diagL_(sol.matrix().diag().size()),
    lowerL_(controlDict_.getOrDefault<label>("storageCoeff", 200) * sol.matrix().diag().size()),
    rowAddrL_(controlDict_.getOrDefault<label>("storageCoeff", 200) * sol.matrix().diag().size()),
    colPtrL_(sol.matrix().diag().size() + 1)
{

    const scalar droptol = controlDict_.getOrDefault<scalar>("droptol", 1.0);
    const scalar pivmin = controlDict_.getOrDefault<scalar>("pivmin", 1e-16);
    label nnz = calcL(diagL_, lowerL_, rowAddrL_, colPtrL_, sol.matrix(), droptol, pivmin);
    Foam::debug::controlDict().set<label>("ICTC_NNZ", nnz);

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::ICTCPreconditioner::calcL
(
    solveScalarField& diagL_,
    solveScalarField& lowerL_,
    labelField& rowAddrL_,
    labelField& colPtrL_,
    const lduMatrix& matrix,
    const scalar droptol,
    const scalar pivmin
)
{

    const label n = diagL_.size();
    const label nl = matrix.lduAddr().lowerAddr().size();

    solveScalar* __restrict__ diagL = diagL_.begin();
    solveScalar* __restrict__ lowerL = lowerL_.begin();
    label* __restrict__ rowAddrL = rowAddrL_.begin();
    label* __restrict__ colPtrL = colPtrL_.begin();

    const scalar* const __restrict__ diagA = matrix.diag().begin();
    const scalar* const __restrict__ lowerA = matrix.lower().begin();
    const label* const __restrict__ rowAddrA = matrix.lduAddr().upperAddr().begin();
    const label* const __restrict__ colAddrA = matrix.lduAddr().lowerAddr().begin();

    solveScalarField scale_(n);
    solveScalar* __restrict__ scale = scale_.begin();
    solveScalarField signScale_(n);
    solveScalar* __restrict__ signScale = signScale_.begin();
    solveScalarField ta_(n);
    solveScalar* __restrict__ ta = ta_.begin();
    labelField itcol(n);
    labelField ifirst_(n);
    label* __restrict__ ifirst = ifirst_.begin();
    labelField list_(n);
    label* __restrict__ list = list_.begin();

    // --- Assumes sign is determined by first diagonal entry
    bool posdef = (diagA[0] > 0.0);
    for (label i = 0; i < n; i++) {
        #ifdef ICTC_DEBUG
        if (posdef != (diagA[i] > 0.0) || diagA[i] == 0.0) {
            FatalErrorIn("Foam::ICTCPreconditioner::calcL") << "matrix determined to be indefinite at column " << i+1 << exit(FatalError);
        }
        #endif
        scale[i] = 1.0 / Foam::sqrt(posdef ? diagA[i] : -diagA[i]);
        signScale[i] = posdef ? scale[i] : -scale[i];
        ifirst[i] = -1;
        list[i] = -1;
        diagL[i] = 0.0;
    }
    colPtrL[0] = 0;
    label ipos = 0;
    label jpos = 0;

    // --- Loop over all columns
    for (label k = 0; k < n; k++) {

        // --- Load column k into ta
        label talen = 0;
    
        // --- Add diagonal element with scaling
        diagL[k] += diagA[k] * signScale[k] * scale[k];
    
        // --- Add off-diagonal elements with scaling
        while (jpos < nl && colAddrA[jpos] == k) {
            label row = rowAddrA[jpos];
            ta[row] = lowerA[jpos] * signScale[k] * scale[row];
            itcol[talen] = row;
            talen++;

            // --- Mark that row has a nonzero value (this is just a flag, not linked list)
            ifirst[row] = 1;
            jpos++;

        }
    
        // --- Update column k using the previous columns
        label j = list[k];
        while (j != -1) { 
            label isj = ifirst[j];
            label iej = colPtrL[j+1] - 1;
            scalar lval = lowerL[isj];
            isj++;
    
            if (isj < iej) {
                ifirst[j] = isj;
                label iptr = j;
                j = list[j];
                list[iptr] = list[rowAddrL[isj]];
                list[rowAddrL[isj]] = iptr;
            } else {
                j = list[j];
            }   
    
            // --- Apply pivot minimum value
            scalar temp = diagL[k];
            if (temp < pivmin) {
                temp = pivmin;
            }   
    
            // --- Apply drop tolerance (alpha = 0)
            if (Foam::mag(lval / Foam::sqrt(temp)) <= droptol) {
                continue;
            }   
    
            for (label i = isj; i <= iej; i++) {
                label row = rowAddrL[i];
                if (ifirst[row] != -1) { 
                    ta[row] = ta[row] - lval * lowerL[i];
                } else {

                    // --- This is a flag, not part of linked list
                    ifirst[row] = 1;
                    itcol[talen] = row;
                    talen++;
                    ta[row] = - lval * lowerL[i];

                }   
            }   
        }   
    
        // --- Drop using droptol
        if (droptol > 0.0) {
            j = 0;
            for (label i = 0; i < talen; i++) {
                label row = itcol[i];
                if (Foam::mag(ta[row]) > droptol) {
                    itcol[j] = row;
                    j++;
                } else {
                    ifirst[row] = -1;
                }   
            }   
            talen = j;
        }   
    
        // --- Sort row indices in ascending order 
        // std::sort(itcol.begin(), itcol.begin() + talen);

        // --- Arrays are small so insertion sort is faster than std::sort
        for (label i = 1; i < talen; i++) {
            label j;
            label key = itcol[i];
            for (j = i; j > 0 && itcol[j-1] > key; j--) {
                itcol[j] = itcol[j-1];
            }
            itcol[j] = key;
        }
    
        // --- Now we can do sqrt
        if (diagL[k] < pivmin) {
            diagL[k] = pivmin;
        }
        #ifdef ICTC_DEBUG
        if (diagL[k] <= 0.0) {
            FatalErrorIn("Foam::ICTCPreconditioner::calcL") << "zero or negative pivot encountered in column " << k+1 << exit(FatalError);
        }
        #endif

        // --- Store reciprocal
        diagL[k] = 1.0 / Foam::sqrt(diagL[k]);
    
        // --- Scale and update remaining diagonals using nondropped elements only
        for (label j = 0; j < talen; j++) {
            label row = itcol[j];
            ta[row] *= diagL[k];
            diagL[row] -= ta[row] * ta[row];
        }   

        #ifdef ICTC_DEBUG
        // --- Check that we have enough storage
        if (ipos + talen > lowerL_.size()) {
            FatalErrorIn("Foam::ICTCPreconditioner::calcL") << "insufficient storage for factor at column " << k+1 << exit(FatalError);
        }
        #endif
    
        // --- Put the largest elements back into the sparse data structure
        for (label icount = 0; icount < talen; icount++) {
            lowerL[ipos] = ta[itcol[icount]];
            rowAddrL[ipos] = itcol[icount];
            ipos++;
        }   
        colPtrL[k+1] = ipos;
    
        // --- Variables ifirst and list keep track of where in column k we are
        if (talen > 0) {
            label isk = colPtrL[k];
            label iptr = rowAddrL[isk];
            list[k] = list[iptr];
            list[iptr] = k;
            ifirst[k] = isk;
        }   
    
        // --- Reset ifirst for next column
        for (label j = 0; j < talen; j++) {
            ifirst[itcol[j]] = -1; 
        }   

    }   

    // --- Scale back the factor
    for (label i = 0; i < n; i++) {
        diagL[i] *= signScale[i];
        for (label j = colPtrL[i]; j < colPtrL[i+1]; j++) {
            lowerL[j] /= signScale[rowAddrL[j]];
        }   
    }   

    return ipos+1;

}


void Foam::ICTCPreconditioner::precondition
(
    solveScalarField& wA,
    const solveScalarField& rA,
    const direction
) const
{
    solveScalar* __restrict__ sol = wA.begin();
    const solveScalar* __restrict__ rhs = rA.begin();

    const solveScalar* __restrict__ diagL = diagL_.begin();
    const solveScalar* __restrict__ lowerL = lowerL_.begin();
    const label* __restrict__ rowAddrL = rowAddrL_.begin();
    const label* __restrict__ colPtrL = colPtrL_.begin();

    const label n = wA.size();

    // --- Copy rhs to sol
    for (label i = 0; i < n; i++) {
        sol[i] = rhs[i];
    }

    // --- Forward solve with L
    for (label i = 0; i < n; i++) {
        sol[i] *= diagL[i];
        scalar s = sol[i];
        for (label k = colPtrL[i]; k < colPtrL[i+1]; k++) {
            sol[rowAddrL[k]] -= lowerL[k] * s;
        }
    }

    // --- Backward solve with L^T
    for (label i = n-1; i >= 0; i--) {
        scalar s = sol[i];
        for (label k = colPtrL[i]; k < colPtrL[i+1]; k++) {
            s -= lowerL[k] * sol[rowAddrL[k]];
        }
        sol[i] = s * diagL[i];
    }

}


// ************************************************************************* //
