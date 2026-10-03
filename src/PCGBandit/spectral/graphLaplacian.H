/*---------------------------------------------------------------------------*\
                    Class graphLaplacian Declaration
\*---------------------------------------------------------------------------*/

#ifndef graphLaplacian_H
#define graphLaplacian_H
#include "SquareMatrix.H"
#include "LLTMatrix.H"
#include "Pair.H"
#include "labelList.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class graphLaplacian Declaration
\*---------------------------------------------------------------------------*/

class graphLaplacian
{
private:

    // Private Data

        const label d_;
        SquareMatrix<scalar> laplacian_;
        scalarField degree_; // degree vector

        labelList componentLabels_; // original row to connected component
        label numComponents_; // number of connected components

        labelList originalRows_; // stored row to original row
        labelList componentOffsets_; // start of each component in originalRows_
        labelList rowOffsets_; // stored row r has entries [rowOffsets_[r], rowOffsets_[r+1])

        labelList columns_; // original column index for each nonzero
        scalarField values_; // coefficient for each nonzero

public:

    // Constructors

        graphLaplacian(const SquareMatrix<scalar>& W);

    // Destructor

        ~graphLaplacian() = default;

    // Member Functions

        //- Apply L to distinct, pre-sized global vectors x and out.
        //  component=-1 uses all rows; otherwise zero out outside that component.
        void apply
        (
            const scalarField& x,
            scalarField& out,
            const label component = -1
        ) const;

        //- Eigenvalue-based upper bound on dTr(mu); computes eigendecomposition on each call.
        scalar dEff(scalar mu) const;
        
        //- Factor mu*L + diag(probs + SMALL) by full Cholesky.
        LLTMatrix<scalar> cholLapPlusProb(const scalarField& probs, const scalar mu) const;
        //- Return p maximizing Tr(diag(p)*A^-1) if trace, else log det(A).
        //  A=mu*L+diag(p+SMALL); optional traceObjective receives the trace
        //  at the returned p, even when maximizing log determinant.
        scalarField optimalDesign
        (
            const scalar mu,
            const bool trace,
            scalar* traceObjective = nullptr
        ) const;
        //--- Maximise Tr(diag(p)*(mu*L + diag(p) + SMALL*I)^-1), return objective value
        scalar dTr(const scalar mu) const;
        //--- Maximise log det(mu*L + diag(p) + SMALL*I) over p, return argmax
        scalarField DOptimalDesign(const scalar mu) const;

        //- Solve (mu*L+diag(probs+SMALL))*hat=e_row by CG, with Cholesky fallback.
        scalarField getHat(const scalarField& probs, const scalar mu, const label row) const;

        //- Return (A^-1*e_row, diag(A^-1)) via Cholesky; A=mu*L+diag(probs+SMALL).
        Pair<scalarField> getHatAndBonus(const scalarField& probs, const scalar mu, const label row) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
