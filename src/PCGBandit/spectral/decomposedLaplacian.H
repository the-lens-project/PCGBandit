/*---------------------------------------------------------------------------*\
                    Class decomposedLaplacian Declaration
\*---------------------------------------------------------------------------*/

#ifndef decomposedLaplacian_H
#define decomposedLaplacian_H
#include "DiagonalMatrix.H"
#include "SquareMatrix.H"
#include "LLTMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class decomposedLaplacian Declaration
\*---------------------------------------------------------------------------*/

class decomposedLaplacian
{
private:

    // Private Data

        const label d_;
        label numZeroEVals_;

        //- Laplacian
        SquareMatrix<scalar> Laplacian_;

        //- Eigenvalues of the Laplacian matrix
        DiagonalMatrix<scalar> Lambda_;
        DiagonalMatrix<scalar> sqrtLambda_;

        //- Eigenrows of the Laplacian matrix
        List<scalarField> X_;

        //- Maximise log determinant or effective dimension on the simplex
        scalarField optimalDesign(const scalar mu, const bool trace) const;

public:

    // Public Data

    // Constructors

        decomposedLaplacian(const SquareMatrix<scalar>& W);

    // Destructor

        ~decomposedLaplacian() = default;

    // Member Functions

        DiagonalMatrix<scalar> EVals() const {
            return Lambda_;
        }

        List<scalarField> ERows() const {
            return X_;
        }

        LLTMatrix<scalar> cholLambdaPlusVPi(const scalarField& Pi, const scalar mu) const;
        LLTMatrix<scalar> cholLaplacianPlusPi(const scalarField& Pi, const scalar mu) const;
        scalar dEff(scalar mu) const;
        //- Maximise log det(mu*L + diag(Pi) + SMALL*I)
        scalarField DOptimalDesign(const scalar mu) const;

        //- Maximise Tr(diag(Pi)*(mu*L + diag(Pi) + SMALL*I)^-1)
        scalarField traceOptimalDesign(const scalar mu = 1) const;
        scalarField getHat(const scalarField& Pi, const scalar mu, const label row) const;
        Pair<scalarField> getHatAndBonus(const scalarField& Pi, const scalar mu, const label row) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
