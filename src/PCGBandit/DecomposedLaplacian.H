/*---------------------------------------------------------------------------*\
                    Class DecomposedLaplacian Declaration
\*---------------------------------------------------------------------------*/

#ifndef DecomposedLaplacian_H
#define DecomposedLaplacian_H
#include "DiagonalMatrix.H"
#include "SquareMatrix.H"
#include "LLTMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class DecomposedLaplacian Declaration
\*---------------------------------------------------------------------------*/

class DecomposedLaplacian
{
private:

    // Private Data

        const label d_;
        label numZeroEVals_;

        const label maxFWIter_ = 1000;
        const scalar tolerance_ = 1e-8;

        //- Laplacian
        SquareMatrix<scalar> Laplacian_;

        //- Eigenvalues of the Laplacian matrix
        DiagonalMatrix<scalar> Lambda_;
        DiagonalMatrix<scalar> sqrtLambda_;

        //- Eigenrows of the Laplacian matrix
        List<scalarField> X_;


public:

    // Public Data

    // Constructors

        DecomposedLaplacian(const SquareMatrix<scalar>& similarityMatrix);

    // Destructor

        ~DecomposedLaplacian() = default;

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
        scalarField DOptimalDesign(const scalar mu) const;
        scalarField getHat(const scalarField& Pi, const scalar mu, const label row) const;
        Pair<scalarField> getHatAndBonus(const scalarField& Pi, const scalar mu, const label row) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
