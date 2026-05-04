/*---------------------------------------------------------------------------*\
                    Class DecomposedLaplacian Implementation
\*---------------------------------------------------------------------------*/

#include "DecomposedLaplacian.H"
#include "EigenMatrix.H"
#include "LLTMatrix.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //
DecomposedLaplacian::DecomposedLaplacian(const SquareMatrix<scalar>& similarityMatrix)
:
    d_(similarityMatrix.n())
{

    // Form the symmetric Laplacian = D - W.
    Laplacian_ = -similarityMatrix;
    for (label i = 0; i < d_; ++i) {
        for (label j = 0; j < d_; ++j) {
            Laplacian_[i][i] += similarityMatrix[i][j];
        }
    }

    // Create a symmetric EigenMatrix object to decompose L
    EigenMatrix em(Laplacian_, true);

    // Extract eigenvalues
    Lambda_ = em.EValsRe();
    Lambda_[0] = 0.0;
    for (numZeroEVals_ = 1; numZeroEVals_ < d_; ++numZeroEVals_) {
        if (Lambda_[numZeroEVals_] >= SMALL) {
            break;
        }
        Lambda_[numZeroEVals_] = 0.0;   
    }
    sqrtLambda_ = DiagonalMatrix<scalar>(d_, 0.0);
    for (label i = numZeroEVals_; i < d_; ++i) {
        if (Lambda_[i] < Lambda_[i-1]) {
            Info<< "Warning: Eigenvalues not sorted. Check the EigenMatrix implementation." << endl;
        }
        sqrtLambda_[i] = sqrt(Lambda_[i]);  
    }

    // Extract eigenvectors
    const SquareMatrix<scalar>& Q = em.EVecs();
    X_.setSize(d_);
    for (label i = 0; i < d_; ++i) {
        X_[i].setSize(d_);
        for (label j = 0; j < d_; ++j) {
            X_[i][j] = Q[i][j];
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar Foam::DecomposedLaplacian::dEff(
    const scalar mu
) const
{

    label omega;
    scalar sumEVals = 0.0;
    scalar sumSqrtEVals = 0.0;
    for (omega = numZeroEVals_; omega < d_; ++omega) {
        sumEVals += Lambda_[omega];
        sumSqrtEVals += sqrtLambda_[omega];
        if (sqrtLambda_[omega] * (1.0 + mu * sumEVals) < mu * Lambda_[omega] * sumSqrtEVals) {
            sumEVals -= Lambda_[omega];
            sumSqrtEVals -= sqrtLambda_[omega];
            break;
        }
    }

    scalar output = scalar(numZeroEVals_);
    for (label i = numZeroEVals_; i < omega; ++i) {
        scalar p = sqrtLambda_[i] * (1.0 + mu * sumEVals) / sumSqrtEVals - mu * Lambda_[i];
        output += p / (mu * Lambda_[i] + p);
    }
    return output;

}

LLTMatrix<scalar> Foam::DecomposedLaplacian::cholLambdaPlusVPi(
    const scalarField& Pi,
    const scalar mu
) const
{

    SquareMatrix<scalar> RegVPi(d_, 0.0);
    for (label i = 0; i < d_; ++i) {
        RegVPi[i][i] = mu * Lambda_[i];
    }

    for (label k = 0; k < d_; ++k) {
        const scalar Pik = Pi[k];
        const scalarField& Xk = X_[k];
        for (label i = 0; i < d_; ++i) {
            const scalar PikXki = Pik * Xk[i];
            for (label j = i; j < d_; ++j) {
                RegVPi[i][j] += PikXki * Xk[j];
            }
        }
    }

    for (label i = 0; i < d_; ++i) {
        for (label j = 0; j < i; ++j) {
            RegVPi[i][j] = RegVPi[j][i];
        }
    }

    return LLTMatrix<scalar>(RegVPi);

}

scalarField Foam::DecomposedLaplacian::DOptimalDesign(
    const scalar mu
) const
{

    const scalar maxFWIter = 1000;
    const scalar tolerance = 1e-6;
    const word dEffinition = "tight";   // "tight" or "dEff"
    const word criterion = "duality";   // "duality" or "KW"
    const word stepSize = "optimal";    // "optimal" or "FW"

    scalar dEff_ = 0.0;
    bool computeSumPigrad = (criterion == "duality");
    if (stepSize == "optimal" || criterion == "KW") {
        if (dEffinition == "dEff") {
            dEff_ = dEff(mu);
        } else if (dEffinition == "tight") {
            computeSumPigrad = true;
        } else {
            FatalErrorInFunction << "Unknown dEffinition option: " << dEffinition << exit(FatalError);
        }
    }

    scalarField Pi(d_, 1.0 / scalar(d_));
    scalarField grad(d_);
    scalarField si(d_);
    scalar sumPigrad = 0.0;
    scalar gap = 0.0;
    scalar gamma = 0.0;

    for (label k = 0; k < maxFWIter; ++k) {

        LLTMatrix<scalar> chol = cholLambdaPlusVPi(Pi, mu);
        for (label i = 0; i < d_; ++i) {
            chol.solve(si, X_[i]);
            grad[i] = sum(X_[i] * si);
        }
        if (computeSumPigrad) {
            sumPigrad = sum(Pi * grad);
            if (dEffinition == "tight") {
                dEff_ = sumPigrad;
            }
        }

        scalar gradMax = max(grad);
        if (criterion == "KW") {
            gap = mag((gradMax - dEff_) / dEff_);
        } else if (criterion == "duality") {
            gap = gradMax - sumPigrad;
        } else {
            FatalErrorInFunction << "Unknown criterion option: " << criterion << exit(FatalError);
        }
        if (gap < tolerance) {
            return Pi;
        }

        if (stepSize == "optimal") {
            gamma = max(0.0, min(1.0, (gradMax / dEff_ - 1.0) / max(gradMax - 1.0, SMALL)));
        } else if (stepSize == "FW") {
            gamma = 2.0 / (k + 2.0); // standard Frank-Wolfe step size guarantees O(1/k) convergence
        } else {
            FatalErrorInFunction << "Unknown stepSize option: " << stepSize << exit(FatalError);
        }

        scalar numMax = 0.0;
        for (label i = 0; i < d_; ++i) {
            if (grad[i] > gradMax - gap) {
                numMax += 1.0;
            }
        }
        Pi *= 1.0 - gamma;
        scalar gammaOverNumMax = gamma / numMax;
        for (label i = 0; i < d_; ++i) {
            if (grad[i] > gradMax - gap) {
                Pi[i] += gammaOverNumMax;
            }
        }
        Pi /= sum(Pi);
    }

    Info<< "Frank-Wolfe solver did not converge: " << gap << endl;
    return Pi;
    
}

LLTMatrix<scalar> Foam::DecomposedLaplacian::cholLaplacianPlusPi(
    const scalarField& Pi,
    const scalar mu
) const
{
    SquareMatrix<scalar> RegPi = mu * Laplacian_;
    for (label i = 0; i < d_; ++i) {
        RegPi[i][i] += Pi[i];
    }
    return LLTMatrix<scalar>(RegPi);
}

scalarField Foam::DecomposedLaplacian::getHat(
    const scalarField& Pi,
    const scalar mu,
    const label row
) const
{

    LLTMatrix<scalar> chol = cholLaplacianPlusPi(Pi, mu);
    scalarField hat(d_);
    scalarField e(d_, 0.0);
    e[row] = 1.0;
    chol.solve(hat, e);
    return hat;
    
}

Pair<scalarField> Foam::DecomposedLaplacian::getHatAndBonus(
    const scalarField& Pi,
    const scalar mu,
    const label row
) const
{

    LLTMatrix<scalar> chol = cholLaplacianPlusPi(Pi, mu);
    Pair<scalarField> output;
    scalarField hat(d_);
    scalarField bonus(d_);
    scalarField e(d_, 0.0);
    for (label i = 0; i < d_; ++i) {
        e[i] = 1.0;
        chol.solve(hat, e);
        bonus[i] = hat[i];
        if (i == row) {
            output[0] = hat;
        }
        e[i] = 0.0;
    }
    output[1] = bonus;
    return output;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
