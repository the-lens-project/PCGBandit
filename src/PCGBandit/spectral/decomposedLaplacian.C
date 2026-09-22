/*---------------------------------------------------------------------------*\
                    Class decomposedLaplacian Implementation
\*---------------------------------------------------------------------------*/

#include "decomposedLaplacian.H"
#include "EigenMatrix.H"
#include "LLTMatrix.H"
#include "SquareMatrix.H"

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //
decomposedLaplacian::decomposedLaplacian(const SquareMatrix<scalar>& W)
:
    d_(W.n())
{
    if (!d_)
    {
        FatalErrorInFunction << "Cannot decompose an empty graph" << exit(FatalError);
    }

    // Form the symmetric Laplacian = D - W.
    Laplacian_ = -W;
    for (label i = 0; i < d_; ++i) {
        for (label j = 0; j < d_; ++j) {
            Laplacian_[i][i] += W[i][j];
        }
    }

    // Create a symmetric EigenMatrix object to decompose L
    EigenMatrix<scalar> em(Laplacian_, true);

    // Each connected component contributes an exact zero eigenvalue.
    // Count components instead of thresholding eigensolver round-off.
    List<bool> visited(d_, false);
    labelList pending(d_);
    numZeroEVals_ = 0;
    for (label root = 0; root < d_; ++root)
    {
        if (visited[root]) continue;
        ++numZeroEVals_;
        visited[root] = true;
        label size = 1;
        pending[0] = root;
        while (size)
        {
            const label i = pending[--size];
            for (label j = 0; j < d_; ++j)
            {
                if (!visited[j] && W(i, j) > 0)
                {
                    visited[j] = true;
                    pending[size++] = j;
                }
            }
        }
    }
    Lambda_ = em.EValsRe();
    for (label i = 0; i < numZeroEVals_; ++i) Lambda_[i] = 0;
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

scalar Foam::decomposedLaplacian::dEff(
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

LLTMatrix<scalar> Foam::decomposedLaplacian::cholLambdaPlusVPi(
    const scalarField& Pi,
    const scalar mu
) const
{

    SquareMatrix<scalar> RegVPi(d_, 0.0);
    for (label i = 0; i < d_; ++i) {
        RegVPi[i][i] = mu * Lambda_[i] + SMALL;
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

scalarField Foam::decomposedLaplacian::DOptimalDesign(const scalar mu) const
{
    return optimalDesign(mu, false);
}


scalarField Foam::decomposedLaplacian::traceOptimalDesign(const scalar mu) const
{
    return optimalDesign(mu, true);
}


scalarField Foam::decomposedLaplacian::optimalDesign
(
    const scalar mu,
    const bool trace
) const
{
    if (!std::isfinite(mu) || mu < 0)
    {
        FatalErrorInFunction << "mu must be finite and non-negative"
            << exit(FatalError);
    }

    const label maxIterations = 100000;
    const scalar tolerance = 1e-8;
    scalarField Pi(d_, 1.0/scalar(d_));
    if (mu == 0 || d_ == 1) return Pi;

    SquareMatrix<scalar> inverse(d_, 0.0);
    SquareMatrix<scalar> metric(trace ? d_ : 0, 0.0);
    scalarField gradient(d_), column(d_), unit(d_, 0.0);
    scalarField hi(d_), hj(d_), gi(d_), gj(d_), u(d_), v(d_);
    scalar gap = GREAT;
    bool refresh = true;

    for (label iteration = 0; iteration <= maxIterations; ++iteration)
    {
        const bool fresh = refresh || iteration == maxIterations
            || iteration % max(label(10), d_) == 0;
        if (fresh)
        {
            const LLTMatrix<scalar> chol = cholLaplacianPlusPi(Pi, mu);
            for (label j = 0; j < d_; ++j)
            {
                unit[j] = 1;
                chol.solve(column, unit);
                for (label i = 0; i < d_; ++i) inverse(i, j) = column[i];
                unit[j] = 0;
            }

            if (trace)
            {
                // metric = H*(mu*L + SMALL*I)*H, with H the inverse.
                // Row differences avoid subtracting nearly equal inverses.
                SquareMatrix<scalar> regularised(d_, 0.0);
                for (label i = 0; i < d_; ++i)
                {
                    for (label j = 0; j < d_; ++j)
                    {
                        scalar value = SMALL*inverse(i, j);
                        for (label k = 0; k < d_; ++k)
                        {
                            value -= mu*Laplacian_(i, k)
                                *(inverse(i, j) - inverse(k, j));
                        }
                        regularised(i, j) = value;
                    }
                }
                for (label i = 0; i < d_; ++i)
                {
                    for (label j = i; j < d_; ++j)
                    {
                        scalar value = 0;
                        for (label k = 0; k < d_; ++k)
                        {
                            value += inverse(i, k)*regularised(k, j);
                        }
                        metric(i, j) = metric(j, i) = value;
                    }
                }
            }
            refresh = false;
        }

        label receiver = 0, donor = -1;
        scalar meanGradient = 0;
        for (label i = 0; i < d_; ++i)
        {
            gradient[i] = trace ? metric(i, i) : inverse(i, i);
            meanGradient += Pi[i]*gradient[i];
            if (gradient[i] > gradient[receiver]) receiver = i;
            if (Pi[i] > 0 && (donor < 0 || gradient[i] < gradient[donor]))
            {
                donor = i;
            }
        }
        // The simplex Frank-Wolfe gap bounds the remaining objective error.
        gap = gradient[receiver] - meanGradient;
        if (gap <= tolerance*max(scalar(1), mag(meanGradient)))
        {
            if (fresh) return Pi;
            refresh = true;
            continue;
        }
        if (iteration == maxIterations) break;

        // Vertex exchange: move mass from the weakest supported arm.
        // D step: Harman et al., arXiv:1801.05661, Appendix A.1.
        const scalar a = inverse(receiver, receiver);
        const scalar b = inverse(donor, donor);
        const scalar c = inverse(receiver, donor);
        const scalar slope = a - b;
        const scalar curvature = max(scalar(0), a*b - c*c);
        scalar step = Pi[donor];
        if (trace)
        {
            const scalar first = gradient[receiver] - gradient[donor];
            const scalar second = b*gradient[receiver] + a*gradient[donor]
                - 2*c*metric(receiver, donor);
            // Along this exchange the objective gain is
            // (first*t - second*t^2)/(1 + slope*t - curvature*t^2).
            scalar lo = 0, hiStep = step;
            for (label k = 0; k < 60; ++k)
            {
                const scalar t = 0.5*(lo + hiStep);
                const scalar det = 1 + slope*t - curvature*t*t;
                const scalar derivative = first - 2*second*t
                    + (first*curvature - second*slope)*t*t;
                if (det > 0 && derivative > 0) lo = t;
                else hiStep = t;
            }
            step = 0.5*(lo + hiStep);
        }
        else if (curvature > 0)
        {
            // Exact maximum of det(A + t*(e_i e_i^T - e_j e_j^T)).
            step = min(step, slope/(2*curvature));
        }

        const scalar det = 1 + slope*step - curvature*step*step;
        if (!(step > 0) || !(det > 0) || !std::isfinite(det))
        {
            if (fresh) break;
            refresh = true;
            continue;
        }

        const scalar kii = step*(1 - step*b)/det;
        const scalar kij = step*step*c/det;
        const scalar kjj = -step*(1 + step*a)/det;
        for (label k = 0; k < d_; ++k)
        {
            hi[k] = inverse(k, receiver);
            hj[k] = inverse(k, donor);
            u[k] = kii*hi[k] + kij*hj[k];
            v[k] = kij*hi[k] + kjj*hj[k];
            if (trace)
            {
                gi[k] = metric(k, receiver);
                gj[k] = metric(k, donor);
            }
        }
        const scalar gii = trace ? metric(receiver, receiver) : 0;
        const scalar gij = trace ? metric(receiver, donor) : 0;
        const scalar gjj = trace ? metric(donor, donor) : 0;
        // Rank-two updates cost O(d^2); periodic Cholesky solves limit drift.
        for (label i = 0; i < d_; ++i)
        {
            for (label j = i; j < d_; ++j)
            {
                const scalar value = inverse(i, j) - u[i]*hi[j] - v[i]*hj[j];
                inverse(i, j) = inverse(j, i) = value;
                if (trace)
                {
                    const scalar value = metric(i, j)
                        - u[i]*gi[j] - v[i]*gj[j] - gi[i]*u[j] - gj[i]*v[j]
                        + u[i]*(gii*u[j] + gij*v[j])
                        + v[i]*(gij*u[j] + gjj*v[j]);
                    metric(i, j) = metric(j, i) = value;
                }
            }
        }
        Pi[receiver] += step;
        Pi[donor] -= step;
    }

    WarningInFunction << (trace ? "Trace" : "D-optimal")
        << " design did not converge; simplex gap = " << gap << endl;
    return Pi;
}

LLTMatrix<scalar> Foam::decomposedLaplacian::cholLaplacianPlusPi(
    const scalarField& Pi,
    const scalar mu
) const
{
    SquareMatrix<scalar> RegPi = mu * Laplacian_;
    for (label i = 0; i < d_; ++i) {
        RegPi[i][i] += Pi[i] + SMALL;
    }
    return LLTMatrix<scalar>(RegPi);
}

scalarField Foam::decomposedLaplacian::getHat(
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

Pair<scalarField> Foam::decomposedLaplacian::getHatAndBonus(
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
