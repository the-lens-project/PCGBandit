/*---------------------------------------------------------------------------*\
                    Class graphLaplacian Implementation
\*---------------------------------------------------------------------------*/

#include "graphLaplacian.H"
#include "conjugateGradient.H"
#include "DiagonalMatrix.H"
#include "EigenMatrix.H"
#include "DynamicList.H"
#include "LLTMatrix.H"
#include "SquareMatrix.H"

#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //
graphLaplacian::graphLaplacian(const SquareMatrix<scalar>& W)
:
    d_(W.n())
{
    if (!d_)
    {
        FatalErrorInFunction << "Cannot construct a Laplacian for an empty graph" << exit(FatalError);
    }

    // Form the symmetric Laplacian = D - W.
    laplacian_ = -W;
    for (label i = 0; i < d_; ++i) {
        for (label j = 0; j < d_; ++j) {
            laplacian_[i][i] += W[i][j];
        }
    }
    degree_ = laplacian_.diag();

    // Discover components and store their rows consecutively as we visit them.
    // Each component contributes one exact zero eigenvalue.
    List<bool> visited(d_, false);
    componentLabels_.setSize(d_);
    originalRows_.setSize(d_);
    componentOffsets_.setSize(d_ + 1);
    labelList pending(d_);
    label row = 0;
    numComponents_ = 0;
    for (label root = 0; root < d_; ++root)
    {
        if (visited[root]) continue;
        const label component = numComponents_++;
        componentOffsets_[component] = row;
        visited[root] = true;
        label size = 1;
        pending[0] = root;
        while (size)
        {
            const label i = pending[--size];
            componentLabels_[i] = component;
            originalRows_[row++] = i;
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
    componentOffsets_[numComponents_] = row;
    componentOffsets_.setSize(numComponents_ + 1);

    // Append each row's nonzeros and record where that row starts.
    DynamicList<label> columns;
    DynamicList<scalar> values;
    rowOffsets_.setSize(d_ + 1);
    forAll(originalRows_, r)
    {
        const label i = originalRows_[r];
        rowOffsets_[r] = values.size();
        for (label j = 0; j < d_; ++j)
        {
            if (laplacian_(i, j) != 0)
            {
                columns.append(j);
                values.append(laplacian_(i, j));
            }
        }
    }
    // The last offset closes the final row, including when that row is empty.
    rowOffsets_[d_] = values.size();
    columns_.transfer(columns);
    values_.transfer(values);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar Foam::graphLaplacian::dEff(
    const scalar mu
) const
{

    EigenMatrix<scalar> em(laplacian_, true);
    DiagonalMatrix<scalar> lambda = em.EValsRe();
    for (label i = 0; i < numComponents_; ++i) lambda[i] = 0;
    DiagonalMatrix<scalar> sqrtLambda(d_, 0.0);
    for (label i = numComponents_; i < d_; ++i) {
        if (lambda[i] < lambda[i-1]) {
            WarningInFunction << "Eigenvalues not sorted. Check the EigenMatrix implementation." << endl;
        }
        sqrtLambda[i] = sqrt(lambda[i]);
    }

    label omega;
    scalar sumEVals = 0.0;
    scalar sumSqrtEVals = 0.0;
    for (omega = numComponents_; omega < d_; ++omega) {
        sumEVals += lambda[omega];
        sumSqrtEVals += sqrtLambda[omega];
        if (sqrtLambda[omega] * (1.0 + mu * sumEVals) < mu * lambda[omega] * sumSqrtEVals) {
            sumEVals -= lambda[omega];
            sumSqrtEVals -= sqrtLambda[omega];
            break;
        }
    }

    if (omega == numComponents_) return scalar(numComponents_);
    const auto eigenvalues = lambda.slice(numComponents_, omega - numComponents_);
    const auto sqrtEigenvalues = sqrtLambda.slice(numComponents_, omega - numComponents_);
    const scalarField p
    (
        sqrtEigenvalues*(1.0 + mu*sumEVals)/sumSqrtEVals - mu*eigenvalues
    );
    return scalar(numComponents_) + sum(p/(mu*eigenvalues + p));

}

scalarField Foam::graphLaplacian::DOptimalDesign(const scalar mu) const
{
    return optimalDesign(mu, false);
}


scalar Foam::graphLaplacian::dTr(const scalar mu) const
{
    scalar dimension;
    optimalDesign(mu, true, &dimension);
    return dimension;
}


scalarField Foam::graphLaplacian::optimalDesign
(
    const scalar mu,
    const bool trace,
    scalar* traceObjective
) const
{

    const label maxFWIter = 100000;
    const scalar tolerance = 1e-8;
    scalarField probs(d_, 1.0/scalar(d_));
    if (traceObjective) *traceObjective = 0;
    if (mu == 0 || d_ == 1)
    {
        if (traceObjective) *traceObjective = scalar(d_)/(1 + scalar(d_)*SMALL);
        return probs;
    }

    SquareMatrix<scalar> inverse(d_, 0.0);
    SquareMatrix<scalar> metric(trace ? d_ : 0, 0.0);
    scalarField gradient(d_), column(d_), unit(d_, 0.0);
    scalarField hi(d_), hj(d_), gi(d_), gj(d_), u(d_), v(d_);
    scalar gap = GREAT;
    bool refresh = true;

    for (label iteration = 0; iteration <= maxFWIter; ++iteration)
    {
        const bool fresh = refresh || iteration == maxFWIter
            || iteration % max(label(10), d_) == 0;
        if (fresh)
        {
            const LLTMatrix<scalar> chol = cholLapPlusProb(probs, mu);
            // Every return below uses a fresh inverse. Reuse its diagonal
            // to obtain the trace objective without any additional solves.
            if (traceObjective) *traceObjective = 0;
            for (label j = 0; j < d_; ++j)
            {
                unit[j] = 1;
                chol.solve(column, unit);
                for (label i = 0; i < d_; ++i) inverse(i, j) = column[i];
                if (traceObjective) *traceObjective += probs[j]*column[j];
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
                            value -= mu*laplacian_(i, k)
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
        for (label i = 0; i < d_; ++i)
        {
            gradient[i] = trace ? metric(i, i) : inverse(i, i);
            if (gradient[i] > gradient[receiver]) receiver = i;
            if (probs[i] > 0 && (donor < 0 || gradient[i] < gradient[donor]))
            {
                donor = i;
            }
        }
        const scalar meanGradient = sumProd(probs, gradient);
        // The simplex Frank-Wolfe gap bounds the remaining objective error.
        gap = gradient[receiver] - meanGradient;
        if (gap <= tolerance*max(scalar(1), mag(meanGradient)))
        {
            if (fresh) return probs;
            refresh = true;
            continue;
        }
        if (iteration == maxFWIter) break;

        // Vertex exchange: move mass from the weakest supported arm.
        // D step: Harman et al., arXiv:1801.05661, Appendix A.1.
        const scalar a = inverse(receiver, receiver);
        const scalar b = inverse(donor, donor);
        const scalar c = inverse(receiver, donor);
        const scalar slope = a - b;
        const scalar curvature = max(scalar(0), a*b - c*c);
        scalar step = probs[donor];
        if (trace)
        {
            const scalar first = gradient[receiver] - gradient[donor];
            const scalar second = b*gradient[receiver] + a*gradient[donor]
                - 2*c*metric(receiver, donor);
            // Along this exchange the objective gain is
            // (first*t - second*t^2)/(1 + slope*t - curvature*t^2).
            scalar lo = 0, hiStep = step;
            for (label k = 0; k < 50; ++k)
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
            if (trace)
            {
                gi[k] = metric(k, receiver);
                gj[k] = metric(k, donor);
            }
        }
        u = kii*hi + kij*hj;
        v = kij*hi + kjj*hj;
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
        probs[receiver] += step;
        probs[donor] -= step;
    }

    WarningInFunction << (trace ? "Trace" : "D-optimal")
        << " design did not converge; simplex gap = " << gap << endl;
    return probs;
}

void Foam::graphLaplacian::apply
(
    const scalarField& x,
    scalarField& y,
    const label component
) const
{
    const label begin = component < 0 ? 0 : componentOffsets_[component];
    const label end = component < 0 ? d_ : componentOffsets_[component + 1];
    if (component >= 0) y = 0;
    for (label r = begin; r < end; ++r)
    {
        scalar value = 0;
        for (label k = rowOffsets_[r]; k < rowOffsets_[r + 1]; ++k)
        {
            value += values_[k]*x[columns_[k]];
        }
        y[originalRows_[r]] = value;
    }
}

LLTMatrix<scalar> Foam::graphLaplacian::cholLapPlusProb(
    const scalarField& probs,
    const scalar mu
) const
{
    SquareMatrix<scalar> LapPlusProb = mu * laplacian_;
    for (label i = 0; i < d_; ++i) {
        LapPlusProb[i][i] += probs[i] + SMALL;
    }
    return LLTMatrix<scalar>(LapPlusProb);
}

scalarField Foam::graphLaplacian::getHat(
    const scalarField& probs,
    const scalar mu,
    const label row
) const
{
    scalarField hat(d_, 0.0);
    if (mu == 0 || degree_[row] == 0)
    {
        hat[row] = 1.0 / (probs[row] + SMALL);
        return hat;
    }

    const label component = componentLabels_[row];
    const scalarField p(probs + SMALL);
    const scalarField invDiag(1.0 / (mu * degree_ + p));

    // Linear system matrix (A = mu*L + diag(p)) on the active component.
    auto A = [&](const scalarField& x, scalarField& out)
    {
        apply(x, out, component);
        out *= mu;
        out += p*x;
    };

    // Target vector (one-hot on the given global row).
    scalarField e(d_, 0.0);
    e[row] = 1.0;

    // Degree preconditioner on the component of the provided row:
    // M = D - mu*degree*degree^T/sum(degree), D = diag(mu*degree + p).
    // Computed via the Sherman-Morrison rank-one update
    scalarField quotient(d_, 0.0);
    scalar denom = 0.0;
    for (label r = componentOffsets_[component]; r < componentOffsets_[component+1]; ++r) {
        label i = originalRows_[r];
        quotient[i] = degree_[i] * invDiag[i];
        denom += quotient[i] * p[i];
    }
    const scalar scale = denom > 0 ? mu / denom : 0;
    auto precondition = [&](const scalarField& r, scalarField& out)
    {
        out = invDiag * r;
        if (scale > 0)
        {
            out += (scale * sumProd(quotient, r)) * quotient;
        }
    };

    // Unit RHS: the squared absolute and relative residuals coincide.
    if (conjugateGradient(A, e, precondition, hat)) return hat;

    // Direct fallback after the bounded iterative solve, on the full system.
    LLTMatrix<scalar> chol = cholLapPlusProb(probs, mu);
    chol.solve(hat, e);
    return hat;
}

Pair<scalarField> Foam::graphLaplacian::getHatAndBonus(
    const scalarField& probs,
    const scalar mu,
    const label row
) const
{

    LLTMatrix<scalar> chol = cholLapPlusProb(probs, mu);
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
