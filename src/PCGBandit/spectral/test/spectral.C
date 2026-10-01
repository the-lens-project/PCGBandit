/*---------------------------------------------------------------------------*\
    Checks similarity rules, Laplacian solves and design optimality.
\*---------------------------------------------------------------------------*/

#include "conjugateGradient.H"
#include "decomposedLaplacian.H"
#include "similarityMatrix.H"

#include <cmath>

using namespace Foam;

namespace
{

label checks = 0;
label failures = 0;

void check(const bool ok, const char* message)
{
    ++checks;
    if (!ok)
    {
        ++failures;
        Info<< "FAIL: " << message << endl;
    }
}

void checkClose(const scalar actual, const scalar expected, const char* message)
{
    check
    (
        std::isfinite(actual)
     && mag(actual - expected) <= 1e-10*max(scalar(1), mag(expected)),
        message
    );
}


dictionary ICArm(const scalar droptol)
{
    dictionary arm;
    arm.set("preconditioner", word(droptol == 1 ? "DIC" : "ICTC"));
    if (droptol != 1) arm.set("droptol", droptol);
    return arm;
}


dictionary GAMGArm(const word& smoother, const label cells = 10, const label merge = 1)
{
    dictionary arm;
    arm.set("preconditioner", word("FGAMG"));
    arm.set("smoother", smoother);
    arm.set("nCellsInCoarsestLevel", cells);
    arm.set("mergeLevels", merge);
    return arm;
}


void checkSimilarityRules()
{
    const List<dictionary> drops({ICArm(1e-2), ICArm(1e-3), ICArm(1e-4)});
    dictionary controls;
    check(similarityMatrix(drops, controls).mode() == similarityMode::path, "default mode");
    controls.set("similarity", word("dist"));
    const similarityMatrix dist(drops, controls);
    check(dist.mode() == similarityMode::dist, "mode from solver controls");
    checkClose(dist(0, 1), 0.5, "one decade of droptol");
    checkClose(dist(0, 2), 1.0/3.0, "two decades of droptol");
    controls.set("similarity", word("path"));
    const similarityMatrix path(drops, controls);
    checkClose(path(0, 1), 1, "adjacent numeric values");
    checkClose(path(1, 2), 1, "second adjacent numeric pair");
    checkClose(path(0, 2), 0, "non-adjacent numeric values");

    const List<dictionary> tinyDrops({ICArm(1e-16), ICArm(1e-17)});
    checkClose(similarityMatrix(tinyDrops, similarityMode::dist)(0, 1), 0.5, "tiny positive droptols remain distinct");
    checkClose(similarityMatrix(tinyDrops, similarityMode::path)(0, 1), 1, "tiny droptols have distinct ranks");

    dictionary explicitDIC = ICArm(1e-2);
    explicitDIC.set("droptol", scalar(1));
    const List<dictionary> equivalent({ICArm(1), explicitDIC});
    checkClose(similarityMatrix(equivalent, similarityMode::dist)(0, 1), 1, "DIC is droptol one");
    checkClose(similarityMatrix(equivalent, similarityMode::path)(0, 1), 0, "identical points differ on no axes");
    checkClose
    (
        similarityMatrix
        (
            List<dictionary>({GAMGArm("GaussSeidel"), GAMGArm("SOR_p1p0")}),
            similarityMode::dist
        )(0, 1),
        1, "GaussSeidel is omega one"
    );

    dictionary gamg;
    gamg.set("preconditioner", word("FGAMG"));
    gamg.set("mergeLevels", label(1));
    const List<dictionary> mixed({ICArm(1e-4), ICArm(1e-2), gamg});
    const similarityMatrix mixedDist(mixed, similarityMode::dist);
    checkClose(mixedDist(0, 2), 1.0/4.5, "missing axis pays its observed span");
    checkClose(mixedDist(0, 2), mixedDist(1, 2), "missing-axis penalty is independent of the value");
    checkClose(similarityMatrix(mixed, similarityMode::path)(0, 2), 0, "path requires the same axes");

    const List<dictionary> smoothers
    ({
        GAMGArm("DICGaussSeidel"), GAMGArm("DICSOR_p1p2"),
        GAMGArm("ICTC_m3"), GAMGArm("symGaussSeidel")
    });
    const similarityMatrix smootherDist(smoothers, similarityMode::dist);
    check(smootherDist.axes().found("smootherFamily"), "smoother family axis");
    check(smootherDist.axes().found("smootherDroptol"), "smoother droptol axis");
    check(smootherDist.axes().found("smootherOmega"), "smoother omega axis");
    check(!smootherDist.axes().found("smoother"), "smoother name is decomposed");
    checkClose(smootherDist(0, 1), 1.0/3.0, "fractional omega distance");
    checkClose(smootherDist(0, 2), 1.0/16.0, "product of family, droptol and missing omega distances");
    const similarityMatrix smootherPath(smoothers, similarityMode::path);
    checkClose(smootherPath(0, 1), 1, "one differing smoother parameter");
    checkClose(smootherPath(0, 2), 0, "different smoother axes");

    dictionary first = GAMGArm("GaussSeidel", 10, 1);
    dictionary second = GAMGArm("GaussSeidel", 1000, 3);
    first.set("nPreSweeps", label(0));
    second.set("nPreSweeps", label(2));
    const List<dictionary> integers({first, second});
    checkClose(similarityMatrix(integers, similarityMode::dist)(0, 1), 1.0/27.0, "coarse cells and integer axes");
    checkClose(similarityMatrix(integers, similarityMode::path)(0, 1), 0, "path excludes multiple differing axes");
    first = ICArm(1);
    second = first;
    first.set("decayRate", scalar(0));
    second.set("decayRate", scalar(0.25));
    checkClose
    (
        similarityMatrix(List<dictionary>({first, second}), similarityMode::dist)(0, 1),
        1.0/3.5, "fractional decay distance"
    );
    first.set("category", word("a"));
    second = first;
    second.set("category", word("b"));
    dictionary third = first;
    third.set("category", word("c"));
    const List<dictionary> categories({first, second, third});
    checkClose(similarityMatrix(categories, similarityMode::dist)(0, 2), 0.75, "categorical distance counts categories");
    checkClose(similarityMatrix(categories, similarityMode::path)(0, 2), 1, "all distinct categories are neighbours");

    for (const bool mixedSpelling : {false, true})
    {
        first = GAMGArm("GaussSeidel");
        second = first;
        first.set("directSolveCoarsest", label(0));
        if (mixedSpelling) second.set("directSolveCoarsest", word("yes"));
        else second.set("directSolveCoarsest", label(1));
        checkClose
        (
            similarityMatrix(List<dictionary>({first, second}), similarityMode::dist)(0, 1),
            2.0/3.0, "boolean values have categorical distance"
        );
    }
    first.set("directSolveCoarsest", word("no"));
    checkClose
    (
        similarityMatrix(List<dictionary>({first, second}), similarityMode::dist)(0, 1),
        2.0/3.0, "word booleans have categorical distance"
    );
    first.set("directSolveCoarsest", label(1));
    checkClose
    (
        similarityMatrix(List<dictionary>({first, second}), similarityMode::dist)(0, 1),
        1, "equivalent boolean spellings agree"
    );
    const similarityMatrix custom
    (
        List<dictionary>({GAMGArm("GaussSeidel"), GAMGArm("SOR_custom")}),
        similarityMode::dist
    );
    check(custom.axes()["smootherFamily"].size() == 2, "unknown smoother suffix keeps its own category");
}


label connectedComponents(const similarityMatrix& matrix)
{
    List<bool> seen(matrix.n(), false);
    DynamicList<label> pending;
    label components = 0;
    for (label root = 0; root < matrix.n(); ++root)
    {
        if (seen[root]) continue;
        ++components;
        seen[root] = true;
        pending.append(root);
        while (pending.size())
        {
            const label i = pending.back();
            pending.resize(pending.size() - 1);
            for (label j = 0; j < matrix.n(); ++j)
            {
                if (!seen[j] && matrix(i, j) > 0)
                {
                    seen[j] = true;
                    pending.append(j);
                }
            }
        }
    }
    return components;
}


struct designMeasures
{
    scalar trace = 0;
    scalarField logDetGradient;
    scalarField traceGradient;
};

// Evaluate both gradients independently in the original (arm) basis.
designMeasures measure
(
    const decomposedLaplacian& lap,
    const scalarField& p,
    const scalar mu
)
{
    designMeasures result;
    result.logDetGradient.setSize(p.size());
    result.traceGradient.setSize(p.size());
    const LLTMatrix<scalar> chol = lap.cholLapPlusProb(p, mu);
    scalarField unit(p.size(), Zero), column(p.size(), Zero);
    forAll(p, i)
    {
        unit[i] = 1;
        chol.solve(column, unit);
        unit[i] = 0;
        result.logDetGradient[i] = column[i];
        result.trace += p[i]*column[i];
        result.traceGradient[i] = column[i] - sum(p*sqr(column));
    }
    return result;
}


scalar relativeGap(const scalarField& p, const scalarField& gradient)
{
    const scalar average = sum(p*gradient);
    return (max(gradient) - average)/max(scalar(1), mag(average));
}


void checkDistribution(const scalarField& p)
{
    bool valid = true;
    for (const scalar probability : p)
    {
        valid = valid && std::isfinite(probability) && probability >= 0;
    }
    check(valid, "finite non-negative probabilities");
    checkClose(sum(p), 1, "probabilities sum to one");
}


void checkConjugateGradient()
{
    SquareMatrix<scalar> matrix(3, Zero);
    matrix(0, 0) = 6;
    matrix(1, 1) = 4;
    matrix(2, 2) = 3;
    matrix(0, 1) = matrix(1, 0) = -1;
    matrix(0, 2) = matrix(2, 0) = 1;
    matrix(1, 2) = matrix(2, 1) = 0.5;
    const scalarField b({2, -3, 0.75});
    const scalar tolerance = 1e-12;
    const scalar toleranceSqr = tolerance*tolerance;
    const auto multiply = [&matrix](const scalarField& input, scalarField& output)
    {
        forAll(input, i)
        {
            output[i] = 0;
            forAll(input, j) output[i] += matrix(i, j)*input[j];
        }
    };
    const auto jacobi = [&matrix](const scalarField& input, scalarField& output)
    {
        forAll(input, i) output[i] = input[i]/matrix(i, i);
    };
    const auto residualSqr = [&matrix, &b](const scalarField& x)
    {
        scalar result = 0;
        forAll(b, i)
        {
            scalar residual = -b[i];
            forAll(b, j) residual += matrix(i, j)*x[j];
            result += sqr(residual);
        }
        return result;
    };
    const LLTMatrix<scalar> chol(matrix);
    scalarField expected(b.size(), Zero);
    chol.solve(expected, b);

    scalarField x(b.size(), scalar(17));
    check
    (
        conjugateGradient(multiply, b, jacobi, x, tolerance, b.size()),
        "generic CG solves a 3x3 SPD system within three iterations"
    );
    check(max(mag(x - expected)) <= 1e-12, "generic CG agrees with Cholesky");
    check(residualSqr(x) <= toleranceSqr, "generic CG meets an independently evaluated residual tolerance");

    const auto exactInverse = [&chol](const scalarField& input, scalarField& output)
    {
        chol.solve(output, input);
    };
    x = 17;
    check
    (
        conjugateGradient(multiply, b, exactInverse, x, tolerance, 1),
        "exact inverse preconditioning converges in one iteration"
    );
    check(residualSqr(x) <= toleranceSqr, "one-iteration exact preconditioning meets the residual tolerance");

    x = 17;
    check
    (
        conjugateGradient(multiply, scalarField(b.size(), Zero), jacobi, x, tolerance),
        "generic CG accepts a zero right-hand side"
    );
    check(sum(sqr(x)) == 0, "generic CG clears the supplied solution for a zero right-hand side");
    x = 17;
    check
    (
        !conjugateGradient(multiply, b, jacobi, x, tolerance, 0),
        "generic CG reports nonconvergence with no iterations for a nonzero right-hand side"
    );
    check(sum(sqr(x)) == 0, "generic CG starts from zero even when no iterations are allowed");
    check
    (
        !conjugateGradient(multiply, b, jacobi, x, tolerance, 1),
        "generic CG respects its iteration limit"
    );

    const auto identity = [](const scalarField& input, scalarField& output)
    {
        output = input;
    };
    // This solve crosses the residual refresh at iteration 25. Restarting
    // there loses conjugacy and prevents convergence within n iterations.
    const label n = 40;
    scalarField longRhs(n, Zero), longExpected(n), product(n);
    longRhs[0] = 1;
    forAll(longExpected, i) longExpected[i] = scalar(n - i)/scalar(n + 1);
    const auto tridiagonal = [](const scalarField& input, scalarField& output)
    {
        forAll(input, i)
        {
            output[i] = 2*input[i];
            if (i > 0) output[i] -= input[i - 1];
            if (i + 1 < input.size()) output[i] -= input[i + 1];
        }
    };
    check
    (
        conjugateGradient(tridiagonal, longRhs, identity, x, tolerance, n),
        "CG converges across a residual refresh within the system dimension"
    );
    check(max(mag(x - longExpected)) <= 1e-11, "CG agrees with the exact tridiagonal solution");
    tridiagonal(x, product);
    check(sumSqr(longRhs - product) <= toleranceSqr, "CG meets the true residual tolerance after a refresh");
}


void checkDisconnectedSolves()
{
    SquareMatrix<scalar> weights(5, Zero);
    weights(0, 3) = weights(3, 0) = 2.5;
    weights(1, 4) = weights(4, 1) = 0.75;
    const decomposedLaplacian lap(weights);
    // Root discovery gives two interleaved pairs followed by the isolated node.
    const labelList componentLabels({0, 1, 2, 0, 1});
    const scalarField x({2, -3, 5, 0.5, 4});
    scalarField expected(x.size(), Zero), actual(x.size(), scalar(99));
    forAll(x, i)
    {
        forAll(x, j) expected[i] += weights(i, j)*(x[i] - x[j]);
    }
    lap.apply(x, actual);
    check(max(mag(actual - expected)) <= 1e-12, "full Laplacian apply matches the dense weighted-graph formula");
    scalarField combined(x.size(), Zero);
    actual = 99;
    for (label component = 0; component < 3; ++component)
    {
        scalarField componentExpected(x.size(), Zero);
        forAll(x, i)
        {
            if (componentLabels[i] == component) componentExpected[i] = expected[i];
        }
        // Reuse the previous output to expose stale entries outside this block.
        lap.apply(x, actual, component);
        check
        (
            max(mag(actual - componentExpected)) <= 1e-12,
            "component apply matches the dense formula at global row indices"
        );
        bool outsideIsZero = true;
        forAll(x, i)
        {
            if (componentLabels[i] != component)
            {
                outsideIsZero = outsideIsZero && actual[i] == 0;
            }
        }
        check(outsideIsZero, "component apply clears every row outside the selected component");
        combined += actual;
    }
    check(max(mag(combined - expected)) <= 1e-12, "component matvecs sum to the full Laplacian matvec");

    const scalarField p({0.08, 0.17, 0.3, 0.2, 0.25});
    for (const scalar mu : {0.0, 0.2})
    {
        const LLTMatrix<scalar> chol = lap.cholLapPlusProb(p, mu);
        forAll(p, row)
        {
            scalarField unit(p.size(), Zero), direct(p.size(), Zero);
            unit[row] = 1;
            chol.solve(direct, unit);
            const scalarField hat = lap.getHat(p, mu, row);
            check
            (
                max(mag(hat - direct)) <= 1e-10,
                "component getHat agrees with the full Cholesky solution"
            );
            scalar residualSqr = 0;
            bool outsideIsZero = true;
            forAll(p, i)
            {
                scalar residual = (p[i] + SMALL)*hat[i] - unit[i];
                forAll(p, j)
                {
                    residual += mu*weights(i, j)*(hat[i] - hat[j]);
                }
                residualSqr += sqr(residual);
                if (componentLabels[i] != componentLabels[row])
                {
                    outsideIsZero = outsideIsZero && hat[i] == 0;
                }
            }
            check(residualSqr <= 1e-12, "component getHat satisfies the full-system residual tolerance");
            check(outsideIsZero, "component getHat is zero outside the selected component");
        }
    }
}


void checkSolves
(
    const similarityMatrix& matrix,
    const decomposedLaplacian& lap,
    const scalarField& p,
    const scalar mu
)
{
    const label row = p.size()/2;
    const scalarField hat = lap.getHat(p, mu, row);
    const Pair<scalarField> both = lap.getHatAndBonus(p, mu, row);
    check(min(both[1]) > 0, "positive inverse diagonal");
    scalar residualSqr = 0, directResidualSqr = 0;
    forAll(p, i)
    {
        scalar value = (p[i] + SMALL)*hat[i] - (i == row ? 1 : 0);
        scalar direct = (p[i] + SMALL)*both[0][i] - (i == row ? 1 : 0);
        forAll(p, j)
        {
            value += mu*matrix(i, j)*(hat[i] - hat[j]);
            direct += mu*matrix(i, j)*(both[0][i] - both[0][j]);
        }
        residualSqr += sqr(value);
        directResidualSqr += sqr(direct);
    }
    // CG uses a unit RHS and a squared L2 residual tolerance of 1e-12.
    // Check each solve's accuracy instead of requiring identical solutions.
    check(residualSqr <= 1e-12*(1 + 1e-8), "getHat meets its squared residual tolerance");
    check(directResidualSqr <= 1e-20, "getHatAndBonus Cholesky residual");
}


// Compare the maximum value against an independently factored distribution
// whose optimality is known by symmetry or from a small reference problem.
scalar checkTraceDimension
(
    const decomposedLaplacian& lap,
    const scalarField& reference,
    const scalar mu
)
{
    const designMeasures expected = measure(lap, reference, mu);
    check
    (
        relativeGap(reference, expected.traceGradient) <= 1e-6,
        "reference distribution has a small trace optimality gap"
    );
    const scalar dimension = lap.dTr(mu);
    checkClose
    (
        dimension, expected.trace,
        "dTr agrees with an independently factored reference optimum"
    );
    return dimension;
}

void checkSymmetricDesigns()
{
    SquareMatrix<scalar> weights(4, scalar(1));
    const decomposedLaplacian lap(weights);
    for (const scalar mu : {0.0, 0.001, 1.0})
    {
        check(max(mag(lap.DOptimalDesign(mu) - scalar(0.25))) < 1e-8, "symmetric D-optimal design is uniform");
        const scalar dimension = checkTraceDimension(lap, scalarField(4, 0.25), mu);
        // Complete graph eigenvalues are 0, 4, 4, 4; uniform Pi is optimal.
        const scalar expected = 1/(1 + 4*SMALL) + 3/(1 + 16*mu + 4*SMALL);
        checkClose(dimension, expected, "complete-graph trace agrees with its spectrum");
    }
    const decomposedLaplacian isolated(SquareMatrix<scalar>(4, Zero));
    checkClose(isolated.dEff(1), 4, "isolated-arm effective dimension has no positive eigenmodes");
    check(max(mag(isolated.DOptimalDesign(1) - scalar(0.25))) < 1e-8, "isolated-arm D-optimal design is uniform");
    checkClose
    (
        checkTraceDimension(isolated, scalarField(4, 0.25), 1),
        4/(1 + 4*SMALL), "isolated-arm trace includes SMALL regularization"
    );

    const decomposedLaplacian single(SquareMatrix<scalar>(1, scalar(1)));
    checkClose(single.dEff(1), 1, "single-arm effective dimension has no positive eigenmodes");
    checkClose
    (
        checkTraceDimension(single, scalarField(1, 1.0), 1),
        1/(1 + SMALL), "single-arm trace includes SMALL regularization"
    );

    SquareMatrix<scalar> path(3, Zero);
    path(0, 1) = path(1, 0) = path(1, 2) = path(2, 1) = 1;
    const decomposedLaplacian pathLap(path);
    const scalarField dDesign = pathLap.DOptimalDesign(0.1);
    check(mag(dDesign[1] - 0.28117800) < 1e-6, "three-node D reference optimum");
    const scalarField traceReference({0.31706111, 0.36587778, 0.31706111});
    const scalar pathDimension = checkTraceDimension(pathLap, traceReference, 0.1);
    check
    (
        pathDimension > measure(pathLap, scalarField(3, 1.0/3.0), 0.1).trace + 1e-4,
        "nonuniform trace optimum improves on uniform probabilities"
    );

    const scalarField p({0.1, 0.2, 0.3, 0.4});
    const scalar mu = 0.1, epsilon = 1e-6;
    scalarField plus(p), minus(p);
    plus[0] += epsilon;
    plus[1] -= epsilon;
    minus[0] -= epsilon;
    minus[1] += epsilon;
    const designMeasures atP = measure(lap, p, mu);
    const scalar derivative =
        (measure(lap, plus, mu).trace - measure(lap, minus, mu).trace)/(2*epsilon);
    check
    (
        mag(derivative - atP.traceGradient[0] + atP.traceGradient[1]) < 1e-8,
        "trace gradient agrees with a simplex finite difference"
    );
}


void runSpace(const char* name, const List<dictionary>& arms)
{
    for (const similarityMode mode : {similarityMode::path, similarityMode::dist})
    {
        Info<< name << "/" << similarityMatrix::modeName(mode) << ": " << arms.size() << " arms";
        const similarityMatrix matrix(arms, mode);
        bool valid = matrix.n() == arms.size();
        for (label i = 0; i < matrix.n(); ++i)
        {
            valid = valid && matrix(i, i) == 1;
            for (label j = 0; j < matrix.n(); ++j)
            {
                const scalar value = matrix(i, j);
                valid = valid && std::isfinite(value) && value == matrix(j, i)
                    && value >= 0 && value <= 1;
                valid = valid && (mode == similarityMode::path ? value == 0 || value == 1 : value > 0);
            }
        }
        check(valid, "similarity matrix invariants");
        const decomposedLaplacian lap(matrix);
        const DiagonalMatrix<scalar> eigenvalues = lap.EVals();
        label zeros = 0;
        bool sorted = true;
        forAll(eigenvalues, i)
        {
            if (eigenvalues[i] == 0) ++zeros;
            sorted = sorted && eigenvalues[i] >= -1e-9;
            if (i) sorted = sorted && eigenvalues[i] >= eigenvalues[i-1] - 1e-9;
        }
        check(sorted, "non-negative sorted Laplacian spectrum");
        check(zeros == connectedComponents(matrix), "nullity equals graph component count");
        const scalar mu = 0.001;
        const scalar dimension = lap.dEff(mu);
        check(dimension >= 1 && dimension <= arms.size(), "effective dimension bounds");
        check(lap.dEff(100*mu) <= dimension + 1e-8, "effective dimension decreases with regularization");

        const scalarField uniform(arms.size(), 1.0/arms.size());
        const scalarField dOptimal = lap.DOptimalDesign(mu);
        const scalar traceDimension = lap.dTr(mu);
        checkDistribution(dOptimal);
        const designMeasures baseline = measure(lap, uniform, mu);
        const designMeasures dDesign = measure(lap, dOptimal, mu);
        const scalar dGap = relativeGap(dOptimal, dDesign.logDetGradient);
        check(dGap <= 1e-6, "D-optimality gap");
        check
        (
            std::isfinite(traceDimension) && traceDimension <= dimension + 1e-8,
            "dTr is finite and bounded by dEff"
        );
        check(traceDimension >= baseline.trace - 1e-8, "dTr is at least the uniform trace");
        check(traceDimension >= dDesign.trace - 1e-8, "dTr is at least the D-optimal design trace");
        checkSolves(matrix, lap, dOptimal, mu);
        checkSolves(matrix, lap, uniform, mu);
        Info<< ", D gap=" << dGap << ", dEff=" << dimension
            << ", trace=" << baseline.trace << " -> " << traceDimension << endl;
    }
}

List<dictionary> smallSpace()
{
    DynamicList<dictionary> arms;
    for (label i = 0; i < 8; ++i)
    {
        arms.append(ICArm(Foam::pow(10.0, -4.0 + 0.5*scalar(i))));
    }
    arms.append(ICArm(1.0));
    return List<dictionary>(arms);
}


List<dictionary> mediumSpace()
{
    DynamicList<dictionary> arms(smallSpace());

    const List<word> smoothers
        ({"GaussSeidel", "DIC", "DICGaussSeidel", "symGaussSeidel"});

    for (const word& smoother : smoothers)
    {
        for (const label nCells : {10, 100, 1000})
        {
            for (const label mergeLevels : {1, 2})
            {
                arms.append(GAMGArm(smoother, nCells, mergeLevels));
            }
        }
    }

    return List<dictionary>(arms);
}


List<dictionary> largeSpace()
{
    DynamicList<dictionary> arms(smallSpace());

    const List<word> smoothers
    ({
        "GaussSeidel", "DIC", "DICGaussSeidel", "symGaussSeidel",
        "SOR_p0p8", "SOR_p1p2", "DICSOR_p0p8", "DICSOR_p1p2",
        "ICTC_m4", "ICTC_m2", "ICTCGaussSeidel_m3"
    });

    for (const word& smoother : smoothers)
    {
        for (const label nCells : {10, 1000})
        {
            for (const label mergeLevels : {1, 2})
            {
                arms.append(GAMGArm(smoother, nCells, mergeLevels));
            }
        }
    }

    DynamicList<dictionary> expanded;
    for (const label probes : {4, 16})
    {
        for (const label history : {0, 8})
        {
            for (const dictionary& arm : arms)
            {
                dictionary copy = arm;
                copy.set("numProbes", probes);
                copy.set("lenHistory", history);
                copy.set("decayRate", scalar(0));
                expanded.append(copy);
            }
        }
    }
    return List<dictionary>(expanded);
}

} // End anonymous namespace


int main(int argc, char** argv)
{
    const bool designsOnly = argc == 2 && string(argv[1]) == "--design-only";
    checkConjugateGradient();
    checkDisconnectedSolves();
    checkSymmetricDesigns();
    if (!designsOnly)
    {
        checkSimilarityRules();
        runSpace("small", smallSpace());
        runSpace("medium", mediumSpace());
        runSpace("large", largeSpace());
    }
    Info<< checks - failures << "/" << checks << " checks passed" << endl;
    return failures ? 1 : 0;
}
