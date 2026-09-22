/*---------------------------------------------------------------------------*\
    Checks similarity rules, Laplacian solves and design optimality.
\*---------------------------------------------------------------------------*/

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
    const LLTMatrix<scalar> chol = lap.cholLaplacianPlusPi(p, mu);
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
    check(max(mag(hat - both[0])) < 1e-12, "getHat agrees with getHatAndBonus");
    check(min(both[1]) > 0, "positive inverse diagonal");
    scalar residual = 0;
    forAll(p, i)
    {
        scalar value = (p[i] + SMALL)*hat[i] - (i == row ? 1 : 0);
        forAll(p, j) value += mu*matrix(i, j)*(hat[i] - hat[j]);
        residual = max(residual, mag(value));
    }
    check(residual < 1e-8, "regularized solve residual");
}


void checkSymmetricDesigns()
{
    SquareMatrix<scalar> weights(4, scalar(1));
    const decomposedLaplacian lap(weights);
    for (const scalar mu : {0.0, 0.001, 1.0})
    {
        check(max(mag(lap.DOptimalDesign(mu) - scalar(0.25))) < 1e-8, "symmetric D-optimal design is uniform");
        check(max(mag(lap.traceOptimalDesign(mu) - scalar(0.25))) < 1e-8, "symmetric trace-optimal design is uniform");
    }
    const decomposedLaplacian isolated(SquareMatrix<scalar>(4, Zero));
    check(max(mag(isolated.DOptimalDesign(1) - scalar(0.25))) < 1e-8, "isolated-arm D-optimal design is uniform");
    check(max(mag(isolated.traceOptimalDesign(1) - scalar(0.25))) < 1e-8, "isolated-arm trace design is uniform");

    SquareMatrix<scalar> path(3, Zero);
    path(0, 1) = path(1, 0) = path(1, 2) = path(2, 1) = 1;
    const decomposedLaplacian pathLap(path);
    const scalarField dDesign = pathLap.DOptimalDesign(0.1);
    const scalarField tDesign = pathLap.traceOptimalDesign(0.1);
    check(mag(dDesign[1] - 0.28117800) < 1e-6, "three-node D reference optimum");
    check(mag(tDesign[1] - 0.36587778) < 1e-6, "three-node trace reference optimum");

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
        const scalarField traceOptimal = lap.traceOptimalDesign(mu);
        checkDistribution(dOptimal);
        checkDistribution(traceOptimal);
        const designMeasures baseline = measure(lap, uniform, mu);
        const designMeasures dDesign = measure(lap, dOptimal, mu);
        const designMeasures traceDesign = measure(lap, traceOptimal, mu);
        const scalar dGap = relativeGap(dOptimal, dDesign.logDetGradient);
        const scalar traceGap = relativeGap(traceOptimal, traceDesign.traceGradient);
        check(dGap <= 1e-6, "D-optimality gap");
        check(traceGap <= 1e-6, "trace optimality gap");
        check(traceDesign.trace >= baseline.trace - 1e-8, "optimized trace is at least the uniform trace");
        checkSolves(matrix, lap, dOptimal, mu);
        checkSolves(matrix, lap, traceOptimal, mu);
        Info<< ", D gap=" << dGap << ", trace gap=" << traceGap
            << ", trace=" << baseline.trace << " -> " << traceDesign.trace << endl;
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


int main()
{
    checkSimilarityRules();
    checkSymmetricDesigns();
    runSpace("small", smallSpace());
    runSpace("medium", mediumSpace());
    runSpace("large", largeSpace());
    Info<< checks - failures << "/" << checks << " checks passed" << endl;
    return failures ? 1 : 0;
}
