/*---------------------------------------------------------------------------*\
                   Class similarityMatrix Implementation
\*---------------------------------------------------------------------------*/

#include "similarityMatrix.H"
#include "char.H"

#include <cmath>

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Allow rounding differences between parsed suffixes and dictionary values.
inline bool sameValue(const scalar x, const scalar y)
{
    return mag(x - y) <= SMALL*max(mag(x), mag(y));
}


// Other numeric axes infer their metric from their values.
axisMetric metricForAxis(const word& axisName)
{
    if
    (
        axisName == "droptol"
     || axisName == "smootherDroptol"
     || axisName == "nCellsInCoarsestLevel"
    )
    {
        return axisMetric::logRatio;
    }

    if (axisName == "smootherOmega" || axisName == "decayRate")
    {
        return axisMetric::fractionalDiff;
    }

    return axisMetric::inferred;
}


struct smootherSplit
{
    word family;
    bool hasDroptol = false;
    bool hasOmega = false;
    scalar droptol = 1.0;
    scalar omega = 1.0;
};


// Decode "m3p5" as -3.5 and "p1p2" as 1.2.
bool parseSuffix(const word& suffix, scalar& value)
{
    if (suffix.size() < 2 || suffix.size() > 12) return false;

    const char sign = suffix[0];
    if (sign != 'm' && sign != 'p') return false;

    string number = suffix.substr(1);
    const string::size_type point = number.find('p');
    if (point == 0 || point == number.size() - 1) return false;

    for (string::size_type i = 0; i < number.size(); ++i)
    {
        if (i == point) number[i] = '.';
        else if (!Foam::isdigit(number[i])) return false;
    }

    if (!readScalar(number, value)) return false;
    if (sign == 'm') value = -value;
    return true;
}


// DIC is droptol=1; GaussSeidel is omega=1.
bool splitSmoother(const word& smoother, smootherSplit& split)
{
    const word::size_type underscore = smoother.rfind('_');
    const bool suffixed = (underscore != word::npos);
    const word stem = smoother.substr(0, underscore);

    scalar parsed = 0.0;
    const bool parsedSuffix =
        suffixed && parseSuffix(smoother.substr(underscore + 1), parsed);

    if
    (
        suffixed
     &&
        (
            !parsedSuffix
         ||
            (
                stem != "ICTC" && stem != "SOR"
             && stem != "ICTCGaussSeidel" && stem != "DICSOR"
            )
        )
    )
    {
        return false;
    }

    if (stem == "symGaussSeidel")
    {
        split.family = "symGS";
        return true;
    }

    if (stem == "ICTC" || stem == "DIC")
    {
        split.family = "ICTC";
        split.hasDroptol = true;
        split.droptol =
            (stem == "ICTC" && parsedSuffix) ? pow(scalar(10), parsed) : 1.0;
        return true;
    }

    if (stem == "SOR" || stem == "GaussSeidel")
    {
        split.family = "SOR";
        split.hasOmega = true;
        split.omega = (stem == "SOR" && parsedSuffix) ? parsed : 1.0;
        return true;
    }

    if
    (
        stem == "ICTCGaussSeidel"
     || stem == "DICGaussSeidel"
     || stem == "DICSOR"
    )
    {
        split.family = "ICTCSOR";
        split.hasDroptol = true;
        split.hasOmega = true;
        split.droptol =
            (stem == "ICTCGaussSeidel" && parsedSuffix)
          ? pow(scalar(10), parsed)
          : 1.0;
        split.omega = (stem == "DICSOR" && parsedSuffix) ? parsed : 1.0;
        return true;
    }

    return false;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * parameterAxis Members * * * * * * * * * * * * //

label parameterAxis::addCategory(const word& category)
{
    forAll(categories_, i)
    {
        if (categories_[i] == category) return i;
    }
    categories_.append(category);
    return categories_.size() - 1;
}


void parameterAxis::addValue(const scalar value)
{
    for (const scalar known : values_)
    {
        if (sameValue(known, value)) return;
    }
    values_.append(value);
}


void parameterAxis::finalise()
{
    Foam::sort(values_);

    if (metric_ != axisMetric::inferred) return;

    // Any fractional value selects the scaled metric for the whole axis.
    metric_ = axisMetric::integerDiff;
    for (const scalar value : values_)
    {
        if (mag(value - std::round(value)) > SMALL)
        {
            metric_ = axisMetric::fractionalDiff;
            return;
        }
    }
}


label parameterAxis::rank(const scalar value) const
{
    forAll(values_, i)
    {
        if (sameValue(values_[i], value)) return i;
    }
    return -1;
}


scalar parameterAxis::distance(const scalar x, const scalar y) const
{
    switch (metric_)
    {
        case axisMetric::logRatio:
            return mag(log10(max(x, VSMALL)) - log10(max(y, VSMALL)));

        case axisMetric::fractionalDiff:
            return 10.0*mag(x - y);

        case axisMetric::categorical:
            return
                sameValue(x, y)
              ? 0.0
              : 1.0/scalar(max(label(1), categories_.size()));

        case axisMetric::integerDiff:
        case axisMetric::inferred:
            break;
    }

    return mag(x - y);
}


scalar parameterAxis::worstCase() const
{
    if (values_.size() < 2) return 0.0;

    if (metric_ == axisMetric::categorical)
    {
        return 1.0/scalar(max(label(1), categories_.size()));
    }

    // Numeric values are sorted; the endpoints give the largest distance.
    return distance(values_[0], values_[values_.size() - 1]);
}


bool parameterAxis::neighbours(const scalar x, const scalar y) const
{
    if (metric_ == axisMetric::categorical) return true;

    const label rankX = rank(x);
    const label rankY = rank(y);
    return (rankX >= 0 && rankY >= 0 && mag(rankX - rankY) == 1);
}


// * * * * * * * * * * * * similarityMatrix Members * * * * * * * * * * * * * //

void similarityMatrix::setValue
(
    HashTable<scalar>& arm,
    const word& axisName,
    const scalar value
)
{
    if (!axes_.found(axisName))
    {
        axes_.insert(axisName, parameterAxis(metricForAxis(axisName)));
    }

    axes_[axisName].addValue(value);
    arm.set(axisName, value);
}


void similarityMatrix::setCategory
(
    HashTable<scalar>& arm,
    const word& axisName,
    const word& category
)
{
    if (!axes_.found(axisName))
    {
        axes_.insert(axisName, parameterAxis(axisMetric::categorical));
    }

    parameterAxis& axis = axes_[axisName];
    const scalar index = scalar(axis.addCategory(category));
    axis.addValue(index);
    arm.set(axisName, index);
}


void similarityMatrix::discoverAxes(const List<dictionary>& armDicts)
{
    arms_.setSize(armDicts.size());

    forAll(armDicts, armi)
    {
        const dictionary& armDict = armDicts[armi];
        HashTable<scalar>& arm = arms_[armi];

        for (const word& key : armDict.toc())
        {
            if (key == "cacheAgglomeration") continue;

            if (key == "preconditioner")
            {
                const word type = armDict.get<word>(key);

                if (type == "ICTC" || type == "DIC")
                {
                    setCategory(arm, key, "IC");
                    if (type == "DIC") setValue(arm, "droptol", 1.0);
                }
                else
                {
                    setCategory(arm, key, type == "FGAMG" ? "GAMG" : type);
                }
                continue;
            }

            if (key == "smoother")
            {
                const word smoother = armDict.get<word>(key);
                smootherSplit split;

                if (splitSmoother(smoother, split))
                {
                    setCategory(arm, "smootherFamily", split.family);
                    if (split.hasDroptol)
                    {
                        setValue(arm, "smootherDroptol", split.droptol);
                    }
                    if (split.hasOmega)
                    {
                        setValue(arm, "smootherOmega", split.omega);
                    }
                }
                else
                {
                    setCategory(arm, "smootherFamily", smoother);
                }
                continue;
            }

            if (key == "directSolveCoarsest")
            {
                setCategory(arm, key, armDict.get<bool>(key) ? "yes" : "no");
                continue;
            }

            if (armDict.lookup(key, keyType::LITERAL).front().isNumber())
            {
                setValue(arm, key, armDict.get<scalar>(key));
            }
            else
            {
                setCategory(arm, key, armDict.get<word>(key));
            }
        }
    }

    forAllIters(axes_, iter)
    {
        iter.val().finalise();
    }
}


void similarityMatrix::buildDist()
{
    const label numArms = arms_.size();

    for (label i = 0; i < numArms; ++i)
    {
        (*this)(i, i) = 1.0;

        for (label j = i + 1; j < numArms; ++j)
        {
            const HashTable<scalar>& armI = arms_[i];
            const HashTable<scalar>& armJ = arms_[j];

            scalar similarity = 1.0;

            forAllConstIters(axes_, iter)
            {
                const word& axisName = iter.key();
                const parameterAxis& axis = iter.val();

                const bool onI = armI.found(axisName);
                const bool onJ = armJ.found(axisName);

                if (!onI && !onJ) continue;

                const scalar d =
                    (onI && onJ)
                  ? axis.distance(armI[axisName], armJ[axisName])
                  : axis.worstCase();

                similarity /= (1.0 + d);
            }

            (*this)(i, j) = similarity;
            (*this)(j, i) = similarity;
        }
    }
}


void similarityMatrix::buildPath()
{
    const label numArms = arms_.size();

    for (label i = 0; i < numArms; ++i)
    {
        (*this)(i, i) = 1.0;

        for (label j = i + 1; j < numArms; ++j)
        {
            const HashTable<scalar>& armI = arms_[i];
            const HashTable<scalar>& armJ = arms_[j];

            if (armI.size() != armJ.size()) continue;

            word differing;
            label numDiffering = 0;
            bool sameAxes = true;

            forAllConstIters(armI, iter)
            {
                if (!armJ.found(iter.key()))
                {
                    sameAxes = false;
                    break;
                }
                // Disabled initialization has no meaningful probe count.
                if
                (
                    iter.key() == "numProbes"
                 && (iter.val() == 0 || armJ[iter.key()] == 0)
                )
                {
                    continue;
                }
                if (!sameValue(iter.val(), armJ[iter.key()]))
                {
                    differing = iter.key();
                    ++numDiffering;
                    if (numDiffering > 1) break;
                }
            }

            if (!sameAxes || numDiffering != 1) continue;

            if (axes_[differing].neighbours(armI[differing], armJ[differing]))
            {
                (*this)(i, j) = 1.0;
                (*this)(j, i) = 1.0;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

similarityMatrix::similarityMatrix
(
    const List<dictionary>& armDicts,
    const similarityMode mode
)
:
    SquareMatrix<scalar>(armDicts.size(), 0.0),
    mode_(mode)
{
    discoverAxes(armDicts);

    if (mode_ == similarityMode::path)
    {
        buildPath();
    }
    else
    {
        buildDist();
    }
}


similarityMatrix::similarityMatrix
(
    const List<dictionary>& armDicts,
    const dictionary& solverControls
)
:
    similarityMatrix(armDicts, readMode(solverControls))
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

similarityMode similarityMatrix::readMode(const dictionary& solverControls)
{
    const word mode = solverControls.getOrDefault<word>("similarity", "path");

    if (mode == "path") return similarityMode::path;
    if (mode == "dist") return similarityMode::dist;

    FatalErrorInFunction
        << "Unknown similarity option " << mode
        << "; expected dist or path" << exit(FatalError);

    return similarityMode::path;
}


word similarityMatrix::modeName(const similarityMode mode)
{
    return (mode == similarityMode::path) ? "path" : "dist";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
