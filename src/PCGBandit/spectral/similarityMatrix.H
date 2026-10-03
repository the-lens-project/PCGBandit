/*---------------------------------------------------------------------------*\
Class
    Foam::similarityMatrix

Description
    Discover parameter axes from the arm dictionaries. The solver's
    "similarity" entry selects:

      - "kernel": product of 1/(1 + distance) over axes. A missing axis incurs
        its largest observed distance.
      - "lattice": connect arms with identical axes that differ at one adjacent
        value. All distinct categorical values are adjacent. Ignore numProbes
        differences when either arm has numProbes=0 (disabled initialization).

    Smoothers have a family axis (ICTC, SOR, ICTCSOR, symGS) and separate
    droptol/omega axes. DIC is droptol=1; GaussSeidel is omega=1.

SourceFiles
    similarityMatrix.C
\*---------------------------------------------------------------------------*/

#ifndef similarityMatrix_H
#define similarityMatrix_H

#include "SquareMatrix.H"
#include "dictionary.H"
#include "HashTable.H"
#include "DynamicList.H"

namespace Foam
{

enum class axisMetric
{
    logRatio,        //!< d(x, y) = |log10(x/y)|, for drop tolerances and cell counts
    integerDiff,     //!< d(x, y) = |x - y|, for sweep and probe counts
    fractionalDiff,  //!< d(x, y) = 10|x - y|, for decay rates and relaxation factors
    categorical,     //!< d(x, y) = (x != y)/nCategories, for unordered names
    inferred         //!< numeric, resolved once every value on the axis is known
};

enum class similarityMode
{
    kernel, //!< product of 1/(1 + d) over the axes
    lattice //!< one iff the arms are a single step apart on one axis
};


/*---------------------------------------------------------------------------*\
                        Class parameterAxis Declaration
\*---------------------------------------------------------------------------*/

// Values are sorted and unique. Categorical values index categories_.
class parameterAxis
{
    // Private Data

        axisMetric metric_;

        DynamicList<scalar> values_;

        DynamicList<word> categories_;


public:

    // Constructors

        parameterAxis()
        :
            metric_(axisMetric::inferred)
        {}

        explicit parameterAxis(const axisMetric metric)
        :
            metric_(metric)
        {}


    // Member Functions

        axisMetric metric() const noexcept { return metric_; }

        const UList<scalar>& values() const noexcept { return values_; }

        const UList<word>& categories() const noexcept { return categories_; }

        label size() const noexcept { return values_.size(); }

        //- Return the category index, adding it if absent
        label addCategory(const word& category);

        void addValue(const scalar value);

        //- Sort values and resolve the metric after discovering all arms
        void finalise();

        //- Index in values(), or -1 if absent
        label rank(const scalar value) const;

        scalar distance(const scalar x, const scalar y) const;

        //- Largest observed distance
        scalar worstCase() const;

        //- Adjacent numeric values, or any categorical pair
        bool neighbours(const scalar x, const scalar y) const;
};


/*---------------------------------------------------------------------------*\
                       Class similarityMatrix Declaration
\*---------------------------------------------------------------------------*/

class similarityMatrix
:
    public SquareMatrix<scalar>
{
    // Private Data

        similarityMode mode_;

        HashTable<parameterAxis> axes_;

        //- Axis name to value for each arm
        List<HashTable<scalar>> arms_;


    // Private Member Functions

        void setValue
        (
            HashTable<scalar>& arm,
            const word& axisName,
            const scalar value
        );

        void setCategory
        (
            HashTable<scalar>& arm,
            const word& axisName,
            const word& category
        );

        void discoverAxes(const List<dictionary>& armDicts);

        void buildKernel(const scalar threshold);

        void buildLattice();


public:

    // Constructors

        similarityMatrix
        (
            const List<dictionary>& armDicts,
            const dictionary& solverControls,
            const scalar threshold = 0.1
        );

        similarityMatrix
        (
            const List<dictionary>& armDicts,
            const similarityMode mode,
            const scalar threshold = 0.1
        );


    // Member Functions

        similarityMode mode() const noexcept { return mode_; }

        const HashTable<parameterAxis>& axes() const noexcept { return axes_; }

        const List<HashTable<scalar>>& arms() const noexcept { return arms_; }


    // Static Member Functions

        //- Read "similarity", defaulting to "lattice"
        static similarityMode readMode(const dictionary& solverControls);

        static word modeName(const similarityMode mode);

};


} // End namespace Foam

#endif

// ************************************************************************* //
