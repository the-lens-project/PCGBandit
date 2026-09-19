/*---------------------------------------------------------------------------*\
Class
    Foam::costEstimation

Description
    Deterministic cost modeling of linear system solves given only the
    preconditioner used and the numbeer of iterations. The goal is to give
    learners loss estimates that do not depend on wallclock and thus are fully
    reproducible (the 'deterministic yes' setting).

SourceFiles
    costEstimation.C
\*---------------------------------------------------------------------------*/

#ifndef costEstimation_H
#define costEstimation_H

#include "lduMatrix.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class costEstimation Declaration
\*---------------------------------------------------------------------------*/

class costEstimation
{
    const lduMatrix& matrix_;
    const label nGeometricD_;
    const bool galerkin_;

    scalar communicationCostEstimate
    (
        const dictionary& subDict,
        const word& preconditioner
    ) const;

    scalar setupCostEstimate
    (
        const dictionary& subDict,
        const word& preconditioner
    ) const;

    scalar subspaceCostEstimate(const dictionary& subDict) const;

public:

    costEstimation
    (
        const lduMatrix& matrix,
        const label nGeometricD,
        const bool galerkin
    )
    :
        matrix_(matrix),
        nGeometricD_(nGeometricD),
        galerkin_(galerkin)
    {}

    scalar perIterationCostEstimate
    (
        const dictionary& subDict,
        const word& preconditioner
    ) const;

    scalar totalCostEstimate
    (
        const dictionary& subDict,
        const label nIterations,
        const label maxIter,
        const label backstop
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
