/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#ifndef FGAMGPreconditioner_H
#define FGAMGPreconditioner_H

#include "../FGAMG/FGAMGSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FGAMGPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class FGAMGPreconditioner
:
    public FGAMGSolver,
    public lduMatrix::preconditioner
{
protected:

    // Protected Data

        //- Number of V-cycles to perform
        label nVcycles_;


	// Create coarse grid correction fields
    	mutable PtrList<solveScalarField> coarseCorrFields;

	// Create coarse grid sources
    	mutable PtrList<solveScalarField> coarseSources;

    	// Create the smoothers for all levels
    	mutable PtrList<lduMatrix::smoother> smoothers;

    	// Scratch fields if processor-agglomerated coarse level meshes
    	// are bigger than original. Usually not needed
    	mutable solveScalarField ApsiScratch;
	mutable solveScalarField finestCorrectionScratch;

        //- Read the control parameters from the controlDict_
        virtual void readControls();

public:

    //- Runtime type information
    TypeName("FGAMG");


    // Constructors

        //- Construct from matrix components and preconditioner solver controls
        FGAMGPreconditioner
        (
            const lduMatrix::solver&,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~FGAMGPreconditioner() = default;


    // Member Functions

        //- Return wA the preconditioned form of residual rA
        virtual void precondition
        (
            solveScalarField& wA,
            const solveScalarField& rA,
            const direction cmpt=0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
