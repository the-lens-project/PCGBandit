/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "FGAMGSolver.H"
#include "SubField.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::FGAMGSolver::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{

    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    solveScalarField& psi = tpsi.ref();

    ConstPrecisionAdaptor<solveScalar, scalar> tsource(source);

    // Setup class containing solver performance data
    solverPerformance solverPerf(typeName, fieldName_);

    // Calculate A.psi used to calculate the initial residual
    solveScalarField Apsi(psi.size());
    matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // Create the storage for the finestCorrection which may be used as a
    // temporary in normFactor
    solveScalarField finestCorrection(psi.size());

    // Calculate normalisation factor
    solveScalar normFactor =
        this->normFactor(psi, tsource(), Apsi, finestCorrection);

    if ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate initial finest-grid residual field
    solveScalarField finestResidual(tsource() - Apsi);

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(finestResidual)(),
        fieldName_,
        true
    );

    // Calculate normalised residual for convergence test
    solverPerf.initialResidual() = gSumMag
    (
        finestResidual,
        matrix().mesh().comm()
    )/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();


    // Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
    )
    {
        // Create coarse grid correction fields
        PtrList<solveScalarField> coarseCorrFields;

        // Create coarse grid sources
        PtrList<solveScalarField> coarseSources;

        // Create the smoothers for all levels
        PtrList<lduMatrix::smoother> smoothers;

        // Scratch fields if processor-agglomerated coarse level meshes
        // are bigger than original. Usually not needed
        solveScalarField scratch1;
        solveScalarField scratch2;

        // Initialise the above data structures
        initVcycle
        (
            coarseCorrFields,
            coarseSources,
            smoothers,
            scratch1,
            scratch2
        );

        do
        {
            Vcycle
            (
                smoothers,
                psi,
                source,
                Apsi,
                finestCorrection,
                finestResidual,

                (scratch1.size() ? scratch1 : Apsi),
                (scratch2.size() ? scratch2 : finestCorrection),

                coarseCorrFields,
                coarseSources,
                cmpt
            );

            // Calculate finest level residual field
            matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);
            finestResidual = tsource();
            finestResidual -= Apsi;

            solverPerf.finalResidual() = gSumMag
            (
                finestResidual,
                matrix().mesh().comm()
            )/normFactor;

            if ((log_ >= 2) || (debug >= 2))
            {
                solverPerf.print(Info.masterStream(matrix().mesh().comm()));
            }
        } while
        (
            (
              ++solverPerf.nIterations() < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_, log_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(finestResidual)(),
        fieldName_,
        false
    );

    return solverPerf;
}


void Foam::FGAMGSolver::Vcycle
(
    const PtrList<lduMatrix::smoother>& smoothers,
    solveScalarField& psi,
    const scalarField& source,
    solveScalarField& Apsi,
    solveScalarField& finestCorrection,
    solveScalarField& finestResidual,

    solveScalarField& scratch1,
    solveScalarField& scratch2,

    PtrList<solveScalarField>& coarseCorrFields,
    PtrList<solveScalarField>& coarseSources,
    const direction cmpt
) const
{
    //debug = 2;

    const label coarsestLevel = matrixLevels_.size() - 1;

    // Restrict finest grid residual for the next level up.
    agglomeration_.restrictField(coarseSources[0], finestResidual, 0, true);

    if (nPreSweeps_ && ((log_ >= 2) || (debug >= 2)))
    {
        Pout<< "Pre-smoothing scaling factors: ";
    }


    // Residual restriction (going to coarser levels)
    for (label leveli = 0; leveli < coarsestLevel; leveli++)
    {
        if (coarseSources.set(leveli + 1))
        {
            // If the optional pre-smoothing sweeps are selected
            // smooth the coarse-grid field for the restricted source
            if (nPreSweeps_)
            {
                coarseCorrFields[leveli] = 0.0;

                smoothers[leveli + 1].scalarSmooth
                (
                    coarseCorrFields[leveli],
                    coarseSources[leveli],  //coarseSource,
                    cmpt,
                    min
                    (
                        nPreSweeps_ +  preSweepsLevelMultiplier_*leveli,
                        maxPreSweeps_
                    )
                );

                // Scale coarse-grid correction field
                // but not on the coarsest level because it evaluates to 1
                if (scaleCorrection_ && leveli < coarsestLevel - 1)
                {
                    solveScalarField::subField ACf
                    (
                        scratch1,
                        coarseCorrFields[leveli].size()
                    );

                    scale
                    (
                        coarseCorrFields[leveli],
                        const_cast<solveScalarField&>
                        (
                            static_cast<const solveScalarField&>(ACf)
                        ),
                        matrixLevels_[leveli],
                        interfaceLevelsBouCoeffs_[leveli],
                        interfaceLevels_[leveli],
                        coarseSources[leveli],
                        cmpt
                    );
                }

                // Correct the residual with the new solution
                // residual can be used by fusing Amul with b-Amul
                matrixLevels_[leveli].residual
                (
                    coarseSources[leveli],
                    coarseCorrFields[leveli],
                    ConstPrecisionAdaptor<scalar, solveScalar>
                    (
                        coarseSources[leveli]
                    )(),
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    cmpt
                );
            }

            // Residual is equal to source
            agglomeration_.restrictField
            (
                coarseSources[leveli + 1],
                coarseSources[leveli],
                leveli + 1,
                true
            );
        }
    }

    if (nPreSweeps_ && ((log_ >= 2) || (debug >= 2)))
    {
        Pout<< endl;
    }


    // Solve Coarsest level with either an iterative or direct solver
    if (coarseCorrFields.set(coarsestLevel))
    {
        solveCoarsestLevel
        (
            coarseCorrFields[coarsestLevel],
            coarseSources[coarsestLevel]
        );
    }

    if ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "Post-smoothing scaling factors: ";
    }

    // Smoothing and prolongation of the coarse correction fields
    // (going to finer levels)

    solveScalarField dummyField(0);

    // Work storage for prolongation
    solveScalarField work;

    for (label leveli = coarsestLevel - 1; leveli >= 0; leveli--)
    {
        if (coarseCorrFields.set(leveli))
        {
            // Create a field for the pre-smoothed correction field
            // as a sub-field of the finestCorrection which is not
            // currently being used
            solveScalarField::subField preSmoothedCoarseCorrField
            (
                scratch2,
                coarseCorrFields[leveli].size()
            );

            // Only store the preSmoothedCoarseCorrField if pre-smoothing is
            // used
            if (nPreSweeps_)
            {
                preSmoothedCoarseCorrField = coarseCorrFields[leveli];
            }


            // Prolong correction to leveli
            const auto& cf = agglomeration_.prolongField
            (
                coarseCorrFields[leveli],   // current level
                work,
                (
                    coarseCorrFields.set(leveli + 1)
                  ? coarseCorrFields[leveli + 1]
                  : dummyField              // dummy value
                ),
                leveli + 1
            );


            // Create A.psi for this coarse level as a sub-field of Apsi
            solveScalarField::subField ACf
            (
                scratch1,
                coarseCorrFields[leveli].size()
            );
            auto& ACfRef = const_cast<solveScalarField&>
            (
                static_cast<const solveScalarField&>(ACf)
            );

            if (interpolateCorrection_)
            {
                // Normal operation : have both coarse level and fine
                // level. No processor agglomeration
                interpolate
                (
                    coarseCorrFields[leveli],
                    ACfRef,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    agglomeration_.restrictAddressing(leveli + 1),
                    cf,
                    cmpt
                );
            }

            // Scale coarse-grid correction field
            // but not on the coarsest level because it evaluates to 1
            if
            (
                scaleCorrection_
             && (interpolateCorrection_ || leveli < coarsestLevel - 1)
            )
            {
                scale
                (
                    coarseCorrFields[leveli],
                    ACfRef,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    coarseSources[leveli],
                    cmpt
                );
            }

            // Only add the preSmoothedCoarseCorrField if pre-smoothing is
            // used
            if (nPreSweeps_)
            {
                coarseCorrFields[leveli] += preSmoothedCoarseCorrField;
            }

            smoothers[leveli + 1].scalarSmooth
            (
                coarseCorrFields[leveli],
                coarseSources[leveli],  //coarseSource,
                cmpt,
                min
                (
                    nPostSweeps_ + postSweepsLevelMultiplier_*leveli,
                    maxPostSweeps_
                )
            );
        }
    }

    // Prolong the finest level correction
    agglomeration_.prolongField
    (
        finestCorrection,
        coarseCorrFields[0],
        0,
        true
    );

    if (interpolateCorrection_)
    {
        interpolate
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            agglomeration_.restrictAddressing(0),
            coarseCorrFields[0],
            cmpt
        );
    }

    if (scaleCorrection_)
    {
        // Scale the finest level correction
        scale
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            finestResidual,
            cmpt
        );
    }

    forAll(psi, i)
    {
        psi[i] += finestCorrection[i];
    }

    smoothers[0].smooth
    (
        psi,
        source,
        cmpt,
        nFinestSweeps_
    );
}


void Foam::FGAMGSolver::initVcycle
(
    PtrList<solveScalarField>& coarseCorrFields,
    PtrList<solveScalarField>& coarseSources,
    PtrList<lduMatrix::smoother>& smoothers,
    solveScalarField& scratch1,
    solveScalarField& scratch2
) const
{

    label maxSize = matrix_.diag().size();

    coarseCorrFields.setSize(matrixLevels_.size());
    coarseSources.setSize(matrixLevels_.size());
    smoothers.setSize(matrixLevels_.size() + 1);

    // Creates copy of the solver (original copy can't be modified directly)
    dictionary modifiedControlDict(controlDict_);

    // Finds range of parameterized smoothers specified in the solver
    List<word> smootherRange = findSmootherRange(modifiedControlDict);

    // Selects initial smoother for the finest level
    if (smootherRange.size() > 0) {
	selectSmootherPerLevel(0, matrixLevels_.size(), smootherRange, modifiedControlDict);
    }

    // Create the smoother for the finest level
    smoothers.set
    (
        0,
        lduMatrix::smoother::New
        (
            fieldName_,
            matrix_,
            interfaceBouCoeffs_,
            interfaceIntCoeffs_,
            interfaces_,
            modifiedControlDict
        )
    );

    forAll(matrixLevels_, leveli)
    {
        if (agglomeration_.nCells(leveli) >= 0)
        {
            label nCoarseCells = agglomeration_.nCells(leveli);

            coarseSources.set(leveli, new solveScalarField(nCoarseCells));
        }

        if (matrixLevels_.set(leveli))
        {
            const lduMatrix& mat = matrixLevels_[leveli];

            label nCoarseCells = mat.diag().size();

            maxSize = max(maxSize, nCoarseCells);

            coarseCorrFields.set(leveli, new solveScalarField(nCoarseCells)); 
		
	    if (smootherRange.size() > 0) {
		selectSmootherPerLevel(leveli+1, matrixLevels_.size(), smootherRange, modifiedControlDict);
	    }

            smoothers.set
            (
                leveli + 1,
                lduMatrix::smoother::New
                (
                    fieldName_,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevelsIntCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    modifiedControlDict
                )
            );
        }
    }

    if (maxSize > matrix_.diag().size())
    {
        // Allocate some scratch storage
        scratch1.setSize(maxSize);
        scratch2.setSize(maxSize);
    }
}

Foam::List<Foam::word> Foam::FGAMGSolver::findSmootherRange
(
    dictionary& controlDict_
) const
{
    
    // Retrieves finest & coarsest level drop tolerances if specified in solution file; otherwise defaults to n/A
    const word finestSmoother(controlDict_.lookup("smoother"));
    const word coarsestSmoother = controlDict_.getOrDefault("coarsestSmoother", finestSmoother);

    // Checks that two different drop tolerance bounds were specified
    if (coarsestSmoother != finestSmoother) {

	// Finds index bounds on list of smoother parameterizations
	label startingIndex = smootherList.find(finestSmoother);
	label endingIndex = smootherList.find(coarsestSmoother);
	bool reversedBounds = false;

	// Throws error if improper drop tolerance bounds were specified
	if (startingIndex == -1 || endingIndex == -1) {
        	FatalErrorInFunction << "Improper drop tolerance specified in solution file. \n"
		<< "Please retry with selections from (ICTC_m0p5, ICTC_m1, ICTC_m1p5, ICTC_m2, ICTC_m2p5, ICTC_m3, ICTC_m3p5, ICTC_m4)" 
		<< " for both finest and coarsest levels." << exit(FatalError);
	} 
	// Flips bounds if specified in unexpected order
	else if (startingIndex > endingIndex) {
		Swap(startingIndex, endingIndex);
		reversedBounds = true;
	}

	// Creates sub-list containing smoothers in the specified bounds
	SubList<word> availableSmoothers(smootherList, endingIndex-startingIndex+1, startingIndex);

	// Reverses the list in the case in which bounds were specified in unexpected order
	if (reversedBounds) {
		reverse(availableSmoothers);
	}

	return List<word>(availableSmoothers);
    }

    // Returns empty list if finest and coarsest drop tolerances not specified
    return List<word>();
}

void Foam::FGAMGSolver::selectSmootherPerLevel
(
    const label level,
    const label maxLevels,
    List<word> availableSmoothers,
    dictionary& controlDict_
) const
{
	// Calculate which smoother to select from list
	scalar levelRatio = scalar(level) / scalar(maxLevels);
	label index = floor(levelRatio*scalar(availableSmoothers.size()-1));
	//Info << "Smoother parameterization selected: " << availableSmoothers[index] << " " << endl;

	// Select smoother
	controlDict_.set("smoother", availableSmoothers[index]);
}


Foam::dictionary Foam::FGAMGSolver::PCGsolverDict
(
    const scalar tol,
    const scalar relTol
) const
{
    dictionary dict(IStringStream("solver PCG; preconditioner DIC;")());
    dict.add("tolerance", tol);
    dict.add("relTol", relTol);

    return dict;
}


Foam::dictionary Foam::FGAMGSolver::PBiCGStabSolverDict
(
    const scalar tol,
    const scalar relTol
) const
{
    dictionary dict(IStringStream("solver PBiCGStab; preconditioner DILU;")());
    dict.add("tolerance", tol);
    dict.add("relTol", relTol);

    return dict;
}


void Foam::FGAMGSolver::solveCoarsestLevel
(
    solveScalarField& coarsestCorrField,
    const solveScalarField& coarsestSource
) const
{
    const label coarsestLevel = matrixLevels_.size() - 1;

    const label coarseComm = matrixLevels_[coarsestLevel].mesh().comm();

    if (directSolveCoarsest_)
    {
        PrecisionAdaptor<scalar, solveScalar> tcorrField(coarsestCorrField);

        coarsestLUMatrixPtr_->solve
        (
            tcorrField.ref(),
            ConstPrecisionAdaptor<scalar, solveScalar>(coarsestSource)()
        );
    }
    //else if
    //(
    //    agglomeration_.processorAgglomerate()
    // && procMatrixLevels_.set(coarsestLevel)
    //)
    //{
    //    //const labelList& agglomProcIDs = agglomeration_.agglomProcIDs
    //    //(
    //    //    coarsestLevel
    //    //);
    //    //
    //    //scalarField allSource;
    //    //
    //    //globalIndex cellOffsets;
    //    //if (Pstream::myProcNo(coarseComm) == agglomProcIDs[0])
    //    //{
    //    //    cellOffsets.offsets() =
    //    //        agglomeration_.cellOffsets(coarsestLevel);
    //    //}
    //    //
    //    //cellOffsets.gather
    //    //(
    //    //    coarseComm,
    //    //    agglomProcIDs,
    //    //    coarsestSource,
    //    //    allSource
    //    //);
    //    //
    //    //scalarField allCorrField;
    //    //solverPerformance coarseSolverPerf;
    //
    //    label solveComm = agglomeration_.procCommunicator(coarsestLevel);
    //
    //    coarsestCorrField = 0;
    //    solverPerformance coarseSolverPerf;
    //
    //    if (Pstream::myProcNo(solveComm) != -1)
    //    {
    //        const lduMatrix& allMatrix = procMatrixLevels_[coarsestLevel];
    //
    //        {
    //            Pout<< "** Master:Solving on comm:" << solveComm
    //                << " with procs:" << UPstream::procID(solveComm) << endl;
    //
    //            if (allMatrix.asymmetric())
    //            {
    //                coarseSolverPerf = PBiCGStab
    //                (
    //                    "coarsestLevelCorr",
    //                    allMatrix,
    //                    procInterfaceLevelsBouCoeffs_[coarsestLevel],
    //                    procInterfaceLevelsIntCoeffs_[coarsestLevel],
    //                    procInterfaceLevels_[coarsestLevel],
    //                    PBiCGStabSolverDict(tolerance_, relTol_)
    //                ).solve
    //                (
    //                    coarsestCorrField,
    //                    coarsestSource
    //                );
    //            }
    //            else
    //            {
    //                coarseSolverPerf = PCG
    //                (
    //                    "coarsestLevelCorr",
    //                    allMatrix,
    //                    procInterfaceLevelsBouCoeffs_[coarsestLevel],
    //                    procInterfaceLevelsIntCoeffs_[coarsestLevel],
    //                    procInterfaceLevels_[coarsestLevel],
    //                    PCGsolverDict(tolerance_, relTol_)
    //                ).solve
    //                (
    //                    coarsestCorrField,
    //                    coarsestSource
    //                );
    //            }
    //        }
    //    }
    //
    //    Pout<< "done master solve." << endl;
    //
    //    //// Scatter to all processors
    //    //coarsestCorrField.setSize(coarsestSource.size());
    //    //cellOffsets.scatter
    //    //(
    //    //    coarseComm,
    //    //    agglomProcIDs,
    //    //    allCorrField,
    //    //    coarsestCorrField
    //    //);
    //
    //    if (debug >= 2)
    //    {
    //        coarseSolverPerf.print(Info.masterStream(coarseComm));
    //    }
    //
    //    Pout<< "procAgglom: coarsestSource   :" << coarsestSource << endl;
    //    Pout<< "procAgglom: coarsestCorrField:" << coarsestCorrField << endl;
    //}
    else
    {
        coarsestCorrField = 0;
        const solverPerformance coarseSolverPerf
        (
            coarsestSolverPtr_->scalarSolve
            (
                coarsestCorrField,
                coarsestSource
            )
        );

        if ((log_ >= 2) || debug)
        {
            coarseSolverPerf.print(Info.masterStream(coarseComm));
        }
    }
}


// ************************************************************************* //
