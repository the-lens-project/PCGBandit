/*---------------------------------------------------------------------------*\
\*---------------------------------------------------------------------------*/

#include "FGAMGSolver.H"
#include "GAMGInterface.H"
#include "PCG.H"
#include "PBiCGStab.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FGAMGSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<FGAMGSolver>
        addFGAMGSolverMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<FGAMGSolver>
        addFGAMGAsymSolverMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FGAMGSolver::FGAMGSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    ),

    // Default values for all controls
    // which may be overridden by those in controlDict
    nPreSweeps_(0),
    preSweepsLevelMultiplier_(1),
    maxPreSweeps_(4),
    nPostSweeps_(2),
    postSweepsLevelMultiplier_(1),
    maxPostSweeps_(4),
    nFinestSweeps_(2),

    cacheAgglomeration_(true),
    interpolateCorrection_(false),
    scaleCorrection_(matrix.symmetric()),
    directSolveCoarsest_(false),

    agglomeration_(GAMGAgglomeration::New(matrix_, controlDict_)),

    matrixLevels_(agglomeration_.size()),
    primitiveInterfaceLevels_(agglomeration_.size()),
    interfaceLevels_(agglomeration_.size()),
    interfaceLevelsBouCoeffs_(agglomeration_.size()),
    interfaceLevelsIntCoeffs_(agglomeration_.size()),

    smootherList
    ({
        "ICTC_m0p5",
        "ICTC_m1",
        "ICTC_m1p5",
        "ICTC_m2",
        "ICTC_m2p5",
        "ICTC_m3",
        "ICTC_m3p5",
        "ICTC_m4"
    })
{
    readControls();

    if (agglomeration_.processorAgglomerate())
    {
        forAll(agglomeration_, fineLevelIndex)
        {
            if (agglomeration_.hasMeshLevel(fineLevelIndex))
            {
                if
                (
                    (fineLevelIndex+1) < agglomeration_.size()
                 && agglomeration_.hasProcMesh(fineLevelIndex+1)
                )
                {
                    // Construct matrix without referencing the coarse mesh so
                    // construct a dummy mesh instead. This will get overwritten
                    // by the call to procAgglomerateMatrix so is only to get
                    // it through agglomerateMatrix


                    const lduInterfacePtrsList& fineMeshInterfaces =
                        agglomeration_.interfaceLevel(fineLevelIndex);

                    PtrList<GAMGInterface> dummyPrimMeshInterfaces
                    (
                        fineMeshInterfaces.size()
                    );
                    lduInterfacePtrsList dummyMeshInterfaces
                    (
                        dummyPrimMeshInterfaces.size()
                    );

                    OCharStream os(IOstreamOption::BINARY);
                    ISpanStream is(IOstreamOption::BINARY);

                    forAll(fineMeshInterfaces, intI)
                    {
                        if (fineMeshInterfaces.set(intI))
                        {
                            os.rewind();

                            refCast<const GAMGInterface>
                            (
                                fineMeshInterfaces[intI]
                            ).write(os);

                            is.reset(os.view());

                            dummyPrimMeshInterfaces.set
                            (
                                intI,
                                GAMGInterface::New
                                (
                                    fineMeshInterfaces[intI].type(),
                                    intI,
                                    dummyMeshInterfaces,
                                    is
                                )
                            );
                        }
                    }

                    forAll(dummyPrimMeshInterfaces, intI)
                    {
                        if (dummyPrimMeshInterfaces.set(intI))
                        {
                            dummyMeshInterfaces.set
                            (
                                intI,
                                &dummyPrimMeshInterfaces[intI]
                            );
                        }
                    }

                    // So:
                    // - pass in incorrect mesh (= fine mesh instead of coarse)
                    // - pass in dummy interfaces
                    agglomerateMatrix
                    (
                        fineLevelIndex,
                        agglomeration_.meshLevel(fineLevelIndex),
                        dummyMeshInterfaces
                    );


                    const labelList& procAgglomMap =
                        agglomeration_.procAgglomMap(fineLevelIndex+1);
                    const List<label>& procIDs =
                        agglomeration_.agglomProcIDs(fineLevelIndex+1);

                    procAgglomerateMatrix
                    (
                        procAgglomMap,
                        procIDs,
                        fineLevelIndex
                    );
                }
                else
                {
                    agglomerateMatrix
                    (
                        fineLevelIndex,
                        agglomeration_.meshLevel(fineLevelIndex + 1),
                        agglomeration_.interfaceLevel(fineLevelIndex + 1)
                    );
                }
            }
            else
            {
                // No mesh. Not involved in calculation anymore
            }
        }
    }
    else
    {
        forAll(agglomeration_, fineLevelIndex)
        {
            // Agglomerate on to coarse level mesh
            agglomerateMatrix
            (
                fineLevelIndex,
                agglomeration_.meshLevel(fineLevelIndex + 1),
                agglomeration_.interfaceLevel(fineLevelIndex + 1)
            );
        }
    }

    if ((log_ >= 2) || (debug & 2))
    {
        for
        (
            label fineLevelIndex = 0;
            fineLevelIndex <= matrixLevels_.size();
            fineLevelIndex++
        )
        {
            if (fineLevelIndex == 0 || matrixLevels_.set(fineLevelIndex-1))
            {
                const lduMatrix& matrix = matrixLevel(fineLevelIndex);
                const lduInterfaceFieldPtrsList& interfaces =
                    interfaceLevel(fineLevelIndex);

                Pout<< "level:" << fineLevelIndex << nl
                    << "    nCells:" << matrix.diag().size() << nl
                    << "    nFaces:" << matrix.lower().size() << nl
                    << "    nInterfaces:" << interfaces.size()
                    << endl;

                forAll(interfaces, i)
                {
                    if (interfaces.set(i))
                    {
                        Pout<< "        " << i
                            << "\ttype:" << interfaces[i].type()
                            << "\tsize:"
                            << interfaces[i].interface().faceCells().size()
                            << endl;
                    }
                }
            }
            else
            {
                Pout<< "level:" << fineLevelIndex << " : no matrix" << endl;
            }
        }
        Pout<< endl;
    }


    if (matrixLevels_.size())
    {
        const label coarsestLevel = matrixLevels_.size() - 1;

        if (matrixLevels_.set(coarsestLevel))
        {
            if (directSolveCoarsest_)
            {
                coarsestLUMatrixPtr_.reset
                (
                    new LUscalarMatrix
                    (
                        matrixLevels_[coarsestLevel],
                        interfaceLevelsBouCoeffs_[coarsestLevel],
                        interfaceLevels_[coarsestLevel]
                    )
                );
            }
            else
            {
                entry* coarseEntry = controlDict_.findEntry
                (
                    "coarsestLevelCorr",
                    keyType::LITERAL_RECURSIVE
                );
                if (coarseEntry && coarseEntry->isDict())
                {
                    coarsestSolverPtr_ = lduMatrix::solver::New
                    (
                        "coarsestLevelCorr",
                        matrixLevels_[coarsestLevel],
                        interfaceLevelsBouCoeffs_[coarsestLevel],
                        interfaceLevelsIntCoeffs_[coarsestLevel],
                        interfaceLevels_[coarsestLevel],
                        coarseEntry->dict()
                    );
                }
                else if (matrixLevels_[coarsestLevel].asymmetric())
                {
                    coarsestSolverPtr_.reset
                    (
                        new PBiCGStab
                        (
                            "coarsestLevelCorr",
                            matrixLevels_[coarsestLevel],
                            interfaceLevelsBouCoeffs_[coarsestLevel],
                            interfaceLevelsIntCoeffs_[coarsestLevel],
                            interfaceLevels_[coarsestLevel],
                            PBiCGStabSolverDict(tolerance_, relTol_)
                        )
                    );
                }
                else
                {
                    coarsestSolverPtr_.reset
                    (
                        new PCG
                        (
                            "coarsestLevelCorr",
                            matrixLevels_[coarsestLevel],
                            interfaceLevelsBouCoeffs_[coarsestLevel],
                            interfaceLevelsIntCoeffs_[coarsestLevel],
                            interfaceLevels_[coarsestLevel],
                            PCGsolverDict(tolerance_, relTol_)
                        )
                    );
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "No coarse levels created, either matrix too small for GAMG"
               " or nCellsInCoarsestLevel too large.\n"
               "    Either choose another solver of reduce "
               "nCellsInCoarsestLevel."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FGAMGSolver::~FGAMGSolver()
{
    if (!cacheAgglomeration_)
    {
        delete &agglomeration_;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::FGAMGSolver::readControls()
{
    lduMatrix::solver::readControls();

    controlDict_.readIfPresent("cacheAgglomeration", cacheAgglomeration_);
    controlDict_.readIfPresent("nPreSweeps", nPreSweeps_);
    controlDict_.readIfPresent
    (
        "preSweepsLevelMultiplier",
        preSweepsLevelMultiplier_
    );
    controlDict_.readIfPresent("maxPreSweeps", maxPreSweeps_);
    controlDict_.readIfPresent("nPostSweeps", nPostSweeps_);
    controlDict_.readIfPresent
    (
        "postSweepsLevelMultiplier",
        postSweepsLevelMultiplier_
    );
    controlDict_.readIfPresent("maxPostSweeps", maxPostSweeps_);
    controlDict_.readIfPresent("nFinestSweeps", nFinestSweeps_);
    controlDict_.readIfPresent("interpolateCorrection", interpolateCorrection_);
    controlDict_.readIfPresent("scaleCorrection", scaleCorrection_);
    controlDict_.readIfPresent("directSolveCoarsest", directSolveCoarsest_);

    if ((log_ >= 2) || debug)
    {
        Info<< "GAMGSolver settings :"
            << " cacheAgglomeration:" << cacheAgglomeration_
            << " nPreSweeps:" << nPreSweeps_
            << " preSweepsLevelMultiplier:" << preSweepsLevelMultiplier_
            << " maxPreSweeps:" << maxPreSweeps_
            << " nPostSweeps:" << nPostSweeps_
            << " postSweepsLevelMultiplier:" << postSweepsLevelMultiplier_
            << " maxPostSweeps:" << maxPostSweeps_
            << " nFinestSweeps:" << nFinestSweeps_
            << " interpolateCorrection:" << interpolateCorrection_
            << " scaleCorrection:" << scaleCorrection_
            << " directSolveCoarsest:" << directSolveCoarsest_
            << endl;
    }
}


const Foam::lduMatrix& Foam::FGAMGSolver::matrixLevel(const label i) const
{
    return i ? matrixLevels_[i-1] : matrix_;
}


const Foam::lduInterfaceFieldPtrsList& Foam::FGAMGSolver::interfaceLevel
(
    const label i
) const
{
    return i ? interfaceLevels_[i-1] : interfaces_;
}


const Foam::FieldField<Foam::Field, Foam::scalar>&
Foam::FGAMGSolver::interfaceBouCoeffsLevel
(
    const label i
) const
{
    return i ? interfaceLevelsBouCoeffs_[i-1] : interfaceBouCoeffs_;
}


const Foam::FieldField<Foam::Field, Foam::scalar>&
Foam::FGAMGSolver::interfaceIntCoeffsLevel
(
    const label i
) const
{
    return i ? interfaceLevelsIntCoeffs_[i-1] : interfaceIntCoeffs_;
}


// ************************************************************************* //
