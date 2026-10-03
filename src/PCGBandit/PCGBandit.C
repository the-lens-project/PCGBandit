/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

//
#include "PCGBandit.H"
#include "PrecisionAdaptor.H"
#include "referenceDictionary.H"
#include "graphLaplacian.H"

#include "clockValue.H"
#include "fvMesh.H"
#include "Pstream.H"
#include "Random.H"

#include "HashTable.H"

//#define PCGB_DEBUG
//#define DUMP_ABSOL

#ifdef DUMP_ABSOL
#include "OStringStream.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    clockValue PCGTime = clockValue();
    scalar PCGCost = 0.0;

    defineTypeNameAndDebug(PCGBandit, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PCGBandit>
        addPCGBanditSymMatrixConstructorToTable_;
    Random rndGen;

    dictionary preconditionerDict;
    dictionary subDict;
    HashTable<configurationSpace> configurationSpaces;
    HashTable<referenceDictionary> learningDicts;

    #ifdef DUMP_ABSOL
    #include "Absol/initializeDumping.H"
    #endif

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PCGBandit::PCGBandit
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
    costs_
    (
        matrix,
        dynamicCast<const fvMesh>(matrix.mesh()).nGeometricD(),
        solverControls.getOrDefault<word>("projection", "galerkin") == "galerkin"
    )
{

    // --- Contextual information specification
    word preconditioner = solverControls.get<word>("preconditioner");
    const fvMesh& mesh = dynamicCast<const fvMesh>(matrix.mesh());
    if (preconditioner == "separate") {
        banditName_ = mesh.name() + "." + fieldName;
    } else if (preconditioner == "joint") {
        banditName_ = "joint";
    } else {
        banditName_ = preconditioner;
    }
    if (relTol_ == 0.0 and Switch(solverControls.getOrDefault<word>("residualContext", "no"))) {
        banditName_ += "Final";
    }

    // --- Learning algorithm specification
    lossEstimator_ = solverControls.getOrDefault<word>("lossEstimator", "RV");
    deterministic_ = Switch(solverControls.getOrDefault<word>("deterministic", "no"));
    backstop_ = solverControls.getOrDefault<label>("backstop", -1);
    static_ = label(solverControls.getOrDefault<label>("static", -1));
    randomUniform_ = Switch(solverControls.getOrDefault<word>("randomUniform", "no"));
    banditAlgorithm_ = solverControls.getOrDefault<word>("banditAlgorithm", "TsallisINF");

    if (!configurationSpaces.found(banditName_))
    {
        if (Pstream::master(matrix.mesh().comm())) {
            rndGen.reset(mesh.time().controlDict().getOrDefault<label>("randomSeed", 0));
        }

        configurationSpaces.set
        (
            banditName_,
            configurationSpace(solverControls, matrix, static_)
        );
        #ifdef PCGB_DEBUG
        Info<< "Preconditioner configurations for " << banditName_ << " : " << configurationSpaces[banditName_].dicts() << endl;
        #endif

    }

    // --- The window ceilings are a property of the configuration space, so
    //     this waits until the space exists.  It is also why subspace_ is a
    //     pointer: banditName_, hence the space, is only known in the body.
    const configurationSpace& configs = configurationSpaces[banditName_];
    subspace_.reset
    (
        new subspaceInitializer
        (
            matrix,
            interfaceBouCoeffs,
            interfaces,
            fieldName,
            solverControls,
            configs.maxLenHistory(),
            configs.maxNumProbes(),
            configs.decayRates()
        )
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PCGBandit::queryLearner
(
    const scalar initialResidual
) const
{
    label i = static_;
    referenceDictionary& learningDict = learningDicts(banditName_);
    const List<dictionary>& preconditionerDicts = configurationSpaces[banditName_].dicts();
    const label comm = matrix_.mesh().comm();

    if (i == -1) {

        if (Pstream::master(comm)) {
            label d = preconditionerDicts.size();
            if (d == 1) {
                i = 0;
            } else if (randomUniform_) {
                i = floor(scalar(d) * rndGen.sample01<scalar>());
            } else if (banditAlgorithm_ == "ThompsonSampling") {
                #include "ThompsonSampling.H"
            } else if (banditAlgorithm_ == "SpectralINF") {
                #include "spectralINF.H"
            } else {
                #include "TsallisINF.H"
            }

            // --- Every learner acts on `loss` only when it is positive and
            //     none of them resets it, so without setting this to zero a
            //     query following a solve that wrote no loss would re-use
            //     the previous round's and attribute it to the wrong arm.
            learningDict.set<scalar>("loss", 0.0);
        }

        Pstream::broadcast(i, comm);

        #ifdef PCGB_DEBUG
        Info<< banditAlgorithm_ << " Selection: ";
    } else {

        Info<< "Static Preconditioner: ";
        #endif
    }

    subDict = preconditionerDicts[i];
    preconditionerDict.set("preconditioner", subDict);

    #ifdef PCGB_DEBUG
    Info<< subDict << endl;
    #endif
}

Foam::solverPerformance Foam::PCGBandit::scalarSolve
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt
) const
{

    #ifdef DUMP_ABSOL
    #include "Absol/startDump.H"
    #endif

    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );
    clockValue initializeTime;
    clockValue preconstructTime;
    clockValue iterationTime;
    clockValue learnerTime;
    clockValue solverTime = clockValue::now();

    bool armDrawn = false;
    label maxIter = maxIter_;
    label backstopIter = maxIter_;
    dictionary backstopDict;
    autoPtr<lduMatrix::preconditioner> preconPtr;

    label nCells = psi.size();
    solveScalar* __restrict__ psiPtr = psi.begin();

    solveScalarField pA(nCells);
    solveScalar* __restrict__ pAPtr = pA.begin();

    solveScalarField wA(nCells);
    solveScalar* __restrict__ wAPtr = wA.begin();

    solveScalar wArA = solverPerf.great_;
    solveScalar wArAold = wArA;

    initializeTime += clockValue::now();
    subspace_->update(psi, cmpt);
    initializeTime -= clockValue::now();

    // --- Calculate A.psi
    matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    solveScalarField rA(source - wA);
    solveScalar* __restrict__ rAPtr = rA.begin();

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        true
    );

    // --- Calculate normalisation factor
    solveScalar normFactor = this->normFactor(psi, source, wA, pA);

    if ((log_ >= 2) || (lduMatrix::debug >= 2))
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())
       /normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    for (label backstop = 0; backstop <= label(backstop_ != 0); backstop++) {

        // --- Set by the subspace correction alone converging the system
        bool converged = false;

        // --- Check convergence, solve if not converged
        if
        (
            minIter_ > 0
         || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
        )
        {

            // --- Select and construct the preconditioner
            if (backstop) {

                // --- Revert to backstopping preconditioner
                preconstructTime += preconstructTime.now();
                if (subDict.get<word>("preconditioner") != "DIC") {
                    preconPtr = lduMatrix::preconditioner::New(*this, backstopDict);
                }
                preconstructTime -= clockValue::now();
                iterationTime += clockValue::now();
            } else {
                // --- Get preconditioner from learning algorithm
                learnerTime = learnerTime.now();
                queryLearner(solverPerf.initialResidual());
                armDrawn = true;

                learnerTime -= clockValue::now();

                // --- Runs subspace initialization of the linear solve
                initializeTime += clockValue::now();
                subspace_->initialize
                (
                    cmpt,
                    subDict.getOrDefault<label>("lenHistory", 0),
                    subDict.getOrDefault<label>("numProbes", 4),
                    subDict.getOrDefault<scalar>("decayRate", 0),
                    normFactor, psi, rA, solverPerf
                );
                initializeTime -= clockValue::now();

                preconstructTime = preconstructTime.now();

                // --- The correction may have converged the system on its own
                converged = minIter_ <= 0 && solverPerf.checkConvergence(tolerance_, relTol_, log_);

                if (!converged) {
                    preconPtr = lduMatrix::preconditioner::New(*this, preconditionerDict);

                    // --- Default backstop iteration computed via a cost estimate ratio
                    if (backstop_ == -1) {
                        backstopIter = label(scalar(maxIter)
                                             * costs_.perIterationCostEstimate(subDict, "DIC")
                                             / costs_.perIterationCostEstimate(subDict, subDict.get<word>("preconditioner")));
                        maxIter = backstopIter;
                    }
                }
                preconstructTime -= clockValue::now();
                iterationTime = iterationTime.now();
            }

            // --- Solver iteration
            if (!converged)
            {
                do
                {

                    // --- Store previous wArA
                    wArAold = wArA;

                    // --- Precondition residual
                    preconPtr->precondition(wA, rA, cmpt);

                    // --- Update search directions:
                    wArA = gSumProd(wA, rA, matrix().mesh().comm());

                    if (solverPerf.nIterations() == 0 || solverPerf.nIterations() == backstopIter)
                    {
                        for (label cell=0; cell<nCells; cell++)
                        {
                            pAPtr[cell] = wAPtr[cell];
                        }
                    }
                    else
                    {
                        solveScalar beta = wArA/wArAold;

                        for (label cell=0; cell<nCells; cell++)
                        {
                            pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                        }
                    }

                    // --- Update preconditioned residual
                    matrix_.Amul(wA, pA, interfaceBouCoeffs_, interfaces_, cmpt);

                    solveScalar wApA = gSumProd(wA, pA, matrix().mesh().comm());

                    // --- Test for singularity
                    if (solverPerf.checkSingularity(mag(wApA)/normFactor)) break;

                    // --- Update solution and residual:

                    solveScalar alpha = wArA/wApA;

                    for (label cell=0; cell<nCells; cell++)
                    {
                        psiPtr[cell] += alpha*pAPtr[cell];
                        rAPtr[cell] -= alpha*wAPtr[cell];
                    }

                    solverPerf.finalResidual() =
                        gSumMag(rA, matrix().mesh().comm())
                       /normFactor;

                } while
                (
                    (
                      ++solverPerf.nIterations() < maxIter
                    && !solverPerf.checkConvergence(tolerance_, relTol_, log_)
                    )
                 || solverPerf.nIterations() < minIter_
                );
            }
            iterationTime -= clockValue::now();
        }

        // --- Exit if converged or if already tried backstopping
        if (backstop == 1 || solverPerf.checkConvergence(tolerance_, relTol_, log_))
        {
            matrix().setResidualField
            (
                ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
                fieldName_,
                false
            );
            break;
        }

        Info << "PCG backstopping at iteration " << backstopIter << endl;
        if (backstop_ == -1) {
            maxIter += maxIter_;
        } else {
            maxIter += backstop_;
        }
        if (subDict.get<word>("preconditioner") == "DIC") {
            backstopIter = maxIter;
        } else {
            backstopDict.set("preconditioner", "DIC");
        }
    }

    solverTime -= clockValue::now();
    learnerTime += clockValue::now();
    scalar costEstimate = 0.0;
    if (armDrawn) {

        // --- Compute solver cost
        if (deterministic_) {
            costEstimate = 1e-9 * costs_.totalCostEstimate(subDict, solverPerf.nIterations(), maxIter_, backstop_);
        } else {
            costEstimate = learnerTime - solverTime;
        }

        // --- Pass cost to learning algorithm
        if (static_ == -1 && !randomUniform_ && Pstream::master(matrix_.mesh().comm())) {
            referenceDictionary& learningDict = learningDicts[banditName_];
            learningDict.set<scalar>("loss", costEstimate);
        }
    }
    learnerTime -= clockValue::now();
    PCGTime -= solverTime;

    Info<< "INFO: banditName=" << banditName_;
    Info<< ", fieldName=" << fieldName_;
    Info<< ", relativeTolerance=" << relTol_;
    Info<< ", tolerance=" << tolerance_;
    Info<< ", initialResidual=" << solverPerf.initialResidual();
    Info<< ", finalResidual=" << solverPerf.finalResidual();
    Info<< ", nIterations=" << solverPerf.nIterations();
    Info<< ", initializeTime=" << -initializeTime;
    Info<< ", preconstructTime=" << -preconstructTime;
    Info<< ", iterationTime=" << -iterationTime;
    Info<< ", learnerTime=" << -learnerTime;
    Info<< ", solverTime=" << -solverTime;
    Info<< ", PCGTime=" << PCGTime;
    if (deterministic_ && armDrawn) {
        Info<< ", costEstimate=" << costEstimate;
        PCGCost += costEstimate;
        Info<< ", PCGCost=" << PCGCost;
    }
    Info<< endl;

    #ifdef DUMP_ABSOL
    #include "Absol/finishDump.H"
    #endif

    return solverPerf;
}

Foam::solverPerformance Foam::PCGBandit::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    return scalarSolve
    (
        tpsi.ref(),
        ConstPrecisionAdaptor<solveScalar, scalar>(source)(),
        cmpt
    );
}

// ************************************************************************* //
