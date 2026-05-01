/*---------------------------------------------------------------------------*\
Application
    spectral

Description
    Build a similarity matrix from eight ICTC preconditioners and one DIC
    preconditioner, decompose the corresponding Laplacian, and print the
    eigenvalues.
\*---------------------------------------------------------------------------*/

#include "DecomposedLaplacian.H"
#include "similarityMatrix.H"

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    using namespace Foam;

    List<dictionary> preconditionerDicts(9);

    for (label i = 0; i < 8; ++i)
    {
        preconditionerDicts[i].set("preconditioner", "ICTC");
        preconditionerDicts[i].set("droptol", Foam::pow(10.0, -4.0 + 0.5*i));
    }

    preconditionerDicts[8].set("preconditioner", "DIC");

    const SquareMatrix<scalar> similarityMatrix =
        preconditionerSimilarityMatrix(preconditionerDicts);

    const DecomposedLaplacian decomposedLaplacian(similarityMatrix);

    Info<< "similarity matrix:" << similarityMatrix << nl;

    Info<< "eigenvalues:" << decomposedLaplacian.EVals() << endl;

    Info<< "eigenrows:" << decomposedLaplacian.ERows() << endl;

    Info << "effective dimension for mu=0.1: " << decomposedLaplacian.dEff(0.1) << endl;

    Info<< "D-optimal design for mu=0.1: " << decomposedLaplacian.DOptimalDesign(0.1) << endl;

    return 0;
}

// ************************************************************************* //