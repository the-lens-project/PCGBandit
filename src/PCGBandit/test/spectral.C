/*---------------------------------------------------------------------------*\
Application
    spectral

Description
    Build a similarity matrix from eight ICTC preconditioners and one DIC
    preconditioner, decompose the corresponding Laplacian, and print the
    eigenvalues.
\*---------------------------------------------------------------------------*/

#include "DecomposedLaplacian.H"
#include "Similarity.H"

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    using namespace Foam;

    label numDroptols = 16;
    List<dictionary> preconditionerDicts(numDroptols + 1);

    for (label i = 0; i < numDroptols; ++i)
    {
        preconditionerDicts[i].set("preconditioner", "ICTC");
        preconditionerDicts[i].set("droptol", Foam::pow(10.0, -4.0 + (4.0 / numDroptols)*i));
    }

    preconditionerDicts[numDroptols].set("preconditioner", "DIC");

    const SquareMatrix<scalar> S_decay = pathMatrix(preconditionerDicts);

    const DecomposedLaplacian decomposedLaplacian(S_decay);

    Info<< "similarity matrix:" << S_decay << nl;

    Info<< "eigenvalues:" << decomposedLaplacian.EVals() << endl;

    Info<< "eigenrows:" << decomposedLaplacian.ERows() << endl;

    for (label i = -8; i <= 4; ++i) {

        const scalar mu = Foam::pow(10.0, scalar(i) / 2.0);
        Info << "effective dimension for mu=" << mu << ": " << decomposedLaplacian.dEff(mu) << endl;
        scalarField Pi = decomposedLaplacian.DOptimalDesign(mu);
        Pi = 1.0 / scalar(numDroptols + 1);
        Info << "fhat=" << decomposedLaplacian.getHat(Pi, mu, numDroptols / 2) << endl;

    }

    return 0;
}

// ************************************************************************* //