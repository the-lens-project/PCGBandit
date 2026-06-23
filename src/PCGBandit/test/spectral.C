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


void small() {

  using namespace Foam;

  label numDroptols = 8;
  List<dictionary> preconditionerDicts(numDroptols + 1);

  for (label i = 0; i < numDroptols; ++i)
  {
      preconditionerDicts[i].set("preconditioner", "ICTC");
      preconditionerDicts[i].set("droptol", Foam::pow(10.0, -4.0 + (4.0 / numDroptols)*i));
  }

  preconditionerDicts[numDroptols].set("preconditioner", "DIC");

  const SquareMatrix<scalar> S = pathMatrix(preconditionerDicts);

  const DecomposedLaplacian decomposedLaplacian(S);

  Info<< "similarity matrix:" << S << nl;

  Info<< "eigenvalues:" << decomposedLaplacian.EVals() << endl;

  Info<< "eigenrows:" << decomposedLaplacian.ERows() << endl;

  for (label i = -8; i <= 4; ++i) {

      const scalar mu = Foam::pow(10.0, scalar(i) / 2.0);
      Info<< "effective dimension for mu=" << mu << ": " << decomposedLaplacian.dEff(mu) << endl;
      scalarField Pi = decomposedLaplacian.DOptimalDesign(mu);
      Info<< "D-optimal design for mu=" << mu << ": " << Pi << endl;
      Pi = 1.0 / scalar(numDroptols + 1);
      Info<< "fhat=" << decomposedLaplacian.getHat(Pi, mu, numDroptols / 2) << endl;

  }

}

void big() {

  using namespace Foam;

  List<word> smoothers = {"GaussSeidel", "DIC", "DICGaussSeidel", "symGaussSeidel"};
  List<label> nCells = {10, 100, 1000};
  List<label> mergeLevels = {1, 2};
  List<word> ICTCSuffixes = {};//{"m4", "m3p5", "m3", "m2p5", "m2", "m1p5", "m1", "m0p5"};

  DynamicList<dictionary> preconditionerDicts;

  for (label i = 0; i < 8; ++i) {
      dictionary dict;
      dict.set("preconditioner", "ICTC");
      dict.set("droptol", Foam::pow(10.0, -4.0 + 0.5*i));
      preconditionerDicts.append(dict);
  }
  dictionary dict;
  dict.set("preconditioner", "DIC");
  preconditionerDicts.append(dict);

  for (label i = 0; i < smoothers.size(); ++i) {
      dictionary dict;
      dict.set("preconditioner", "FGAMG");
      dict.set("smoother", smoothers[i]);
      for (label j = 0; j < nCells.size(); ++j) {
          dict.set("nCellsInCoarsestLevel", nCells[j]);
          for (label k = 0; k < mergeLevels.size(); ++k) {
              dict.set("mergeLevels", mergeLevels[k]);
              preconditionerDicts.append(dict);
          }
      }
  }
  for (label i = 0; i < ICTCSuffixes.size(); ++i) {
      dictionary dict;
      dict.set("preconditioner", "FGAMG");
      dict.set("smoother", "ICTC_" + ICTCSuffixes[i]);
      for (label j = 0; j < ICTCSuffixes.size(); ++j) {
          dict.set("coarsestSmoother", "ICTC_" + ICTCSuffixes[j]);
          for (label k = 0; k < nCells.size(); ++k) {
              dict.set("nCellsInCoarsestLevel", nCells[k]);
              for (label l = 0; l < mergeLevels.size(); ++l) {
                  dict.set("mergeLevels", mergeLevels[l]);
                  preconditionerDicts.append(dict);
              }
          }
      }
  }


  SquareMatrix<scalar> S = pathMatrix(preconditionerDicts);
  label n = preconditionerDicts.size();
  label row = n / 2;

  const DecomposedLaplacian decomposedLaplacian(S);
  scalar mu = 0.1;
  scalarField Pi = decomposedLaplacian.DOptimalDesign(mu);
  Info<< "D-optimal design for mu=" << mu << ": " << Pi << endl;
  Pi = 1.0 / scalar(n);
  Info<< "fhat=" << decomposedLaplacian.getHat(Pi, mu, row) << endl;

  Info<< preconditionerDicts[row] << endl;
  Info<< "HAS NEIGHBORS:" << endl;
  for (label i = 0; i < n; ++i) {
      if (S(i, row) == 1.0) {
          Info<< preconditionerDicts[i] << endl;
      }
  }

  if (n < 50) { 
      for (label i = 0; i < n; ++i) {
          for (label j = 0; j < n; ++j) {
              Info<< S(i, j) << ",";
          }
          Info<< endl;
      }
  }

}


int main(int argc, char *argv[]) {

    (void)argc;
    (void)argv;

    small();
    big();

    return 0;
}

// ************************************************************************* //
