/*---------------------------------------------------------------------------*\
Application
    spectral

Description
    Builds path/similarity matrices over three config spaces of increasing
    size -- small() (eight ICTC preconditioners and one DIC preconditioner),
    medium() (small() plus GaussSeidel/DIC/DICGaussSeidel/symGaussSeidel GAMG
    smoothers over a grid of nCellsInCoarsestLevel/mergeLevels), and large()
    (medium() with the smoother axis further extended by the SOR/DICSOR
    omega family) -- decomposes the corresponding Laplacian, and prints
    eigenvalues/D-optimal designs/neighbor lists for inspection.
\*---------------------------------------------------------------------------*/

#include "decomposedLaplacian.H"
#include "similarityMatrix.H"


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

  const decomposedLaplacian decomposedLaplacian(S);

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

void medium() {

  using namespace Foam;

  List<word> smoothers = {"GaussSeidel", "DIC", "DICGaussSeidel", "symGaussSeidel"};
  List<label> nCells = {10, 100, 1000};
  List<label> mergeLevels = {1, 2};

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


  SquareMatrix<scalar> S = pathMatrix(preconditionerDicts);
  label n = preconditionerDicts.size();
  label row = n / 2;

  const decomposedLaplacian decomposedLaplacian(S);
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

  for (label i = 0; i < n; ++i) {
      for (label j = 0; j < n; ++j) {
          Info<< S(i, j) << ",";
      }
      Info<< endl;
  }

}

void large() {

  using namespace Foam;

  List<word> smoothers = {
      "GaussSeidel", "DIC", "DICGaussSeidel", "symGaussSeidel",
      "SOR_p0p8", "SOR_p1p2", "DICSOR_p0p8", "DICSOR_p1p2"
  };
  List<label> nCells = {10, 100, 1000};
  List<label> mergeLevels = {1, 2};

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

  SquareMatrix<scalar> S = pathMatrix(preconditionerDicts);
  label n = preconditionerDicts.size();
  label row = n / 2;

  const decomposedLaplacian decomposedLaplacian(S);
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

}


int main(int argc, char *argv[]) {

  (void)argc;
  (void)argv;

  small();
  medium();
  large();

  return 0;
}

// ************************************************************************* //
