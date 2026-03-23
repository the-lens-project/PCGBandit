#!/bin/sh

DEFAULT="solver PCGBandit; preconditioner separate; backstop 10000;"

NAME=${1:-pitzDaily}
SEED=0
if [ $NAME = "closedPipe" ] || [ $NAME = "fringingBField" ] ; then
  NPROC=16
  cd FreeMHD && wmake solvers/epotMultiRegionInterFoam && cd ..
else
  NPROC=1
fi
if [ "$2" = "debug" ] ; then
  DEBUG="true"
  DEFAULT="$DEFAULT deterministic yes;"
fi

mkdir -p $NAME
cd $NAME


################################################################################
##################################   cases   ###################################
################################################################################

if [ $NAME == "boxTurb32" ] ; then

  cp -r $FOAM_TUTORIALS/DNS/dnsFoam/boxTurb16 boxTurb32
  cd boxTurb32

  sed -i '16a\libs ( libICTC.so libFGAMG.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '18a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  sed -i 's/writeInterval   0\.25/writeInterval   100/' system/controlDict
  if [ $DEBUG ] ; then
    sed -i 's/endTime         10/endTime         0.1/' system/controlDict
  fi

  sed -i 's/16/32/g' system/blockMeshDict
  sed -i 's/nonuniform List<vector>/uniform (0 0 0);/' 0.orig/U
  sed -i '20,4119d' 0.orig/U
  sed -i 's/deltaT          0\.025/deltaT          0.005/' system/controlDict

  sed -i '21d' system/fvSolution
  sed -i 's/preconditioner  DIC;/'"$DEFAULT"'/' system/fvSolution

  runSimulation() {
    sed -i 's/'"$DEFAULT"'/'"$1"'/' system/fvSolution
    cp -r 0.orig 0
    blockMesh > log.blockMesh
    boxTurb > log.boxTurb
    echo $2 $1
    dnsFoam > log.dnsFoam
    mv log.dnsFoam $2
    foamCleanTutorials
    sed -i 's/'"$1"'/'"$DEFAULT"'/' system/fvSolution
  }

fi

################################################################################

if [ $NAME == "pitzDaily" ] ; then

  cp -r $FOAM_TUTORIALS/incompressible/pimpleFoam/RAS/pitzDaily .
  cd pitzDaily

  sed -i '16a\libs ( libICTC.so libFGAMG.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '18a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  sed -i 's/writeInterval   0\.01/writeInterval   1/' system/controlDict
  if [ $DEBUG ] ; then
    sed -i 's/endTime         0\.3/endTime         0.0003/' system/controlDict
  fi

  sed -i 's/(1 3 0\.3)/(1 3 0.15)/' system/blockMeshDict
  sed -i 's/(2 4 0\.25)/(2 4 0.125)/' system/blockMeshDict
  sed -i 's/(1 1 0\.25)/(1 1 0.125)/' system/blockMeshDict
  sed -i 's/(18 30 1)/(36 60 1)/' system/blockMeshDict
  sed -i 's/(180 27 1)/(360 54 1)/' system/blockMeshDict
  sed -i 's/(180 30 1)/(360 60 1)/' system/blockMeshDict
  sed -i 's/(25 27 1)/(50 54 1)/' system/blockMeshDict
  sed -i 's/(25 30 1)/(50 60 1)/' system/blockMeshDict

  sed -i '21d' system/fvSolution
  sed -i 's/smoother         DICGaussSeidel;/'"$DEFAULT"'/' system/fvSolution

  runSimulation() {
    sed -i 's/'"$DEFAULT"'/'"$1"'/' system/fvSolution
    blockMesh > log.blockMesh
    echo $2 $1
    pimpleFoam > log.pimpleFoam
    mv log.pimpleFoam $2
    foamCleanTutorials
    sed -i 's/'"$1"'/'"$DEFAULT"'/' system/fvSolution
  }

fi

################################################################################

if [ $NAME == "interStefanProblem" ] ; then

  mkdir interStefanProblem
  cd interStefanProblem
  for SUBFOLDER in 0.orig constant system ; do
    mkdir $SUBFOLDER
    cp -r $FOAM_TUTORIALS/verificationAndValidation/multiphase/StefanProblem/setups.orig/common/$SUBFOLDER/* $SUBFOLDER/.
    cp -r $FOAM_TUTORIALS/verificationAndValidation/multiphase/StefanProblem/setups.orig/interCondensatingEvaporatingFoam/$SUBFOLDER/* $SUBFOLDER/.
  done

  sed -i '17a\libs ( libICTC.so libFGAMG.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '19a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  sed -i 's/writeInterval   5/writeInterval   100/' system/controlDict
  if [ $DEBUG ] ; then
    sed -i 's/endTime         50/endTime         1.4/' system/controlDict
  fi

  sed -i 's/(400 2 1)/(800 4 1)/' system/blockMeshDict

  sed -i '45d' system/fvSolution
  sed -i '/solver *PCG;/d' system/fvSolution
  sed -i '/smoother *DIC;/d' system/fvSolution
  sed -i 's/preconditioner *DIC;/'"$DEFAULT"'/' system/fvSolution

  runSimulation() {
    sed -i 's/'"$DEFAULT"'/'"$1"'/' system/fvSolution
    cp -r 0.orig 1.36
    blockMesh > log.blockMesh
    renumberMesh -overwrite -constant > log.renumberMesh
    checkMesh > log.checkMesh
    setAlphaField > log.setAlphaField
    echo $2 $1
    interCondensatingEvaporatingFoam > log.interCondensatingEvaporatingFoam
    mv log.interCondensatingEvaporatingFoam $2
    foamCleanTutorials
    sed -i 's/'"$1"'/'"$DEFAULT"'/' system/fvSolution
  }

fi

################################################################################

if [ $NAME == "porousDamBreak" ] ; then

  cp -r $FOAM_TUTORIALS/verificationAndValidation/multiphase/interIsoFoam/porousDamBreak porousDamBreak
  cd porousDamBreak

  sed -i '16a\libs ( libICTC.so libFGAMG.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '18a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  sed -i 's/writeInterval   0\.05/writeInterval   10/' system/controlDict
  if [ $DEBUG ] ; then
    sed -i 's/endTime         4/endTime         0.004/' system/controlDict
  fi

  sed -i 's/(75 93 1)/(150 186 1)/' system/blockMeshDict
  sed -i 's/(73 93 1)/(146 186 1)/' system/blockMeshDict
  sed -i 's/(76 93 1)/(152 186 1)/' system/blockMeshDict
  sed -i 's/(75 53 1)/(150 106 1)/' system/blockMeshDict
  sed -i 's/(73 53 1)/(146 106 1)/' system/blockMeshDict
  sed -i 's/(76 53 1)/(152 106 1)/' system/blockMeshDict

  sed -i '/solver *PCG;/d' system/fvSolution
  sed -i 's/preconditioner *DIC;/'"$DEFAULT"'/' system/fvSolution

  runSimulation() {
    sed -i 's/'"$DEFAULT"'/'"$1"'/' system/fvSolution
    cp -r 0.orig 0
    blockMesh > log.blockMesh
    setFields > log.setFields
    echo $2 $1
    interIsoFoam > log.interIsoFoam
    mv log.interIsoFoam $2
    foamCleanTutorials
    sed -i 's/'"$1"'/'"$DEFAULT"'/' system/fvSolution
  }

fi

################################################################################

if [ $NAME == "closedPipe" ] ; then

  cp -r ../FreeMHD/closedPipe .
  cd closedPipe
  wmake libso dynamicCode/outletUxB

  sed -i '19a\libs ( libICTC.so libFGAMG.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '21a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  sed -i 's/writeInterval   1e-5/writeInterval   1e-1/' system/controlDict
  if [ $DEBUG ] ; then
    sed -i 's/endTime         0.025/endTime         2.5e-6/' system/controlDict
  fi

  for region in $(foamListRegions) ; do
    sed -i '/solver *PCG;/d' system/$region/fvSolution
    sed -i 's/preconditioner *DIC;/'"$DEFAULT"'/' system/$region/fvSolution
  done

  runSimulation() {
    cp -r 0.orig 0
    blockMesh > log.blockMesh
    topoSet > log.topoSet
    splitMeshRegions -cellZonesOnly -overwrite -fileHandler collated > log.splitMeshRegions
    for region in $(foamListRegions) ; do
      sed -i 's/'"$DEFAULT"'/'"$1"'/' system/$region/fvSolution
      changeDictionary -region $region -fileHandler collated > log.changeDictionary.$region
      setExprFields -region $region -fileHandler collated > log.setExprFields.$region
    done

    echo $2 $1
    sed -i 's/numberOfSubdomains.*;/numberOfSubdomains '"$NPROC"';/' system/decomposeParDict
    for region in $(foamListRegions) ; do
      sed -i 's/numberOfSubdomains.*;/numberOfSubdomains '"$NPROC"';/' system/$region/decomposeParDict
    done
    decomposePar -allRegions -force -fileHandler collated > log.decomposePar
    mpirun -n $NPROC epotMultiRegionInterFoam -parallel > log.epotMultiRegionInterFoam
    mv log.epotMultiRegionInterFoam $2

    rm -rf log.*
    rm -rf 0
    rm -rf constant/cellToRegion
    rm -rf postProcessing
    for region in $(foamListRegions) ; do
      rm -rf constant/$region/polyMesh
      sed -i 's/'"$1"'/'"$DEFAULT"'/' system/$region/fvSolution
    done
  }

fi

################################################################################

if [ $NAME == 'fringingBField' ] ; then

  cp -r ../FreeMHD/fringingBField .
  cd fringingBField
  wmake libso dynamicCode/outletUxB
  
  sed -i '19a\libs ( libICTC.so libFGAMG.so libPCGBandit.so );\
  ' system/controlDict
  sed -i '21a\randomSeed\t\t\t'"$SEED"';\
  ' system/controlDict
  sed -i 's/writeInterval   0.1/writeInterval   1/' system/controlDict
  if [ $DEBUG ] ; then
    sed -i 's/endTime         0.5/endTime         5e-6/' system/controlDict
  fi

  for region in $(foamListRegions) ; do
    sed -i '/solver *PCG;/d' system/$region/fvSolution
    sed -i 's/preconditioner *DIC;/'"$DEFAULT"'/' system/$region/fvSolution
  done

  runSimulation() {
    cp -r 0.orig 0
    blockMesh > log.blockMesh
    splitMeshRegions -cellZonesOnly -overwrite -fileHandler collated > log.splitMeshRegions
    for region in $(foamListRegions) ; do
      sed -i 's/'"$DEFAULT"'/'"$1"'/' system/$region/fvSolution
      changeDictionary -region $region -fileHandler collated > log.changeDictionary.$region
      setExprFields -region $region -fileHandler collated > log.setExprFields.$region
    done

    echo $2 $1
    sed -i 's/numberOfSubdomains.*;/numberOfSubdomains '"$NPROC"';/' system/decomposeParDict
    for region in $(foamListRegions) ; do
      sed -i 's/numberOfSubdomains.*;/numberOfSubdomains '"$NPROC"';/' system/$region/decomposeParDict
    done
    decomposePar -allRegions -force -fileHandler collated > log.decomposePar
    mpirun -n $NPROC epotMultiRegionInterFoam -parallel > log.epotMultiRegionInterFoam
    mv log.epotMultiRegionInterFoam $2

    rm -rf log.*
    rm -rf 0
    rm -rf constant/cellToRegion
    rm -rf postProcessing
    for region in $(foamListRegions) ; do
      rm -rf constant/$region/polyMesh
      sed -i 's/'"$1"'/'"$DEFAULT"'/' system/$region/fvSolution
    done
  }

fi

################################################################################
###################################   runs   ###################################
################################################################################

# PCGBandit
SOLVER="$DEFAULT smootherTune yes; nCellsInCoarsestLevelTune yes; mergeLevelsTune yes; numDroptols 8;"
runSimulation "$SOLVER" ../PCGBandit

# DIC
SOLVER=$DEFAULT
runSimulation "$SOLVER" ../DIC

# GAMG with DICGaussSeidel smoother
SOLVER="$DEFAULT smootherTune (DICGaussSeidel); DICTune no;"
runSimulation "$SOLVER" ../GAMG

################################################################################

cd ..
rm -rf $NAME
