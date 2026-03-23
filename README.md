# PCGBandit

Package for tuning PCG preconditioners in OpenFOAM simulations on-the-fly (PCGBandit).
In addition, this repository implements the following:

* An implementation of a thresholded incomplete Cholesky preconditioner (`ICTC`) with a customizable drop tolerance parameter (`droptol`).
* An implementation of GAMG (`FGAMG`) that precomputes ICTC (and DIC) smoothing factors before solving the system rather than at each iteration; it also allows the ICTC drop tolerance to vary across multigrid levels.

## Setup

### Docker setup

Run `sh launch.sh` from the repository root to pull the OpenFOAM Docker image, build all three libraries (`PCGBandit`, `ICTC`, `FGAMG`), and drop into an interactive container shell.
The repository is mounted at `/home/openfoam` inside the container.

### Building manually

From an OpenFOAM-sourced environment:
```
cd src/PCGBandit && wmake libso && cd ../..
cd src/ICTC && wmake libso && cd ../..
cd src/FGAMG && wmake libso && cd ../..
```

## Configuring your case

Add the following line to your OpenFOAM case directory's `system` subfolder:
```
libs ( libICTC.so libFGAMG.so libPCGBandit.so );
```

Then, for any PCG solver to be tuned, replace its `solver` and `preconditioner` specifications in `fvSolution` with PCGBandit settings.
A minimal configuration is:
```
solver          PCGBandit;
preconditioner  separate;
```
This tunes over the `DIC` preconditioner only (equivalent to standard PCG with DIC). The `separate` setting for `preconditioner` specifies tuning each linear system in the simulation separately; the alternative option is `joint`, with custom groupings possible using other names.
Set `randomSeed` in `controlDict` (not `fvSolution`) to control the random seed.

### Tuning the ICTC preconditioner drop tolerance

Tuning the ICTC preconditioner drop tolerance is controlled by the following parameters:
| Keyword | Default | Description |
|---------|---------|-------------|
| `numDroptols` | `0` | number of drop tolerances to tune over |
| `minLogDroptol` | `-4` | smallest drop tolerance (10<sup>-4</sup>) |
| `maxLogDroptol` | `-0.5` | largest drop tolerance (10<sup>-0.5</sup>) |

The solver will consider `numDroptols` different drop tolerances, evenly-spaced in log-space between `minLogDroptol` and `maxLogDroptol`.
For example, with `numDroptols 8; minLogDroptol -4; maxLogDroptol -0.5`, the drop tolerances will be approximately: `(0.0001 0.000316 0.001 0.00316 0.01 0.0316 0.1 0.316)`.
The `ICTC` preconditioner may also be used on its own as a preconditioner for PCG by replacing (e.g.) `preconditioner DIC;` with `preconditioner ICTC;` in the relevant `fvSolution` file.

### Tuning multigrid parameters

To tune over GAMG configurations, add one or more `<param>Tune` keywords, where `<param>` is a valid GAMG parameter.
Each `<param>Tune` keyword accepts either `yes` (use built-in defaults), `no` (default), or an explicit list of values:
```
solver          PCGBandit;
preconditioner  separate;
smootherTune    yes;            // default: (GaussSeidel DIC DICGaussSeidel symGaussSeidel)
nCellsInCoarsestLevelTune yes;  // default: (10 100 1000)
mergeLevelsTune yes;            // default: (1 2)
numDroptols     8;              // number of ICTC drop tolerances to tune
```

The available `<param>Tune` keywords and their defaults are:
| Keyword | Default values |
|---------|---------------|
| `smootherTune` | `GaussSeidel DIC DICGaussSeidel symGaussSeidel` |
| `agglomeratorTune` | `faceAreaPair algebraicPair` |
| `directSolveCoarsestTune` | `no yes` |
| `nCellsInCoarsestLevelTune` | `10 100 1000` |
| `mergeLevelsTune` | `1 2` |
| `nPreSweepsTune` | `0 2` |
| `nPostSweepsTune` | `1 2` |
| `nFinestSweepsTune` | `2` |
| `nVcyclesTune` | `1 2` |

Additional optional keywords (in addition to regular PCG keywords):
| Keyword | Default | Description |
|---------|---------|-------------|
| `residualContext` | `no` | whether to tune the `Final` solve in a correction loop separately |
| `DICTune` | `yes` | include `DIC` as a candidate preconditioner |
| `cacheAgglomeration` | `yes` (auto `no` if tuning agglomerator/nCellsInCoarsestLevel/mergeLevels) | GAMG agglomeration caching |
| `banditAlgorithm` | `TsallisINF` | bandit algorithm (`TsallisINF` or `ThompsonSampling`) |
| `lossEstimator` | `RV` | loss estimator for the bandit algorithm (`RV` or `IW`) |
| `backstop` | `-1` | backstop iteration limit (`-1` = auto) |
| `static` | `-1` | index of static preconditioner schedule (`-1` = off) |
| `deterministic` | `no` | deterministic mode for reproducible runs |

When tuning GAMG smoothers, you may also specify `ICTC` or `ICTCGaussSeidel` as options in `smootherTune` or `coarsestSmootherTune`; this approach is best-used with `FGAMG` loaded.
These will automatically expand to a family of options with different drop tolerances.
Drop tolerances are specified using logarithmic suffixes and controlled by:

| Keyword | Default | Description |
|---------|---------|-------------|
| `numSmootherDroptols` | (auto to the number of available drop tolerance (`8`)) | number of drop tolerances to tune over |
| `minSmootherLogDroptol` | `m4` | smallest drop tolerance for ICTC smoothers (10<sup>-4</sup>) |
| `maxSmootherLogDroptol` | `m0p5` | largest drop tolerance for ICTC smoothers (10<sup>-0.5</sup>) |

The available suffix values and their corresponding drop tolerances are:
| Suffix | Drop tolerance |
|--------|----------------|
| `m4` | 10<sup>-4</sup> = 0.0001 |
| `m3p5` | 10<sup>-3.5</sup> ≈ 0.000316 |
| `m3` | 10<sup>-3</sup> = 0.001 |
| `m2p5` | 10<sup>-2.5</sup> ≈ 0.00316 |
| `m2` | 10<sup>-2</sup> = 0.01 |
| `m1p5` | 10<sup>-1.5</sup> ≈ 0.0316 |
| `m1` | 10<sup>-1</sup> = 0.1 |
| `m0p5` | 10<sup>-0.5</sup> ≈ 0.316 |

You can also independently tune the smoother on the finest level and the coarse levels (requires `FGAMG`):
```
smootherTune       yes;                         // tune finest-level smoother
coarsestSmootherTune yes;                       // tune coarse-level smoother independently
```
With `coarsestSmootherTune yes` and ICTC/ICTCGaussSeidel in `smootherTune`, each combination of finest-level and coarse-level drop tolerances becomes a candidate configuration.

## Examples

The script `examples/run.sh` runs preconfigured simulations inside the Docker container.
Its first argument specifies the case name and its optional second argument `debug` runs with a shorter end time and deterministic cost estimation.

The two FreeMHD cases (`closedPipe`, `fringingBField`) require the case files shipped in `examples/FreeMHD.zip`.
Unzip it before running:
```
cd examples
unzip FreeMHD.zip
```

Example commands (run from within the container at `/home/openfoam/examples`):
```
bash run.sh boxTurb32              # OpenFOAM tutorial DNS/dnsFoam/boxTurb16 (at 2x resolution)
bash run.sh boxTurb32 debug        # same, short run for testing
bash run.sh pitzDaily              # OpenFOAM tutorial incompressible/pimpleFoam/RAS/pitzDaily (at 2x resolution)
bash run.sh interStefanProblem     # OpenFOAM tutorial verificationAndValidation/multiphase/StefanProblem (at 2x resolution)
bash run.sh porousDamBreak         # OpenFOAM tutorial verificationAndValidation/multiphase/interIsoFoam/porousDamBreak (at 2x resolution)
bash run.sh closedPipe             # FreeMHD case (requires FreeMHD.zip, 16 MPI ranks)
bash run.sh fringingBField         # FreeMHD case (requires FreeMHD.zip, 16 MPI ranks)
```

Each invocation runs three solver configurations and saves their logs alongside the case directory:
* `PCGBandit` — full bandit tuning over `ICTC`, `DIC`, and `GAMG` configurations
* `DIC` — baseline `DIC`-only preconditioner
* `GAMG` — `GAMG` with `DICGaussSeidel` smoother only

## References

1. Khodak, Jung, Wynne, Chow, Kolemen. *PCGBandit: One-shot acceleration of transient PDE solvers via online-learned preconditioners.* 2025.
2. Khodak, Chow, Balcan, Talwalkar. *Learning to relax: Setting solver parameters across a sequence of linear system instances.* ICLR 2024.
3. Wynne, Saenz, Al-Salami, Xu, Sun, Hu, Hanada, Kolemen. *FreeMHD: Validation and verification of the open-source, multi-domain, multi-phase solver for electrically conductive flows.* Phys. Plasmas **32** (1) 2025.
