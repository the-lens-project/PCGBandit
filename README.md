# PCGBandit

PCGBandit tunes PCG preconditioners during OpenFOAM simulations. The package includes:

- `ICTC`: incomplete Cholesky with a configurable drop tolerance.
- `FGAMG`: multigrid with cached ICTC/DIC smoothing factors and smoother parameters that can vary across levels.
- `SOR` and `DICSOR`: parameterized SOR smoothers, with optional DIC smoothing before SOR.
- `subspaceInitialization`: initial-guess corrections from previous iterates, including the standalone `siPCG` solver.
- `SpectralINF`: a bandit algorithm that shares feedback between similar configurations.

## Setup

### Docker

Run `bash launch.sh` from the repository root. It pulls the OpenFOAM image, builds all five libraries, and opens a shell with the repository mounted at `/home/openfoam`.

### Manual build

From the repository root in an OpenFOAM-sourced environment:

```sh
cd src/ICTC && wmake libso && cd ../..
cd src/SOR && wmake libso && cd ../..
cd src/subspaceInitialization && wmake libso && cd ../..
cd src/PCGBandit && wmake libso && cd ../..
cd src/FGAMG && wmake libso && cd ../..
```

## Configuring a case

To enable all components, add this to `system/controlDict`:

```foam
libs ( libICTC.so libSOR.so libsubspaceInitialization.so libFGAMG.so libPCGBandit.so );
```

`libPCGBandit` automatically loads its required `libsubspaceInitialization` dependency. Load `libICTC` and `libSOR` when selecting their preconditioners or smoothers. Load `libFGAMG` before `libPCGBandit` to use FGAMG for multigrid candidates; otherwise PCGBandit uses built-in GAMG.

In the relevant `system/fvSolution` solver entry, keep the existing tolerances and replace `solver` and `preconditioner` with:

```foam
solver          PCGBandit;
preconditioner  separate;
```

This selects only DIC until more candidates are configured. `separate` maintains a bandit per mesh region and field; `joint` shares one bandit, and any other name defines a custom group. The first solver dictionary in a group defines its candidate configurations, so use consistent tuning settings within each group.

Set `randomSeed` in `controlDict` to control random sampling.

### Solver controls

These supplement the usual PCG controls:

| Keyword | Default | Description |
|---------|---------|-------------|
| `DICTune` | `yes` | Include DIC as a candidate. |
| `residualContext` | `no` | Use a separate bandit for solves with `relTol 0`. |
| `banditAlgorithm` | `TsallisINF` | `TsallisINF`, `SpectralINF`, or `ThompsonSampling`. |
| `lossEstimator` | `RV` | `RV` or `IW` for TsallisINF; unused by the other algorithms. |
| `deterministic` | `no` | Use estimated operation cost instead of elapsed solver time as feedback. |
| `randomUniform` | `no` | Choose candidates uniformly without learning. |
| `static` | `-1` | Fix a zero-based candidate index; `-1` enables selection. |
| `backstop` | `-1` | DIC fallback: `-1` chooses a threshold from estimated costs; `0` disables fallback; a positive value sets the additional iteration budget. |
| `cacheAgglomeration` | `yes` | Cache multigrid agglomeration. Disabled when tuning multiple agglomeration settings, except with a fixed `static` candidate. |

### ICTC drop tolerance

| Keyword | Default | Description |
|---------|---------|-------------|
| `numDroptols` | `0` | Number of ICTC preconditioner candidates. |
| `minLogDroptol` | `-4` | Base-10 logarithm of the smallest drop tolerance. |
| `maxLogDroptol` | `-0.5` | Base-10 logarithm of the largest drop tolerance. |

Tolerances are evenly spaced in log space. With eight candidates and the default bounds, they are approximately `(0.0001 0.000316 0.001 0.00316 0.01 0.0316 0.1 0.316)`. A single candidate uses `minLogDroptol`.

To use ICTC directly with PCG:

```foam
solver PCG;
preconditioner
{
    preconditioner ICTC;
    droptol        1e-3;
}
```

### Multigrid tuning

Each keyword below accepts `yes` to use its default grid, `no` to leave that axis untuned, or an explicit list. The chosen axes form a Cartesian product of multigrid candidates, added alongside DIC and ICTC candidates.

```foam
smootherTune                 yes;
nCellsInCoarsestLevelTune    yes;
mergeLevelsTune             yes;
numDroptols                 8;
```

| Keyword | Grid used by `yes` |
|---------|--------------------|
| `smootherTune` | `(GaussSeidel DIC DICGaussSeidel symGaussSeidel)` |
| `agglomeratorTune` | `(faceAreaPair algebraicPair)` |
| `directSolveCoarsestTune` | `(no yes)` |
| `nCellsInCoarsestLevelTune` | `(10 100 1000)` |
| `mergeLevelsTune` | `(1 2)` |
| `nPreSweepsTune` | `(0 2)` |
| `nPostSweepsTune` | `(1 2)` |
| `nFinestSweepsTune` | `(2)` |
| `nVcyclesTune` | `(1 2)` |

Without smoother tuning, multigrid candidates use `DICGaussSeidel`. Explicit `nCellsInCoarsestLevelTune` list values are clamped to the minimum local cell count across ranks, with a lower bound of one.

#### ICTC smoothers

Include `ICTC` or `ICTCGaussSeidel` in an explicit `smootherTune` list to expand that family over drop tolerances. FGAMG caches their factors between smoothing calls.

| Keyword | Default | Description |
|---------|---------|-------------|
| `numSmootherDroptols` | `4` | Target number of suffixes to select. |
| `minSmootherLogDroptol` | `m4` | Lower drop-tolerance bound, encoded as a suffix. |
| `maxSmootherLogDroptol` | `m0p5` | Upper drop-tolerance bound. |

Suffixes range from `m5` to `m0p5` in half-decade steps: `m5`, `m4p5`, `m4`, ..., `m1`, `m0p5`. For example, `ICTC_m3p5` uses a drop tolerance of 10^-3.5. The default selection is `(m4 m3 m2 m1)`; `numSmootherDroptols 8` selects all suffixes within the default bounds.

The ICTC and SOR smoother grids use an integer stride through their suffix lists. The requested count is a target: the actual count can differ, and the upper bound can be omitted.

#### SOR smoothers

Include `SOR` or `DICSOR` in `smootherTune` to expand over relaxation factors:

```foam
smootherTune       (SOR DICSOR);
minSmootherOmega   p0p6;
maxSmootherOmega   p1p4;
numSmootherOmegas  5;
```

| Keyword | Default | Description |
|---------|---------|-------------|
| `minSmootherOmega` | `p0p8` | Lower relaxation-factor bound. |
| `maxSmootherOmega` | `p1p2` | Upper relaxation-factor bound. |
| `numSmootherOmegas` | `5` | Target number of suffixes to select. |

Available factors run from `p0p1` to `p1p9` in steps of 0.1. The example selects `(p0p6 p0p8 p1p0 p1p2 p1p4)`; the default grid is `(p0p8 p0p9 p1p0 p1p1 p1p2)`.

`SOR_p1p0` and `DICSOR_p1p0` use the SOR implementations at omega=1. `GaussSeidel` and `DICGaussSeidel` retain the built-in implementations and remain separate candidates if explicitly included. FGAMG and PCGBandit do not link `libSOR` automatically.

#### Smoother parameters across levels

When configuring FGAMG directly:

```foam
solver             FGAMG;
smoother           ICTC_m4;
coarsestSmoother    ICTC_m1;
```

FGAMG interpolates through the available values within an ICTC, ICTCGaussSeidel, SOR, or DICSOR family. GaussSeidel and DICGaussSeidel can serve as omega=1 endpoints in the corresponding SOR families while retaining their built-in implementations. Omitting `coarsestSmoother` uses the same smoother on all levels.

PCGBandit does not implement `coarsestSmootherTune`; each candidate uses one smoother setting across levels.

### Subspace initialization

Subspace initialization corrects the incoming PCG guess using a randomized subspace of previous iterates. It is disabled by default. Choose a window with `lenHistory`, or an exponentially weighted moving-average sketch with `decayRate`:

```foam
solver          PCGBandit;
preconditioner  separate;
lenHistory      8;              // alternatively: decayRate 0.5;
numProbes       4;
```

The same plain settings work with `solver siPCG; preconditioner DIC;`, provided `libsubspaceInitialization` is loaded. Other solver integrations can use the initializer with an `fvMesh` registry.

| Keyword | Default | Description |
|---------|---------|-------------|
| `lenHistory` | `0` | Maximum stored iterates; `0` disables the window. |
| `decayRate` | `0.0` | Weight multiplier per stored solve. Use `0 < decayRate < 1` for exponential forgetting; `0` disables the EWMA sketch. |
| `numProbes` | `4` | Number of random directions; `0` disables the correction. |
| `projection` | `galerkin` | Energy projection for a symmetric definite matrix; `leastSquares` instead minimizes the residual 2-norm. |
| `truncTol` | `1e-14` | Discard reduced eigenmodes at or below this fraction of the largest eigenvalue magnitude. |
| `persistState` | `no` | Write window or EWMA state at simulation write times and resume on restart. |

PCG's matrix requirements still apply with either projection. `siPCG` rejects positive `lenHistory` and `decayRate` together; PCGBandit treats them as alternative candidate configurations.

History counts solves that reach initialization, not timesteps. Each sample is the incoming iterate before correction. Correction starts as soon as usable history exists, with the probe count capped by the available history and configured capacity. Every maintained window and sketch is updated, including when an off arm is selected.

A solve that exits at the initial convergence check adds no sample and counts as no arm pull. If initialization itself converges, it still counts as a pull and skips preconditioner construction and CG iterations when `minIter` permits.

State is shared per mesh, field, and component, including between `p` and `pFinal`. The first caller fixes capacity, seed, and persistence; use matching settings and tuning grids for dictionaries sharing a field. Later requests exceeding capacity are clamped with a warning.

A window allocates one local solution vector per slot in the maximum configured `lenHistory`. EWMA allocates the maximum configured `numProbes` vectors per rate. A correction additionally allocates about twice its probe count in local vectors. Larger probe counts increase projection work quadratically through the reduced Gram matrix; they do not necessarily reduce total solve time.

#### Tuning initialization

`lenHistoryTune`, `decayRateTune`, and `numProbesTune` accept `yes`, `no`, or a list. `yes` uses the grid below; `no` keeps the plain setting or default. Lists override plain settings, must contain non-negative values, and require integers for lengths and probe counts.

```foam
lenHistoryTune  yes;            // (0 4 8 12 16)
decayRateTune   yes;            // (0 0.25 0.5 0.75)
numProbesTune   yes;            // (4 8 12 16)
```

An empty history or rate list excludes that form; an empty probe list disables initialization. Supply unique decay rates with distinct registry-key representations: rate lists are not deduplicated.

Probe counts are paired with each history length and each rate. Window and EWMA configurations form a union, and each is combined with every preconditioner candidate. Window probe counts are clamped to the window length; duplicate configurations collapse. The three default grids produce 23 subspace configurations, or 115 arms with `numDroptols 4` and the default DIC candidate.

Zero-valued choices add an off configuration. Tuning does not otherwise add one; if no subspace configurations remain, the ordinary preconditioner candidates are retained. Window configurations precede EWMA configurations, following input order with probe count varying inside each. Consequently, off is not necessarily arm zero.

#### Restarting

Without persistence, subspace state starts empty after a restart. With `persistState yes`, binary state files are written as `contextWindow:<field>:<component>` and `EWMASketch:<field>:<component>:<rate>` in the time directory.

Reuse state only with the same mesh ordering, decomposition, capacities, and seed. Size mismatches are rejected, but equal-sized changes in cell ordering are not detected. Restart fields and state must come from the same checkpoint; binary field output avoids additional rounding. Persistence restores subspace state, not the bandit's learning history, and does not by itself guarantee an identical simulation trajectory.

The [design document](SUBSPACE_INITIALIZATION.md) contains the equations, implementation rationale, and historical benchmark results. The methods build on the [windowed subspace approach](https://arxiv.org/abs/2309.02156) and its [exponentially weighted extension](https://arxiv.org/abs/2511.18071).

## Examples

Run `examples/run.sh` with Bash inside the container. The optional `debug` argument shortens the run and enables deterministic cost estimation.

The FreeMHD cases require the bundled archive to be extracted first:

```sh
cd /home/openfoam/examples
unzip FreeMHD.zip
```

From `/home/openfoam/examples`:

```sh
bash run.sh boxTurb32              # dnsFoam/boxTurb16, twice the resolution
bash run.sh boxTurb32 debug        # short test run
bash run.sh pitzDaily              # pimpleFoam/RAS/pitzDaily, twice the resolution
bash run.sh interStefanProblem     # StefanProblem, twice the resolution
bash run.sh porousDamBreak         # interIsoFoam/porousDamBreak, twice the resolution
bash run.sh closedPipe             # FreeMHD, 16 MPI ranks
bash run.sh fringingBField          # FreeMHD, 16 MPI ranks
```

By default, each invocation runs three configurations through PCGBandit and saves their logs under `examples/<case>/`:

- `PCGBandit`: tuning over ICTC, DIC, and multigrid candidates.
- `DIC`: one DIC candidate.
- `GAMG`: one multigrid candidate with DICGaussSeidel smoothing. Despite the log name, these examples load FGAMG and use it.

## References

1. Khodak, Jung, Wynne, Chow, Kolemen. [One-shot acceleration of transient PDE solvers via online-learned preconditioners](https://arxiv.org/abs/2509.08765). 2025 preprint.
2. Khodak, Chow, Balcan, Talwalkar. [Learning to relax: Setting solver parameters across a sequence of linear system instances](https://arxiv.org/abs/2310.02246). ICLR 2024.
3. Wynne, Saenz, Al-Salami, Xu, Sun, Hu, Hanada, Kolemen. [FreeMHD: Validation and verification of the open-source, multi-domain, multi-phase solver for electrically conductive flows](https://arxiv.org/abs/2409.08950). Physics of Plasmas 32 (1), 2025.
