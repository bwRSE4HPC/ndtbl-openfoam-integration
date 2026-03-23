# bwRSE4HPCFoam Project

[![Build](https://github.com/bwRSE4HPC/bwRSE4HPCFoam/actions/workflows/ci.yml/badge.svg)](https://github.com/bwRSE4HPC/bwRSE4HPCFoam/actions)


## What this repo is

This repository is a small OpenFOAM extension focused on benchmarking
memory use and runtime for tabulated combustion lookups, not on running a
full CFD time-marching solver.

At a high level it provides:

- A custom OpenFOAM application: `bwRSE4HPCFoam`
- A custom user library: `libbwRSE4HPCcombustionModels`
- A tutorial case: `tutorials/D4DummySetup`

The codebase is compact: about 66 files total, with roughly 2,073 lines of
`.C`/`.H` source under `src/` and `applications/`.

## Core idea

The solver initializes an OpenFOAM case, constructs the thermodynamics and
transport models, instantiates a custom combustion model, performs exactly one
`reaction->correct()` call, and reports:

- combustion model initialization time
- table lookup time
- total CPU time
- resident memory usage (summed across MPI ranks)

So this is essentially a synthetic benchmark harness for table interpolation in
an OpenFOAM-style combustion model.

## Main components

### 1. Solver application

`applications/solver/bwRSE4HPCFoam/bwRSE4HPCFoam.C`

- Entry point for the benchmark executable.
- Uses standard OpenFOAM setup snippets (`createTime.H`, `createMesh.H`,
  `createFields.H`, etc.).
- Creates the selected combustion model via `combustionModel::New(...)`.
- Calls `reaction->correct()` once.
- Reads `/proc/self/status` to estimate RSS memory and reduces it across MPI
  processes with `Foam::reduce`.

`applications/solver/bwRSE4HPCFoam/createFields.H`

- Builds thermodynamics, velocity, pressure, turbulence, and transport fields.
- Instantiates the custom combustion model.
- Records initialization timing around that construction.

### 2. Combustion model framework

`src/combustionModels/combustionModel/*`

- `combustionModel` is the abstract base class.
- It stores references to mesh, thermo, turbulence, and transport models.
- Selection is done through OpenFOAM's runtime selection table.
- `combustionModel::New(...)` reads `constant/combustionProperties` and chooses
  the implementation.

Important detail:

- If the dictionary selects `fgmModel`, the factory rewrites that to
  `D4DummyModel`.
- That means `fgmModel` acts more like a user-facing alias than a directly
  instantiated class in this repo.

### 3. Dummy FGM model

`src/combustionModels/fgmModel/D4DummyModel/*`

- `D4DummyModel` is the only real active tabulated model here.
- It reads four scalar fields from the case: `Param1`, `Param2`, `Param3`,
  `Param4`.
- It computes four output fields: `Table1`, `Table2`, `Table3`, `Table4`.
- In `correct()`, it performs a 4D interpolation for every internal cell and
  every boundary face.

This is the heart of the benchmark.

### 4. Table lookup support

`src/combustionModels/tableSolver/*`

- `tableSolver` owns the loaded tables and the interpolation helpers.
- It converts normalized input coordinates into upper-bound indices and local
  interpolation positions.

### 5. Fallback model

`src/combustionModels/noCombustion/*`

- Minimal no-op model used as the default if no combustion dictionary is found.
- `correct()` does nothing.

## How configuration flows

1. The tutorial case sets `combustionModel fgmModel;` in
   `constant/combustionProperties`.
2. The factory maps `fgmModel` to `D4DummyModel`.
3. `D4DummyModel` expects the case to provide:
   - `0/Param1` through `0/Param4`
   - `constant/Table1_table` through `constant/Table4_table`
4. `tableSolver` loads those table dictionaries.
5. `correct()` interpolates the tabulated values into `Table1` through `Table4`.

## Build and run

Top-level helper scripts:

- `Allwmake`: builds the library and solver with `wmake`
- `Allclean`: cleans the library and solver artifacts

Tutorial helper scripts:

- `tutorials/D4DummySetup/Allrun`
  - symlinks a chosen table resolution into `constant/`
  - copies `0_orig` to `0`
  - runs `blockMesh`
  - optionally decomposes and runs `bwRSE4HPCFoam` under MPI
- `tutorials/D4DummySetup/Allclean`
  - removes generated times, logs, mesh, and linked tables

The tutorial is clearly the intended way to exercise the benchmark.

## Practical observations

- This code is tightly coupled to OpenFOAM conventions and build tooling.
- The solver is intentionally minimal: it is measuring setup and lookup cost,
  not solving a full transient problem.
- The interpolation path is specialized for exactly four dimensions in the
  provided implementation.
- Some tutorial metadata looks inherited from older cases (for example
  `system/controlDict` still says `application fgmFoam`), but `Allrun` actually
  launches `bwRSE4HPCFoam`.

## What to look at first

If you want to understand the repo quickly, start here:

1. `applications/solver/bwRSE4HPCFoam/bwRSE4HPCFoam.C`
2. `applications/solver/bwRSE4HPCFoam/createFields.H`
3. `src/combustionModels/combustionModel/combustionModelNew.C`
4. `src/combustionModels/fgmModel/D4DummyModel/D4DummyModel.C`
5. `src/combustionModels/tableSolver/tableSolver/tableSolver.C`
