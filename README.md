# ndtblIntegrationFoam

[![License](https://img.shields.io/github/license/bwRSE4HPC/ndtbl-openfoam-integration)](https://opensource.org/license/gpl-3.0)
[![Build](https://github.com/bwRSE4HPC/ndtbl-openfoam-integration/actions/workflows/ci.yml/badge.svg)](https://github.com/bwRSE4HPC/ndtbl-openfoam-integration/actions)
[![pre-commit.ci](https://results.pre-commit.ci/badge/github/bwRSE4HPC/ndtbl-openfoam-integration/main.svg)](https://results.pre-commit.ci/latest/github/bwRSE4HPC/ndtbl-openfoam-integration/main)

OpenFOAM benchmark application for multidimensional table lookup with [`ndtbl`](https://github.com/bwRSE4HPC/ndtbl).

The application loads synthetic FGM tables, evaluates them in an OpenFOAM time loop, and reports initialization time, lookup time, memory use, and memory-map residency. It is a benchmark and integration example, not a general-purpose CFD solver.

## Requirements

- Linux
- OpenFOAM Foundation v10
- Python 3.11 or newer
- MPI for parallel runs (normally provided with OpenFOAM)

## Build

```bash
git clone --recurse-submodules https://github.com/bwrse4hpc/ndtbl-openfoam-integration.git
cd ndtbl-openfoam-integration

# Load your OpenFOAM v10 environment first.
python3 -m pip install ./ext/ndtbl/python/ndtbl
./Allwmake
```

Build products are written to `FOAM_USER_LIBBIN` and `FOAM_USER_APPBIN`.

## Run the example

```bash
cd tutorials/D3TransCube
TABLERESOLUTION=21 PARALLEL=false ./Allrun
```

`Allrun` creates the mesh and synthetic table, runs `ndtblIntegrationFoam`, and writes `log.ndtblIntegrationFoam`. Use `./Allclean` to remove generated case data.

Common settings are passed as environment variables:

| Variable | Default | Purpose |
| --- | --- | --- |
| `TABLERESOLUTION` | `101` | Points per table axis |
| `TABLE_SUFFIX` | `f` | Values stored as `f` (float32) or `d` (float64) |
| `TABLE_LAYOUT` | `combined` | One combined file or `split` files |
| `PARALLEL` | `true` | Enable an MPI run |
| `NPROCS` | `4` | Number of MPI ranks |
| `TABLE_SOURCE` | unset | Use a pre-generated table instead of generating one |
| `TMPDIR` | unset | Stage tables on node-local storage |

Example parallel run:

```bash
TABLERESOLUTION=101 TABLE_LAYOUT=split NPROCS=4 ./Allrun
```

## Benchmark sweep

```bash
cd tutorials/D3TransCube
./run_ndtbl_benchmarks.sh \
    --resolution 21 \
    --resolution 101 \
    --table-layout combined \
    --table-layout split \
    --mpi 0 \
    --mpi 4
```

Results and solver logs are written to `tutorials/D3TransCube/benchmark-results/`. Run the script with `--help` for all options.

## Memory diagnostics

Set `ndtblDiagnostics true` in `constant/combustionProperties` to write rank-level diagnostics to:

```text
postProcessing/ndtblResidency/<start-time>/residency.tsv
postProcessing/ndtblResidency/<start-time>/numa.tsv
```

`residency.tsv` retains its version-1 format. `numa.tsv` has its own `# ndtbl_numa_tsv_version=1` header and these columns:

```text
sample time label group rank host numa_available policy kernel_page_size_bytes numa_node pages
```

Both files use the same residency snapshot at each diagnostic sample. The master rank writes gathered records from all ranks. NUMA records contain one row per reported NUMA node, with page counts for the whole mapping containing the table payload. Join the files on `(sample, rank, group)` within a run; interpret NUMA node IDs together with `host`.

If NUMA information is unavailable, a placeholder row has `numa_available=0` and `nan` for policy, page size, node, and page count. An available mapping with no node counters has `numa_available=1`, retains its policy and any known page size, and uses `nan` for node and page count. An unknown kernel page size is also `nan`. Zero page counts remain zero. Raw proc lines are not written.

## Reproducing the paper example

The `D3TransCube` tutorial is the synthetic unit-cube transport case used for the memory-residency evaluation in the accompanying ndtbl paper. The published configuration uses a combined, single-precision table with `501^3` points and three fields (approximately 1.44 GiB), five MPI ranks, and a decomposition of `(5 1 1)`. The figure in the paper places ranks 0--2 on one node and ranks 3--4 on a second node.

After building the application, run the paper configuration from an allocation whose MPI placement gives the required three-rank/two-rank node split:

```bash
cd tutorials/D3TransCube
TABLERESOLUTION=501 \
TABLE_SUFFIX=f \
TABLE_LAYOUT=combined \
PARALLEL=true \
NPROCS=5 \
./Allrun
```

On a cluster, set `TMPDIR` to node-local storage if the shared table should be copied once per allocated node. In a Slurm allocation, the staging helper uses `srun` for those per-node copies. The script advises the kernel to evict the generated or staged table from the page cache before the solver starts; a completely cold cache and the resulting timing still depend on the operating system and other activity on the nodes.

The diagnostics used by the paper are written to:

```text
tutorials/D3TransCube/postProcessing/ndtblResidency/0/residency.tsv
```

The compile-time mmap configuration is recorded in `src/combustionModels/Make/options`. For the paper configuration, mmap and diagnostics are enabled, population is disabled. Generating the full table requires several GiB of free memory and disk space; the memory-lock limit must permit approximately 1.44 GiB per rank.

## Repository layout

- `applications/solver/ndtblIntegrationFoam/`: benchmark application
- `src/combustionModels/`: OpenFOAM models and `ndtbl` integration
- `tutorials/D3TransCube/`: runnable example and benchmark driver
- `scripts/`: table-staging helpers for local and Slurm runs
- `ext/ndtbl/`: pinned `ndtbl` submodule

## Citation and license

Citation metadata is provided in [`CITATION.cff`](CITATION.cff). For reproducible scholarly use, cite the archived release rather than the moving default branch.

This project is licensed under [GPL-3.0-or-later](LICENSE.md). See [NOTICE.md](NOTICE.md) for copyright information.

The [`ndtbl`](https://github.com/bwRSE4HPC/ndtbl) submodule is licensed separately under the [MIT License](ext/ndtbl/LICENSE.md). Its copyright holders are listed in [ext/ndtbl/COPYING.md](ext/ndtbl/COPYING.md).
