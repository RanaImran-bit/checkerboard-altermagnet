# Development environment (`qmc-platform`)

Reproducible conda env providing the Fortran/MPI/BLAS toolchain for the QMC
code plus the Python stack for the ED reference and the validation platform.

## Why micromamba + pip

On this Apple-Silicon (arm64) machine the original plan hit two snags:

1. **System has no MPI**, and the QMC code hard-requires it
   (`include 'mpif.h'`). The env supplies `openmpi` + a matching `gfortran`.
2. **QuSpin is not on conda-forge for osx-arm64** — the `weinbe58` channel
   stops at py3.10/osx-64. PyPI *does* ship native osx-arm64/cp311 wheels
   (`quspin` 1.0.1 → compiled `quspin-extensions`), so QuSpin is installed
   via the `pip:` section of `environment.yml`.
3. **conda 4.8.5's classic solver was unusably slow** on the conda-forge
   solve. We use a bundled standalone **micromamba** binary
   (`.tools/bin/micromamba`, gitignored) whose libsolv resolver solves in
   seconds.

## Create / recreate

```bash
# one-time: fetch the micromamba binary (single file, no solve)
mkdir -p .tools && curl -Ls https://micro.mamba.pm/api/micromamba/osx-arm64/latest \
  | tar -xj -C .tools bin/micromamba

# build the env
.tools/bin/micromamba create -y -r "$HOME/anaconda3" -n qmc-platform -f environment.yml
```

## Activate (required before any build)

```bash
source tools/env.sh
```

Always use `conda activate` (via `tools/env.sh`), **not** a manual `PATH`
edit: the conda-forge compiler activation hooks set `SDKROOT` to the macOS
SDK. Without it the conda `gfortran` cannot link the C runtime
(`ld: library not found for -lm`).

## Verified toolchain (2026-06-23)

- Python 3.11.15 (arm64), numpy 2.4.6, scipy, numba 0.65.1
- `mpif90` / `mpirun` (openmpi), `gfortran` 14.3.0, openblas
- QuSpin 1.0.1 native arm64
- MPI Fortran compiles + runs on 2 ranks; OpenBLAS `dgemm` from Fortran OK
