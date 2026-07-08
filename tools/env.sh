#!/usr/bin/env bash
# Activate the qmc-platform conda env for building/running the Fortran QMC
# code and the Python ED/platform.
#
# IMPORTANT: source this (don't just put the env bin/ on PATH). Proper
# `conda activate` runs the conda-forge compiler activation hooks that set
# SDKROOT / sysroot; without them the conda gfortran fails to link against
# the macOS SDK (e.g. "library not found for -lm").
#
#   source tools/env.sh
#
# The env was created with the bundled micromamba (.tools/bin/micromamba)
# from environment.yml. To (re)create it:
#   .tools/bin/micromamba create -y -r "$HOME/anaconda3" -n qmc-platform -f environment.yml

# NB: do NOT `set -u` here — the conda-forge compiler activation scripts
# reference unset vars (e.g. GFORTRAN) and would abort under nounset.
_ANACONDA="${ANACONDA_HOME:-$HOME/anaconda3}"
# shellcheck disable=SC1091
source "${_ANACONDA}/etc/profile.d/conda.sh"
conda activate qmc-platform

echo "qmc-platform active: python=$(python --version 2>&1 | awk '{print $2}')" \
     "mpif90=$(command -v mpif90)" "CONDA_PREFIX=${CONDA_PREFIX}"
