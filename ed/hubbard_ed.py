#!/usr/bin/env python3
"""Exact-diagonalization reference for the single-band Hubbard model on an
Lx x Ly square lattice with periodic boundaries, using QuSpin.

    H = -t sum_<ij>,s (c^+_{is} c_{js} + h.c.) + U sum_i n_{i up} n_{i dn}

Emits the SAME common schema as tools/parse_out.py so ED and QMC runs are
directly comparable by the platform/UI. ED values are exact -> error = 0.

This is the ground-truth reference for small clusters (the CPQMC benchmark
runs Nup=Ndn=1 on 4x4 = 256 states, trivially exact here).

Usage:
    python ed/hubbard_ed.py --lx 4 --ly 4 --nup 1 --ndn 1 --U 3 --t 1
    python ed/hubbard_ed.py ... -o ed.json
"""
from __future__ import annotations
import os
# must precede numpy/quspin import: conda openblas (libgomp) + quspin core
# (llvm libomp) clash; force single-thread (deterministic) and allow load.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import argparse, json, sys
import numpy as np
from quspin.basis import spinful_fermion_basis_general
from quspin.operators import hamiltonian


def square_bonds(lx: int, ly: int):
    """Nearest-neighbor bonds (+x and +y, each counted once) with PBC.
    Site index i = x + lx*y."""
    bonds = []
    for y in range(ly):
        for x in range(lx):
            i = x + lx * y
            jx = ((x + 1) % lx) + lx * y          # +x neighbor
            jy = x + lx * ((y + 1) % ly)          # +y neighbor
            bonds.append((i, jx))
            if ly > 1:
                bonds.append((i, jy))
    return bonds


def build(lx, ly, nup, ndn, t, U):
    N = lx * ly
    basis = spinful_fermion_basis_general(N, Nf=(nup, ndn))
    bonds = square_bonds(lx, ly)

    # hopping -t (c^+_i c_j + h.c.); QuSpin "+-|" = up species, "|+-" = down
    hop_pm = [[-t, i, j] for (i, j) in bonds]      # c^+_i c_j
    hop_mp = [[+t, i, j] for (i, j) in bonds]      # c_i c^+_j  (h.c. partner)
    interaction = [[U, i, i] for i in range(N)]    # U n_up n_dn

    no_check = dict(check_pcon=False, check_symm=False, check_herm=False)
    static_K = [["+-|", hop_pm], ["-+|", hop_mp],
                ["|+-", hop_pm], ["|-+", hop_mp]]
    static_V = [["n|n", interaction]]
    H = hamiltonian(static_K + static_V, [], basis=basis, dtype=np.float64, **no_check)
    K = hamiltonian(static_K, [], basis=basis, dtype=np.float64, **no_check)
    V = hamiltonian(static_V, [], basis=basis, dtype=np.float64, **no_check)
    return basis, H, K, V


def ground_state(lx, ly, nup, ndn, t, U):
    basis, H, K, V = build(lx, ly, nup, ndn, t, U)
    Hd = H.toarray()
    w, v = np.linalg.eigh(Hd)
    psi = v[:, 0]
    e0 = float(w[0])
    ke = float(psi @ (K.toarray() @ psi))
    pe = float(psi @ (V.toarray() @ psi))
    return dict(e0=e0, ke=ke, pe=pe, ns=basis.Ns)


def record(lx, ly, nup, ndn, t, U):
    g = ground_state(lx, ly, nup, ndn, t, U)
    nsites = lx * ly
    return {
        "code": "ed-quspin",
        "model": {"nsites": nsites, "ne": nup + ndn, "nup": nup, "ndn": ndn,
                  "U": U, "t0": t, "lx": lx, "ly": ly},
        "observables": {
            "energy_total": {"value": g["e0"], "error": 0.0},
            "energy_kinetic": {"value": g["ke"], "error": 0.0},
            "energy_potential": {"value": g["pe"], "error": 0.0},
            "energy_per_site": {"value": g["e0"] / nsites, "error": 0.0},
        },
        "ed_basis_dim": g["ns"],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4)
    ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=1)
    ap.add_argument("--ndn", type=int, default=1)
    ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--U", type=float, default=3.0)
    ap.add_argument("-o", "--out")
    a = ap.parse_args()
    rec = record(a.lx, a.ly, a.nup, a.ndn, a.t, a.U)
    js = json.dumps(rec, indent=2)
    if a.out:
        from pathlib import Path
        Path(a.out).write_text(js); print(f"wrote {a.out}")
    else:
        print(js)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
