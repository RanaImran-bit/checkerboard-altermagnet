#!/usr/bin/env python3
"""Reference for the two-orbital (d_xz, d_yz) altermagnetic Hubbard model that
code/src implements (see code/src/t.f90 :: Initlatt).

Hopping (spin-independent, identical to the Fortran real-space construction):
  eps_x (k) = -2 t1 cos kx - 2 t2 cos ky - 4 t3 cos kx cos ky      (orbital d_xz)
  eps_y (k) = -2 t2 cos kx - 2 t1 cos ky - 4 t3 cos kx cos ky      (orbital d_yz)
  eps_xy(k) = -4 t4 sin kx sin ky                                  (inter-orbital)
Anisotropy t1 != t2 is the altermagnetic d-wave splitting.

Two references:
  * free_fermion_energy:  U=0 exact GS energy = sum of lowest NE single-particle
    levels (per spin). No QuSpin -> fast, no OpenMP issues. The unambiguous
    target for debugging the CPQMC propagation/measurement at zero interaction.
  * ed_energy (optional, --U/--uxy/--v): full interacting ED via QuSpin on a
    small cluster.

Emits the common results schema (tools/parse_out.py / compare.py compatible).

Usage:
  python ed/altermagnet_ed.py --lx 4 --ly 4 --nup 16 --ndn 16 \
      --t1 -1 --t2 -1 --t3 -1 --t4 -1            # free fermion (default U=0)
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import argparse, json, sys
import numpy as np


def build_hopping(lx, ly, t1, t2, t3, t4):
    """nsites x nsites one-body matrix, nsites = 2*lx*ly, orbital-blocked exactly
    as code/src/t.f90 (orbital 1 = sites 1..lxy, orbital 2 = lxy+1..2*lxy)."""
    lxy = lx * ly
    nsites = 2 * lxy
    # site index i(ix,iy,orb): orbital 1 first (all ix,iy), then orbital 2
    def idx(ix, iy, orb):
        ix %= lx; iy %= ly
        base = (orb - 1) * lxy
        return base + ix * ly + iy            # matches Fortran loop ix outer, iy inner

    tk = np.zeros((nsites, nsites))
    for orb in (1, 2):
        for ix in range(lx):
            for iy in range(ly):
                i = idx(ix, iy, orb)
                # intra-orbital NN
                if orb == 1:
                    tk[i, idx(ix + 1, iy, 1)] -= t1
                    tk[i, idx(ix - 1, iy, 1)] -= t1
                    tk[i, idx(ix, iy + 1, 1)] -= t2
                    tk[i, idx(ix, iy - 1, 1)] -= t2
                else:
                    tk[i, idx(ix + 1, iy, 2)] -= t2
                    tk[i, idx(ix - 1, iy, 2)] -= t2
                    tk[i, idx(ix, iy + 1, 2)] -= t1
                    tk[i, idx(ix, iy - 1, 2)] -= t1
                # intra-orbital NNN (diagonals)
                tk[i, idx(ix + 1, iy + 1, orb)] -= t3
                tk[i, idx(ix - 1, iy - 1, orb)] -= t3
                tk[i, idx(ix + 1, iy - 1, orb)] -= t3
                tk[i, idx(ix - 1, iy + 1, orb)] -= t3
                # inter-orbital hybridization on diagonals
                other = 2 if orb == 1 else 1
                tk[i, idx(ix + 1, iy + 1, other)] += t4
                tk[i, idx(ix - 1, iy - 1, other)] += t4
                tk[i, idx(ix + 1, iy - 1, other)] -= t4
                tk[i, idx(ix - 1, iy + 1, other)] -= t4
    return tk


def free_fermion_energy(lx, ly, nup, ndn, t1, t2, t3, t4):
    tk = build_hopping(lx, ly, t1, t2, t3, t4)
    # symmetrize defensively (construction is already symmetric for these terms)
    w = np.linalg.eigvalsh(0.5 * (tk + tk.T))
    e_up = float(np.sort(w)[:nup].sum())
    e_dn = float(np.sort(w)[:ndn].sum())
    return e_up + e_dn, tk, w


def _idx(lx, ly, ix, iy, orb):
    lxy = lx * ly
    return (orb - 1) * lxy + (ix % lx) * ly + (iy % ly)


def ed_energy(lx, ly, nup, ndn, t1, t2, t3, t4, uxx, uxy, v):
    """Interacting ED of the two-orbital model via QuSpin (small clusters).

    Interaction Hamiltonian (matching code/src InitV / Vee, with n = n_up+n_dn):
      uxx * sum_{orbital-site a} n_{a,up} n_{a,dn}                  (on-site Hubbard)
      uxy * sum_{real site r}   n_{r,orb1} n_{r,orb2}              (inter-orbital)
      v   * sum_{+x,+y bonds}   sum_{oi,oj} s(oi,oj) n_{i,oi} n_{j,oj}
            with s = +1 if oi==oj else -1.
    """
    from quspin.operators import hamiltonian
    static, basis, nsites = _build_static(lx, ly, nup, ndn, t1, t2, t3, t4, uxx, uxy, v)
    no_check = dict(check_pcon=False, check_symm=False, check_herm=False)
    H = hamiltonian(static, [], basis=basis, dtype=np.float64, **no_check)
    w = np.linalg.eigvalsh(H.toarray())
    return float(w[0]), basis.Ns


def _build_static(lx, ly, nup, ndn, t1, t2, t3, t4, uxx, uxy, v):
    """Common-schema static operator list + basis for the two-orbital model."""
    from quspin.basis import spinful_fermion_basis_general
    lxy = lx * ly
    nsites = 2 * lxy
    tk = build_hopping(lx, ly, t1, t2, t3, t4)
    hop_pm = [[tk[i, j], i, j] for i in range(nsites) for j in range(nsites)
              if abs(tk[i, j]) > 1e-15]

    def dens_dens(coef, a, b):
        return ([["nn|", [[coef, a, b]]], ["|nn", [[coef, a, b]]],
                 ["n|n", [[coef, a, b]]], ["n|n", [[coef, b, a]]]])

    static = [["+-|", hop_pm], ["|+-", hop_pm]]
    blocks = {}
    def add(terms):
        for op, lst in terms:
            blocks.setdefault(op, []).extend(lst)

    if uxx != 0:
        add([["n|n", [[uxx, a, a] for a in range(nsites)]]])
    if uxy != 0:
        for ix in range(lx):
            for iy in range(ly):
                add(dens_dens(uxy, _idx(lx, ly, ix, iy, 1), _idx(lx, ly, ix, iy, 2)))
    if v != 0:
        for ix in range(lx):
            for iy in range(ly):
                for (jx, jy) in (((ix + 1) % lx, iy), (ix, (iy + 1) % ly)):
                    for oi in (1, 2):
                        for oj in (1, 2):
                            s = v if oi == oj else -v
                            add(dens_dens(s, _idx(lx, ly, ix, iy, oi), _idx(lx, ly, jx, jy, oj)))
    static += [[op, lst] for op, lst in blocks.items()]
    basis = spinful_fermion_basis_general(nsites, Nf=(nup, ndn))
    return static, basis, nsites


def ed_correlations(lx, ly, nup, ndn, t1, t2, t3, t4, uxx, uxy, v):
    """Ground-state equal-time observables for the two-orbital model:
      green_up/green_dn : G^s_ij = <c^+_{i,s} c_{j,s}>
      nn                : charge correlation <n_i n_j>   (n = n_up + n_dn)
      szsz              : spin correlation <S^z_i S^z_j>  (S^z = (n_up - n_dn)/2)
    Diagonal/density correlations come from the (diagonal) number operators;
    off-diagonal G from the hopping operators (small clusters only)."""
    from quspin.operators import hamiltonian
    static, basis, nsites = _build_static(lx, ly, nup, ndn, t1, t2, t3, t4, uxx, uxy, v)
    nc = dict(check_pcon=False, check_symm=False, check_herm=False)
    H = hamiltonian(static, [], basis=basis, dtype=np.float64, **nc)
    w, V = np.linalg.eigh(H.toarray())
    psi = V[:, 0]; p2 = psi ** 2
    # diagonal occupation vectors n_{i,s}(state) for fast density correlations
    du = [np.real(hamiltonian([["n|", [[1.0, i]]]], [], basis=basis, dtype=np.float64, **nc).diagonal())
          for i in range(nsites)]
    dd = [np.real(hamiltonian([["|n", [[1.0, i]]]], [], basis=basis, dtype=np.float64, **nc).diagonal())
          for i in range(nsites)]
    nn = np.zeros((nsites, nsites)); szsz = np.zeros((nsites, nsites))
    Gu = np.zeros((nsites, nsites)); Gd = np.zeros((nsites, nsites))
    for i in range(nsites):
        ni = du[i] + dd[i]; szi = 0.5 * (du[i] - dd[i])
        Gu[i, i] = float((p2 * du[i]).sum()); Gd[i, i] = float((p2 * dd[i]).sum())
        for j in range(nsites):
            nj = du[j] + dd[j]; szj = 0.5 * (du[j] - dd[j])
            nn[i, j] = float((p2 * ni * nj).sum())
            szsz[i, j] = float((p2 * szi * szj).sum())
            if i != j:
                Hu = hamiltonian([["+-|", [[1.0, i, j]]]], [], basis=basis, dtype=np.float64, **nc)
                Hd = hamiltonian([["|+-", [[1.0, i, j]]]], [], basis=basis, dtype=np.float64, **nc)
                Gu[i, j] = float(psi @ Hu.dot(psi)); Gd[i, j] = float(psi @ Hd.dot(psi))
    return {"e0": float(w[0]), "nsites": nsites,
            "green_up": Gu.tolist(), "green_dn": Gd.tolist(),
            "nn": nn.tolist(), "szsz": szsz.tolist()}


def record(lx, ly, nup, ndn, t1, t2, t3, t4):
    nsites = 2 * lx * ly
    e0, tk, w = free_fermion_energy(lx, ly, nup, ndn, t1, t2, t3, t4)
    return {
        "code": "ed-altermagnet-free",
        "model": {"nsites": nsites, "ne": nup + ndn, "nup": nup, "ndn": ndn,
                  "lx": lx, "ly": ly, "t1": t1, "t2": t2, "t3": t3, "t4": t4,
                  "U": 0.0},
        "observables": {
            "energy_total": {"value": e0, "error": 0.0},
            "energy_kinetic": {"value": e0, "error": 0.0},     # U=0 => all kinetic
            "energy_potential": {"value": 0.0, "error": 0.0},
            "energy_per_site": {"value": e0 / nsites, "error": 0.0},
        },
        "single_particle_levels": [round(x, 8) for x in np.sort(w).tolist()],
    }


def record_interacting(lx, ly, nup, ndn, t1, t2, t3, t4, uxx, uxy, v):
    nsites = 2 * lx * ly
    e0, ns = ed_energy(lx, ly, nup, ndn, t1, t2, t3, t4, uxx, uxy, v)
    return {
        "code": "ed-altermagnet",
        "model": {"nsites": nsites, "ne": nup + ndn, "nup": nup, "ndn": ndn,
                  "lx": lx, "ly": ly, "t1": t1, "t2": t2, "t3": t3, "t4": t4,
                  "U": uxx, "uxy": uxy, "v": v},
        "observables": {
            "energy_total": {"value": e0, "error": 0.0},
            "energy_per_site": {"value": e0 / nsites, "error": 0.0},
        },
        "ed_basis_dim": ns,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4)
    ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=16)
    ap.add_argument("--ndn", type=int, default=16)
    ap.add_argument("--t1", type=float, default=-1.0)
    ap.add_argument("--t2", type=float, default=-1.0)
    ap.add_argument("--t3", type=float, default=-1.0)
    ap.add_argument("--t4", type=float, default=-1.0)
    ap.add_argument("--uxx", type=float, default=0.0)
    ap.add_argument("--uxy", type=float, default=0.0)
    ap.add_argument("--v", type=float, default=0.0)
    ap.add_argument("--corr", action="store_true", help="output equal-time GF + spin/charge correlations")
    ap.add_argument("-o", "--out")
    a = ap.parse_args()
    if a.corr:
        rec = ed_correlations(a.lx, a.ly, a.nup, a.ndn, a.t1, a.t2, a.t3, a.t4,
                              a.uxx, a.uxy, a.v)
    elif a.uxx or a.uxy or a.v:
        rec = record_interacting(a.lx, a.ly, a.nup, a.ndn, a.t1, a.t2, a.t3, a.t4,
                                 a.uxx, a.uxy, a.v)
    else:
        rec = record(a.lx, a.ly, a.nup, a.ndn, a.t1, a.t2, a.t3, a.t4)
    js = json.dumps(rec, indent=2)
    if a.out:
        from pathlib import Path
        Path(a.out).write_text(js); print(f"wrote {a.out}")
    else:
        print(js)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
