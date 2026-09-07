#!/usr/bin/env python -u
"""Exact diagonalization of the 12-site Lieb Hubbard cluster.

Settles whether the altermagnetic order the 48-site CPQMC reports is physical or
an artifact of the constrained-path trial. The scan gave m_B - m_C ~ 1.9 with an
altermagnetic trial and 0.12 to 0.27 with paramagnetic or free-electron trials,
a factor of 7 to 16, so neither number can be trusted on its own. ED on the
cluster the reference work itself used (2x2 cells, 12 sites) is exact and settles
it.

IMPORTANT, and the reason this is not a moment comparison: on a finite cluster
the exact ground state is a spin symmetry eigenstate, so <n_i,up - n_i,dn> = 0
identically. The order parameter must be read from its FLUCTUATION. Define

    O = (1/N_B) sum_{i in B} (n_i,up - n_i,dn)  -  (1/N_C) sum_{i in C} (...)

In a symmetry-broken state <O> = m_B - m_C, which is what CPQMC reports. In the
exact eigenstate <O> = 0 while sqrt(<O^2>) carries the same magnitude. O is
diagonal in the occupation basis, so both are cheap once the ground state is in
hand.

The energy is the second, cleaner benchmark: it is a variational upper bound for
CPQMC and needs no symmetry argument at all.

    python -u ed_lieb.py <t2> [U] [L]
"""
import sys, os, time, itertools
import numpy as np
from scipy.sparse import csr_matrix, kron, identity, diags
from scipy.sparse.linalg import eigsh
import importlib.util

t2 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3
U  = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
L  = int(sys.argv[3])   if len(sys.argv) > 3 else 4          # fine grid, cells = L/2
T0 = -1.0

here = os.path.dirname(os.path.abspath(__file__))
srcs = [os.path.join(here, "fortran", "wf_lieb.py"), os.path.join(here, "wf_lieb.py")]
src  = next(p for p in srcs if os.path.exists(p))
body = open(src).read().split("# ================================================================\n# Main")[0]
open("/tmp/_wfl_ed.py", "w").write("EPS_A = 0.0\n" + body)
spec = importlib.util.spec_from_file_location("_wfl_ed", "/tmp/_wfl_ed.py")
w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)

sites = w.lieb_sites(L, L)
NS = len(sites)
lab = ["A" if (ix % 2 and iy % 2) else ("B" if iy % 2 else "C") for ix, iy in sites]
K = w.GetK(L, L, 1, T0, t2)
NE = 2 * (L // 2) ** 2                      # n = 4 electrons per cell, compensated
print(f"Lieb ED: L={L} ({NS} sites: {lab.count('A')}A {lab.count('B')}B {lab.count('C')}C)")
print(f"  U={U} on B and C only, t0={T0}, t2={t2}, NUP=NDN={NE}")

# ---- single-spin sector -----------------------------------------------------
occs = [sum(1 << p for p in c) for c in itertools.combinations(range(NS), NE)]
occs.sort()
index = {o: i for i, o in enumerate(occs)}
D = len(occs)
print(f"  sector dim {D}, full space {D*D}")

rows, cols, vals = [], [], []
for a, o in enumerate(occs):
    for i in range(NS):
        if not (o >> i) & 1:
            continue
        for j in range(NS):
            if K[i, j] == 0.0 or (o >> j) & 1:
                continue
            # c_j^dag c_i : fermion sign from the occupations strictly between
            o1 = o & ~(1 << i)
            lo, hi = (i, j) if i < j else (j, i)
            mask = ((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1)
            sgn = -1.0 if bin(o1 & mask).count("1") % 2 else 1.0
            rows.append(index[o1 | (1 << j)]); cols.append(a); vals.append(sgn * K[i, j])
T = csr_matrix((vals, (rows, cols)), shape=(D, D))
print(f"  hopping nnz {T.nnz}")

# ---- interaction: U on B and C only, diagonal in the occupation basis --------
bc = np.array([i for i, t in enumerate(lab) if t != "A"])
nbc = np.array([[(o >> int(i)) & 1 for i in bc] for o in occs], dtype=np.int8)
t0 = time.time()
Hd = U * (nbc[:, None, :] * nbc[None, :, :]).sum(-1).ravel().astype(np.float64)
print(f"  interaction diagonal built in {time.time()-t0:.1f}s")

I = identity(D, format="csr")
H = kron(T, I, format="csr") + kron(I, T, format="csr") + diags(Hd)
print(f"  H: {H.shape[0]} x {H.shape[0]}, nnz {H.nnz}")

t0 = time.time()
E, V = eigsh(H, k=1, which="SA", maxiter=20000, tol=1e-10)
psi = V[:, 0]; E0 = float(E[0])
print(f"  ground state in {time.time()-t0:.1f}s   E0 = {E0:.8f}   E0/site = {E0/NS:.8f}")

# ---- the altermagnetic order parameter, from its fluctuation ----------------
B = np.array([i for i, t in enumerate(lab) if t == "B"])
C = np.array([i for i, t in enumerate(lab) if t == "C"])
def sub_occ(idx):
    return np.array([[(o >> int(i)) & 1 for i in idx] for o in occs], dtype=np.int8).sum(1)
nB, nC = sub_occ(B), sub_occ(C)
# O = <n_up - n_dn> averaged over B, minus the same over C
O = (nB[:, None] - nB[None, :]) / len(B) - (nC[:, None] - nC[None, :]) / len(C)
O = O.ravel().astype(np.float64)
p = psi ** 2
mean_O = float((p * O).sum()); rms_O = float(np.sqrt((p * O ** 2).sum()))
print(f"\n  <O>        = {mean_O:+.6f}   (must be ~0: the exact state is a symmetry eigenstate)")
print(f"  sqrt(<O^2>) = {rms_O:.6f}   <-- compare with CPQMC m_B - m_C")
# a uniform paramagnet still has a nonzero rms from pure spin fluctuation; quote
# the U = 0 value of the same quantity as the floor when interpreting this.
print(f"\n  cluster: {NS} sites, {NE}+{NE} electrons, n = {2*NE/NS:.4f} per site "
      f"= {2*NE/((L//2)**2):.0f} per unit cell")
