#!/usr/bin/env python3
"""Single-step reference for the Fortran Vee two-site HS update (code/src/Vee.f90).

For a walker (phi_up, phi_dn) and trial (psiT_up, psiT_dn), one Vee call applies
an interaction channel V*n_{a,si}*n_{b,sj} by scaling row a of the si-determinant
by expV1 and row b of the sj-determinant by expV2, for a chosen HS field ising.

This script computes, for each (channel spins, ising):
  - EXACT overlap ratio  det(psiT^T phi')/det(psiT^T phi)  by literally scaling rows
  - Vee's fast formula    rdet1*rdet2  (Sherman-Morrison; sequential if same-spin)
and reports the max discrepancy. If they disagree, the Vee GF-ratio update is the
bug for i!=j (uxx uses i==j, which is a special case that already validates).

Run:  python pyqmc/vee_check.py
"""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
from altermagnet_ed import build_hopping

np.random.seed(0)
lx = ly = 2
lxy = lx * ly
n = 2 * lxy              # orbital-sites
nup = ndn = 4
dt = 0.01
V = 2.0                  # uxy coupling (repulsive)

# trial = lowest orbitals of the two-orbital hopping
K = build_hopping(lx, ly, -1, -1, -1, -1)
w, vecs = np.linalg.eigh(K)
psiT_up = vecs[:, :nup].copy()
psiT_dn = vecs[:, :ndn].copy()
# walker: trial propagated a bit (random unitary mix in occupied space) so phi != psiT
phi_up = (np.eye(n) + 0.1 * np.random.randn(n, n)) @ psiT_up
phi_dn = (np.eye(n) + 0.1 * np.random.randn(n, n)) @ psiT_dn

# HS constants (repulsive branch, MkExpV): cosh(alpha) = exp(dt*V/2)
alpha = np.arccosh(np.exp(0.5 * dt * V))
def expV1(s):  # site a (operator role 1)
    return np.exp(alpha * s - 0.5 * dt * V)
def expV2(s):  # site b (operator role 2)
    return np.exp(-alpha * s - 0.5 * dt * V)

def overlap(psiT, phi):
    return np.linalg.det(psiT.T @ phi)

def green_diag(psiT, phi, k):
    # G_kk = [phi (psiT^T phi)^{-1} psiT^T]_kk
    A = psiT.T @ phi
    G = phi @ np.linalg.solve(A, psiT.T)
    return G[k, k], G

def scaled_overlap(psiT, phi, k, b):
    phi2 = phi.copy(); phi2[k, :] *= b
    return overlap(psiT, phi2), phi2

# interaction site pair: inter-orbital on-site (a, b) = (0, lxy)  [0-indexed]
a, b = 0, lxy
print(f"two-orbital 2x2, uxy on-site pair (a={a}, b={b}), V={V}, dt={dt}, alpha={alpha:.5f}\n")
print(f"{'channel(si,sj)':16s}{'ising':>6}{'exact ratio':>16}{'Vee rdet1*rdet2':>18}{'|diff|':>12}")

maxdiff = 0.0
for (si, sj) in [(1, 1), (1, 2), (2, 1), (2, 2)]:   # spinlsi, spinlsj (1=up,2=dn)
    pu0, pd0 = phi_up.copy(), phi_dn.copy()
    ov0 = overlap(psiT_up, pu0) * overlap(psiT_dn, pd0)
    same_spin = (si == sj)
    for s in (1, -1):
        # ---- EXACT: scale row a of si-det by expV1, row b of sj-det by expV2 ----
        pu, pd = pu0.copy(), pd0.copy()
        if si == 1: pu[a, :] *= expV1(s)
        else:       pd[a, :] *= expV1(s)
        if sj == 1: pu[b, :] *= expV2(s)
        else:       pd[b, :] *= expV2(s)
        ovN = overlap(psiT_up, pu) * overlap(psiT_dn, pd)
        exact = ovN / ov0

        # ---- Vee fast formula ----
        psiT_i = psiT_up if si == 1 else psiT_dn
        phi_i = pu0 if si == 1 else pd0
        psiT_j = psiT_up if sj == 1 else psiT_dn
        phi_j = pu0 if sj == 1 else pd0
        d1 = expV1(s) - 1.0
        d2 = expV2(s) - 1.0
        # Fortran builds g = inverse overlap (psiT^T phi)^{-1}; G_kk = phi(k,:) g psiT(k,:)
        gi = np.linalg.inv(psiT_i.T @ phi_i)                  # N x N
        Gaa = phi_i[a, :] @ gi @ psiT_i[a, :]
        rdet1 = 1.0 + d1 * Gaa
        if not same_spin:
            gj = np.linalg.inv(psiT_j.T @ phi_j)
            Gbb = phi_j[b, :] @ gj @ psiT_j[b, :]
            rdet2 = 1.0 + d2 * Gbb
        else:
            # Fortran same-spin: Sherman-Morrison on the inverse-overlap g_i
            gl_a = gi @ psiT_i[a, :]        # g . phiT(a,:)
            gr_a = phi_i[a, :] @ gi        # phi(a,:) . g
            g_temp = gi - (d1 / rdet1) * np.outer(gl_a, gr_a)
            Gbb_p = phi_i[b, :] @ g_temp @ psiT_i[b, :]
            rdet2 = 1.0 + d2 * Gbb_p
        vee = rdet1 * rdet2
        diff = abs(exact - vee); maxdiff = max(maxdiff, diff)
        print(f"({si},{sj}){'':10s}{s:>6}{exact:>16.6f}{vee:>18.6f}{diff:>12.2e}")

print(f"\nMAX |exact - Vee| = {maxdiff:.3e}  -> {'OK (formula matches)' if maxdiff < 1e-8 else 'BUG in Vee rdet/GF update'}")
