# Does the trial selection now pick the variationally best state?
#
# Iteration() used to return the SUM OF MEAN-FIELD EIGENVALUES and Run() picked
# the seed minimising it. That sum double counts the interaction, so it ranks
# states differently from the real Hartree-Fock energy, and on the Lieb lattice
# the two disagree about the winner. This recomputes the variational energy
# straight from the converged orbitals,
#
#       E = Tr(K rho_up) + Tr(K rho_dn) + sum_i U_i <n_i,up> <n_i,dn>
#
# independently of whatever Iteration() reports, and checks that the value the
# code now returns matches it and that the lowest one is the state we expect.
import sys, numpy as np
sys.path.insert(0, "/tmp/ck")
from wfmod import Iteration, GetK, hubbard_u, lieb_sites, Getnup, Getndn

L, U, t2, t0 = 4, 10.0, 0.3, -1.0
NUP = NDN = 8
K  = GetK(L, L, 1, t0, t2)
Uv = hubbard_u(L, L, U)
sites = lieb_sites(L, L)
lab = ["A" if (ix % 2 and iy % 2) else ("B" if iy % 2 else "C") for ix, iy in sites]

print(f"L={L} (12 sites)  U={U} on B and C only  t2={t2}  NUP=NDN={NUP}")
print(f"{'seed':>8}{'returned E':>13}{'E_var':>13}{'match':>7}"
      f"{'m_A':>9}{'m_B':>9}{'m_C':>9}{'m_B-m_C':>10}")
rows = []
for seed in ["am", "neel", "pm", "random", "random"]:
    nup, ndn, E, vu, vd = Iteration(L, L, U, t2, t0, NUP=NUP, NDN=NDN, seed=seed)
    ru = vu[:, :NUP] @ vu[:, :NUP].T
    rd = vd[:, :NDN] @ vd[:, :NDN].T
    nu, nd = np.diag(ru).copy(), np.diag(rd).copy()
    Evar = np.sum(K * ru) + np.sum(K * rd) + float(np.sum(Uv * nu * nd))
    m = nu - nd
    mm = {s: float(np.mean([m[i] for i, t in enumerate(lab) if t == s])) for s in "ABC"}
    ok = "yes" if abs(E - Evar) < 1e-6 else "NO"
    print(f"{seed:>8}{E:>13.6f}{Evar:>13.6f}{ok:>7}"
          f"{mm['A']:>9.4f}{mm['B']:>9.4f}{mm['C']:>9.4f}{mm['B']-mm['C']:>10.4f}")
    rows.append((seed, E, Evar, mm))

win = min(rows, key=lambda r: r[2])
print(f"\n  variationally lowest: {win[0]}   E_var = {win[2]:.6f}   "
      f"m_B-m_C = {win[3]['B']-win[3]['C']:+.4f}")
print(f"  code would choose   : {min(rows, key=lambda r: r[1])[0]}")
print("\n  the code's choice IS the variational minimum"
      if min(rows, key=lambda r: r[1])[0] == win[0] else
      "\n  *** MISMATCH: the code does not pick the variational minimum ***")
