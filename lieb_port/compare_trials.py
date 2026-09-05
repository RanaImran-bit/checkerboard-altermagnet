# Compare the two trial wave functions by their actual variational energy.
#
# 251 and 113new ran the same parameter point and measured different physics.
# CPQMC is constrained path, so the trial decides the answer, and the two runs
# ended up with different trials because Run() picks the seed with the smallest
# SUM OF EIGENVALUES rather than the smallest Hartree-Fock energy. Those rank
# states differently, so the criterion can prefer the worse state.
#
# For a Slater determinant the mean-field energy is
#       E = sum_sigma Tr(K rho_sigma) + sum_i U_i <n_i,up> <n_i,dn>
# which is what this compares. Lower wins.
import sys, numpy as np
sys.path.insert(0, "/tmp/ck")
from wfmod import GetK, hubbard_u, lieb_sites

L, U, t2, t0 = 4, 10.0, 0.3, -1.0
K  = GetK(L, L, 1, t0, t2)
Uv = hubbard_u(L, L, U)
sites = lieb_sites(L, L)
lab = ["A" if (ix % 2 and iy % 2) else ("B" if iy % 2 else "C") for ix, iy in sites]

def report(tag, fu, fd):
    # wfup.txt holds phiT_up transposed: one row per occupied orbital
    pu = np.loadtxt(fu); pd = np.loadtxt(fd)
    if pu.ndim == 1: pu = pu[None, :]
    if pd.ndim == 1: pd = pd[None, :]
    pu, pd = pu.T, pd.T                       # -> (nsites, norb)
    rho_u = pu @ pu.T; rho_d = pd @ pd.T      # orthonormal columns
    nu, nd = np.diag(rho_u).copy(), np.diag(rho_d).copy()
    E = np.sum(K * rho_u) + np.sum(K * rho_d) + float(np.sum(Uv * nu * nd))
    m = nu - nd
    mm = {s: float(np.mean([m[i] for i, t in enumerate(lab) if t == s])) for s in "ABC"}
    print(f"{tag:>10}{E:>14.6f}{mm['A']:>10.4f}{mm['B']:>10.4f}{mm['C']:>10.4f}"
          f"{mm['B']-mm['C']:>12.4f}{nu.sum()+nd.sum():>9.1f}")
    return E

print(f"{'trial':>10}{'E_var':>14}{'m_A':>10}{'m_B':>10}{'m_C':>10}{'m_B-m_C':>12}{'Nelec':>9}")
e251 = report("251", "wf251_up.txt", "wf251_dn.txt")
e113 = report("113new", "wf113_up.txt", "wf113_dn.txt")
print(f"\n  difference E(113new) - E(251) = {e113-e251:+.6f}")
print("  -> 251's trial is variationally BETTER" if e251 < e113 else
      "  -> 113new's trial is variationally BETTER")
