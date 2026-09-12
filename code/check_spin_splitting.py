"""Does our model have momentum-domain spin polarisation at all?

The PRL defines the altermagnetic order parameter as

    Delta_tot = sum_k |Delta n(k)| = sum_k |n_up(k) - n_dn(k)|          (Eq. 2)

with the ABSOLUTE VALUE inside the sum, so it is nonzero even though the NET
magnetisation sum_k [n_up(k) - n_dn(k)] vanishes. That is exactly the signature
of an altermagnet: compensated, but spin-split in momentum space.

Their Hamiltonian (Eq. 1) uses SPIN-DEPENDENT nearest-neighbour hopping,
    t_x,up = t + tA,  t_y,up = t - tA,   t_x,dn = t - tA,  t_y,dn = t + tA
so the two spins see different dispersions and Delta_tot is nonzero by
construction.

OUR checkerboard hopping is SPIN-INDEPENDENT (checkerboard.py: "nsites x nsites
checkerboard one-body matrix (spin-independent)"), and every run uses
nup = ndn. This script computes Delta_tot for both.
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import checkerboard as cb

L, NUP = 12, 72          # half filling
N = L * L


def occ_from_K(K, nup):
    """Non-interacting momentum occupations: fill the lowest nup levels and read
    off the occupation of each k-state."""
    e, v = np.linalg.eigh(K)
    occ = np.zeros(N)
    for m in range(nup):
        occ += np.abs(v[:, m]) ** 2      # site-basis density of the filled levels
    return occ


def delta_tot_ours(delta):
    """Our checkerboard: one spin-independent K, so K_up = K_dn exactly."""
    K = cb.checkerboard_hopping(L, L, -1.0, 0.3, -delta)
    nk_up = occ_from_K(K, NUP)
    nk_dn = occ_from_K(K, NUP)           # identical matrix -> identical result
    return np.abs(nk_up - nk_dn).sum()


def delta_tot_prl(tA):
    """The PRL model: spin-dependent NN hopping, evaluated in momentum space
    where it is diagonal."""
    k = 2 * np.pi * np.arange(L) / L
    kx, ky = np.meshgrid(k, k, indexing="ij")
    e_up = -2 * ((1 + tA) * np.cos(kx) + (1 - tA) * np.cos(ky))
    e_dn = -2 * ((1 - tA) * np.cos(kx) + (1 + tA) * np.cos(ky))
    # fill the lowest NUP states for each spin independently
    def fill(e):
        f = np.zeros(e.size)
        f[np.argsort(e.ravel())[:NUP]] = 1.0
        return f
    return np.abs(fill(e_up) - fill(e_dn)).sum(), np.abs(e_up - e_dn).max()


print(f"L={L}, half filling (nup=ndn={NUP}), U=0\n")
print("OUR checkerboard model  (delta on DIAGONAL bonds, spin-INdependent):")
print(f"{'delta':>7}{'Delta_tot':>14}{'max|e_up-e_dn|':>18}")
for d in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:
    K = cb.checkerboard_hopping(L, L, -1.0, 0.3, -d)
    # K_up and K_dn are literally the same array -> spin splitting is zero
    print(f"{d:7.1f}{delta_tot_ours(d):14.2e}{0.0:18.2e}")

print("\nPRL model  (tA on NEAREST-NEIGHBOUR bonds, spin-DEPENDENT):")
print(f"{'tA':>7}{'Delta_tot':>14}{'max|e_up-e_dn|':>18}")
for tA in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:
    dt, split = delta_tot_prl(tA)
    print(f"{tA:7.1f}{dt:14.1f}{split:18.3f}")

print("\ne_up(k) - e_dn(k) = -4*tA*(cos kx - cos ky)  for the PRL model")
print("                  =  0                        for ours, at every k and every delta")
