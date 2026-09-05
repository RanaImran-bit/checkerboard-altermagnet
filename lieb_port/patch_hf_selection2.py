# Compute the trial's variational energy from the ORBITALS, not from the
# eigenvalue sum.
#
# The first attempt used  E_HF = sum(eps) - sum_i U_i n_up n_dn.  That identity
# only holds at true self consistency. The SCF mixes densities
# (nup <- (1-a)*nup + a*new) and stops when the eigenvalue sum settles to 1e-4,
# so the density used to BUILD H is not the density carried by the eigenvectors,
# and for states that never really converge the two disagree wildly. Measured at
# 12 sites, U=10, t2=0.3, the paramagnetic seed gave -11.13 from that identity
# against a true +31.87, so it still won the comparison and still poisoned the
# choice.
#
# Taking the energy straight from the orbitals removes the assumption entirely:
#
#   E = Tr(K rho_up) + Tr(K rho_dn) + sum_i U_i <n_i,up> <n_i,dn>
#   rho_sigma = V[:, :N_sigma] V[:, :N_sigma]^T
#
# K_up already carries eps_A on the A diagonal, so it is the full one-body part.
import sys

wf = open("wf_lieb.py").read()
old = """    # Hartree-Fock energy: undo the double counting of the interaction, using the
    # self-consistent densities rather than the mixed ones.
    n_u = Getnup(eig_vup, NUP)
    n_d = Getndn(eig_vdn, NDN)
    E_hf = E - float(np.sum(U * n_u * n_d))
    return nup, ndn, E_hf, eig_vup, eig_vdn"""
new = """    # Variational energy of the determinant the orbitals actually describe.
    # Deliberately NOT E - sum(U n_up n_dn): that identity assumes exact self
    # consistency, which the mixed-density SCF does not reach for states that
    # only converge loosely, and it then flatters exactly those states.
    r_u = eig_vup[:, :NUP] @ eig_vup[:, :NUP].T
    r_d = eig_vdn[:, :NDN] @ eig_vdn[:, :NDN].T
    n_u = Getnup(eig_vup, NUP)
    n_d = Getndn(eig_vdn, NDN)
    E_hf = float(np.sum(K_up * r_u) + np.sum(K_up * r_d) + np.sum(U * n_u * n_d))
    return nup, ndn, E_hf, eig_vup, eig_vdn"""
if wf.count(old) != 1:
    sys.exit(f"expected 1 occurrence, found {wf.count(old)}")
open("wf_lieb.py", "w").write(wf.replace(old, new))
print("trial energy now computed from the orbitals")
