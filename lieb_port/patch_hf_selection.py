# Select the trial wave function by the Hartree-Fock energy, not by the sum of
# mean-field eigenvalues.
#
# Run() compared seeds using  E = sum(eps_up[:NUP]) + sum(eps_dn[:NDN]).  Each
# eigenvalue already carries the mean field of the opposite spin, so that sum
# double counts the interaction:
#
#   sum(eps_up) + sum(eps_dn) = Tr(K rho_up) + Tr(K rho_dn) + 2 sum_i U_i nup nd
#   E_HF                      = Tr(K rho_up) + Tr(K rho_dn) +   sum_i U_i nup nd
#   => E_HF = sum(eps) - sum_i U_i <n_i,up> <n_i,dn>
#
# The two rank states differently, and on the Lieb lattice they disagree about
# which state wins. Measured at 12 sites, U=10 on B and C, t2=0.3, n=4 per cell:
#
#            sum(eps)      E_HF     state
#   am        +0.0562    -2.0209    altermagnet, m_B - m_C = 1.94
#   (113new)  -0.0434    -1.7404    paramagnet,  m_B - m_C = 0.00
#
# The paramagnet has the LOWER eigenvalue sum and the HIGHER true energy, so the
# old criterion preferred it. CPQMC is constrained path, so that choice decides
# the answer: 251 happened to pick the altermagnet and measured m_B - m_C = 1.92,
# while 113new picked the paramagnet and measured 0.014 on identical parameters.
#
# NOTE: the same criterion appears in the checkerboard wf_unified.py. That source
# is deliberately left untouched here, but it is worth checking separately
# whether any checkerboard run picked its trial the same way.
import sys

wf = open("wf_lieb.py").read()

def sub1(old, new, tag):
    global wf
    if wf.count(old) != 1:
        sys.exit(f"[{tag}] expected 1 occurrence, found {wf.count(old)}")
    wf = wf.replace(old, new)

sub1("""        E  = sum(eig_eup[:NUP]) + sum(eig_edn[:NDN])
        E1 = E2
        E2 = E""",
"""        # Convergence is still tracked on the eigenvalue sum, which is what the
        # SCF actually stationarises. E_hf below is the variational energy and is
        # what Run() compares across seeds.
        E  = sum(eig_eup[:NUP]) + sum(eig_edn[:NDN])
        E1 = E2
        E2 = E""", "energy")

sub1("""        if abs(E2 - E1) < 0.0001:
            break
        times += 1

    return nup, ndn, E, eig_vup, eig_vdn""",
"""        if abs(E2 - E1) < 0.0001:
            break
        times += 1

    # Hartree-Fock energy: undo the double counting of the interaction, using the
    # self-consistent densities rather than the mixed ones.
    n_u = Getnup(eig_vup, NUP)
    n_d = Getndn(eig_vdn, NDN)
    E_hf = E - float(np.sum(U * n_u * n_d))
    return nup, ndn, E_hf, eig_vup, eig_vdn""", "return")

sub1("""        print(f"attempt {i+1}/{nitr} (seed={seeds[i]}): E = {E:.6f}", flush=True)""",
     """        print(f"attempt {i+1}/{nitr} (seed={seeds[i]}): E_HF = {E:.6f}", flush=True)""",
     "print")
sub1('print(f"HF ground state energy: {E_gs:.6f}")',
     'print(f"HF ground state energy: {E_gs:.6f}   (variational, not the eigenvalue sum)")',
     "gs print")
open("wf_lieb.py", "w").write(wf)
print("trial now selected by the Hartree-Fock energy")
