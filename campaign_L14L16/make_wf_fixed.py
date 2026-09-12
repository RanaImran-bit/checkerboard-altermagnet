# Build wf_unified_fixed.py from the original, correcting BOTH defects the audit
# exposed. The original wf_unified.py is NOT modified.
#
# Defect 1 - wrong selection criterion.
#   Run() picked the seed minimising sum(eps_up[:NUP]) + sum(eps_dn[:NDN]). Each
#   eigenvalue already carries the opposite spin's mean field, so that sum double
#   counts the interaction and does not rank determinants by energy. Across 304
#   audited runs it chose a variationally worse trial 109 times (36%), and in all
#   28 affected half-filling cells the chosen trial OVERSTATED the staggered
#   moment. Replaced by the energy computed from the orbitals themselves:
#       E = Tr(K rho_up) + Tr(K rho_dn) + U sum_i <n_up> <n_dn>
#   Deliberately not E - sum(U n_up n_dn): that identity holds only at exact self
#   consistency, which the mixed-density SCF does not reach for loosely converged
#   states, i.e. exactly where the ranking matters.
#
# Defect 2 - the SCF stops too early.
#   Convergence was |dE| < 1e-4 on the eigenvalue sum with 3 resets of 1000
#   iterations. In several audited runs a FRESH neel solve beat the saved trial,
#   so the original had not converged. Tightened to 1e-8, 4000 iterations, and
#   more restarts, and the seed list is widened so the landscape is sampled
#   properly rather than relying on one lucky draw.
import re, sys

s = open("/tmp/ck/wf_unified_orig.py").read()

def sub1(old, new, tag):
    global s
    if s.count(old) != 1:
        sys.exit(f"[{tag}] expected 1 occurrence, found {s.count(old)}")
    s = s.replace(old, new)

sub1("""        if abs(E2 - E1) < 0.0001:
            break
        times += 1

    return nup, ndn, E, eig_vup, eig_vdn""",
"""        if abs(E2 - E1) < 1e-8:
            break
        times += 1

    # Variational energy of the determinant these orbitals describe. This, not
    # the eigenvalue sum, is what Run() must compare across seeds.
    r_u = eig_vup[:, :NUP] @ eig_vup[:, :NUP].T
    r_d = eig_vdn[:, :NDN] @ eig_vdn[:, :NDN].T
    n_u = Getnup(eig_vup, NUP)
    n_d = Getndn(eig_vdn, NDN)
    E_hf = float(np.sum(K_up * r_u) + np.sum(K_down * r_d) + U * np.sum(n_u * n_d))
    return nup, ndn, E_hf, eig_vup, eig_vdn""", "energy+tol")

sub1("    MAX_RESETS = 3", "    MAX_RESETS = 6", "resets")
sub1("        if times == 1000:", "        if times == 4000:", "iters")
sub1("    seeds = ['neel', 'pm'] + ['random'] * 5",
     "    seeds = ['neel', 'pm'] + ['random'] * 10", "seeds")
sub1('print(f"attempt {i+1}/{nitr} (seed={seeds[i]}): E = {E:.6f}", flush=True)',
     'print(f"attempt {i+1}/{nitr} (seed={seeds[i]}): E_HF = {E:.6f}", flush=True)', "print")

open("/tmp/ck/wf_unified_fixed.py", "w").write(s)
print("wrote /tmp/ck/wf_unified_fixed.py")
