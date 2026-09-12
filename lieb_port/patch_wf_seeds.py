# Trim the Hartree-Fock seed schedule, which dominates the job runtime.
#
# Measured on 113new: the 48-site trial stage takes about 50 minutes, longer
# than the QMC startup, and nearly all of it is wasted. The schedule was
# inherited from the checkerboard code, where frustration made random restarts
# worth running. On Lieb the physical seeds win by two orders of magnitude and
# the random ones never converge at all:
#
#   12 sites, U=10, t2=0.3:   am 0.056192   neel 0.056430   pm 26.60   random 6.69-6.75
#   48 sites, U=10, t2=0.0:   random seeds exhaust their resets at E = 59.75
#
# So: two random controls instead of five, and one reset instead of three. The
# physical seeds and the paramagnetic control are untouched, so the schedule
# still cannot silently miss a better state without at least reporting it.
#
# Also fixes a cosmetic bug of mine: U became a per-site array in
# patch_kaushal_franz.py, so the non-convergence warning printed all 48 entries.
import sys

wf = open("wf_lieb.py").read()

def sub1(old, new, tag):
    global wf
    if wf.count(old) != 1:
        sys.exit(f"[{tag}] expected 1 occurrence, found {wf.count(old)}")
    wf = wf.replace(old, new)

sub1("    MAX_RESETS = 3", "    MAX_RESETS = 1", "resets")
sub1("    seeds = ['am', 'neel', 'pm'] + ['random'] * 5",
     "    seeds = ['am', 'neel', 'pm'] + ['random'] * 2", "seeds")

# U is a vector now; report the scalar that names the run, not all N entries
sub1('print(f"Warning: U={U} seed={seed}: {MAX_RESETS} resets exhausted, "',
     'print(f"Warning: U={Uscalar} seed={seed}: {MAX_RESETS} resets exhausted, "', "warn1")
sub1('print(f"Warning: U={U} seed={seed} did not converge in 1000 steps. "',
     'print(f"Warning: U={Uscalar} seed={seed} did not converge in 1000 steps. "', "warn2")
sub1("    U = hubbard_u(Lx, Ly, U)             # per site: 0 on A, U on B and C",
     "    Uscalar = U                          # keep the scalar for messages\n"
     "    U = hubbard_u(Lx, Ly, U)             # per site: 0 on A, U on B and C", "uscalar")

open("wf_lieb.py", "w").write(wf)
print("seed schedule trimmed to am/neel/pm + 2 random, MAX_RESETS=1; U print fixed")
