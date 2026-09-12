# Was the checkerboard trial wave function variationally the best one found?
#
# Run() picked the seed with the smallest SUM OF MEAN-FIELD EIGENVALUES. That sum
# double counts the interaction and does not rank Slater determinants by energy,
# so the selected trial can be worse than another one the same SCF already found.
# CPQMC is constrained path, so the trial decides the physics.
#
# For every run this
#   1. loads the trial actually used (wfup.txt, wfdn.txt) and computes its true
#          E = Tr(K rho_up) + Tr(K rho_dn) + U sum_i <n_up> <n_dn>
#      together with its staggered moment,
#   2. re-runs the Hartree-Fock from the same seed schedule, scoring every seed
#      by that same variational energy, and takes the best,
#   3. reports the gap. A positive gap means a better trial was reachable and the
#      run was constrained to a worse state.
#
# Random seeds are stochastic, so step 2 does not reproduce the original run
# exactly. That is fine: the question is whether a better trial exists, not
# whether we can replay the old draw.
#
# READ ONLY. wf_unified.py is imported from a stripped copy in /tmp and the
# original is never written to.
import os, re, sys, csv, numpy as np

SRC = sys.argv[1]                    # a wf_unified.py to import
DIRS = sys.argv[2:]                  # run directories to audit

strip = "/tmp/_wfu_audit.py"
s = open(SRC).read().split("# ================================================================\n# Run —")[0]
open(strip, "w").write(s)
sys.path.insert(0, "/tmp")
import importlib.util
spec = importlib.util.spec_from_file_location("_wfu_audit", strip)
wfu = importlib.util.module_from_spec(spec); spec.loader.exec_module(wfu)

def parse(folder):
    L  = int(re.search(r"L(\d+)", folder).group(1))
    u  = float(re.search(r"u(-?\d+(?:\.\d+)?)", folder).group(1))
    tA = float(re.search(r"tA(-?\d+(?:\.\d+)?)", folder).group(1))
    tt = float(re.search(r"tt(-?\d+(?:\.\d+)?)", folder).group(1))
    return L, u, tA, tt

def occ(path):
    a = np.loadtxt(path)
    if a.ndim == 1: a = a[None, :]
    return a.T                                   # (nsites, norb)

def energy(K, pu, pd, U):
    ru = pu @ pu.T; rd = pd @ pd.T
    nu, nd = np.diag(ru).copy(), np.diag(rd).copy()
    return float(np.sum(K*ru) + np.sum(K*rd) + U*np.sum(nu*nd)), nu, nd

def stagger(L, nu, nd):
    m = nu - nd
    sgn = np.array([(-1.0)**(((i % L)+1) + ((i // L)+1)) for i in range(len(m))])
    return float(np.abs(np.sum(sgn*m))/len(m))

w = csv.writer(sys.stdout)
w.writerow(["folder","L","U","tA","tt","NUP","NDN","E_used","ms_used",
            "E_best","ms_best","gap","best_seed","verdict"])
for d in DIRS:
    try:
        folder = os.path.basename(d.rstrip("/"))
        L, U, tA, tt = parse(folder)
        pu, pd = occ(os.path.join(d,"wfup.txt")), occ(os.path.join(d,"wfdn.txt"))
        NUP, NDN = pu.shape[1], pd.shape[1]
        K = wfu.GetK(L, L, 1, tt, tA)
        E_used, nu, nd = energy(K, pu, pd, U)
        ms_used = stagger(L, nu, nd)

        best = None
        for seed in ["neel","pm"] + ["random"]*5:
            try:
                a, b, _E, vu, vd = wfu.Iteration(L, L, U, tA, tt, NUP=NUP, NDN=NDN, seed=seed)
            except Exception:
                continue
            e, n1, n2 = energy(K, vu[:, :NUP], vd[:, :NDN], U)
            if best is None or e < best[0]:
                best = (e, stagger(L, n1, n2), seed)
        if best is None:
            continue
        gap = E_used - best[0]
        v = ("OK" if gap < 1e-4 else
             ("WORSE_TRIAL_same_order" if abs(ms_used-best[1]) < 0.02 else "WORSE_TRIAL_DIFFERENT_ORDER"))
        w.writerow([folder,L,U,tA,tt,NUP,NDN,f"{E_used:.6f}",f"{ms_used:.4f}",
                    f"{best[0]:.6f}",f"{best[1]:.4f}",f"{gap:.6f}",best[2],v])
        sys.stdout.flush()
    except Exception as e:
        print(f"# skip {d}: {e}", file=sys.stderr)
