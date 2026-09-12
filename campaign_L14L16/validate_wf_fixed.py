# Does the corrected generator actually find the better trial?
#
# Re-solves cells the audit flagged and compares against both numbers it
# recorded: the trial the original run used, and the best the audit could find.
# The corrected generator should land at or below audit ms_best in energy.
import sys, numpy as np, importlib.util
src = "/tmp/ck/wf_unified_fixed.py"
s = open(src).read().split("# ================================================================\n# Run —")[0]
open("/tmp/_wff.py", "w").write(s)
spec = importlib.util.spec_from_file_location("_wff", "/tmp/_wff.py")
w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)

def stagger(L, nu, nd):
    m = nu - nd
    sgn = np.array([(-1.0)**(((i % L)+1) + ((i // L)+1)) for i in range(len(m))])
    return float(abs(np.sum(sgn*m))/len(m))

# L,  U,  delta,  N,   ms the original used, ms the audit found
CASES = [
    (8,  3.5, 0.3, 32, 0.4617, 0.1424),
    (8,  4.0, 0.4, 32, 0.4865, 0.1751),
    (10, 3.0, 0.3, 50, 0.3092, 0.0845),
    (10, 4.5, 0.5, 50, 0.4010, 0.1772),
    (12, 4.0, 0.5, 72, 0.2760, 0.1145),
]
print(f"{'L':>3}{'U':>6}{'delta':>7}{'ms_orig':>10}{'ms_audit':>10}{'ms_fixed':>10}"
      f"{'E_fixed':>12}{'seed':>9}")
for L, U, de, N, ms_o, ms_a in CASES:
    tA, tt = -de, 0.3
    K = w.GetK(L, L, 1, tt, tA)
    best = None
    for seed in ["neel", "pm"] + ["random"]*10:
        try:
            _, _, E, vu, vd = w.Iteration(L, L, U, tA, tt, NUP=N, NDN=N, seed=seed)
        except Exception:
            continue
        if best is None or E < best[0]:
            ru = vu[:, :N] @ vu[:, :N].T; rd = vd[:, :N] @ vd[:, :N].T
            best = (E, stagger(L, np.diag(ru), np.diag(rd)), seed)
    print(f"{L:>3}{U:>6.1f}{de:>7.1f}{ms_o:>10.4f}{ms_a:>10.4f}"
          f"{best[1]:>10.4f}{best[0]:>12.5f}{best[2]:>9}")
