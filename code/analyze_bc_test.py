"""Does the delta ~ 0.3 slope reversal survive antiperiodic boundary conditions?

At half filling with PERIODIC boundaries the U = 0 gap of the checkerboard model is
EXACTLY zero at every L and delta (2 to 12 states at E_F), so which states are
occupied is decided by the trial wavefunction rather than by the physics.
Antiperiodic boundaries shift the k-mesh off the high-symmetry points and close the
shell. If the SIGN of d(Delta_tot)/dU agrees between the two, the reversal is
physical; if it flips, it was partly a filling artifact.

Magnitudes are NOT expected to agree, since the two runs sample different k-meshes.
Only the sign of the slope is a controlled comparison, and both sides are run at
identical settings (beta = 3, N_w = 160) so that the boundary phase is the only
difference. These are therefore NOT production numbers and must never share an axis
with the Fortran data.

  python analyze_bc_test.py pbc.csv apbc.csv
"""
import sys
import numpy as np
import pandas as pd


def agg(df):
    g = df.groupby(["delta", "U"])["delta_tot"]
    return g.mean(), g.std(ddof=1) / np.sqrt(g.count())


def main(fp_pbc, fp_apbc, ulo=2.0, uhi=5.0):
    p, a = pd.read_csv(fp_pbc), pd.read_csv(fp_apbc)
    print(f"PBC rows={len(p)} apx={p.apx.unique()} | APBC rows={len(a)} apx={a.apx.unique()}\n")
    print(f"{'delta':>6}  {'BC':>5} | {'Dtot(U=%g)' % ulo:>18} {'Dtot(U=%g)' % uhi:>18} | "
          f"{'slope':>18}  sign")
    res = {}
    for d in sorted(p.delta.unique()):
        for nm, df in (("PBC", p), ("APBC", a)):
            m, e = agg(df)
            v0, e0 = m[(d, ulo)], e[(d, ulo)]
            v1, e1 = m[(d, uhi)], e[(d, uhi)]
            sl, es = v1 - v0, np.hypot(e0, e1)
            res[(d, nm)] = (sl, es, sl / es if es > 0 else np.nan)
            print(f"{d:6.1f}  {nm:>5} | {v0:11.4f} +/-{e0:6.4f} {v1:11.4f} +/-{e1:6.4f} | "
                  f"{sl:+11.4f} +/-{es:6.4f}  {'+' if sl > 0 else '-'} "
                  f"({sl/es:+.1f} sigma)")
        print()
    print("VERDICT per delta:")
    ok = 0
    for d in sorted(p.delta.unique()):
        sp, _, gp = res[(d, "PBC")]
        sa, _, ga = res[(d, "APBC")]
        same = (sp > 0) == (sa > 0)
        ok += same
        weak = min(abs(gp), abs(ga)) < 2.0
        print(f"  delta={d}: PBC {'+' if sp > 0 else '-'} ({gp:+.1f}s)   "
              f"APBC {'+' if sa > 0 else '-'} ({ga:+.1f}s)   -> "
              f"{'AGREE' if same else 'FLIP'}"
              f"{'   (one side < 2 sigma: near the zero crossing, do not quote)' if weak else ''}")
    print(f"\nsigns agree at {ok}/{len(p.delta.unique())} deltas")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0] if a else "../data/bc_test/pbc.csv",
         a[1] if len(a) > 1 else "../data/bc_test/apbc.csv")
