"""Build vertex summaries for L=8, 10, 12 from the b32 runs collected off 251.

The L=12 summary existed already. The L=8 and L=10 runs turned out to cover the
SAME full (U, delta) grid at the SAME settings (beta=32, tau window 0.8, N_w=500,
six seeds), sitting unused in 251:~/collected/eqt_b32/. Only the U=4 slices had
ever been pulled down, as chi_fss_L{8,10}_U4.csv, which is why the vertex was
believed to be a single-size result.

Cross-checks performed:
  - settings identical across all three sizes
  - the U=4 rows reproduce the existing chi_fss_L{8,10}_U4.csv
  - the rebuilt L=12 reproduces chi_L12_vertex_summary.csv
"""
import os, glob, numpy as np, pandas as pd

D  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
CH = ["son", "sext", "d", "dxy"]

def summarize(df):
    rows = []
    for (L, U, dl), s in df.groupby(["L", "U", "delta"]):
        r = dict(L=int(L), U=U, delta=dl, nseed=len(s))
        for c in CH:
            v = s[f"chi_{c}_vertex"].to_numpy(float)
            r[c] = v.mean()
            r[f"{c}_err"] = v.std(ddof=1)/np.sqrt(len(v)) if len(v) > 1 else np.nan
        rows.append(r)
    return pd.DataFrame(rows).sort_values(["L","U","delta"]).reset_index(drop=True)

src = sorted(glob.glob(os.path.join(D, "chi_fss_b32", "*.csv")))
src += sorted(glob.glob(os.path.join(D, "chi_L12", "eqtime_L12_U*.csv")))
raw = pd.concat([pd.read_csv(f) for f in src], ignore_index=True)
raw = raw.drop_duplicates(subset=["L","U","delta","seed"])

print("=== settings, must be identical across sizes ===")
for L, g in raw.groupby("L"):
    print(f"  L={int(L):>3}: beta={sorted(g.beta_proj.unique())} tau={sorted(g.tau_max.unique())} "
          f"nw={sorted(g.nw.unique())} seeds={g.seed.nunique()} U={sorted(g.U.unique())}")
assert raw.beta_proj.nunique() == 1 and raw.tau_max.nunique() == 1 and raw.nw.nunique() == 1, \
    "settings differ between sizes, do not merge these"

S = summarize(raw)
print("\n=== grid completeness ===")
for L, g in S.groupby("L"):
    p = g.pivot_table(index="U", columns="delta", values="d")
    print(f"  L={int(L):>3}: {len(g)} cells, {p.shape[0]}x{p.shape[1]} grid, "
          f"missing {int(p.isna().sum().sum())}")

print("\n=== cross-check: U=4 rows against the previously collected chi_fss files ===")
for L in (8, 10):
    f = os.path.join(D, f"chi_fss_L{L}_U4.csv")
    if not os.path.exists(f): print(f"  L={L}: no reference file"); continue
    ref = summarize(pd.read_csv(f))
    m = S[(S.L == L) & (S.U == 4.0)].merge(ref, on=["U","delta"], suffixes=("","_ref"))
    worst = max(float((m[c] - m[f"{c}_ref"]).abs().max()) for c in CH)
    print(f"  L={L}: {len(m)} cells, max |difference| {worst:.2e}")

f12 = os.path.join(D, "chi_L12_vertex_summary.csv")
ref = pd.read_csv(f12)
m = S[S.L == 12].merge(ref, on=["U","delta"], suffixes=("","_ref"))
worst = max(float((m[c] - m[f"{c}_ref"]).abs().max()) for c in CH)
print(f"  L=12 vs chi_L12_vertex_summary.csv: {len(m)} cells, max |difference| {worst:.2e}")

out = os.path.join(D, "chi_vertex_summary_3L.csv")
S.to_csv(out, index=False)
print(f"\nwrote {out}: {len(S)} cells, L = {sorted(S.L.unique())}")
