# --- rebuild chi_L12_vertex_summary.csv from the raw runs. Paste-and-run on 251. ---
# The shipped summary stopped at U=4.5 even though eqtime_L12_U5.csv had been sitting
# in chi_L12/ since 23 Aug with all six seeds and all seven delta. Any figure using
# UMAX = min(fortran, vertex) therefore silently cut off at 4.5.
import os, glob, numpy as np, pandas as pd

CH  = ["son", "sext", "d", "dxy"]
SRC = next((p for p in ("chi_L12", "../data/chi_L12",
                        os.path.expanduser("~/results/chi_L12"),
                        os.path.expanduser("~/checkerboard-altermagnet/data/chi_L12"))
            if os.path.isdir(p)), None)
if SRC is None: raise SystemExit("cannot find chi_L12/. Try:  find ~ -name 'eqtime_L12_U*.csv'")
OUT = os.path.join(os.path.dirname(SRC.rstrip("/")) or ".", "chi_L12_vertex_summary.csv")

fs = sorted(glob.glob(os.path.join(SRC, "eqtime_L12_U*.csv")))
df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
print(f"{len(fs)} files, {len(df)} rows, U = {sorted(df.U.unique())}")

rows = []
for (U, dl), s in df.groupby(["U", "delta"]):
    r = dict(U=U, delta=dl, nseed=len(s))
    for c in CH:
        v = s[f"chi_{c}_vertex"].to_numpy(float)
        r[c] = v.mean()
        r[f"{c}_err"] = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan
    rows.append(r)
S = pd.DataFrame(rows).sort_values(["U", "delta"]).reset_index(drop=True)

if os.path.exists(OUT):
    old = pd.read_csv(OUT)
    m = old.merge(S, on=["U", "delta"], suffixes=("_old", "_new"))
    worst = max(float((m[f"{c}_old"] - m[f"{c}_new"]).abs().max()) for c in CH)
    print(f"shared cells with the existing file: {len(m)}, max |difference| {worst:.2e}")
    os.replace(OUT, OUT + ".bak")
    print(f"old file kept as {os.path.basename(OUT)}.bak")
S.to_csv(OUT, index=False)
print(f"wrote {OUT}: {len(S)} cells, U = {sorted(S.U.unique())}")
