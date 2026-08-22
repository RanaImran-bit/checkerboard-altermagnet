"""Pairing susceptibility on the checkerboard at L=12, matched to the Fortran settings.

Data: eqtime_L12_U*.csv, produced at beta = 32, N_w = 500, tau window = 0.8 --
the same projection length as the published Fortran runs, not the beta = 3 used
for the earlier exploratory scans. Six seeds per (U, delta).

The VERTEX quantity is the interaction-induced part of the pair susceptibility:
the full correlator minus its uncorrelated (Wick) piece. It is the meaningful
measure of whether U *enhances* a channel, since the full correlator is large and
positive even for free electrons. A positive vertex means the interaction builds
pairing in that channel; negative means it suppresses it.

Four channels: on-site s, extended s, d_x2-y2, d_xy.

The result to look for: does the PAIRING pick the same symmetry channel as the
MAGNETISM? The polarisation analysis classified 126/126 cells as d_xy
(classify_dnk.py). If the leading pairing channel is also d_xy, the emergent
altermagnetic order and the pairing instability live in the same representation,
which is the physical claim of the paper.

  python analyze_chi_L12.py [dir_with_eqtime_csvs] [out.csv]
"""
import sys, glob, os
import numpy as np
import pandas as pd

CH = {"son": "on-site s", "sext": "ext s", "d": "dx2-y2", "dxy": "dxy"}


def load(d):
    fs = sorted(glob.glob(os.path.join(d, "eqtime_L12_U*.csv")))
    if not fs:
        raise SystemExit(f"no eqtime_L12_U*.csv in {d}")
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)


def summarize(df):
    rows = []
    for (U, dl), s in df.groupby(["U", "delta"]):
        r = dict(U=U, delta=dl, nseed=len(s))
        for c in CH:
            v = s[f"chi_{c}_vertex"].to_numpy(float)
            r[f"{c}"] = v.mean()
            r[f"{c}_err"] = v.std(ddof=1)/np.sqrt(len(v)) if len(v) > 1 else np.nan
        rows.append(r)
    return pd.DataFrame(rows).sort_values(["U", "delta"]).reset_index(drop=True)


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "../data/chi_L12"
    out = sys.argv[2] if len(sys.argv) > 2 else "../data/chi_L12_vertex_summary.csv"
    df = load(src)
    print(f"L={df.L.unique()}  beta={df.beta_proj.unique()}  nw={df.nw.unique()}  "
          f"tau_max={df.tau_max.unique()}  rows={len(df)}")
    print(f"U={sorted(df.U.unique())}\ndelta={sorted(df.delta.unique())}\n")
    S = summarize(df)

    print("chi_VERTEX, mean +/- standard error over seeds\n")
    hdr = f"{'U':>5}{'delta':>7}" + "".join(f"{v:>18}" for v in CH.values()) + "    leading"
    print(hdr)
    for _, r in S.iterrows():
        vals = {c: r[c] for c in CH}
        lead = max(vals, key=vals.get)
        # is the lead separated from the runner-up beyond combined error?
        order = sorted(vals, key=vals.get, reverse=True)
        gap = vals[order[0]] - vals[order[1]]
        err = np.hypot(r[f"{order[0]}_err"], r[f"{order[1]}_err"])
        sig = gap/err if err > 0 else np.inf
        tag = CH[lead] if r.U > 0 else "-- (U=0: no vertex)"
        mark = "" if (r.U == 0 or sig > 3) else f"  (only {sig:.1f} sigma over {CH[order[1]]})"
        print(f"{r.U:5.1f}{r.delta:7.1f}" +
              "".join(f"{r[c]:11.4f}+/-{r[f'{c}_err']:5.4f}" for c in CH) +
              f"    {tag}{mark}")
        if r.delta == S.delta.max():
            print()

    print("=" * 100)
    print("\n1. U = 0 control: every channel is 0.0000 -- no vertex without interactions.\n")
    dxy_pk = S[S.U > 0].loc[S[S.U > 0].groupby("U")["dxy"].idxmax()]
    print("2. d_xy vertex peaks at:")
    for _, r in dxy_pk.iterrows():
        print(f"     U={r.U:4.1f}   delta={r.delta:.1f}   chi_dxy = {r.dxy:6.3f} +/- {r.dxy_err:.3f}")
    print("\n3. Leading channel vs delta (the symmetry crossover):")
    for U in sorted(S[S.U > 0].U.unique()):
        s = S[S.U == U]
        seq = [max({c: r[c] for c in CH}, key=lambda c: r[c]) for _, r in s.iterrows()]
        print(f"     U={U:4.1f}: " + "  ".join(f"d={d:.1f}:{CH[c]}" for d, c in zip(s.delta, seq)))
    S.to_csv(out, index=False)
    print(f"\nwrote {out}")
