"""(U, delta) heatmaps at half filling: AFM structure factor and the AM order parameter.

Two different quantities, deliberately side by side, because they move in
OPPOSITE directions:

  S_AFM(q*)  how strong the spin correlations are      -> falls with delta
  Psi_dxy    whether they break the diagonal symmetry  -> rises with delta

S is taken at the TRUE peak of S(q), not at a fixed (pi,pi): the ordering vector
moves with delta (spiral at small delta, Neel near half filling for
delta = 0.2 to 0.5, stripe above).

Neither is a polarisation. A polarisation is the one-body <S^z_i>, identically
zero here because every run has n_up = n_dn (S_z = 0 sector).
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "axes.labelsize": 16, "xtick.labelsize": 13, "ytick.labelsize": 13,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
pk = pd.read_csv(f"{D}/peak_master.csv"); pk = pk[(pk.L == 12) & np.isclose(pk.n, 1.0)]
mg = pd.read_csv(f"{D}/magnetic_master.csv")
mg = mg[(mg.L == 12) & np.isclose(mg.n, 1.0)]

A = pk.groupby(["U", "delta"]).S_peak.mean().unstack("delta")
B = mg.groupby(["U", "delta"]).Psi_dxy.mean().abs().unstack("delta") * 1e3
B = B.drop(index=0.0, errors="ignore")          # U=0 has no interaction-driven signal

fig, ax = plt.subplots(1, 2, figsize=(14.5, 5.4))
for k, (M, ttl, cl) in enumerate([
        (A, r"(a) AFM structure factor  $S_{\rm AFM}(\mathbf{q}^*)$", "viridis"),
        (B, r"(b) AM order parameter  $|\Psi_{d_{xy}}|\times 10^{3}$", "magma")]):
    a = ax[k]
    im = a.pcolormesh(M.columns.values, M.index.values, M.values,
                      cmap=cl, shading="gouraud", rasterized=True)
    a.set_xlabel(r"anisotropy  $\delta$"); a.set_ylabel(r"interaction  $U$")
    a.set_title(ttl, fontsize=16)
    a.set_yticks(M.index.values); a.set_xticks(M.columns.values)
    for i, U in enumerate(M.index.values):
        for j, dl in enumerate(M.columns.values):
            v = M.values[i, j]
            a.text(dl, U, f"{v:.2f}", ha="center", va="center", fontsize=9,
                   color="white" if v < np.nanpercentile(M.values, 55) else "black")
    # gouraud cells are centred on the data points, so the outermost labels sit
    # half outside the axes without a margin
    dx = (M.columns.values[-1] - M.columns.values[0]) * 0.06
    dy = (M.index.values[-1] - M.index.values[0]) * 0.09
    a.set_xlim(M.columns.values[0] - dx, M.columns.values[-1] + dx)
    a.set_ylim(M.index.values[0] - dy, M.index.values[-1] + dy)
    fig.colorbar(im, ax=a, pad=0.02)
fig.tight_layout()
fig.savefig(f"{D}/../figs/fig_heatmap_U_delta.png", dpi=600, bbox_inches="tight",
            facecolor="white")
print("S_AFM(q*) at half filling, rows U, cols delta:")
print(A.round(3).to_string())
print("\n|Psi_dxy| x1e3:")
print(B.round(2).to_string())
