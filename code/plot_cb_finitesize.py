#!/usr/bin/env python3
"""Fig: finite-size stability of the channel swap (meeting point #12).
chi_dxy and chi_dx2-y2 vs delta at a fixed doped filling (n ~ 0.78), overlaid for L = 6, 8, 10, 12.
L=6 from the committed U=4 grid (chi_grid_all.csv); L=8,10,12 from the fss_L*.csv scans.

  python plot_cb_finitesize.py           # auto-picks whichever fss_L*.csv exist
Reads: chi_grid_all.csv (L=6, U=4) + fss_L8.csv / fss_L10.csv / fss_L12.csv
Writes fig_finitesize.png.
"""
import os, glob
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

NTARGET = 0.78
frames = {}
# L=6 from the U=4 grid
if os.path.exists("chi_grid_all.csv"):
    c = pd.read_csv("chi_grid_all.csv"); c = c[np.isclose(c.U, 4.0)]; frames[6] = c
# L=8,10,12 from the finite-size scans
for f in sorted(glob.glob("fss_L*.csv")):
    L = int(f.split("fss_L")[1].split(".")[0]); frames[L] = pd.read_csv(f)

def near(df):
    ns = df.n.unique(); return df[np.isclose(df.n, ns[np.argmin(np.abs(ns - NTARGET))])]

colors = {6: '#1f77b4', 8: '#2ca02c', 10: '#d62728', 12: '#9467bd'}
fig, ax = plt.subplots(1, 2, figsize=(13, 5.2), facecolor='white')
for a, (col, lab) in zip(ax, [('chi_dxy', r'$\chi_{d_{xy}}$'), ('chi_d', r'$\chi_{d_{x^2-y^2}}$')]):
    for L in sorted(frames):
        g = near(frames[L]).groupby('delta')[col].mean().reset_index().sort_values('delta')
        a.plot(g.delta, g[col], '-o', lw=2, ms=7, color=colors.get(L, 'k'), label=f'L={L}')
    a.axhline(0, color='gray', ls='--', lw=0.8); a.set_xlabel(r'anisotropy $\delta$', fontsize=14)
    a.set_ylabel(lab, fontsize=14); a.set_title(fr'{lab} vs $\delta$  (n$\approx${NTARGET})', fontsize=14)
    a.legend(fontsize=12, frameon=False)
fig.suptitle('Finite-size stability of the channel swap (U=4, doped)', fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig("fig_finitesize.png", dpi=200, bbox_inches='tight', facecolor='white')
print("saved fig_finitesize.png  (sizes: %s)" % sorted(frames))
