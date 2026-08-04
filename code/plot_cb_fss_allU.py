#!/usr/bin/env python3
"""Finite-size figures (#12) across the FULL (n, delta, U) cube at L = 8, 10, 12.
Reads fss_full_all.csv (cols: L,nup,delta,U,seed,chi_son,chi_sext,chi_d,chi_dxy,n)
plus chi_grid_all.csv for the L=6 baseline, and writes, for EVERY U:

  fig_fss_maps_U{U}.png    - (n,delta) phase maps: rows = L(8,10,12), cols = 4 channels
  fig_fss_lines_U{U}.png   - 4 channels vs delta at n~0.78 (doped), overlaid L = 6,8,10,12
  fig_fss_swap_U{U}.png    - the channel swap: chi_dxy and chi_dx2-y2 vs delta, overlaid across L

Usage:
  python plot_cb_fss_allU.py [fss_full_all.csv] [chi_grid_all.csv]
"""
import sys, os
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

FSS = sys.argv[1] if len(sys.argv) > 1 else "fss_full_all.csv"
L6  = sys.argv[2] if len(sys.argv) > 2 else "chi_grid_all.csv"
NDOPED = 0.78                      # doped filling for the line/swap cuts
NMIN   = 0.45                      # keep everything (grid starts at 0.5)
CH  = [('chi_son', 'on-site s'), ('chi_sext', 'extended-s'),
       ('chi_d', r'$d_{x^2-y^2}$'), ('chi_dxy', r'$d_{xy}$')]
COL = ['#7f7f7f', '#2ca02c', '#1f77b4', '#d62728']   # on-site s, ext-s, dx2-y2, dxy
LCOL = {6: '#1f77b4', 8: '#2ca02c', 10: '#d62728', 12: '#9467bd'}

d = pd.read_csv(FSS); d = d[d.n >= NMIN]
US = sorted(d.U.unique())
SIZES = sorted(d.L.unique())

# L=6 baseline (has no L column; add it), aligned columns
l6 = None
if os.path.exists(L6):
    l6 = pd.read_csv(L6)
    if 'L' not in l6.columns:
        l6['L'] = 6


def near_filling(df, ntarget):
    ns = np.sort(df.n.unique())
    return df[np.isclose(df.n, ns[np.argmin(np.abs(ns - ntarget))])]


for U in US:
    dU = d[np.isclose(d.U, U)]
    gU = dU.groupby(['L', 'n', 'delta'])[[c for c, _ in CH]].mean().reset_index()

    # ---- Fig 1: (n,delta) maps, rows = L, cols = channel ----
    fig, axes = plt.subplots(len(SIZES), 4, figsize=(19, 4.4 * len(SIZES)), facecolor='white')
    for r, L in enumerate(SIZES):
        gL = gU[gU.L == L]
        for c, (col, lab) in enumerate(CH):
            a = axes[r, c]
            v = gL[col].values
            vmin, vmax = np.nanpercentile(v, 5), np.nanpercentile(v, 95)
            sc = a.scatter(gL.n, gL.delta, c=v, cmap='viridis', vmin=vmin, vmax=vmax,
                           s=340, marker='s', edgecolors='k', linewidths=0.4)
            if r == 0:
                a.set_title(lab, fontsize=17)
            if c == 0:
                a.set_ylabel(f'L={L}\n' + r'$\delta$', fontsize=14)
            if r == len(SIZES) - 1:
                a.set_xlabel('Filling $n$', fontsize=13)
            a.set_xlim(0.45, 1.05); a.set_ylim(-0.05, 0.45); a.tick_params(labelsize=11)
            fig.colorbar(sc, ax=a, fraction=0.046, pad=0.04)
    fig.suptitle(fr'Finite-size pairing susceptibility over $(n,\delta)$ at $U={U:g}$',
                 fontsize=18, y=1.0)
    plt.tight_layout()
    plt.savefig(f'fig_fss_maps_U{U:g}.png', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()

    # ---- Fig 2: 4-channel line cuts vs delta at doped filling, overlaid across L ----
    fig, ax = plt.subplots(1, 4, figsize=(20, 4.7), facecolor='white')
    for a, (col, lab) in zip(ax, CH):
        for L in ([6] + SIZES if l6 is not None else SIZES):
            src = l6[np.isclose(l6.U, U)] if L == 6 else dU[dU.L == L]
            if len(src) == 0 or col not in src.columns:
                continue
            sub = near_filling(src, NDOPED)
            s = sub.groupby('delta')[col].agg(['mean', 'sem']).reset_index().sort_values('delta')
            a.errorbar(s.delta, s['mean'], yerr=s['sem'].fillna(0), marker='o', ms=6, lw=1.8,
                       capsize=3, color=LCOL.get(L, 'k'), label=f'L={L}')
        a.axhline(0, color='gray', ls='--', lw=0.8)
        a.set_title(lab, fontsize=16); a.set_xlabel(r'anisotropy $\delta$', fontsize=13)
        a.set_ylabel(r'$\chi^{\rm vertex}$', fontsize=13); a.legend(fontsize=10, frameon=False)
    fig.suptitle(fr'Channel cuts vs $\delta$ at $n\approx{NDOPED:g}$, $U={U:g}$ (finite-size overlay)',
                 fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(f'fig_fss_lines_U{U:g}.png', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()

    # ---- Fig 3: the channel swap (dxy vs dx2-y2) vs delta, overlaid across L ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2), facecolor='white')
    for a, (col, lab) in zip(ax, [('chi_dxy', r'$\chi_{d_{xy}}$'), ('chi_d', r'$\chi_{d_{x^2-y^2}}$')]):
        for L in ([6] + SIZES if l6 is not None else SIZES):
            src = l6[np.isclose(l6.U, U)] if L == 6 else dU[dU.L == L]
            if len(src) == 0 or col not in src.columns:
                continue
            sub = near_filling(src, NDOPED)
            s = sub.groupby('delta')[col].agg(['mean', 'sem']).reset_index().sort_values('delta')
            a.errorbar(s.delta, s['mean'], yerr=s['sem'].fillna(0), marker='o', ms=7, lw=2,
                       capsize=3, color=LCOL.get(L, 'k'), label=f'L={L}')
        a.axhline(0, color='gray', ls='--', lw=0.8); a.set_xlabel(r'anisotropy $\delta$', fontsize=14)
        a.set_ylabel(lab, fontsize=14); a.set_title(fr'{lab} vs $\delta$ (n$\approx${NDOPED:g})', fontsize=14)
        a.legend(fontsize=12, frameon=False)
    fig.suptitle(fr'Finite-size stability of the channel swap ($U={U:g}$, doped)', fontsize=15, y=1.01)
    plt.tight_layout()
    plt.savefig(f'fig_fss_swap_U{U:g}.png', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"U={U:g}: saved fig_fss_maps_U{U:g}.png, fig_fss_lines_U{U:g}.png, fig_fss_swap_U{U:g}.png")

print("done. Us =", [f'{u:g}' for u in US])
