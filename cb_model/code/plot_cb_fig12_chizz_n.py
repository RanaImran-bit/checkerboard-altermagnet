#!/usr/bin/env python3
"""Fig 12: staggered magnetic susceptibility chi_zz(pi,pi) vs filling n.
Reads chizz_vs_n.csv (columns: n, chizz_pipi, err) produced from checkerboard_chispin.py.
  python plot_cb_fig12_chizz_n.py ../data/chizz_vs_n.csv
"""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- EDIT HERE: style ----
FS_LABEL = 18       # axis label font size
FS_TICK  = 13       # tick font size
FS_TITLE = 16       # title font size
MS       = 10       # marker size
COLOR    = '#8c2d8c'
TITLE    = r'$6\times6$,  $\delta=0.4$,  $U=4$'
YLABEL   = r'$\chi_{zz}(\pi,\pi)$ (staggered magnetic)'
OUT      = "Fig12_chizz_vs_n.png"
# --------------------------

CSV = sys.argv[1] if len(sys.argv) > 1 else "../data/chizz_vs_n.csv"
d = pd.read_csv(CSV).sort_values('n')

fig, ax = plt.subplots(figsize=(7, 5.5), facecolor='white')
ax.errorbar(d.n, d.chizz_pipi, yerr=d.err, marker='s', ms=MS, lw=2, capsize=4, color=COLOR)
ax.set_xlabel(r'Filling $n$', fontsize=FS_LABEL)
ax.set_ylabel(YLABEL, fontsize=FS_LABEL - 2)
ax.set_title(TITLE, fontsize=FS_TITLE)
ax.tick_params(labelsize=FS_TICK)
plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
print("saved", OUT)
