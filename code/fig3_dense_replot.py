"""FIGURE 3 -- dominant d channel over (n, delta), replotted to delta = 0.7.

Same code as the Pairing_and_AM_L12.ipynb cell. Nothing about the method
changed. The cell plots whatever deltas the file contains, and the file already
holds delta up to 0.7 -- it was the SAVED OUTPUT IMAGE in the notebook that was
stale, generated before the delta = 0.5, 0.6, 0.7 runs landed. Re-executing the
cell is enough; no edit was needed to reach 0.7.

Coverage is now the full 19 fillings x 8 anisotropies at 6 seeds each, with no
empty cells, so shading='gouraud' is smoothing between measured points rather
than across gaps.

This is the notebook cell verbatim, only repointed at the current data file.
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet"
dn = pd.read_csv(f"{D}/data/dense_L12_U4_all.csv")
g = dn.groupby(['n', 'delta'])[['chi_d', 'chi_dxy']].mean().reset_index()
g['diff'] = g.chi_dxy - g.chi_d
NS, DS = np.sort(g.n.unique()), np.sort(g.delta.unique())
Z = g.pivot(index='delta', columns='n', values='diff').reindex(index=DS, columns=NS).values
v = np.nanpercentile(np.abs(Z), 95)

fig, a = plt.subplots(figsize=(9, 5.4), facecolor='white')
im = a.pcolormesh(NS, DS, Z, cmap='coolwarm', shading='gouraud', vmin=-v, vmax=v)
a.contour(NS, DS, Z, levels=[0], colors='k', linewidths=2.4)
a.set_xlabel('filling $n$', fontsize=16); a.set_ylabel(r'anisotropy $\delta$', fontsize=17)
a.tick_params(labelsize=14)
cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.03)
cb.set_label(r'$\chi_{d_{xy}} - \chi_{d_{x^2-y^2}}$', fontsize=15)
a.set_title(rf'$L=12$, $U=4$, {len(NS)} fillings ($\Delta n={np.diff(NS)[0]:.3f}$)', fontsize=16)
plt.tight_layout()
plt.savefig(f"{D}/figs/fig3_phasediagram_dense_L12_U4_to07.png", dpi=600,
            bbox_inches='tight', facecolor='white')

print(f"delta range plotted: {DS.min():g} to {DS.max():g}   ({len(DS)} values)")
print(f"colour scale +/- {v:.2f}")
print("\nchi_dxy - chi_d at half filling, by delta:")
hf = g[np.isclose(g.n, g.n.max())].sort_values('delta')
for _, r in hf.iterrows():
    print(f"   delta={r.delta:g}: {r['diff']:+.2f}")
print("\nfraction of the (n,delta) plane where dxy wins:", f"{(Z > 0).mean()*100:.0f}%")
# plt.show()  # headless: savefig only, so the script does not block
