import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt

CSV  = '/home/phd25imran/analysis/data/magnetic_master.csv'
L, NTARGET = 12, 1.0

d = pd.read_csv(CSV); d = d[d.L == L]
g = d.groupby(['U','n','delta']).m.agg(['mean','sem']).reset_index().rename(columns={'mean':'m','sem':'e'})
u0 = g[g.U == 0][['n','delta','m','e']].rename(columns={'m':'m0','e':'e0'})
g = g.merge(u0, on=['n','delta'])
g['dm']  = g.m - g.m0
g['err'] = np.hypot(g.e, g.e0)
g['mAM'] = g.dm * g.delta
g['mAM_err'] = g.err * g.delta

ns = np.sort(g.n.unique()); N = ns[np.argmin(abs(ns - NTARGET))]
s  = g[np.isclose(g.n, N)]
US = sorted(u for u in s.U.unique() if u > 0); DS = np.sort(s.delta.unique())

fig, ax = plt.subplots(2, 2, figsize=(13, 10), facecolor='white')

P = s.pivot(index='delta', columns='U', values='mAM')[US]
im = ax[0,0].pcolormesh(np.arange(len(US)+1), np.r_[DS - 0.05, DS[-1]+0.05], P.values,
                        cmap='magma', shading='flat')
ax[0,0].set_xticks(np.arange(len(US))+0.5); ax[0,0].set_xticklabels([f'{u:g}' for u in US])
ax[0,0].set_yticks(DS)
ax[0,0].set_xlabel('$U$', fontsize=16); ax[0,0].set_ylabel(r'$\delta$', fontsize=17)
ax[0,0].set_title(r'(a)  $m_{\rm AM}=\Delta m\,\delta$', fontsize=16)
fig.colorbar(im, ax=ax[0,0])

for dd, c in zip(DS, plt.cm.viridis(np.linspace(0,.9,len(DS)))):
    t = s[np.isclose(s.delta, dd)].sort_values('U')
    ax[0,1].errorbar(t.U, t.dm, t.err, marker='o', color=c, lw=1.8, ms=5,
                     capsize=2, label=rf'$\delta={dd:g}$')
ax[0,1].set_xlabel('$U$', fontsize=16); ax[0,1].set_ylabel(r'$\Delta m$', fontsize=16)
ax[0,1].set_title(r'(b)  interaction builds the moment', fontsize=16)
ax[0,1].legend(fontsize=9, ncol=2)

for u, c in zip(US, plt.cm.plasma(np.linspace(0,.8,len(US)))):
    t = s[np.isclose(s.U, u)].sort_values('delta')
    ax[1,0].errorbar(t.delta, t.dm, t.err, marker='s', color=c, lw=1.8, ms=5,
                     capsize=2, label=f'$U={u:g}$')
ax[1,0].set_xlabel(r'$\delta$', fontsize=17); ax[1,0].set_ylabel(r'$\Delta m$', fontsize=16)
ax[1,0].set_title(r'(c)  anisotropy suppresses it', fontsize=16)
ax[1,0].legend(fontsize=11)

for u, c in zip(US, plt.cm.plasma(np.linspace(0,.8,len(US)))):
    t = s[np.isclose(s.U, u)].sort_values('delta')
    ax[1,1].errorbar(t.delta, t.mAM, t.mAM_err, marker='^', color=c, lw=1.8, ms=6,
                     capsize=2, label=f'$U={u:g}$')
ax[1,1].set_xlabel(r'$\delta$', fontsize=17); ax[1,1].set_ylabel(r'$m_{\rm AM}$', fontsize=16)
ax[1,1].set_title(r'(d)  the product', fontsize=16)
ax[1,1].legend(fontsize=11)
for a in ax.ravel(): a.tick_params(labelsize=12)
fig.suptitle(rf'Altermagnetic order parameter, $L={L}$, $n={N:.3f}$', fontsize=18)
plt.tight_layout(rect=[0,0,1,0.96])
plt.savefig('/home/phd25imran/analysis/fig_mam.png', dpi=105, facecolor='white')
print('n =', N)
print(s[s.U>0].groupby('U').apply(lambda t: f"dm({DS[0]:g})={t[np.isclose(t.delta,DS[0])].dm.values[0]:.4f} -> dm(0.7)={t[np.isclose(t.delta,0.7)].dm.values[0]:.4f}", include_groups=False).to_string())
