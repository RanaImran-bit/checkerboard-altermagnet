import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
OUT='/home/phd25imran/analysis'
CSV='/home/phd25imran/analysis/data/pairing_master.csv'
L=12
CH=[('chi_son',r'on-site $s$'),('chi_sext',r'extended $s$'),
    ('chi_d',r'$d_{x^2-y^2}$'),('chi_dxy',r'$d_{xy}$')]

d=pd.read_csv(CSV); d=d[d.L==L]; d['N']=d.L**2
for c,_ in CH: d[c+'_n']=d[c]/d.N
US=sorted(d.U.unique())
g=d.groupby(['U','n','delta'])[[c+'_n' for c,_ in CH]].mean().reset_index()

# ONE shared symmetric scale across every panel. 2nd/98th percentile over all
# channels and all U, so the on-site-s extreme does not flatten everything else.
allv=np.concatenate([g[c+'_n'].values for c,_ in CH])
lim=float(np.max(np.abs(allv)))
# ONE shared scale, but symmetric-log: the on-site-s extreme reaches -0.59 while the
# physics in the d channels lives at +/-0.04. A linear shared scale flattens the latter
# to white. linthresh sets where the scale turns from linear to logarithmic.
NORM=SymLogNorm(linthresh=0.01, vmin=-lim, vmax=lim, base=10)
print(f'shared range +/-{lim:.4f}, linear below 0.01  (raw {allv.min():+.3f} to {allv.max():+.3f})')

ns=np.sort(g.n.unique()); ds=np.sort(g.delta.unique())
def edges(v):
    v=np.asarray(v,float); m=(v[1:]+v[:-1])/2
    return np.r_[v[0]-(m[0]-v[0]), m, v[-1]+(v[-1]-m[-1])]
NE,DE=edges(ns),edges(ds)

fig,ax=plt.subplots(len(US),4,figsize=(15.5,3.15*len(US)),facecolor='white',squeeze=False)
for i,U in enumerate(US):
    for j,(c,lab) in enumerate(CH):
        a=ax[i][j]
        piv=g[np.isclose(g.U,U)].pivot(index='delta',columns='n',values=c+'_n')
        im=a.pcolormesh(NE,DE,piv.values,cmap='RdBu_r',norm=NORM,shading='flat')
        a.set_xticks([0.5,0.7,0.9]); a.set_yticks([0,0.2,0.4,0.6])
        a.tick_params(labelsize=10)
        if i==0: a.set_title(lab,fontsize=16)
        if j==0: a.set_ylabel(f'$U={U:g}$\n'+r'$\delta$',fontsize=13)
        else: a.set_yticklabels([])
        if i==len(US)-1: a.set_xlabel('Filling $n$',fontsize=13)
        else: a.set_xticklabels([])
cb=fig.colorbar(im,ax=ax,fraction=0.016,pad=0.015)
cb.set_label(r'pairing vertex  $\chi_\alpha/N$   (symlog)   red = attractive,  blue = repulsive',fontsize=12)
fig.suptitle(rf'$(n,\delta)$ phase diagram, all channels and interactions, $L={L}$  —  one shared colour scale',
             fontsize=17)
plt.savefig(f'{OUT}/fig_pd.png',dpi=105,facecolor='white',bbox_inches='tight')
print('nothing clipped; full range shown')
