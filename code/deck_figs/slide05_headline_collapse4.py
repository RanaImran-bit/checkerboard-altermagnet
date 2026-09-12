import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt

CSV='/home/phd25imran/analysis/data/pairing_master.csv'
NTARGET, LS, U = 1.0, [8,10,12], 4.0
CH=[('chi_son',r'on-site $s$'),('chi_sext',r'extended $s$'),
    ('chi_d',r'$d_{x^2-y^2}$'),('chi_dxy',r'$d_{xy}$')]
MK={8:'o',10:'s',12:'^'}

d=pd.read_csv(CSV); d=d[d.L.isin(LS)]; d['N']=d.L**2
for c,_ in CH: d[c+'_n']=d[c]/d.N
ns=np.sort(d.n.unique()); N=ns[np.argmin(abs(ns-NTARGET))]
d=d[np.isclose(d.n,N)&np.isclose(d.U,U)]
g=d.groupby(['L','delta'])[[c+'_n' for c,_ in CH]].agg(['mean','sem'])

fig,ax=plt.subplots(1,4,figsize=(18,4.4),facecolor='white')
for a,(c,lab) in zip(ax,CH):
    for L in LS:
        s=g.loc[L]
        a.errorbar(s.index.values,s[(c+'_n','mean')],s[(c+'_n','sem')],
                   marker=MK[L],lw=1.7,ms=6,capsize=2,label=f'$L={L}$')
    a.axhline(0,color='k',lw=0.9,ls=':')
    a.set_title(lab,fontsize=17); a.set_xlabel(r'$\delta$',fontsize=16)
    a.tick_params(labelsize=12)
ax[0].set_ylabel(r'$\chi_\alpha/N$',fontsize=16); ax[0].legend(fontsize=12)
fig.suptitle(rf'Per-site pairing vertex, all four channels, $U={U:g}$, $n={N:.3f}$',fontsize=17)
plt.tight_layout(rect=[0,0,1,0.92]); plt.savefig('/home/phd25imran/analysis/fig_collapse4.png',dpi=105,facecolor='white')

print(f'U={U:g} n={N:.3f}   chi/N   (delta=0 -> peak -> delta=0.7),  collapse spread across L')
for c,lab in CH:
    print(f'\n {lab}')
    for L in LS:
        v=g.loc[L][(c+'_n','mean')].values; dd=g.loc[L].index.values
        print(f'   L={L:2d}: '+'  '.join(f'{x:+.4f}' for x in v))
    v12=g.loc[12][(c+'_n','mean')].values
    sp=[g.loc[L][(c+'_n','mean')].values for L in LS]
    spread=100*np.max([abs(a-b) for a in sp for b in sp])/max(abs(v12).max(),1e-9)
    sign='no sign change' if (v12>0).all() or (v12<0).all() else 'CHANGES SIGN'
    print(f'   delta= '+'  '.join(f'{x:+.4f}'.replace('+','').rjust(7) for x in dd))
    print(f'   -> {sign};  max spread across L = {spread:.0f}% of peak')
