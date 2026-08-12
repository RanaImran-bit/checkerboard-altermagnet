import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt

CSV='/home/phd25imran/analysis/data/pairfull_delta0.csv'
L, UFIX, NFIX = 12, 4.0, 1.0
CH=[('son',r'on-site $s$'),('sext',r'extended $s$'),
    ('d',r'$d_{x^2-y^2}$'),('dxy',r'$d_{xy}$')]

d=pd.read_csv(CSV); d=d[d.L==L]
for c,_ in CH:
    d[f'bub_{c}']=d[f'chi_{c}_full']-d[f'chi_{c}_vertex']

fig,ax=plt.subplots(2,4,figsize=(18,8.4),facecolor='white')

# row 1: vs U at half filling
s=d[np.isclose(d.n,NFIX)]
g=s.groupby('U')[[f'chi_{c}_full' for c,_ in CH]+[f'bub_{c}' for c,_ in CH]
                 +[f'chi_{c}_vertex' for c,_ in CH]].agg(['mean','sem'])
for j,(c,lab) in enumerate(CH):
    a=ax[0,j]; U=g.index.values
    a.errorbar(U,g[(f'chi_{c}_full','mean')],g[(f'chi_{c}_full','sem')],
               marker='o',ms=8,lw=2,color='C3',capsize=3,label=r'full  $P$')
    a.plot(U,g[(f'bub_{c}','mean')],marker='o',ms=8,lw=2,color='C0',
           mfc='white',mew=2,ls='--',label=r'bubble  $\bar P$')
    a.set_title(lab,fontsize=17); a.set_xlabel('$U$',fontsize=15)
    a.tick_params(labelsize=12)
    v=g[(f'chi_{c}_vertex','mean')].values
    a.text(0.05,0.05,'vertex '+('>0 attractive' if v[1:].mean()>0 else '<0 repulsive'),
           transform=a.transAxes,fontsize=12,
           color='darkred' if v[1:].mean()>0 else 'navy')
ax[0,0].set_ylabel(rf'$\chi$   ($n={NFIX:g}$)',fontsize=15); ax[0,0].legend(fontsize=12)

# row 2: vs filling at fixed U
s2=d[np.isclose(d.U,UFIX)]
g2=s2.groupby('n')[[f'chi_{c}_full' for c,_ in CH]+[f'bub_{c}' for c,_ in CH]].agg(['mean','sem'])
for j,(c,lab) in enumerate(CH):
    a=ax[1,j]; nn=g2.index.values
    a.errorbar(nn,g2[(f'chi_{c}_full','mean')],g2[(f'chi_{c}_full','sem')],
               marker='s',ms=7,lw=2,color='C3',capsize=3,label=r'full  $P$')
    a.plot(nn,g2[(f'bub_{c}','mean')],marker='s',ms=7,lw=2,color='C0',
           mfc='white',mew=2,ls='--',label=r'bubble  $\bar P$')
    a.set_xlabel('Filling $n$',fontsize=15); a.tick_params(labelsize=12)
ax[1,0].set_ylabel(rf'$\chi$   ($U={UFIX:g}$)',fontsize=15); ax[1,0].legend(fontsize=12)
fig.suptitle(rf'Full vs uncorrelated pair-field susceptibility, $L={L}$, $\delta=0$'
             '\n'+r'$P>\bar P$ means the interaction is attractive in that channel',fontsize=17)
plt.tight_layout(rect=[0,0,1,0.91]); plt.savefig('/home/phd25imran/analysis/fig_fullsus.png',dpi=105,facecolor='white')

print(f'L={L}, delta=0, n={NFIX:g}:   full / bubble / vertex')
for c, lab in CH:
    for u in g.index:
        if u == 0: continue
        f = g.loc[u, (f"chi_{c}_full", "mean")]
        b = g.loc[u, (f"bub_{c}", "mean")]
        v = g.loc[u, (f"chi_{c}_vertex", "mean")]
        print(f'  {c:5s} U={u:g}:  full {f:7.2f}   bubble {b:7.2f}   vertex {v:+7.2f}')
