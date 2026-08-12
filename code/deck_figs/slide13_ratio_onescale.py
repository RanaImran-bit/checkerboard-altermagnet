import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt

CSV='/home/phd25imran/analysis/data/pairfull_delta0.csv'
L, UFIX, NFIX = 12, 4.0, 1.0
CH=[('son',r'on-site $s$','#7f7f7f','o'),('sext',r'extended $s$','#2ca02c','s'),
    ('d',r'$d_{x^2-y^2}$','#1f77b4','^'),('dxy',r'$d_{xy}$','#d62728','D')]

d=pd.read_csv(CSV); d=d[d.L==L]
for c,_,_,_ in CH:
    d[f'bub_{c}']=d[f'chi_{c}_full']-d[f'chi_{c}_vertex']
    d[f'rat_{c}']=d[f'chi_{c}_full']/d[f'bub_{c}']

fig,ax=plt.subplots(1,2,figsize=(14,5.6),facecolor='white')

s=d[np.isclose(d.n,NFIX)]
g=s.groupby('U')[[f'rat_{c}' for c,_,_,_ in CH]].agg(['mean','sem'])
for c,lab,col,mk in CH:
    a=g.index.values
    ax[0].errorbar(a,g[(f'rat_{c}','mean')],g[(f'rat_{c}','sem')],
                   marker=mk,ms=8,lw=2,color=col,capsize=3,label=lab)
ax[0].set_xlabel('$U$',fontsize=17); ax[0].set_title(rf'at half filling $n={NFIX:g}$',fontsize=16)

s2=d[np.isclose(d.U,UFIX)]
g2=s2.groupby('n')[[f'rat_{c}' for c,_,_,_ in CH]].agg(['mean','sem'])
for c,lab,col,mk in CH:
    nn=g2.index.values
    ax[1].errorbar(nn,g2[(f'rat_{c}','mean')],g2[(f'rat_{c}','sem')],
                   marker=mk,ms=8,lw=2,color=col,capsize=3,label=lab)
ax[1].set_xlabel('Filling $n$',fontsize=17); ax[1].set_title(rf'at $U={UFIX:g}$',fontsize=16)

for a in ax:
    a.axhline(1,color='k',lw=1.4,ls='--')
    a.tick_params(labelsize=13); a.set_ylabel(r'$P/\bar{P}$',fontsize=17)
    a.legend(fontsize=12, loc='center right', framealpha=0.95)
    a.text(0.02, 0.97, 'above 1: interaction ATTRACTIVE', transform=a.transAxes,
           fontsize=11, va='top', color='darkred')
    a.text(0.02, 0.03, 'below 1: repulsive', transform=a.transAxes,
           fontsize=11, va='bottom', color='navy')
fig.suptitle(rf'All four pairing channels on one scale, $L={L}$, $\delta=0$',fontsize=18)
plt.tight_layout(rect=[0,0,1,0.93]); plt.savefig('/home/phd25imran/analysis/fig_ratio.png',dpi=105,facecolor='white')

print('P/Pbar at half filling:   (>1 attractive, <1 repulsive)')
hdr = ''.join(f'{c:>14}' for c, _, _, _ in CH)
print(f'{"U":>4}{hdr}')
for u in g.index:
    print(f'{u:4g}' + ''.join(f'{g.loc[u, (f"rat_{c}", "mean")]:14.3f}' for c, _, _, _ in CH))
print()
print('P/Pbar vs filling at U=4:')
print(f'{"n":>6}{hdr}')
for n in g2.index:
    print(f'{n:6.3f}' + ''.join(f'{g2.loc[n, (f"rat_{c}", "mean")]:14.3f}' for c, _, _, _ in CH))
