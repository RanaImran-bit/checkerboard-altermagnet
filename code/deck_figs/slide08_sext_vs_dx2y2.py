import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
OUT='/home/phd25imran/analysis'
P='/home/phd25imran/analysis/data/'
GRN,BLU='#2ca02c','#1f77b4'

d=pd.read_csv(P+'pairing_master.csv'); d=d[d.L==12]; d['N']=d.L**2
for c in ['chi_sext','chi_d']: d[c+'_n']=d[c]/d.N
f=pd.read_csv(P+'pairfull_delta0.csv'); f=f[f.L==12]
for c in ['sext','d']:
    f['bub_'+c]=f[f'chi_{c}_full']-f[f'chi_{c}_vertex']
    f['rat_'+c]=f[f'chi_{c}_full']/f['bub_'+c]

fig,ax=plt.subplots(1,4,figsize=(19,4.3),facecolor='white')

# (a) vs delta at half filling, all U
s=d[np.isclose(d.n,1.0)&(d.U>0)]
for U,ls in zip([2.,4.,6.,8.],['-',':','--','-.']):
    t=s[np.isclose(s.U,U)].groupby('delta')[['chi_sext_n','chi_d_n']].mean()
    ax[0].plot(t.index,t.chi_sext_n,ls,color=GRN,lw=2,marker='s',ms=5)
    ax[0].plot(t.index,t.chi_d_n,ls,color=BLU,lw=2,marker='^',ms=5)
ax[0].set_xlabel(r'$\delta$',fontsize=15); ax[0].set_ylabel(r'vertex $\chi/N$',fontsize=14)
ax[0].set_title('(a)  both decay with anisotropy'+'\n'+r'(line style: $U=2,4,6,8$)',fontsize=13)
ax[0].plot([],[],color=GRN,lw=2,marker='s',label='extended $s$')
ax[0].plot([],[],color=BLU,lw=2,marker='^',label=r'$d_{x^2-y^2}$')
ax[0].legend(fontsize=11); ax[0].axhline(0,color='k',lw=.9,ls=':')

# (b) vs filling at U=4, several delta
s2=d[np.isclose(d.U,4)]
for dd,al in zip([0.0,0.2,0.4,0.7],[1.0,.75,.5,.3]):
    t=s2[np.isclose(s2.delta,dd)].groupby('n')[['chi_sext_n','chi_d_n']].mean()
    ax[1].plot(t.index,t.chi_sext_n,color=GRN,alpha=al,lw=2,marker='s',ms=5)
    ax[1].plot(t.index,t.chi_d_n,color=BLU,alpha=al,lw=2,marker='^',ms=5)
ax[1].axhline(0,color='k',lw=1.1,ls='--')
ax[1].set_xlabel('Filling $n$',fontsize=15)
ax[1].set_title(r'(b)  both turn attractive near $n=1$'+'\n'+r'(fading: $\delta=0,0.2,0.4,0.7$)',fontsize=13)

# (c) which leads: difference
for dd,al in zip([0.0,0.2,0.4,0.7],[1.0,.75,.5,.3]):
    t=s2[np.isclose(s2.delta,dd)].groupby('n')[['chi_sext_n','chi_d_n']].mean()
    ax[2].plot(t.index,t.chi_sext_n-t.chi_d_n,color='#6D2E46',alpha=al,lw=2,marker='o',ms=5,
               label=rf'$\delta={dd:g}$')
ax[2].axhline(0,color='k',lw=1.1,ls='--')
ax[2].set_xlabel('Filling $n$',fontsize=15)
ax[2].set_ylabel(r'$\chi_{s\text{-}ext}-\chi_{d_{x^2-y^2}}$',fontsize=13)
ax[2].set_title(r'(c)  extended $s$ leads only near $n=1$',fontsize=14); ax[2].legend(fontsize=10,loc='lower right')

# (d) P/Pbar at delta=0
g=f[np.isclose(f.n,1.0)].groupby('U')[['rat_sext','rat_d']].agg(['mean','sem'])
ax[3].errorbar(g.index,g[('rat_sext','mean')],g[('rat_sext','sem')],color=GRN,lw=2,
               marker='s',ms=7,capsize=3,label='extended $s$')
ax[3].errorbar(g.index,g[('rat_d','mean')],g[('rat_d','sem')],color=BLU,lw=2,
               marker='^',ms=7,capsize=3,label=r'$d_{x^2-y^2}$')
ax[3].axhline(1,color='k',lw=1.1,ls='--')
ax[3].set_xlabel('$U$',fontsize=15); ax[3].set_ylabel(r'$P/\bar{P}$',fontsize=14)
ax[3].set_title(r'(d)  $\delta=0$: they cross near $U=6$',fontsize=14); ax[3].legend(fontsize=11)
for a in ax: a.tick_params(labelsize=11)
plt.tight_layout(); plt.savefig(f'{OUT}/fig_sd.png',dpi=110,facecolor='white')

print('crossover: extended s minus dx2-y2 at half filling, per delta (U=4)')
for dd in sorted(s2.delta.unique()):
    t=s2[np.isclose(s2.delta,dd)].groupby('n')[['chi_sext_n','chi_d_n']].mean()
    print(f'  d={dd:g}: sext={t.chi_sext_n.loc[1.0]:+.4f}  d={t.chi_d_n.loc[1.0]:+.4f}  diff={t.chi_sext_n.loc[1.0]-t.chi_d_n.loc[1.0]:+.4f}')
print()
print('P/Pbar at half filling, delta=0:')
for U in g.index: print(f'  U={U:g}: sext={g.loc[U,("rat_sext","mean")]:.3f}  d={g.loc[U,("rat_d","mean")]:.3f}')
