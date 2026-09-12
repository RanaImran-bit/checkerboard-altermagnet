import matplotlib; matplotlib.use("Agg")
import numpy as np, glob, os, pandas as pd, matplotlib.pyplot as plt

SQDIR='/home/phd25imran/analysis/data/sq'
L, NTARGET = 12, 1.0

rows=[]
for f in sorted(glob.glob(os.path.join(SQDIR,'sq_*.npz'))):
    z=np.load(f,allow_pickle=True); cols=[str(c) for c in z['cols']]
    ix={c:i for i,c in enumerate(cols)}; meta=z['meta']; Sq=z['Sq']
    for k in range(len(meta)):
        r=meta[k]; LL=int(r[ix['L']]); S=Sq[k]; h=LL//2
        i,j=np.unravel_index(np.argmax(S),S.shape)
        rows.append(dict(L=LL,U=float(r[ix['U']]),delta=float(r[ix['delta']]),
                         n=float(r[ix['n']]),seed=int(r[ix['seed']]),
                         Spipi=float(S[h,h]), Sstar=float(S[i,j]),
                         qx=2*i/LL, qy=2*j/LL))
D=pd.DataFrame(rows)
D['m_pipi']=np.sqrt(np.maximum(D.Spipi,0)/D.L**2)
D['m_star']=np.sqrt(np.maximum(D.Sstar,0)/D.L**2)

def build(col):
    g=D.groupby(['L','U','n','delta'])[col].agg(['mean','sem']).reset_index() \
       .rename(columns={'mean':'m','sem':'e'})
    u0=g[g.U==0][['L','n','delta','m','e']].rename(columns={'m':'m0','e':'e0'})
    g=g.merge(u0,on=['L','n','delta']); g['dm']=g.m-g.m0; g['err']=np.hypot(g.e,g.e0)
    return g
Gp, Gs = build('m_pipi'), build('m_star')

ns=np.sort(D.n.unique()); N=ns[np.argmin(abs(ns-NTARGET))]
sp=Gp[(Gp.L==L)&np.isclose(Gp.n,N)]; ss=Gs[(Gs.L==L)&np.isclose(Gs.n,N)]
US=[2.,4.,6.,8.]; DS=np.sort(sp.delta.unique())

fig,ax=plt.subplots(1,3,figsize=(16.5,4.8),facecolor='white')
for u,c in zip(US, plt.cm.plasma(np.linspace(0,.8,len(US)))):
    t=sp[np.isclose(sp.U,u)].sort_values('delta')
    ax[0].errorbar(t.delta,t.dm,t.err,marker='s',color=c,lw=1.7,ms=5,capsize=2,label=f'$U={u:g}$')
    t=ss[np.isclose(ss.U,u)].sort_values('delta')
    ax[1].errorbar(t.delta,t.dm,t.err,marker='o',color=c,lw=1.7,ms=5,capsize=2,label=f'$U={u:g}$')
ax[0].set_title(r'(a)  $\Delta m$ at fixed $(\pi,\pi)$',fontsize=15)
ax[1].set_title(r'(b)  $\Delta m$ at the true peak $q^*$',fontsize=15)
for a in ax[:2]:
    a.set_xlabel(r'$\delta$',fontsize=16); a.set_ylabel(r'$\Delta m$',fontsize=15)
    a.legend(fontsize=10); a.tick_params(labelsize=11)

sub=D[(D.L==L)&np.isclose(D.n,N)&(D.U>0)]
q=sub.groupby('delta').apply(lambda t:(t.Spipi/t.Sstar).mean(),include_groups=False)
ax[2].plot(q.index,q.values,marker='D',color='crimson',lw=2,ms=7)
ax[2].axhline(1,color='k',lw=0.9,ls=':')
ax[2].axhline(0.95,color='grey',lw=0.9,ls='--')
ax[2].set_xlabel(r'$\delta$',fontsize=16)
ax[2].set_ylabel(r'$S(\pi,\pi)\,/\,S(q^*)$',fontsize=15)
ax[2].set_title(r'(c)  is $(\pi,\pi)$ still the peak?',fontsize=15)
ax[2].tick_params(labelsize=11)
fig.suptitle(rf'Ordering wavevector, $L={L}$, $n={N:.3f}$',fontsize=17)
plt.tight_layout(rect=[0,0,1,0.93]); plt.savefig('/home/phd25imran/analysis/fig_qstar.png',dpi=105,facecolor='white')

print(f'L={L} n={N:.3f}   fraction of magnetic weight at (pi,pi), and modal q*')
for d in DS:
    t=sub[np.isclose(sub.delta,d)]
    mode=t.groupby(['qx','qy']).size().idxmax()
    print(f'  delta={d:g}:  S(pi,pi)/S(q*) = {(t.Spipi/t.Sstar).mean():.3f}   modal q* = ({mode[0]:.2f}pi,{mode[1]:.2f}pi)')
print()
for u in US:
    a=sp[np.isclose(sp.U,u)].sort_values('delta').dm.values
    b=ss[np.isclose(ss.U,u)].sort_values('delta').dm.values
    print(f'U={u:g}: dm(pi,pi) {a[0]:.4f}->{a[-1]:.4f} (x{a[0]/a[-1]:.1f})   dm(q*) {b[0]:.4f}->{b[-1]:.4f} (x{b[0]/b[-1]:.1f})')
