import matplotlib; matplotlib.use("Agg")
import numpy as np, glob, os, matplotlib.pyplot as plt

OUT='/home/phd25imran/analysis'
SQDIRS=[f'{D}/lo', f'{D}/hi2']
L, U, NTARGET = 10, 6.0, 1.0

recs=[]
for Dd in SQDIRS:
    for f in sorted(glob.glob(os.path.join(Dd,'sq_*.npz'))):
        z=np.load(f,allow_pickle=True)
        cols=[str(c) for c in z['cols']]; meta=z['meta']; Sq=z['Sq']
        ix={c:i for i,c in enumerate(cols)}
        for k in range(len(meta)):
            recs.append(dict(L=int(meta[k][ix['L']]),U=float(meta[k][ix['U']]),
                             delta=float(meta[k][ix['delta']]),
                             n=float(meta[k][ix['n']]),Sq=Sq[k]))
inv={}
for r in recs: inv.setdefault((r['L'],r['U']),set()).add(r['delta'])
print('AVAILABLE (L,U) -> deltas:')
for k in sorted(inv): print(f'   L={k[0]:2d} U={k[1]:g}: {sorted(inv[k])}')

sel=[r for r in recs if r['L']==L and np.isclose(r['U'],U)]
if not sel: raise SystemExit(f'no grids for L={L} U={U:g}')
ns=sorted({r['n'] for r in sel}); N=min(ns,key=lambda x:abs(x-NTARGET))
sel=[r for r in sel if np.isclose(r['n'],N)]
DS=sorted({r['delta'] for r in sel})
avg={d:np.mean([r['Sq'] for r in sel if np.isclose(r['delta'],d)],axis=0) for d in DS}

vmax=max(v.max() for v in avg.values())
step=2.0/L                                   # cell width in units of pi
ext=[-step/2, 2-step/2, -step/2, 2-step/2]   # q/pi from 0 to 2, cell-centred

fig,ax=plt.subplots(1,len(DS),figsize=(2.45*len(DS),3.3),facecolor='white',squeeze=False)
for a,d in zip(ax[0],DS):
    S=avg[d]
    im=a.imshow(S.T,origin='lower',extent=ext,cmap='inferno',vmin=0,vmax=vmax,
                interpolation='nearest',aspect='equal')
    i,j=np.unravel_index(np.argmax(S),S.shape)
    a.plot(2*i/L,2*j/L,marker='x',color='cyan',ms=10,mew=2.5)
    a.plot(1,1,marker='+',color='w',ms=9,mew=1.4)      # (pi,pi) for reference
    a.set_title(rf'$\delta={d:g}$',fontsize=14)
    a.set_xticks([0,1,2]); a.set_yticks([0,1,2])
    a.set_xticklabels(['0',r'$\pi$',r'$2\pi$'],fontsize=10)
    a.set_yticklabels(['0',r'$\pi$',r'$2\pi$'] if a is ax[0,0] else [],fontsize=10)
    a.set_xlabel(r'$q_x$',fontsize=12)
ax[0,0].set_ylabel(r'$q_y$',fontsize=12)
fig.colorbar(im,ax=ax[0].tolist(),fraction=0.015,pad=0.015,label=r'$S^{zz}(q)$')
fig.suptitle(rf'Spin structure factor, $L={L}$, $U={U:g}$, $n={N:.3f}$   '
             r'($\times$ peak,  $+$ is $(\pi,\pi)$)',fontsize=14)
plt.savefig('/home/phd25imran/analysis/fig_sqmap3.png',dpi=110,facecolor='white',bbox_inches='tight')

h=L//2
print('\ndelta   peak q             S_max    S(pi,pi)  ratio')
for d in DS:
    S=avg[d]; i,j=np.unravel_index(np.argmax(S),S.shape)
    print(f' {d:4.1f}   ({2*i/L:.2f}pi,{2*j/L:.2f}pi)    {S[i,j]:.4f}   {S[h,h]:.4f}   {S[h,h]/S[i,j]:.3f}')
