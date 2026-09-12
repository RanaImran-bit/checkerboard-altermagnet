"""Delta n(k) = n_up(k) - n_dn(k) at selected points in the (U, delta) plane.

The direct check: does the momentum-space spin polarisation have d-wave
symmetry, and which d? Each panel is one (L, U, delta) cell. The altermagnetic signature to look for is C4-ODDNESS: rotating k by 90 degrees
must FLIP the sign of Delta n(k). A ferromagnet would be C4-even (one sign
everywhere) and a conventional Neel antiferromagnet would have Delta n(k) = 0 at
every k, since its two sublattices are related by a translation.

The projections P_dxy and P_d are printed to stdout rather than written on the
panels, so the figure stays clean.

NO WHITE EDGES. Deduplicating the k-grid drops the +pi points (they repeat -pi),
so the data now stops at pi - 2pi/L and griddata returns NaN beyond it -- a
white band on the top and right. Delta n(k) is periodic, so the fix is to tile
the measured points into the eight neighbouring zone images before
interpolating; the [-pi, pi] box is then covered by real data everywhere.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy.interpolate import griddata

mpl.rcParams.update({
    'font.family':'serif','font.serif':['DejaVu Serif'],
    'mathtext.fontset':'dejavuserif','axes.linewidth':1.2,
    'xtick.direction':'out','ytick.direction':'out',
    'axes.labelsize':15,'xtick.labelsize':12,'ytick.labelsize':12})

D='/Users/liujiaxin/Desktop/checkerboard-altermagnet/data'
z=np.load(f'{D}/dnk_maps_dedup.npz',allow_pickle=True)
maps,meta=z['maps'],z['meta']
# spread across the plane: the Delta_tot ridge, a high-U cell, mid, and a dead corner
SEL=[(12,3,0.2),(12,5,0.4),(12,4,0.3),(12,4.5,0.7),(12,2,0.7),(12,3,0.1)]

fig,ax=plt.subplots(2,3,figsize=(16.5,10.2),facecolor='white')
vmax=max(np.abs(maps[np.where((meta[:,0]==L)&np.isclose(meta[:,1],U)&
        np.isclose(meta[:,2],dl))[0][0]][:,2]).max() for L,U,dl in SEL)
for n,(L,U,dl) in enumerate(SEL):
    a=ax.flat[n]
    i=np.where((meta[:,0]==L)&np.isclose(meta[:,1],U)&np.isclose(meta[:,2],dl))[0][0]
    kx,ky,dn=maps[i][:,0],maps[i][:,1],maps[i][:,2]
    # tile into the neighbouring zone images so the interpolation has data right
    # out to the boundary; Delta n(k) is periodic with period 2pi in each axis
    tx=np.concatenate([kx+a*2*np.pi for a in (-1,0,1) for _ in (-1,0,1)])
    ty=np.concatenate([ky+b*2*np.pi for _ in (-1,0,1) for b in (-1,0,1)])
    td=np.tile(dn,9)
    gx,gy=np.meshgrid(np.linspace(-np.pi,np.pi,240),np.linspace(-np.pi,np.pi,240))
    gz=griddata((tx,ty),td,(gx,gy),method='linear')
    assert not np.isnan(gz).any(), 'still NaN after tiling'

    im=a.contourf(gx,gy,gz,levels=60,cmap='RdBu_r',vmin=-vmax,vmax=vmax)
    gd=np.cos(kx)-np.cos(ky); gxy=np.sin(kx)*np.sin(ky)
    Pd=(dn*gd).sum()/np.abs(gd).sum(); Pxy=(dn*gxy).sum()/np.abs(gxy).sum()
    tot=np.abs(dn).sum()
    a.set_title(rf'$L={L}$, $U={U:g}$, $\delta={dl:g}$', fontsize=15)
    a.set_aspect('equal'); a.set_xticks([-np.pi,0,np.pi]); a.set_yticks([-np.pi,0,np.pi])
    a.set_xticklabels([r'$-\pi$','0',r'$\pi$']); a.set_yticklabels([r'$-\pi$','0',r'$\pi$'])
    if n>=3: a.set_xlabel(r'$k_x$')
    if n%3==0: a.set_ylabel(r'$k_y$')
    print(f'L={L} U={U:g} d={dl:g}: sum|dn|={tot:7.3f}  P_dxy={Pxy:+.5f} '
          f'({100*abs((dn*gxy).sum())/tot:4.1f}%)  P_d={Pd:+.5f} ({100*abs((dn*gd).sum())/tot:4.1f}%)')
cb=fig.colorbar(im,ax=ax,pad=0.02,fraction=0.022)
cb.set_label(r'$\Delta n(\mathbf{k}) = n_\uparrow(\mathbf{k})-n_\downarrow(\mathbf{k})$',fontsize=15)
plt.savefig(f'{D}/../figs/fig_dnk_points.png',dpi=600,bbox_inches='tight',facecolor='white')
