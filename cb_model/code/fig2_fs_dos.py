import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.ndimage import minimum_filter, gaussian_filter1d
import os

t  = 1.0
tp = -0.3
deltas = [0.0, 0.4]
n_fill = 0.8

def bands(kx, ky, delta):
    e0 = -4*tp*np.cos(kx)*np.cos(ky)
    g  =  2*t *(np.cos(kx)+np.cos(ky))
    m  =  4*delta*np.sin(kx)*np.sin(ky)
    root = np.sqrt(g**2 + m**2)
    return e0+root, e0-root, root

EMPTY = '#bcd4ea'; FS = '#c0392b'
dos_colors = ['#1f77b4', '#c0392b']

# --- k grids ---
Lf = 500
kf = np.linspace(-np.pi, np.pi, Lf, endpoint=False)
KX, KY = np.meshgrid(kf, kf, indexing='ij')

fig = plt.figure(figsize=(17, 5.7), facecolor='white')
gs  = GridSpec(1, 3, width_ratios=[1, 1, 1.35], wspace=0.30, figure=fig)

# ---------- panels (a),(b): Fermi surface + van Hove ----------
for i, delta in enumerate(deltas):
    ax = fig.add_subplot(gs[0, i]); ax.set_facecolor('white')
    Ep, Em, root = bands(KX, KY, delta)
    allE = np.sort(np.concatenate([Ep.ravel(), Em.ravel()]))
    mu = allE[int(round((n_fill/2.0)*allE.size))-1]

    ax.contourf(KX, KY, Em, levels=[mu, Em.max()+1], colors=[EMPTY], zorder=1)
    ax.contour (KX, KY, Em, levels=[mu], colors=[FS], linewidths=3.0, linestyles='solid', zorder=6)

    gx, gy = np.gradient(Em, kf, kf); gmag = np.hypot(gx, gy)
    Exx, _   = np.gradient(gx, kf, kf); Exy, Eyy = np.gradient(gy, kf, kf)
    detH = Exx*Eyy - Exy*Exy
    cand = (gmag == minimum_filter(gmag, size=11)) & (gmag < 0.05) & (root > 0.20) & (detH < 0)
    ys, xs = np.where(cand); kept = []
    for yy, xx in zip(ys, xs):
        x, y = KX[yy, xx], KY[yy, xx]
        if all(np.hypot(x-px, y-py) > 0.4 for px, py in kept):
            kept.append((x, y))
    if kept:
        vx, vy = zip(*kept)
        ax.scatter(vx, vy, s=200, c='gold', edgecolors='black', marker='*', linewidth=1.5, zorder=10)

    ax.set_xlabel(r'$k_x$', fontsize=22)
    if i == 0: ax.set_ylabel(r'$k_y$', fontsize=22)
    ax.set_xticks([-np.pi,0,np.pi]); ax.set_yticks([-np.pi,0,np.pi])
    ax.set_xticklabels([r'$-\pi$','0',r'$\pi$'], fontsize=18)
    ax.set_yticklabels([r'$-\pi$','0',r'$\pi$'], fontsize=18)
    ax.set_xlim(-np.pi,np.pi); ax.set_ylim(-np.pi,np.pi); ax.set_box_aspect(1)
    ax.set_title(rf"$\delta = {delta:.1f}$", fontsize=22, pad=8)
    ax.text(0.04, 0.97, f"({chr(97+i)})", transform=ax.transAxes, fontsize=18, fontweight='bold', va='top', ha='left')
    for sp in ax.spines.values(): sp.set_color('black'); sp.set_linewidth(1.6)

# ---------- panel (c): DOS ----------
axc = fig.add_subplot(gs[0, 2]); axc.set_facecolor('white')
N = 900
kk = np.linspace(-np.pi, np.pi, N, endpoint=False)
KXd, KYd = np.meshgrid(kk, kk)
def bands_dos(delta):
    e0 = -4*tp*np.cos(KXd)*np.cos(KYd)
    g  =  2*t *(np.cos(KXd)+np.cos(KYd))
    m  =  4*delta*np.sin(KXd)*np.sin(KYd)
    root = np.sqrt(g**2 + m**2)
    return e0+root, e0-root
bins = 600
for delta, c in zip(deltas, dos_colors):
    Ep, Em = bands_dos(delta)
    allE = np.concatenate([Ep.ravel(), Em.ravel()])
    mu = np.sort(allE)[int(round((n_fill/2.0)*allE.size))-1]
    dos, edges = np.histogram(allE, bins=bins, density=True); ec = 0.5*(edges[:-1]+edges[1:])
    dos = gaussian_filter1d(dos, 1.5)
    axc.fill_between(ec, dos, color=c, alpha=0.12)
    axc.plot(ec, dos, color=c, lw=2.6, label=rf'$\delta = {delta:.1f}$')
    axc.axvline(mu, color=c, ls='--', lw=1.6, alpha=0.8)
axc.set_xlim(-3.0, 5.4); axc.set_ylim(bottom=0)
axc.set_box_aspect(1/1.35)      # same panel HEIGHT as the square (a),(b) panels
axc.set_xlabel(r'Energy  $E/t$', fontsize=22)
axc.set_ylabel(r'DOS', fontsize=22)
axc.tick_params(labelsize=18)
from matplotlib.lines import Line2D
handles = [Line2D([0],[0], color=dos_colors[0], lw=2.6,           label=r'$\delta=0.0$'),
           Line2D([0],[0], color=dos_colors[1], lw=2.6,           label=r'$\delta=0.4$'),
           Line2D([0],[0], color='0.4',         lw=1.6, ls='--',  label=r'$E_F$')]
axc.legend(handles=handles, fontsize=15, frameon=False, loc='upper right')
axc.text(0.04, 0.97, '(c)', transform=axc.transAxes, fontsize=18, fontweight='bold', va='top', ha='left')
for sp in axc.spines.values(): sp.set_color('black'); sp.set_linewidth(1.6)

out=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","figures","fig2_fs_dos.png")
plt.savefig(out, dpi=110, bbox_inches='tight', facecolor='white'); print("saved", out)
