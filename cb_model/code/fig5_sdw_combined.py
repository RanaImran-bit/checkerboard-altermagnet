"""Fig 5: SDW momentum space (a)-(c) + real space (d), one row. RUN ON 251.
Three representative points from the phase diagram (Fig. 3), all at delta = 0.20,
so decreasing filling walks Neel -> spiral -> diagonal:
    Neel      n=0.929  delta=0.20  Q=(pi,pi)
    Spiral    n=0.847  delta=0.20  Q=(0.86pi, pi)
    Diagonal  n=0.745  delta=0.20  Q=(0.57pi, 0.57pi)
(a)-(c): S^z(k) maps.   (d): C_SDW(i) along x for the three points.
Layout: 5 gridspec columns; column 3 is an INVISIBLE spacer so the (c)-(d) gap is
wider than the (a)-(b)-(c) gaps and panel (d) can be widened independently.
Folders in ~/L14_archive have tA = -delta (negative) and tt = 0.3."""
import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter

plt.rcParams.update({'font.family':'DejaVu Sans','figure.facecolor':'white','axes.facecolor':'white'})

base_path  = "/home/phd25imran/L14_archive"
output_dir = "/home/phd25imran/Checkerboard_Model"
L, U, TT   = 14, 4.0, 0.3
pi = np.pi
n_show = L // 2 + 1                       # real-space distances i = 0 .. L/2  (= 8 columns 0..7)
os.makedirs(output_dir, exist_ok=True)

# (label, filling n, anisotropy delta, marker, color)
POINTS = [
    ("Néel",     0.929, 0.20, 'o', 'blue'),
    ("Spiral",   0.847, 0.20, 's', 'red'),
    ("Diagonal", 0.745, 0.20, '^', 'green'),
]

# tA in the folder name is NEGATIVE (tA = -delta); the '-?' allows the minus sign.
FOLDER_PATTERN = re.compile(r'L(\d+)n([\d.]+)u([\d.]+)tA(-?[\d.]+)tt([\d.-]+)N(\d+)')

def find_folder(base, L, U, n, delta, tt, tol=0.02):
    """Best folder matching (n, delta, tt) with tA = -delta. Returns folder name or None."""
    target_tA = -delta
    best, best_score = None, 1e9
    for f in os.listdir(base):
        m = FOLDER_PATTERN.match(f)
        if not m:
            continue
        if int(m.group(1)) != L or abs(float(m.group(3)) - U) > 1e-6:
            continue
        nv, tAv, ttv = float(m.group(2)), float(m.group(4)), float(m.group(5))
        if abs(nv - n) > tol or abs(tAv - target_tA) > tol or abs(ttv - tt) > tol:
            continue
        score = abs(nv - n) + abs(tAv - target_tA) + abs(ttv - tt)
        if score < best_score and os.path.isfile(os.path.join(base, f, 'dir-kVals', 'sdwz.dat')):
            best_score, best = score, f
    return best

def read_grid(path):
    """S^z(k) on the (kx, ky) grid. Returns (kx_unique, ky_unique, Z[nx, ny])."""
    d = pd.read_csv(path, sep=r'\s+', header=None, skiprows=1)
    xu, yu = np.unique(d[0].values), np.unique(d[1].values)
    xi = {v: i for i, v in enumerate(xu)}
    yi = {v: i for i, v in enumerate(yu)}
    z = np.full((len(xu), len(yu)), np.nan)
    for k in range(len(d)):
        z[xi[d.iloc[k, 0]], yi[d.iloc[k, 1]]] = d.iloc[k, 2]
    return xu, yu, z

def deal_sdw_function(folderpath):
    """Real-space correlation C_SDW(i,j) = (1/N) sum_k S^z(k) cos(kx*i + ky*j).
    Same routine as the previous paper (drops the duplicate (pi,-pi) corner)."""
    fp = os.path.join(folderpath, 'sdwz.dat')
    data = pd.read_csv(fp, header=None, delimiter=r"\s+", skiprows=1)
    fdat = data[(data[0] != 3.141592741) | (data[1] != -3.141592741)]
    kx, ky, sk = np.array(fdat.iloc[:, 0]), np.array(fdat.iloc[:, 1]), np.array(fdat.iloc[:, 2])
    n = int(np.round(np.sqrt(len(fdat))))
    F = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            F[i, j] = np.sum(sk * np.cos(kx * i + ky * j)) / (n * n)
    return F

# ---- resolve each point's folder ONCE (used by both the map and the curve) ----
folders = []
for (order, n_val, delta, mk, col) in POINTS:
    f = find_folder(base_path, L, U, n_val, delta, TT)
    if f is None:
        print(f"MISSING {order}: n={n_val}, delta={delta}")
    else:
        print(f"{order}: {f}")
    folders.append(f)

# ---- figure: 5 columns (3 maps, invisible spacer, wide panel d) ----
fig, axes = plt.subplots(1, 5, figsize=(20, 6.0),
                         gridspec_kw={'width_ratios': [1, 1, 1, 0.15, 1.15]})  # (d) column wider
plt.subplots_adjust(left=0.05, right=0.99, top=0.90, bottom=0.13, wspace=0.42)
axes[3].axis('off')                       # invisible spacer between (c) and (d)
lab = ['(a)', '(b)', '(c)', '(d)']
x = np.arange(n_show)

# (a)-(c): momentum-space S^z(k) maps
for c, (order, n_val, delta, mk, col) in enumerate(POINTS):
    ax, folder = axes[c], folders[c]
    if folder is None:
        ax.text(0.5, 0.5, 'missing', transform=ax.transAxes, ha='center')
        continue
    sp = os.path.join(base_path, folder, 'dir-kVals', 'sdwz.dat')
    xu, yu, z = read_grid(sp)
    cf = ax.contourf(xu, yu, z.T, levels=50, cmap='RdBu_r', extend='both')   # z.T -> (ny, nx)
    cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(cf, cax=cax); cb.ax.tick_params(labelsize=16)
    cb.formatter = FormatStrFormatter('%.2f'); cb.update_ticks()
    ax.set_aspect('equal')
    ax.set_xticks([-pi, 0, pi]); ax.set_yticks([-pi, 0, pi])
    ax.set_xticklabels([r'$-\pi$', '0', r'$\pi$'], fontsize=20)
    ax.set_yticklabels([r'$-\pi$', '0', r'$\pi$'], fontsize=20)
    ax.set_xlabel(r'$k_x$', fontsize=22)
    if c == 0:
        ax.set_ylabel(r'$k_y$', fontsize=22)
    ax.set_title(rf"{lab[c]} $n={n_val:.3f},\ \delta={delta:.2f}$", fontsize=20)

# (d): real-space C_SDW(i) along x for all three points
axd = axes[4]
axd.set_box_aspect(0.816)                 # (d) same HEIGHT as (a)-(c) maps, ~23% wider
for (order, n_val, delta, mk, col), folder in zip(POINTS, folders):
    if folder is None:
        continue
    kdir = os.path.join(base_path, folder, 'dir-kVals')
    Cx = np.real(deal_sdw_function(kdir))[0:n_show, 0]                       # j = 0  -> along x
    axd.plot(x, Cx, marker=mk, ms=9, lw=1.7, color=col, markeredgecolor='black',
             markeredgewidth=0.5, label=rf"$n={n_val:.3f}$", zorder=5)
axd.axhline(0, color='gray', ls='--', lw=1.0, alpha=0.6)
axd.set_xlim(-0.3, n_show - 0.7); axd.set_xticks(range(n_show))
axd.tick_params(direction='in', labelsize=15, top=True, right=True)
axd.set_xlabel(r'$i$', fontsize=18); axd.set_ylabel(r'$C_{\mathrm{SDW}}(i)$', fontsize=18)
axd.set_title(rf"{lab[3]} $\delta=0.20$", fontsize=20)
axd.legend(frameon=False, fontsize=13, loc='upper right')

out = os.path.join(output_dir, f"5_sdw_combined_L{L}.png")
fig.savefig(out, dpi=600, bbox_inches='tight', facecolor='white')
print("saved", out)
