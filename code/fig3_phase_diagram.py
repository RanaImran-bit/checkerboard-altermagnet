"""Fig 3: magnetic phase diagram in the (n, delta) plane from CPQMC S^z(k).
RUN ON 251 (needs raw sdwz.dat in ~/L14_archive).
Every point is classified by its DOMINANT peak (Q = argmax S^z(q)) into one of the four
correlations: Neel (pi,pi), Stripe (pi,0)/(0,pi), Spiral (pi,q)/(q,pi), Diagonal (q,q)
(nearest category). Long- vs short-range is decided separately by R_p / C_SDW(r).
y-axis is the anisotropy delta = -tA."""
import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

def extract_extrema(sdwz_path):
    try:
        data = pd.read_csv(sdwz_path, sep=r'\s+', header=None, skiprows=1)
        if data.shape[1] < 3: return []
        kx = data.iloc[:,0].to_numpy(); ky = data.iloc[:,1].to_numpy(); sz = data.iloc[:,2].to_numpy()
        order = np.argsort(-np.abs(sz))
        return [(float(kx[i]), float(ky[i]), float(sz[i])) for i in order]
    except Exception as e:
        print(f"ERROR reading {sdwz_path}: {e}"); return []

def find_folders(base_path, L=14, U=4.0):
    pat = re.compile(r'L(?P<L>\d+)n(?P<n>[\d.]+)u(?P<u>[\d.]+)tA(?P<tA>-?[\d.]+)tt(?P<tt>[\d.]+)N(?P<N>\d+)')
    out = []
    for folder in os.listdir(base_path):
        m = pat.match(folder)
        if not m: continue
        if int(m.group('L')) != L or abs(float(m.group('u')) - U) > 1e-6: continue
        sdwz = os.path.join(base_path, folder, "dir-kVals", "sdwz.dat")
        if not os.path.isfile(sdwz): continue
        out.append({'folder': folder, 'n': float(m.group('n')),
                    'delta': -float(m.group('tA')), 'sdwz_path': sdwz})
    out.sort(key=lambda x: (x['n'], x['delta']))
    return out

base_path  = "/home/phd25imran/L14_archive"
output_dir = "/home/phd25imran/Checkerboard_Model"
L, U = 14, 4.0
os.makedirs(output_dir, exist_ok=True)

folders = find_folders(base_path, L=L, U=U)
print(f"found {len(folders)} folders")

rows = []
for info in folders:
    ex = extract_extrema(info['sdwz_path'])
    if ex: rows.append({'n': info['n'], 'delta': info['delta'], 'peak': ex[0]})  # global-max peak
df = pd.DataFrame(rows)

pi = np.pi
# Neel and Stripe get fixed colors; Spiral and Diagonal SHARE the q/pi colorbar
# (distinguished only by marker: triangle = spiral, star = diagonal).
style  = {"M": ("o","red"), "X": ("s","blue")}
labels = {"M": "Néel "   + r"$(\pi,\pi)$",
          "X": "Stripe " + r"$(\pi,0)/(0,\pi)$"}
LEG_GREEN = "green"                                                # representative swatch for q-colored markers
norm = Normalize(0.0, 1.0); cmap = plt.cm.viridis

def classify(kx, ky, tol=0.25):
    """Assign the dominant peak to the NEAREST of the four correlations.
    Returns (class, q) where q (in units of pi) is defined for Spiral and Diagonal."""
    qx = np.arccos(np.cos(kx)); qy = np.arccos(np.cos(ky))          # fold to [0, pi]
    if qx > pi-tol and qy > pi-tol:                    return "M", None   # Neel (pi,pi)
    if (qx > pi-tol and qy < tol) or (qx < tol and qy > pi-tol): return "X", None  # Stripe
    d_spiral = min(pi-qx, pi-qy)                                    # distance to the zone edge
    d_diag   = abs(qx-qy)/np.sqrt(2.0)                              # distance to the diagonal
    if d_spiral <= d_diag: return "N", min(qx, qy)/pi              # Spiral (pi,q), q = interior comp.
    return "Q", (qx+qy)/2.0/pi                                     # Diagonal (q,q), q = diagonal proj.

fig, ax = plt.subplots(figsize=(11, 6.5), facecolor='white'); ax.set_facecolor('white')
plotted = set(); qs = []
for _, row in df.iterrows():
    kx, ky, val = row['peak']
    cls, q = classify(kx, ky)
    if cls == "N":                                                 # spiral -> triangle, q-colored
        ax.scatter(row['n'], row['delta'], marker='^', color=cmap(norm(q)), s=130,
                   edgecolors='black', linewidths=0.8, zorder=3); qs.append(q); plotted.add("N")
    elif cls == "Q":                                               # diagonal -> star, SAME q-colorbar
        ax.scatter(row['n'], row['delta'], marker='*', color=cmap(norm(q)), s=200,
                   edgecolors='black', linewidths=0.8, zorder=3); qs.append(q); plotted.add("Q")
    else:
        mk, col = style[cls]
        ax.scatter(row['n'], row['delta'], marker=mk, color=col, s=115,
                   edgecolors='black', linewidths=0.6, zorder=3); plotted.add(cls)

ax.set_xlim(0.48, 1.0); ax.set_ylim(-0.02, 0.42)
ax.set_xlabel(r'Filling $n$', fontsize=16); ax.set_ylabel(r'Anisotropy $\delta$', fontsize=16)
ax.tick_params(labelsize=13)
for sp in ax.spines.values(): sp.set_color('black')

handles = []
for c in ["M", "X"]:
    if c in plotted:
        handles.append(plt.Line2D([0],[0], marker=style[c][0], color='w', markerfacecolor=style[c][1],
                                   markeredgecolor='black', markersize=11, label=labels[c]))
if "Q" in plotted:
    handles.append(plt.Line2D([0],[0], marker='*', color='w', markerfacecolor=LEG_GREEN,
                              markeredgecolor='black', markersize=14, label="Diagonal " + r"$(q,q)$"))
if "N" in plotted:
    handles.append(plt.Line2D([0],[0], marker='^', color='w', markerfacecolor=LEG_GREEN,
                              markeredgecolor='black', markersize=11, label="Spiral " + r"$(\pi,q)$"))
ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.85),
          title="Dominant SDW peak", frameon=True, fontsize=10)

if qs:
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.55, pad=0.02)
    cbar.set_label(r'ordering wavevector $q/\pi$', fontsize=12); cbar.ax.tick_params(labelsize=11)

plt.tight_layout()
out = os.path.join(output_dir, f"phase_diagram_L{L}_U{U}.png")
plt.savefig(out, dpi=600, bbox_inches='tight', facecolor='white'); print("saved", out)
