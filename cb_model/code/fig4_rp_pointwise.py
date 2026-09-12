"""Fig 4a: sharpness index R_p, POINT-WISE, in the (n, delta) plane.
RUN ON 251 (needs raw sdwz.dat in ~/L14_archive). Same R_p definition as the
group's published code; single panel because t' is fixed and delta = -tA is the knob."""
import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA_ROOT  = "/home/phd25imran/L14_archive"
OUTPUT_DIR = "/home/phd25imran/Checkerboard_Model"
TARGET_FILE = "sdwz.dat"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# NOTE: tA is NEGATIVE in our folders (tA = -delta), so allow the minus sign.
FOLDER_PATTERN = re.compile(
    r'L(?P<L>\d+)n(?P<n>[\d.]+)u(?P<u>[\d.]+)tA(?P<tA>-?[\d.]+)tt(?P<tt>[\d.-]+)N(?P<N>\d+)')

def find_nearest_indices(data, kx, ky, max_idx, num_nearest):
    distances = np.sqrt((data.iloc[:,0]-kx)**2 + (data.iloc[:,1]-ky)**2)
    distances.iloc[max_idx] = np.inf
    return distances.nsmallest(num_nearest).index

def calculate_rp(filepath):
    try:
        data = pd.read_csv(filepath, skiprows=1, header=None, sep=r'\s+')
        if data.empty or data.shape[1] < 3: return None
        val_col = 2
        max_idx = data[val_col].idxmax()
        kx, ky = data.iloc[max_idx,0], data.iloc[max_idx,1]
        max_val = data.iloc[max_idx, val_col]
        max_dist_x = data.iloc[:,0].abs().max(); max_dist_y = data.iloc[:,1].abs().max()
        kx_edge = np.isclose(np.abs(kx), max_dist_x); ky_edge = np.isclose(np.abs(ky), max_dist_y)
        num_nearest = 1 if (kx_edge and ky_edge) else (2 if (kx_edge or ky_edge) else 4)
        nearest_idx = find_nearest_indices(data, kx, ky, max_idx, num_nearest)
        avg_neighbor = data.iloc[nearest_idx, val_col].mean()
        return (max_val - avg_neighbor)/max_val if max_val != 0 else 0
    except Exception:
        return None

results = []
for folder in os.listdir(DATA_ROOT):
    m = FOLDER_PATTERN.match(folder)
    if not m or m.group('u') != "4.0": continue
    fp = os.path.join(DATA_ROOT, folder, "dir-kVals", TARGET_FILE)
    if not os.path.exists(fp): continue
    rp = calculate_rp(fp)
    if rp is not None:
        results.append({'n': float(m.group('n')), 'delta': -float(m.group('tA')), 'rp': rp})

df = pd.DataFrame(results)
print(f"{len(df)} points, R_p range [{df.rp.min():.2f}, {df.rp.max():.2f}]")

vmin, vmax = df.rp.min(), df.rp.max()
fig, ax = plt.subplots(figsize=(9, 6.2), facecolor='white')
sc = ax.scatter(df.n, df.delta, c=df.rp, s=320, marker='*', cmap='magma',
                edgecolors='black', linewidths=1.0, vmin=vmin, vmax=vmax)
ax.set_xlabel(r'Filling $n$', fontsize=18); ax.set_ylabel(r'Anisotropy $\delta$', fontsize=18)
ax.set_xlim(0.47, 1.02); ax.set_ylim(-0.03, 0.43); ax.tick_params(labelsize=14)
cb = fig.colorbar(sc, ax=ax, pad=0.02); cb.set_label(r'Sharpness $R_p$', fontsize=16); cb.ax.tick_params(labelsize=13)
for sp in ax.spines.values(): sp.set_linewidth(1.5)
plt.tight_layout()
out = os.path.join(OUTPUT_DIR, "Rp_pointwise_L14_U4.0.png")
plt.savefig(out, dpi=600, bbox_inches='tight', facecolor='white'); print("saved", out)
