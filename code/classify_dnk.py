"""Symmetry classification of the momentum-resolved spin polarisation Dn(k).

Two dimensionless, basis-independent ratios:

    A_odd = sum_k |Dn(k) + Dn(R90 k)| / sum_k |Dn(k)|      0 => C4-ODD  => d-wave
    M_odd = sum_k |Dn(k) + Dn(M k)|   / sum_k |Dn(k)|      0 => d_x2-y2, 2 => d_xy
                                                           (M: (kx,ky) -> (ky,kx))

These replace projecting Dn onto a single lattice harmonic, which understates the
d-wave content badly: most of the weight sits in HIGHER d_xy harmonics, so a small
P_dxy does not mean "not d-wave". The ratios above use only the symmetry, so they
cannot be fooled that way.

Two traps this code exists to avoid:

  1. Duplicate k-points. n_up.dat has (L+1)^2 rows because +pi is repeated as -pi.
     Worse, float32 pi = 3.141592741 > np.pi, so neither `> np.pi` nor plain
     modular folding removes it. Points within 1e-4 of |pi| are snapped to -pi.

  2. Array-rotation shortcuts. np.rot90 is NOT the C4 partner map on an FFT-indexed
     mesh: it is off by one grid step, and returns 0.42-0.63 for pure d-wave
     harmonics that must return exactly 0. (It happens to be right on the
     half-shifted antiperiodic mesh, which makes it more dangerous, not less.)
     Partners are found here by explicit lookup on folded coordinates.

Usage:  python classify_dnk.py <maps.npz> [out.csv]
"""
import sys
import numpy as np


def snap(a):
    """Fold to (-pi, pi] with the float32 +pi -> -pi duplicate removed."""
    a = np.asarray(a, float)
    a = np.where(np.abs(np.abs(a) - np.pi) < 1e-4, -np.pi, a)
    return (a + np.pi) % (2 * np.pi) - np.pi


def classify(arr):
    """arr: (nk, 3) columns kx, ky, Dn. Returns (n_unique_k, A_odd, M_odd)."""
    kx, ky, v = snap(arr[:, 0]), snap(arr[:, 1]), arr[:, 2].astype(float)
    seen = {}
    for a, b, val in zip(kx, ky, v):
        seen[(round(a, 6), round(b, 6))] = val
    keys = list(seen.keys())
    V = np.array([seen[k] for k in keys])
    look = {k: i for i, k in enumerate(keys)}

    def partner(op):
        idx = np.empty(len(keys), dtype=int)
        for i, (a, b) in enumerate(keys):
            p = op(a, b)
            k = (round(float(snap(p[0])), 6), round(float(snap(p[1])), 6))
            if k not in look:
                return None
            idx[i] = look[k]
        return idx

    den = np.abs(V).sum()
    if den < 1e-12:                      # Dn == 0 identically (e.g. U = 0)
        return len(keys), np.nan, np.nan
    i4, im = partner(lambda a, b: (-b, a)), partner(lambda a, b: (b, a))
    ao = np.abs(V + V[i4]).sum() / den if i4 is not None else np.nan
    mo = np.abs(V + V[im]).sum() / den if im is not None else np.nan
    return len(keys), ao, mo


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "../data/fortran_L12_maps.npz"
    out = sys.argv[2] if len(sys.argv) > 2 else None
    d = np.load(src, allow_pickle=True)
    maps, meta, cols = d["maps"], d["meta"], [str(c) for c in d["cols"]]
    rows = []
    for arr, mt in zip(maps, meta):
        nk, ao, mo = classify(np.asarray(arr, float))
        rows.append(list(mt) + [nk, ao, mo])
    R = np.array(rows, float)
    hdr = cols + ["nk", "A_odd", "M_odd"]
    print("  ".join(f"{h:>8}" for h in hdr))
    for r in R:
        print("  ".join(f"{x:8.4f}" for x in r))
    ok = R[~np.isnan(R[:, -2])]
    if len(ok):
        print(f"\ninteracting cells: {len(ok)}")
        print(f"  A_odd {ok[:,-2].min():.4f} - {ok[:,-2].max():.4f}"
              f"   C4-odd (< 0.5): {int((ok[:,-2] < 0.5).sum())}/{len(ok)}")
        print(f"  M_odd {ok[:,-1].min():.4f} - {ok[:,-1].max():.4f}"
              f"   d_xy (> 1.0):   {int((ok[:,-1] > 1.0).sum())}/{len(ok)}")
    if out:
        np.savetxt(out, R, delimiter=",", header=",".join(hdr), comments="")
        print(f"\nwrote {out}")
