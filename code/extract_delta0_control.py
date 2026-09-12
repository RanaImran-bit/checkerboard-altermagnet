"""Delta_tot at delta = 0: the sampling noise floor at production settings.

At delta = 0 both checkerboard diagonals equal t1 on every site, the A/B distinction
disappears, and the model is the plain t-t' square lattice. Its Neel sublattices are
related by TRANSLATION, so n_up(k) = n_dn(k) identically and Delta_tot must vanish by
symmetry. Whatever is measured is therefore the statistical floor.

The U = 0 column cannot serve as this control: there the Green function is exact and
the zero is arithmetic, not statistical. This is the only genuine null we have.

Duplicate k-points: n_up.dat has (L+1)^2 rows because +pi repeats as -pi, and
float32 pi = 3.141592741 defeats both a naive `> np.pi` test and modular folding.
Points within 1e-4 of |pi| are snapped to -pi before deduplication.

  python extract_delta0_control.py <dir with U*/n_up.dat,n_dn.dat> [out.csv]
"""
import sys, os, glob
import numpy as np
import pandas as pd


def snap(a):
    a = np.asarray(a, float)
    a = np.where(np.abs(np.abs(a) - np.pi) < 1e-4, -np.pi, a)
    return (a + np.pi) % (2 * np.pi) - np.pi


def cell(d):
    up = np.loadtxt(os.path.join(d, "n_up.dat"), skiprows=1)
    dn = np.loadtxt(os.path.join(d, "n_dn.dat"), skiprows=1)
    seen = {}
    for a, b, x, y, e1, e2 in zip(snap(up[:, 0]), snap(up[:, 1]),
                                  up[:, 2], dn[:, 2], up[:, 3], dn[:, 3]):
        seen[(round(a, 6), round(b, 6))] = (x, y, e1, e2)
    V = np.array(list(seen.values()))
    dtot = float(np.abs(V[:, 0] - V[:, 1]).sum())
    err = float(np.sqrt((np.hypot(V[:, 2], V[:, 3]) ** 2).sum()))
    return len(seen), dtot, err


def main(src, out=None):
    rows = []
    for d in sorted(glob.glob(os.path.join(src, "U*")),
                    key=lambda s: float(os.path.basename(s)[1:])):
        U = float(os.path.basename(d)[1:])
        try:
            nk, dtot, err = cell(d)
        except Exception as e:
            print(f"U={U}: {e}"); continue
        rows.append(dict(U=U, nk=nk, dtot=dtot, dtot_err=err,
                         dtot_N=dtot / nk, dtot_N_err=err / nk))
    df = pd.DataFrame(rows)
    print("delta = 0 (plain t-t' square lattice: altermagnetism FORBIDDEN by symmetry)")
    print("L = 12, half filling, Fortran production settings\n")
    print(f"{'U':>5}{'n_k':>5}{'Delta_tot/N':>14}{'+/-':>11}")
    for _, r in df.iterrows():
        print(f"{r.U:5.1f}{int(r.nk):5d}{r.dtot_N:14.6f}{r.dtot_N_err:11.6f}")
    print(f"\nNOISE FLOOR = {df.dtot_N.mean():.6f} +/- {df.dtot_N.std(ddof=1):.6f} "
          f"(n={len(df)} cells)")
    if out:
        df.to_csv(out, index=False); print(f"wrote {out}")
    return df


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0] if a else "../data/delta0", a[1] if len(a) > 1 else None)
