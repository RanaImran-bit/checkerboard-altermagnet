"""Delta_tot and S(Q) for the two new L=14 cells at U=2, delta = 0.1 and 0.3.

These close the gap that made the largest-lattice check of the small-U peak
possible at delta=0.2 and nowhere else. Provenance verified before extraction:
in.dat line 5 reads t0=-1, t1=+0.3, t2=-delta, tam=0 and mc2duph.f90 carries the
checkerboard parity branch, so these are our model rather than the sibling one.

Conventions match extract_dnk_fortran.py: Delta_tot = sum_k |n_up - n_dn| over the
(L+1)^2 output grid, and Delta_tot/N with N = L^2.
"""
import os, numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
L, N = 14, 196

print(f"{'delta':>7}{'Delta_tot':>12}{'err':>9}{'Dtot/N':>10}{'net':>10}{'snr':>8}"
      f"{'corr_dxy':>10}{'S(Q)':>10}{'S(Q)/N':>10}")
rows = []
for d in ("0.1", "0.3"):
    u = pd.read_csv(f"{D}/L14_U2_new/d{d}_n_up.dat", sep=r"\s+", skiprows=1,
                    header=None, names=["kx","ky","v","e"])
    w = pd.read_csv(f"{D}/L14_U2_new/d{d}_n_dn.dat", sep=r"\s+", skiprows=1,
                    header=None, names=["kx","ky","v","e"])
    dn  = u.v.values - w.v.values
    err = np.hypot(u.e.values, w.e.values)
    f   = np.sin(u.kx.values) * np.sin(u.ky.values)
    dtot = float(np.abs(dn).sum())
    s = np.loadtxt(f"{D}/L14_U2_new/d{d}_sdwz.dat", skiprows=1)
    i = int(np.argmax(s[:,2]))
    r = dict(L=L, n=1.0, U=2.0, delta=float(d), delta_tot=dtot,
             delta_tot_err=float(np.sqrt((err**2).sum())), dtot_N=dtot/N,
             net=float(dn.sum()), snr=float(np.abs(dn/err).mean()),
             corr_dxy=float(np.corrcoef(dn, f)[0,1]),
             S=s[i,2], S_over_N=s[i,2]/N, Qx=s[i,0]/np.pi, Qy=s[i,1]/np.pi)
    rows.append(r)
    print(f"{r['delta']:>7.1f}{r['delta_tot']:>12.4f}{r['delta_tot_err']:>9.4f}"
          f"{r['dtot_N']:>10.4f}{r['net']:>10.2e}{r['snr']:>8.1f}"
          f"{r['corr_dxy']:>10.3f}{r['S']:>10.4f}{r['S_over_N']:>10.5f}")
pd.DataFrame(rows).to_csv(f"{D}/L14_U2_new/summary.csv", index=False)

print("\n=== the small-U peak, now at four sizes (Delta_tot/N at U=2) ===")
f8 = pd.read_csv(f"{D}/fortran_L8_L10_L12_dedup.csv")
f8 = f8[np.isclose(f8.n,1.0) & (f8.U==2.0)]
dn14 = pd.read_csv(f"{D}/fortran_dnk.csv"); dn14 = dn14[np.isclose(dn14.n,1.0) & (dn14.L==14) & (dn14.U==2.0)]
print(f"{'delta':>7}{'L=8':>10}{'L=10':>10}{'L=12':>10}{'L=14':>10}   trend")
for dl in (0.1, 0.2, 0.3):
    row = [float(f8[(f8.L==L) & np.isclose(f8.delta,dl)].dtot_N.iloc[0]) for L in (8,10,12)]
    if abs(dl-0.2) < 1e-9:
        v14 = float(dn14[np.isclose(dn14.delta,dl)].dtot_N.iloc[0])
    else:
        v14 = [r for r in rows if abs(r['delta']-dl)<1e-9][0]['dtot_N']
    row.append(v14)
    tr = "falls at L=14" if row[3] < row[2] else "rises at L=14"
    print(f"{dl:>7.1f}" + "".join(f"{v:>10.4f}" for v in row) + f"   {tr}")
