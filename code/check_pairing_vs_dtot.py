"""Do the d channels actually track the magnetic polarisation, or just track delta?

Task 1 asked for the two d channels plotted "against the polarization Delta_tot".
Using Delta_tot as a background colour only answers that by eye. Both the pairing
vertex and Delta_tot are functions of (U, delta), so a raw correlation between
them is largely the shared dependence on delta and proves nothing.

The test that means something is a PARTIAL correlation: strip the linear
dependence on delta (and then on U as well) out of both quantities and correlate
what is left. If pairing tracks the magnetism for a physical reason, the
correlation survives. If it was only the common dependence on geometry, it dies.
"""
import os, numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def resid(y, X):
    """Residual of y after least-squares removal of the columns of X (plus a constant)."""
    A = np.column_stack([np.ones(len(y))] + list(X))
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ beta


def report(name, df, cols):
    d_, x_, m_ = df[cols[0]].values, df[cols[1]].values, df["dtot_N"].values
    dl, U = df["delta"].values, df["U"].values
    print(f"\n--- {name}  ({len(df)} cells) ---")
    print(f"{'':>26}{'d_x2-y2':>12}{'d_xy':>12}")
    r0 = [np.corrcoef(v, m_)[0,1] for v in (d_, x_)]
    print(f"{'raw corr with Delta_tot':>26}{r0[0]:>12.3f}{r0[1]:>12.3f}")
    r1 = [np.corrcoef(resid(v, [dl]), resid(m_, [dl]))[0,1] for v in (d_, x_)]
    print(f"{'partial, delta removed':>26}{r1[0]:>12.3f}{r1[1]:>12.3f}")
    r2 = [np.corrcoef(resid(v, [dl, U]), resid(m_, [dl, U]))[0,1] for v in (d_, x_)]
    print(f"{'partial, delta and U out':>26}{r2[0]:>12.3f}{r2[1]:>12.3f}")
    return r0, r1, r2


f = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n, 1.0) & (f.U >= 2)]
v = pd.read_csv(os.path.join(D, "chi_L12_vertex_summary.csv"))
v = v[v.U >= 2].merge(f[f.L == 12][["U", "delta", "dtot_N"]], on=["U", "delta"])

print("=" * 62)
print("Does pairing track the magnetism, beyond their shared dependence on delta?")
print("=" * 62)
for L in sorted(f.L.unique()):
    report(f"equal time, L={L}", f[f.L == L], ["d", "dxy"])
report("time-integrated vertex, L=12", v, ["d", "dxy"])

print("\n\n=== at FIXED delta, does d_xy rise with Delta_tot across U? ===")
print("(within a column both are driven by U alone, so this isolates the U axis)")
print(f"{'delta':>7}{'eq-time r':>12}{'vertex r':>12}   n")
for dl in sorted(v.delta.unique()):
    a = f[(f.L == 12) & np.isclose(f.delta, dl)].sort_values("U")
    b = v[np.isclose(v.delta, dl)].sort_values("U")
    ra = np.corrcoef(a.dxy, a.dtot_N)[0,1]
    rb = np.corrcoef(b.dxy, b.dtot_N)[0,1]
    print(f"{dl:>7.1f}{ra:>12.3f}{rb:>12.3f}   {len(b)}")
