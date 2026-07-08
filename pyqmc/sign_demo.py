#!/usr/bin/env python3
"""Sign-problem demo: free projection (exact, exponential) vs constrained path
(biased, polynomial), and how a better trial slows the sign decay.

One CPMC step samples a discrete HS field per interaction term with overlap ratio
R. The ONLY difference between the two methods is how a candidate field is scored:
  - constrained path: q = max(R, 0)  -> walkers that cross the node are removed;
    the weight stays positive (no sign problem) at the cost of a trial-dependent bias.
  - free projection : q = |R|, and sign(R) is carried on a SIGNED weight -> exact,
    but the average sign <s> = sum w / sum|w| decays ~ exp(-gamma * tau * N), so the
    energy estimator E = sum(w E_L)/sum(w) has an error that BLOWS UP with tau.

We project e^{-tau H}|psi_T>, measuring E(tau) and <s>(tau) across many independent
populations so the free-projection error explosion is visible, and the CP curve
stays flat. With --release we equilibrate under CP then RELEASE the constraint and
watch the sign decay rate vs trial quality.

    python pyqmc/sign_demo.py --lx 4 --ly 4 --nup 5 --ndn 5 --U 6 --t2 0.4
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
from cpqmc import LatticeModel, TrialWF, Estimators, square_hopping


def hopping(lx, ly, t, t2):
    """NN square-lattice hopping, optionally + frustrating NNN t2 (breaks
    bipartiteness -> a real sign problem even near half filling)."""
    K = square_hopping(lx, ly, t)
    if t2:
        def idx(x, y): return (x % lx) * ly + (y % ly)
        for x in range(lx):
            for y in range(ly):
                i = idx(x, y)
                for (dx, dy) in ((1, 1), (1, -1)):
                    K[i, idx(x + dx, y + dy)] -= t2
                    K[idx(x + dx, y + dy), i] -= t2
    return K


def project(model, trial, est, nw, nsteps, constrained, rng, meas_every=2):
    """Project one population for nsteps slices; return per-measured-step
    (tau_index, sum_w_EL, sum_w, sum_abs_w)."""
    n = model.n
    pu = np.stack([trial.up.copy() for _ in range(nw)])
    pd = np.stack([trial.dn.copy() for _ in range(nw)])
    w = np.ones(nw)                                   # SIGNED weights
    out = []
    for it in range(nsteps):
        for i in range(nw):
            if w[i] == 0:
                continue
            u = model.expK @ pu[i]; d = model.expK @ pd[i]
            O = trial.overlap(u, d)
            if O == 0:
                w[i] = 0.0; continue
            dead = False
            for (a, sa, b, sb, e1, e2, cf) in model.terms:
                R = np.empty(2); cand = []
                for k in (0, 1):
                    u2 = u.copy(); d2 = d.copy()
                    (u2 if sa == 0 else d2)[a, :] *= e1[k]
                    (u2 if sb == 0 else d2)[b, :] *= e2[k]
                    O2 = trial.overlap(u2, d2)
                    R[k] = cf[k] * O2 / O
                    cand.append((u2, d2, O2))
                q = np.maximum(R, 0.0) if constrained else np.abs(R)
                tot = q.sum()
                if tot <= 0:
                    w[i] = 0.0; dead = True; break
                k = rng.choice(2, p=q / tot)
                u, d, O = cand[k]
                if constrained:
                    w[i] *= 0.5 * tot
                else:                                  # importance-sampled signed update
                    w[i] *= 0.5 * np.sign(R[k]) * np.abs(R).sum()
            if dead:
                continue
            pu[i] = model.expK @ u; pd[i] = model.expK @ d
        # rescale magnitudes to avoid over/underflow (cancels in all ratios)
        mag = np.mean(np.abs(w[w != 0])) if np.any(w != 0) else 1.0
        if mag > 0:
            w = w / mag
        if (it + 1) % meas_every == 0:
            sw = se = sa = 0.0
            for i in range(nw):
                if w[i] == 0:
                    continue
                try:
                    el = est._local_energy(pu[i], pd[i])
                except np.linalg.LinAlgError:
                    continue
                if not np.isfinite(el):
                    continue
                se += w[i] * el; sw += w[i]; sa += abs(w[i])
            out.append(((it + 1), se, sw, sa))
    return out


def _worker(args):
    """One independent population (both free + CP), for multiprocessing."""
    P, seed = args
    K = hopping(P["lx"], P["ly"], P["t"], P["t2"])
    model = LatticeModel(K, P["lx"], P["ly"], P["dt"], P["U"], 0.0, 0.0)
    vecs = model.eigvecs
    trial = TrialWF(vecs[:, :P["nup"]].copy(), vecs[:, :P["ndn"]].copy(), mode="fixed")
    est = Estimators(model, trial)
    rf = project(model, trial, est, P["nw"], P["nsteps"], False, np.random.default_rng(seed))
    rc = project(model, trial, est, P["nw"], P["nsteps"], True, np.random.default_rng(10_000 + seed))
    return rf, rc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=5); ap.add_argument("--ndn", type=int, default=5)
    ap.add_argument("--U", type=float, default=6.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--t2", type=float, default=0.4, help="frustrating NNN hopping (sign problem)")
    ap.add_argument("--dt", type=float, default=0.05)
    ap.add_argument("--nsteps", type=int, default=40); ap.add_argument("--nw", type=int, default=40)
    ap.add_argument("--npop", type=int, default=24, help="independent populations (for error bars)")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    K = hopping(a.lx, a.ly, a.t, a.t2)
    model = LatticeModel(K, a.lx, a.ly, a.dt, a.U, 0.0, 0.0)
    vecs = model.eigvecs
    trial = TrialWF(vecs[:, :a.nup].copy(), vecs[:, :a.ndn].copy(), mode="fixed")
    est = Estimators(model, trial)

    # ED reference
    try:
        from hubbard_ed import build
        from scipy.sparse.linalg import eigsh
        # ED uses pure NN; with t2 rebuild a custom ED is overkill -> only when t2=0
        if a.t2 == 0:
            _, H, _, _ = build(a.lx, a.ly, a.nup, a.ndn, a.t, a.U)
            E0 = float(eigsh(H.tocsr(), k=1, which="SA", return_eigenvectors=False)[0])
        else:
            E0 = float("nan")
    except Exception:
        E0 = float("nan")

    print(f"# {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} U={a.U} t2={a.t2}  dt={a.dt}  "
          f"npop={a.npop} nw={a.nw}  ED={E0}")
    print(f"{'tau':>6} {'<sign>_free':>12} {'E_free':>10} {'err_free':>9} {'E_CP':>10} {'err_CP':>8}")
    nm = a.nsteps // 2
    acc = {m: {"fe": [], "fs": [], "cp": []} for m in range(nm)}
    P = dict(lx=a.lx, ly=a.ly, nup=a.nup, ndn=a.ndn, U=a.U, t=a.t, t2=a.t2, dt=a.dt,
             nsteps=a.nsteps, nw=a.nw)
    import multiprocessing as mp
    nproc = min(a.npop, mp.cpu_count())
    with mp.Pool(nproc) as pool:
        for (rf, rc) in pool.map(_worker, [(P, a.seed + p) for p in range(a.npop)]):
            for m, (it, se, sw, sa) in enumerate(rf):
                if sa > 0:
                    acc[m]["fe"].append(se / sw if sw != 0 else np.nan)
                    acc[m]["fs"].append(sw / sa)
            for m, (it, se, sw, sa) in enumerate(rc):
                if sw > 0:
                    acc[m]["cp"].append(se / sw)
    for m in range(nm):
        tau = (m + 1) * 2 * a.dt
        fe = np.array(acc[m]["fe"]); fs = np.array(acc[m]["fs"]); cp = np.array(acc[m]["cp"])
        fe = fe[np.isfinite(fe)]
        ef = fe.mean() if len(fe) else np.nan
        errf = fe.std() / np.sqrt(len(fe)) if len(fe) > 1 else np.nan
        ec = cp.mean() if len(cp) else np.nan
        errc = cp.std() / np.sqrt(len(cp)) if len(cp) > 1 else np.nan
        print(f"{tau:>6.2f} {np.nanmean(fs):>12.4f} {ef:>10.4f} {errf:>9.4f} {ec:>10.4f} {errc:>8.4f}")


if __name__ == "__main__":
    main()
