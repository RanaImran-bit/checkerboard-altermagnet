# qmc_algo — new QMC method development

Clean workspace for developing a **new/improved QMC algorithm**, separate from the
validated platform code so experiments can move fast without touching the gated
production tree.

## Layout
```
src/        the new method (modules)
tests/      ED-gated unit + regression tests (the contract: match ED where exact)
notes/      design notes, derivations, the running lab notebook
results/    run outputs, scans, figures (gitignore heavy artifacts)
```

## Reuse, don't reinvent
The repo already has validated infrastructure — lean on it as ground truth:
- `pyqmc/cpqmc.py` — validated CP-AFQMC (LatticeModel / TrialWF / WalkerEnsemble /
  Propagator / Estimators / CPMC). Reference implementation + reusable pieces.
- `pyqmc/regression.py` — 5 fixed-seed bit-for-bit gate. Any shared code we touch
  must keep this green.
- ED references: `ed/hubbard_ed.py`, `ed/altermagnet_ed.py`, and the `pyqmc/validate_*.py`
  Lehmann/Slater-Condon checks. These are the exact gates for the new method.
- `docs/theory/cpqmc_theory.tex` — the derivation the current code maps to; extend
  it with the new method's math.

## Working rule
Every new estimator/propagator is validated on a small cluster against ED **first**
(exact at U=0, within known bias interacting) before any scaling run. Dynamic
estimators must be checked on an *interacting* (asymmetric) case — the U=0 pass is
necessary but not sufficient (see HANDOFF, WS2 lesson).

## Method
_TBD — see notes/00_design.md once the specific idea is pinned._
