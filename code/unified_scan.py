"""Minimal deploy shim for `unified_scan`.

`cpqmc.CPMC.run_bp_chi_spin` does `from unified_scan import _shift_index`, and that is the
ONLY thing it needs from this module (verified: cpqmc.py references unified_scan exactly once,
at the import, and uses only _shift_index). The real pyqmc/unified_scan.py additionally imports
agp_dwave -> agp -> agp_bp_vertex -> pair_rspace, none of which are used by run_bp_chi_spin;
shipping it to a compute node therefore fails with

    ModuleNotFoundError: No module named 'agp_dwave'

This shim provides the one helper verbatim so the magnetic driver stays self-contained
(cpqmc.py + checkerboard.py + unified_scan.py) with no further dependencies.

Function copied unchanged from qmc-platform-master/pyqmc/unified_scan.py:44.
"""
import numpy as np


def _shift_index(lx, ly):
    """shift[dx,dy] = array n(m) mapping site m -> site (x+dx, y+dy) (PBC). Site i=x*ly+y."""
    n = lx * ly; xm = np.arange(n) // ly; ym = np.arange(n) % ly
    S = np.empty((lx, ly, n), dtype=int)
    for dx in range(lx):
        for dy in range(ly):
            S[dx, dy] = ((xm + dx) % lx) * ly + (ym + dy) % ly
    return S
