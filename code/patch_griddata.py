"""One-off: switch the parameter-map figures from griddata to reggrid.

griddata triangulates a rectangular lattice, and the triangulation is not stable
under adding rows (see gridinterp.py). Every (U, delta) and (n, delta) map in
these scripts was affected, drifting up to 64% of the panel range when a row was
dropped. This rewrites those call sites.

The k-space maps (kx, ky) are NOT touched here: several tile and snap their points
first, so whether they form a clean rectangle needs checking case by case.

Run once from code/. Idempotent: already-patched files are reported as skipped.
"""
import os, re, sys

IMPORT_REPO = "from gridinterp import reggrid\n"

# repo scripts: explicit (old, new) pairs, matched verbatim
REPO = {
 "fig7_pairing_phase_diagram.py": [(
  "G = griddata((df.n, df.delta), df['diff'], (NI, DI), method='linear')",
  "G = reggrid(df, 'n', 'delta', 'diff', NI, DI)")],
 "fig_4ch_3L.py": [(
  "Z = griddata((s.delta.values, s.U.values), s[key].values, (gx, gy), method='linear')",
  "Z = reggrid(s, 'delta', 'U', key, gx, gy)")],
 "fig_4ch_3L_shared.py": [(
  "Z = griddata((s.delta.values, s.U.values), s[key].values, (gx, gy), method='linear')",
  "Z = reggrid(s, 'delta', 'U', key, gx, gy)")],
 "fig_4ch_shared_L12.py": [(
  "Z = griddata((d.delta.values, d.U.values), d[key].values, (gx, gy), method='linear')",
  "Z = reggrid(d, 'delta', 'U', key, gx, gy)")],
 "fig_dtot_L8_L10_L12.py": [(
  "Z = griddata((s.delta.values, s.U.values), s.dtot_N.values, (gx, gy), method='linear')",
  "Z = reggrid(s, 'delta', 'U', 'dtot_N', gx, gy)")],
 "fig_dtot_bg_pair_contours.py": [(
  "BG = griddata((d.delta.values, d.U.values), d.dtot_N.values, (gx, gy), method='linear')",
  "BG = reggrid(d, 'delta', 'U', 'dtot_N', gx, gy)"), (
  "F = griddata((d.delta.values, d.U.values), d[key].values, (gx, gy), method='linear')",
  "F = reggrid(d, 'delta', 'U', key, gx, gy)")],
 "fig_dtot_colour_ch_contours_3L.py": [(
  "BG = griddata((s.delta.values, s.U.values), s.dtot_N.values, (gx, gy), method='linear')",
  "BG = reggrid(s, 'delta', 'U', 'dtot_N', gx, gy)"), (
  "F = griddata((s.delta.values, s.U.values), s[key].values, (gx, gy), method='linear')",
  "F = reggrid(s, 'delta', 'U', key, gx, gy)")],
 "fig_dtot_contours_fillings.py": [(
  "BG = griddata((t.delta.values, t.U.values), t.dtot_N.values, (gx, gy), method='linear')",
  "BG = reggrid(t, 'delta', 'U', 'dtot_N', gx, gy)"), (
  "F = griddata((t.delta.values, t.U.values), t[key].values, (gx, gy), method='linear')",
  "F = reggrid(t, 'delta', 'U', key, gx, gy)")],
 "fig_dtot_map_L12.py": [(
  "Z = griddata((d.delta.values, d.U.values), d.dtot_N.values, (gx, gy), method='linear')",
  "Z = reggrid(d, 'delta', 'U', 'dtot_N', gx, gy)")],
 "fig_dtot_pair_L12.py": [(
  "gi = griddata((d.delta.values, d.U.values), d[key].values, (gx, gy),\n"
  "                  method='linear')",
  "gi = reggrid(d, 'delta', 'U', key, gx, gy)")],
 "fig_dtot_pairing_same_axes.py": [(
  'gi = griddata((T[xk].values, T[yk].values), T[fld].values, (gx, gy),\n'
  '                      method="linear")',
  "gi = reggrid(T, xk, yk, fld, gx, gy)")],
 "fig_dtot_prl.py": [(
  'gi = griddata((h.delta.values, h.U.values), h.dtot_N.values, (gx, gy), method="linear")',
  "gi = reggrid(h, 'delta', 'U', 'dtot_N', gx, gy)"), (
  'gi2 = griddata((g10.n.values, g10.delta.values), g10.dtot_N.values,\n'
  '               (gx2, gy2), method="linear")',
  "gi2 = reggrid(g10, 'n', 'delta', 'dtot_N', gx2, gy2)")],
 "fig_dnk_fortran.py": [(
  'gi = griddata((h.delta.values, h.U.values), h.delta_tot.values,\n'
  '              tuple(np.meshgrid(np.linspace(0, 0.4, 300), np.linspace(0, 5, 300))),\n'
  '              method="linear")',
  "_gx, _gy = np.meshgrid(np.linspace(0, 0.4, 300), np.linspace(0, 5, 300))\n"
  "gi = reggrid(h, 'delta', 'U', 'delta_tot', _gx, _gy)"), (
  'gi2 = griddata((g10.n.values, g10.delta.values), g10.delta_tot.values,\n'
  '               tuple(np.meshgrid(np.linspace(0.5, 0.98, 300),\n'
  '                                 np.linspace(0, 0.4, 300))), method="linear")',
  "_gx2, _gy2 = np.meshgrid(np.linspace(0.5, 0.98, 300), np.linspace(0, 0.4, 300))\n"
  "gi2 = reggrid(g10, 'n', 'delta', 'delta_tot', _gx2, _gy2)")],
 "fig_eq_vs_uneq_phase.py": [(
  '        griddata(pts, T[c].values, (NG, DG), method="linear",\n'
  '                 fill_value=np.nan), 4.0) for c in KEYS])',
  '        reggrid(T, "n", "delta", c, NG, DG,\n'
  '                fill_value=np.nan), 4.0) for c in KEYS])')],
 "make_fig_dchannels.py": [(
  'Z = griddata((s.delta.values, s.U.values), (s.dxy - s.d).values, (gx, gy),\n'
  '                 method="linear")',
  'Z = reggrid(s.assign(_diff=s.dxy - s.d), "delta", "U", "_diff", gx, gy)')],
}

# my scripts that keep an inline grid(): swap the cubic/linear pair for reggrid
INLINE_OLD = '''    P = (s.delta.values, s.U.values)
    a = griddata(P, v, (gx, gy), method="cubic")
    b = griddata(P, v, (gx, gy), method="linear")
    a[np.isnan(a)] = b[np.isnan(a)]
    return a'''
INLINE_NEW = '''    return reggrid(s, "delta", "U", v, gx, gy)'''

done, skipped, missed = [], [], []
for fn, pairs in REPO.items():
    if not os.path.exists(fn): missed.append((fn, "file not found")); continue
    src = open(fn).read()
    if "reggrid" in src: skipped.append(fn); continue
    ok = True
    for old, new in pairs:
        if old not in src: missed.append((fn, old.split("\n")[0][:60])); ok = False; break
        src = src.replace(old, new)
    if not ok: continue
    # add the import after the last existing import line
    lines = src.split("\n")
    last = max(i for i, l in enumerate(lines)
               if l.startswith("import ") or l.startswith("from "))
    lines.insert(last + 1, IMPORT_REPO.rstrip())
    open(fn, "w").write("\n".join(lines))
    done.append(fn)

print(f"patched  ({len(done)}):")
for f in done: print("   ", f)
if skipped: print(f"skipped, already using reggrid ({len(skipped)}):");  [print("   ", f) for f in skipped]
if missed:
    print(f"NOT patched ({len(missed)}) -- pattern not found, handle by hand:")
    for f, why in missed: print(f"    {f}: {why}")
