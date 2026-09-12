"""Second pass: the d-channel scripts written this session.

Repo scripts import gridinterp. The *_251.py scripts are meant to be pasted into a
notebook on node 251, so they cannot import anything local and carry the fixed
interpolator inline instead.
"""
import os

INLINE_DOC = '''    """Interpolate a column onto the display mesh.

    NOT scipy.griddata. That treats the points as scattered and triangulates, but
    this data is a rectangular lattice where every cell can split along either
    diagonal, both valid Delaunay. Qhull's tie-break is unstable: dropping one row
    of the scan moved contours elsewhere in the panel by up to 64% of its range.
    The figure was then not a function of the data alone.

    RegularGridInterpolator with method="linear" is bilinear per cell: unique,
    local, reproducible, bounded by the measurements."""
    p = s.pivot_table(index="U", columns="delta", values=col)
    assert not p.isna().any().any(), "grid has holes, pivot_table left NaN"
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method="linear", bounds_error=False, fill_value=None)
    return fn(np.stack([gy.ravel(), gx.ravel()], -1)).reshape(gx.shape)'''

# (file, is_251, old grid() body, call-site rewrites)
JOBS = [
 ("make_fig_dchannels_map.py", False,
  '''    P = (s.delta.values, s.U.values)
    a = griddata(P, v, (gx, gy), method="cubic")
    b = griddata(P, v, (gx, gy), method="linear")
    a[np.isnan(a)] = b[np.isnan(a)]
    return a''',
  '''    return reggrid(s, "delta", "U", col, gx, gy)'''),
 ("make_fig_dchannels_eqtime.py", False,
  '''    P = (s.delta.values, s.U.values)
    a = griddata(P, v, (gx, gy), method="cubic")
    b = griddata(P, v, (gx, gy), method="linear")
    a[np.isnan(a)] = b[np.isnan(a)]
    return a''',
  '''    return reggrid(s, "delta", "U", col, gx, gy)'''),
 ("make_fig_dchannels_eqtime_vs_vertex.py", False,
  '''    P = (s[col_x].values, s[col_y].values)
    a = griddata(P, v_, (gx, gy), method="cubic")
    b = griddata(P, v_, (gx, gy), method="linear")
    a[np.isnan(a)] = b[np.isnan(a)]
    return a''',
  '''    return reggrid(s, col_x, col_y, col, gx, gy)'''),
 ("fig_dchannels_251.py", True,
  '''    P = (s.delta.values, s.U.values)
    a = griddata(P, v, (gx,gy), "cubic"); b = griddata(P, v, (gx,gy), "linear")
    a[np.isnan(a)] = b[np.isnan(a)];  return a''', None),
 ("fig_dchannels_eqtime_251.py", True,
  '''    P = (s.delta.values, s.U.values)
    a = griddata(P, v, (gx,gy), "cubic"); b = griddata(P, v, (gx,gy), "linear")
    a[np.isnan(a)] = b[np.isnan(a)];  return a''', None),
 ("fig_dchannels_eqtime_vs_vertex_251.py", True,
  '''    P = (s.delta.values, s.U.values)
    a = griddata(P, val, (gx,gy), "cubic"); b = griddata(P, val, (gx,gy), "linear")
    a[np.isnan(a)] = b[np.isnan(a)];  return a''', None),
]

CALLS = [("grid(s, s.dtot_N.values)",       'grid(s, "dtot_N")'),
         ("grid(bg, bg.dtot_N.values)",     'grid(bg, "dtot_N")'),
         ("grid(s, s.d.values)",            'grid(s, "d")'),
         ("grid(s, s.dxy.values)",          'grid(s, "dxy")'),
         ("grid(s, (s.dxy - s.d).values)",  'grid(s.assign(_diff=s.dxy - s.d), "_diff")'),
         ("grid(s, (s.dxy-s.d).values)",    'grid(s.assign(_diff=s.dxy-s.d), "_diff")')]

done, skipped, missed = [], [], []
for fn, is251, old, new in JOBS:
    if not os.path.exists(fn): missed.append((fn, "not found")); continue
    src = open(fn).read()
    if "RegularGridInterpolator" in src or "reggrid" in src:
        skipped.append(fn); continue
    if old not in src: missed.append((fn, "grid() body not matched")); continue
    src = src.replace(old, INLINE_DOC if is251 else new)
    if is251:
        src = src.replace("from scipy.interpolate import griddata",
                          "from scipy.interpolate import RegularGridInterpolator")
        # inline version takes a column name; rename the parameter
        for a, b in (("def grid(s, v):", "def grid(s, col):"),
                     ("def grid(s, val):", "def grid(s, col):")):
            src = src.replace(a, b)
    else:
        lines = src.split("\n")
        last = max(i for i, l in enumerate(lines)
                   if l.startswith("import ") or l.startswith("from "))
        lines.insert(last + 1, "from gridinterp import reggrid")
        src = "\n".join(lines)
        for a, b in (("def grid(s, v):", "def grid(s, col):"),
                     ('def grid(s, v_, col_x="delta", col_y="U"):',
                      'def grid(s, col, col_x="delta", col_y="U"):')):
            src = src.replace(a, b)
    for a, b in CALLS: src = src.replace(a, b)
    open(fn, "w").write(src)
    done.append(fn)

print(f"patched ({len(done)}):");  [print("   ", f) for f in done]
if skipped: print(f"already fixed ({len(skipped)}):"); [print("   ", f) for f in skipped]
if missed:  print(f"NOT patched ({len(missed)}):");    [print(f"    {f}: {w}") for f, w in missed]
