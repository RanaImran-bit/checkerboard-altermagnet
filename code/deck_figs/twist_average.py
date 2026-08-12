"""1:2:1 twist average of the pairing vertex. Run once the (-1,-1) leg lands."""
import numpy as np, pandas as pd, glob, os, sys

PER = '/home/phd25imran/analysis/data/pairing_master.csv'   # periodic, weight 1
T21 = sys.argv[1] if len(sys.argv) > 1 else 'tw_m1p1.csv'   # (-1,+1), weight 2
T11 = sys.argv[2] if len(sys.argv) > 2 else 'tw_m1m1.csv'   # (-1,-1), weight 1
L, U, NUP = 12, 4.0, 72
CH = ['chi_son', 'chi_sext', 'chi_d', 'chi_dxy']
N = L * L

p = pd.read_csv(PER)
p = p[(p.L == L) & np.isclose(p.U, U) & (p.nup == NUP)]
per = p.groupby('delta')[CH].agg(['mean', 'sem']) / N       # these columns ARE the vertex

def leg(path):
    d = pd.concat([pd.read_csv(f) for f in glob.glob(path)], ignore_index=True) \
        if '*' in path else pd.read_csv(path)
    d = d[(d.L == L) & np.isclose(d.U, U) & (d.nup == NUP)]
    d = d.rename(columns={c + '_vertex': c for c in CH})
    return d.groupby('delta')[CH].agg(['mean', 'sem']) / N

a, b = leg(T21), leg(T11)
ds = sorted(set(per.index) & set(a.index) & set(b.index))
print(f'twist-averaged pairing vertex  chi/N   L={L}, U={U:g}, n=1.0')
print(f'weights 1 : 2 : 1  for  (+1,+1) : (-1,+1) : (-1,-1)\n')
print(f'{"delta":>6}{"periodic":>11}{"(-1,+1)":>11}{"(-1,-1)":>11}{"AVERAGE":>12}{"+/-":>9}   swap?')
for x in ds:
    v = [per.loc[x, ('chi_dxy', 'mean')], a.loc[x, ('chi_dxy', 'mean')], b.loc[x, ('chi_dxy', 'mean')]]
    e = [per.loc[x, ('chi_dxy', 'sem')], a.loc[x, ('chi_dxy', 'sem')], b.loc[x, ('chi_dxy', 'sem')]]
    avg = (v[0] + 2 * v[1] + v[2]) / 4
    err = np.sqrt(e[0]**2 + (2 * e[1])**2 + e[2]**2) / 4
    print(f'{x:6.1f}{v[0]:+11.4f}{v[1]:+11.4f}{v[2]:+11.4f}{avg:+12.4f}{err:9.4f}'
          f'   {"ATTRACTIVE" if avg > 2*err else ("repulsive" if avg < -2*err else "unresolved")}')

print('\nall four channels, twist-averaged:')
print(f'{"delta":>6}' + ''.join(f'{c.replace("chi_",""):>12}' for c in CH))
for x in ds:
    row = ''
    for c in CH:
        avg = (per.loc[x, (c, 'mean')] + 2 * a.loc[x, (c, 'mean')] + b.loc[x, (c, 'mean')]) / 4
        row += f'{avg:+12.4f}'
    print(f'{x:6.1f}{row}')
