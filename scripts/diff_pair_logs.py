import csv
from math import isclose

PY = 'PAIRPY.CSV'
C = 'PAIRDIAG.CSV'

def read_csv(path):
    rows = []
    with open(path, 'r', newline='') as f:
        for line in f:
            rows.append(line.strip())
    return rows

p = read_csv(PY)
c = read_csv(C)

print('PY lines:', len(p), 'C lines:', len(c))

# Compare face headers
py_faces = [l for l in p if l.startswith('face,')]
c_faces = [l for l in c if l.startswith('face,')]
print('\nFaces in PY:', len(py_faces), 'Faces in C:', len(c_faces))

# Parse face lines into numeric fields
import re
float_re = re.compile(r"-?\d+\.\d+")
for i,(pl,cl) in enumerate(zip(py_faces, c_faces)):
    print('\nFace', i)
    print('PY:', pl)
    print('C :', cl)
    pnums = float_re.findall(pl)
    cnums = float_re.findall(cl)
    print('PY numeric count:', len(pnums), 'C numeric count:', len(cnums))
    # show first few numeric fields
    print(' PY nums:', pnums[:8])
    print(' C nums:', cnums[:8])
    # Show relative diffs for first 4 comparable fields if available
    for k in range(min(4, len(pnums), len(cnums))):
        pyv = float(pnums[k]); cv = float(cnums[k]);
        if pyv == 0:
            diff = float('inf') if cv != 0 else 0.0
        else:
            diff = (cv - pyv) / abs(pyv)
        print(f'  Field {k}: PY={pyv:.6g} C={cv:.6g} rel_diff={diff:.6g}')

# Find pair_should_swap lines
py_pair = [l for l in p if l.startswith('pair_should_swap') or l.startswith('pair_should_swap_')]
c_pair = [l for l in c if l.startswith('pair_should_swap')]
print('\nPY pair lines:')
for l in py_pair: print(' ', l)
print('\nC pair lines:')
for l in c_pair: print(' ', l)

# Show any TRI lines in C
tri_lines = [l for l in c if l.startswith('TRI')]
print('\nTRI lines in C:', len(tri_lines))
for l in tri_lines:
    print(' ', l)

# Heuristic: find any corrupted fields (non-numeric tokens) in eval lines
print('\nChecking eval lines for non-numeric tokens...')
bad = False
for l in c:
    if l.startswith('eval,'):
        parts = l.split(',')
        for idx,part in enumerate(parts[3:7], start=3):
            try:
                float(part)
            except ValueError:
                print(' Bad numeric field in C eval line:', l)
                bad = True
                break
if not bad:
    print(' No corrupted eval fields found in C')

print('\nDone')
