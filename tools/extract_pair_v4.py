import csv
PVF = r"Windows\build_vs\bin\Release\pair_values_v4.csv"
TARGETS = [(1,0),(0,1)]
rows = []
with open(PVF,newline='') as f:
    rdr = csv.reader(f)
    header = next(rdr)
    h = {n:i for i,n in enumerate(header)}
    for r in rdr:
        try:
            f1 = int(r[h['f1']]); f2 = int(r[h['f2']])
        except Exception:
            continue
        if (f1,f2) in TARGETS:
            rows.append((f1,f2,r))

print('Found',len(rows),'rows for pairs',TARGETS)
for f1,f2,r in rows:
    stage = r[h.get('stage','stage')]
    idx = r[h.get('idx','idx')] if 'idx' in h else ''
    xo = r[h.get('xo','xo')] if 'xo' in h else ''
    yo = r[h.get('yo','yo')] if 'yo' in h else ''
    zo = r[h.get('zo','zo')] if 'zo' in h else ''
    print(f'pair=({f1},{f2})', stage, idx, 'xo=',xo,'yo=',yo,'zo=',zo)
rows = []
with open(PVF,newline='') as f:
    rdr = csv.reader(f)
    header = next(rdr)
    h = {n:i for i,n in enumerate(header)}
    for r in rdr:
        try:
            f1 = int(r[h['f1']]); f2 = int(r[h['f2']])
        except Exception:
            continue
        if (f1,f2) == TARGET:
            rows.append(r)

print('Found',len(rows),'rows for pair',TARGET)
for r in rows:
    stage = r[h.get('stage','stage')]
    idx = r[h.get('idx','idx')] if 'idx' in h else ''
    xo = r[h.get('xo','xo')] if 'xo' in h else ''
    yo = r[h.get('yo','yo')] if 'yo' in h else ''
    zo = r[h.get('zo','zo')] if 'zo' in h else ''
    print(stage, idx, 'xo=',xo,'yo=',yo,'zo=',zo)
