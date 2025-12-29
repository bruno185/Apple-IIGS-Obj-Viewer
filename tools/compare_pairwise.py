#!/usr/bin/env python3
"""Compare C pair_debug CSV (V4) with Pascal simulation decisions for Devant tests."""
import csv
import sys
import math
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float

import os
import sys
C_CSV = os.environ.get('C_CSV', "Windows\\build_vs\\bin\\Release\\pair_debug_v4_onload.csv")
OBJ_PATH = os.environ.get('OBJ_PATH', "Windows\\3d objects\\q2.obj")
# Allow command-line overrides: python compare_pairwise.py OBJ_PATH [C_CSV]
if len(sys.argv) > 1:
    OBJ_PATH = sys.argv[1]
if len(sys.argv) > 2:
    C_CSV = sys.argv[2]

EPS = 0.01

def arr(r):
    return 0.0 if abs(r) < EPS else r

# load OBJ and faces
verts, faces = parse_obj(OBJ_PATH)
# Use the viewer observer params matching the log (H=30, V=20, D=39.214455)
H = 30; V = 20; W = 0; D = 39.214455
xo, yo, zo, x2d, y2d = transform_and_project(verts, H, V, W, D)
faces_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)

# Dump Pascal per-value trace to mirror C's pair_values_v4.csv for line-by-line diff
pascal_vals = os.environ.get('PASCAL_VALS', "Windows\\build_vs\\bin\\Release\\pair_values_pascal.csv")
cos_h = math.cos(math.radians(H)); sin_h = math.sin(math.radians(H))
cos_v = math.cos(math.radians(V)); sin_v = math.sin(math.radians(V))
xobs = D * cos_h * cos_v
yobs = D * sin_h * cos_v
zobs = D * sin_v
with open(pascal_vals, 'w', newline='') as pf:
    pf.write('pass,i,j,f1,f2,stage,idx,xo,yo,zo,wx,wy,wz,tempo_obs,tempo_world,valobs,Ap,Bp,Cp,Dp,Aq,Bq,Cq,Dq\n')
    for i in range(len(faces_info)-1):
        for j in range(i+1, len(faces_info)):
            f1=i; f2=j
            Ap = arr(to_float(faces_info[f1]['a'])); Bp = arr(to_float(faces_info[f1]['b'])); Cp = arr(to_float(faces_info[f1]['c'])); Dp = arr(to_float(faces_info[f1]['d']))
            Aq = arr(to_float(faces_info[f2]['a'])); Bq = arr(to_float(faces_info[f2]['b'])); Cq = arr(to_float(faces_info[f2]['c'])); Dq = arr(to_float(faces_info[f2]['d']))
            pf.write(f"0,{i},{j},{f1},{f2},planes,{Ap:.6f},{Bp:.6f},{Cp:.6f},{Dp:.6f},{Aq:.6f},{Bq:.6f},{Cq:.6f},{Dq:.6f}\n")
            for vid in faces_info[f2]['idxs']:
                xo_v = to_float(xo[vid]); yo_v = to_float(yo[vid]); zo_v = to_float(zo[vid])
                wx,wy,wz = verts[vid]
                tempo_obs = arr(Ap * xo_v + Bp * yo_v + Cp * zo_v + Dp)
                tempo_world = arr(Ap * wx + Bp * wy + Cp * wz + Dp)
                pf.write(f"0,{i},{j},{f1},{f2},vertex_pq,{vid},{xo_v:.6f},{yo_v:.6f},{zo_v:.6f},{wx:.6f},{wy:.6f},{wz:.6f},{tempo_obs:.6f},{tempo_world:.6f}\n")
            for vid in faces_info[f1]['idxs']:
                xo_v = to_float(xo[vid]); yo_v = to_float(yo[vid]); zo_v = to_float(zo[vid])
                wx,wy,wz = verts[vid]
                tempo_obs = arr(Aq * xo_v + Bq * yo_v + Cq * zo_v + Dq)
                tempo_world = arr(Aq * wx + Bq * wy + Cq * wz + Dq)
                pf.write(f"0,{i},{j},{f1},{f2},vertex_qp,{vid},{xo_v:.6f},{yo_v:.6f},{zo_v:.6f},{wx:.6f},{wy:.6f},{wz:.6f},{tempo_obs:.6f},{tempo_world:.6f}\n")
            valobs_p = arr(Ap * xobs + Bp * yobs + Cp * zobs + Dp)
            valobs_q = arr(Aq * xobs + Bq * yobs + Cq * zobs + Dq)
            pf.write(f"0,{i},{j},{f1},{f2},valobs,{valobs_p:.6f},{valobs_q:.6f}\n")

# read C CSV
rows = []
with open(C_CSV, newline='') as f:
    rdr = csv.DictReader(x for x in f if not x.startswith('#'))
    for r in rdr:
        rows.append(r)

mismatches = []
for r in rows:
    f1 = int(r['f1']); f2 = int(r['f2'])
    # parse C-reported values
    c_dev_pq = int(r['devant_pq']); c_dev_qp = int(r['devant_qp'])
    # use C's recorded observer coordinates if present
    try:
        xobs = arr(float(r['xobs']))
        yobs = arr(float(r['yobs']))
        zobs = arr(float(r['zobs']))
    except Exception:
        # fallback: compute from faces_info metadata (not ideal)
        xobs = yobs = zobs = 0.0

    # Pascal-plane coeffs (arrondi applied to coefficients per Pascal)
    Ap = arr(to_float(faces_info[f1]['a'])); Bp = arr(to_float(faces_info[f1]['b'])); Cp = arr(to_float(faces_info[f1]['c'])); Dp = arr(to_float(faces_info[f1]['d']))
    Aq = arr(to_float(faces_info[f2]['a'])); Bq = arr(to_float(faces_info[f2]['b'])); Cq = arr(to_float(faces_info[f2]['c'])); Dq = arr(to_float(faces_info[f2]['d']))

    # compute dev_pq (Pascal logic)
    pos=neg=0
    for vid in faces_info[f2]['idxs']:
        tempo = arr(Ap * verts[vid][0] + Bp * verts[vid][1] + Cp * verts[vid][2] + Dp)
        if tempo > 0: pos += 1
        elif tempo < 0: neg += 1
    if pos > 0 and neg > 0:
        p_dev = 0
    else:
        cotepoint = 1 if pos > 0 else -1
        position = 1 if arr(Ap * xobs + Bp * yobs + Cp * zobs + Dp) > 0 else -1
        p_dev = 1 if position == cotepoint else -1

    # compute dev_qp
    pos2=neg2=0
    for vid in faces_info[f1]['idxs']:
        tempo = arr(Aq * verts[vid][0] + Bq * verts[vid][1] + Cq * verts[vid][2] + Dq)
        if tempo > 0: pos2 += 1
        elif tempo < 0: neg2 += 1
    if pos2 > 0 and neg2 > 0:
        q_dev = 0
    else:
        cotepoint2 = 1 if pos2 > 0 else -1
        position2 = 1 if arr(Aq * xobs + Bq * yobs + Cq * zobs + Dq) > 0 else -1
        q_dev = 1 if position2 == cotepoint2 else -1

    if p_dev != c_dev_pq or q_dev != c_dev_qp:
        mismatches.append({
            'pair': (f1, f2),
            'c_dev_pq': c_dev_pq, 'pas_dev_pq': p_dev,
            'c_dev_qp': c_dev_qp, 'pas_dev_qp': q_dev,
            'Ap': Ap, 'Bp': Bp, 'Cp': Cp, 'Dp': Dp,
            'Aq': Aq, 'Bq': Bq, 'Cq': Cq, 'Dq': Dq,
            'xobs': xobs, 'yobs': yobs, 'zobs': zobs,
            'pos_pq': pos, 'neg_pq': neg, 'pos_qp': pos2, 'neg_qp': neg2
        })

# report
if not mismatches:
    print('Aucune divergence trouvée entre C (pair_debug) et Pascal pour les paires listées.')
    print('Compare complete: no mismatches')
    sys.exit(0)

print('Divergences trouvées pour les paires suivantes:')
for m in mismatches:
    a,b = m['pair']
    print(f"Pair {a},{b}: C dev_pq={m['c_dev_pq']} vs Pascal dev_pq={m['pas_dev_pq']}; C dev_qp={m['c_dev_qp']} vs Pascal dev_qp={m['pas_dev_qp']}")
    print(f"  Ap,Bp,Cp,Dp = {m['Ap']:.6f},{m['Bp']:.6f},{m['Cp']:.6f},{m['Dp']:.6f}")
    print(f"  Aq,Bq,Cq,Dq = {m['Aq']:.6f},{m['Bq']:.6f},{m['Cq']:.6f},{m['Dq']:.6f}")
    print(f"  xobs,yobs,zobs = {m['xobs']:.6f},{m['yobs']:.6f},{m['zobs']:.6f}")
    print(f"  pos/neg counts: p->q pos/neg = {m['pos_pq']}/{m['neg_pq']}, q->p pos/neg = {m['pos_qp']}/{m['neg_qp']}")
    print('')

sys.exit(0)
