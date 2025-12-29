#!/usr/bin/env python3
import csv
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float

C_CSV = "Windows\\build_vs\\bin\\Release\\pair_debug_v4_onload.csv"
OBJ_PATH = "Windows\\3d objects\\q2.obj"
H=30; V=20; W=0; D=39.214455

verts, faces = parse_obj(OBJ_PATH)
xo, yo, zo, x2d, y2d = transform_and_project(verts, H, V, W, D)
faces_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)

# load C rows
rows = []
with open(C_CSV, newline='') as f:
    # skip comments
    for line in f:
        if line.startswith('#'):
            continue
        rows.append(line.strip())
# parse header
hdr = rows[0].split(',')
entries = [dict(zip(hdr, r.split(','))) for r in rows[1:]]

pairs_to_check = [(0,1),(0,2),(1,2)]
for e in entries:
    f1 = int(e['f1']); f2 = int(e['f2'])
    if (f1,f2) not in pairs_to_check:
        continue
    print(f'PAIR {f1} vs {f2}')
    Ap=float(e['Ap']); Bp=float(e['Bp']); Cp=float(e['Cp']); Dp=float(e['Dp'])
    Aq=float(e['Aq']); Bq=float(e['Bq']); Cq=float(e['Cq']); Dq=float(e['Dq'])
    xobs=float(e['xobs']); yobs=float(e['yobs']); zobs=float(e['zobs'])
    print(' C planes: P',Ap,Bp,Cp,Dp,' Q',Aq,Bq,Cq,Dq)
    # compute per-vertex tv using Pascal xo/yo/zo (to_float converts fixed to float)
    print('\n Per-vertex tv for P plane applied to Q vertices (Pascal coords):')
    for idx in faces_info[f2]['idxs']:
        tv = Ap * to_float(xo[idx]) + Bp * to_float(yo[idx]) + Cp * to_float(zo[idx]) + Dp
        arr = 0.0 if abs(tv) < 0.01 else tv
        print(f' v{idx}: tv={tv:.6f} arr={arr:.6f} sign={"pos" if arr>0 else ("neg" if arr<0 else "zero")}')
    print('\n Per-vertex tv for Q plane applied to P vertices (Pascal coords):')
    for idx in faces_info[f1]['idxs']:
        tv = Aq * to_float(xo[idx]) + Bq * to_float(yo[idx]) + Cq * to_float(zo[idx]) + Dq
        arr = 0.0 if abs(tv) < 0.01 else tv
        print(f' v{idx}: tv={tv:.6f} arr={arr:.6f} sign={"pos" if arr>0 else ("neg" if arr<0 else "zero")}')
    print('\n Pascal plane coeffs (float from fixed):')
    Ap_p = to_float(faces_info[f1]['a']); Bp_p = to_float(faces_info[f1]['b']); Cp_p = to_float(faces_info[f1]['c']); Dp_p = to_float(faces_info[f1]['d'])
    Aq_p = to_float(faces_info[f2]['a']); Bq_p = to_float(faces_info[f2]['b']); Cq_p = to_float(faces_info[f2]['c']); Dq_p = to_float(faces_info[f2]['d'])
    print(' Pascal planes: P',Ap_p,Bp_p,Cp_p,Dp_p,' Q',Aq_p,Bq_p,Cq_p,Dq_p)
    print('\n Per-vertex tv for Pascal P plane on Q vertices:')
    for idx in faces_info[f2]['idxs']:
        tv = Ap_p * to_float(xo[idx]) + Bp_p * to_float(yo[idx]) + Cp_p * to_float(zo[idx]) + Dp_p
        arr = 0.0 if abs(tv) < 0.01 else tv
        print(f' v{idx}: tv={tv:.6f} arr={arr:.6f} sign={"pos" if arr>0 else ("neg" if arr<0 else "zero")}')
    print('\n---\n')
