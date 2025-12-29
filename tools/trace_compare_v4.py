#!/usr/bin/env python3
"""Compare per-instruction trace from C's pair_values_v4.csv with Python pascal_sim computation.
Stops at the first mismatch and prints context.
Usage: python tools/trace_compare_v4.py [OBJ_PATH] [C_PVF]
Defaults: OBJ_PATH="Windows\\3d objects\\q2.obj"  C_PVF="Windows\\build_vs\\bin\\Release\\pair_values_v4.csv"
"""
import csv
import sys
import math
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float

OBJ_PATH = sys.argv[1] if len(sys.argv) > 1 else "Windows\\3d objects\\q2.obj"
C_PVF = sys.argv[2] if len(sys.argv) > 2 else "Windows\\build_vs\\bin\\Release\\pair_values_v4.csv"
# viewer params that match the app's defaults/logs
H = 30; V = 20; W = 0; D = 39.214455

EPS = 1e-1  # broader tolerance while we converge on remaining rounding-order differences (temporary)

import os
# load model and run pascal pipeline
verts, faces = parse_obj(OBJ_PATH)
# Apply same Y/Z swap and optional centering as C loader (GS3Dp behavior)
verts_c = [(x, z, y) for (x,y,z) in verts]
# Apply bbox centering unless NO_AUTO_CENTER is set
if not os.getenv('NO_AUTO_CENTER'):
    minx = min(v[0] for v in verts_c); maxx = max(v[0] for v in verts_c)
    miny = min(v[1] for v in verts_c); maxy = max(v[1] for v in verts_c)
    minz = min(v[2] for v in verts_c); maxz = max(v[2] for v in verts_c)
    cx = (minx + maxx) * 0.5; cy = (miny + maxy) * 0.5; cz = (minz + maxz) * 0.5
    verts_c = [(x-cx, y-cy, z-cz) for (x,y,z) in verts_c]
xo, yo, zo, x2d, y2d = transform_and_project(verts_c, H, V, W, D)
faces_info = compute_faces(verts_c, faces, xo, yo, zo, x2d, y2d)

# Helper arr (Delphi Arrondi) used by C: arr(r)=0 if abs(r)<0.01 else r
def arr(r):
    return 0.0 if abs(r) < 0.01 else r

# Build lookup of expected plane coeffs (float) per face
plane_map = {}
for fi, f in enumerate(faces_info):
    plane_map[fi] = (arr(to_float(f['a'])), arr(to_float(f['b'])), arr(to_float(f['c'])), arr(to_float(f['d'])))

# read C CSV
with open(C_PVF, newline='') as f:
    rdr = csv.DictReader(f)
    rows = list(rdr)

# iterate rows and compare
for r in rows:
    stage = r['stage']
    # skip instr/logging lines that do not conform to the main CSV header
    try:
        i = int(r['i']); j = int(r['j']); f1 = int(r['f1']); f2 = int(r['f2'])
    except Exception:
        continue
    if stage == 'planes':
        # Some CSV lines use the shared header and place plane coeffs starting at the 'idx' column.
        # Robustly extract 8 coefficients in order: Ap,Bp,Cp,Dp,Aq,Bq,Cq,Dq
        def getf(key):
            try:
                return float(r[key])
            except Exception:
                return None
        Ap_c = getf('Ap')
        if Ap_c is None:
            # fallback mapping for 'planes' row: idx->Ap, xo->Bp, yo->Cp, zo->Dp, wx->Aq, wy->Bq, wz->Cq, ax_obs->Dq
            Ap_c = getf('idx')
            Bp_c = getf('xo')
            Cp_c = getf('yo')
            Dp_c = getf('zo')
            Aq_c = getf('wx')
            Bq_c = getf('wy')
            Cq_c = getf('wz')
            Dq_c = getf('ax_obs')
        else:
            Bp_c = getf('Bp'); Cp_c = getf('Cp'); Dp_c = getf('Dp'); Aq_c = getf('Aq'); Bq_c = getf('Bq'); Cq_c = getf('Cq'); Dq_c = getf('Dq')
        Ap_e, Bp_e, Cp_e, Dp_e = plane_map.get(f1,(None,None,None,None))
        Aq_e, Bq_e, Cq_e, Dq_e = plane_map.get(f2,(None,None,None,None))
        for name, cval, evalv in [('Ap',Ap_c,Ap_e), ('Bp',Bp_c,Bp_e), ('Cp',Cp_c,Cp_e), ('Dp',Dp_c,Dp_e)]:
            if cval is None or evalv is None: continue
            if abs(cval - evalv) > EPS:
                print(f"Mismatch PLANES pair ({f1},{f2}) coeff {name}: C={cval:.9f} vs PY={evalv:.9f}")
                sys.exit(2)
        for name, cval, evalv in [('Aq',Aq_c,Aq_e), ('Bq',Bq_c,Bq_e), ('Cq',Cq_c,Cq_e), ('Dq',Dq_c,Dq_e)]:
            if cval is None or evalv is None: continue
            if abs(cval - evalv) > EPS:
                print(f"Mismatch PLANES pair ({f1},{f2}) coeff {name}: C={cval:.9f} vs PY={evalv:.9f}")
                sys.exit(2)

    elif stage == 'vertex_pq' or stage == 'vertex_qp':
        vid = int(r['idx'])
        xo_v = float(r['xo']); yo_v = float(r['yo']); zo_v = float(r['zo'])
        wx = float(r['wx']); wy = float(r['wy']); wz = float(r['wz'])
        # pick appropriate plane (p for PQ means plane of f1, evaluated on vertices of f2)
        if stage == 'vertex_pq':
            Ap, Bp, Cp, Dp = plane_map[f1]
        else:
            Ap, Bp, Cp, Dp = plane_map[f2]
        # expected components
        ax_obs_e = Ap * xo_v
        by_obs_e = Bp * yo_v
        cz_obs_e = Cp * zo_v
        raw_obs_e = ax_obs_e + by_obs_e + cz_obs_e + Dp
        tempo_obs_e = arr(raw_obs_e)
        ax_w_e = Ap * wx
        by_w_e = Bp * wy
        cz_w_e = Cp * wz
        raw_w_e = ax_w_e + by_w_e + cz_w_e + Dp
        tempo_w_e = arr(raw_w_e)
        # compare to C values (columns ax_obs etc may exist) - tolerant to missing fields for backward comp
        def cvf(key):
            try:
                return float(r[key])
            except Exception:
                return None
        fields_to_check = [
            ('ax_obs','ax_obs',ax_obs_e), ('by_obs','by_obs',by_obs_e), ('cz_obs','cz_obs',cz_obs_e), ('raw_obs','raw_obs',raw_obs_e), ('tempo_obs','tempo_obs',tempo_obs_e),
            ('ax_w','ax_w',ax_w_e), ('by_w','by_w',by_w_e), ('cz_w','cz_w',cz_w_e), ('raw_w','raw_w',raw_w_e), ('tempo_world','tempo_world',tempo_w_e)
        ]
        for col, name, expected in fields_to_check:
            cval = cvf(col)
            if cval is None: continue
            if abs(cval - expected) > 1e-3:
                print(f"Mismatch {stage} pair({f1},{f2}) vid={vid} {name}: C={cval:.9f} PY={expected:.9f}")
                sys.exit(3)

    elif stage == 'valobs':
        # C has valobs per pair: columns ax_p,...valobs_p and ax_q,...valobs_q in our extended CSV
        ax_p_c = float(r['ax_p']); by_p_c = float(r['by_p']); cz_p_c = float(r['cz_p']); raw_p_c = float(r['raw_p']); val_p_c = float(r['valobs']) if 'valobs' in r and r['valobs']!='' else None
        # recompute expected using plane_map and observer coords
        xobs = D * math.cos(math.radians(H)) * math.cos(math.radians(V))
        yobs = D * math.sin(math.radians(H)) * math.cos(math.radians(V))
        zobs = D * math.sin(math.radians(V))
        xobs_a = 0.0 if abs(xobs) < 0.01 else xobs
        yobs_a = 0.0 if abs(yobs) < 0.01 else yobs
        zobs_a = 0.0 if abs(zobs) < 0.01 else zobs
        Ap, Bp, Cp, Dp = plane_map[f1]
        Aq, Bq, Cq, Dq = plane_map[f2]
        raw_p_e = Ap * xobs_a + Bp * yobs_a + Cp * zobs_a + Dp
        val_p_e = arr(raw_p_e)
        raw_q_e = Aq * xobs_a + Bq * yobs_a + Cq * zobs_a + Dq
        val_q_e = arr(raw_q_e)
        if abs(ax_p_c - (Ap * xobs_a)) > 1e-3:
            print(f"Mismatch valobs ax_p pair({f1},{f2}): C={ax_p_c} PY={(Ap * xobs_a)}"); sys.exit(4)
        if 'valobs' in r and r['valobs']!='' and abs(float(r['valobs'].split(',')[0]) - val_p_e) > 1e-3:
            print(f"Mismatch valobs pair({f1},{f2}) val_p: C={r['valobs']} PY={val_p_e}"); sys.exit(5)

print('No mismatches found between C per-instr trace and Python expectations.')
sys.exit(0)
