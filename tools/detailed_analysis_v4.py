#!/usr/bin/env python3
"""Detailed analysis of pair_values_v4.csv vs pascal_sim expectations.
Focus on pair (1,0) first; also iterates all pairs and reports mismatches counts.
"""
import csv, math, sys
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float

OBJ = "Windows\\3d objects\\q2.obj"
PVF = "Windows\\build_vs\\bin\\Release\\pair_values_v4.csv"
H=30; V=20; W=0; D=39.214455

EPS_PLANES = 1e-3
EPS_COMPONENT = 1e-3

# Delphi Arrondi semantics
def arr(r):
    return 0.0 if abs(r) < 0.01 else r

verts, faces = parse_obj(OBJ)
xo, yo, zo, x2d, y2d = transform_and_project(verts, H, V, W, D)
faces_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)

# Build expected plane coeffs map
plane_map = {fi: (arr(to_float(f['a'])), arr(to_float(f['b'])), arr(to_float(f['c'])), arr(to_float(f['d']))) for fi,f in enumerate(faces_info)}

# Read CSV rows into a list grouped by pairs
rows = []
with open(PVF, newline='') as f:
    rdr = csv.reader(f)
    header = next(rdr)
    for r in rdr:
        rows.append(r)

# build mapping for header names to indices
hmap = {name:idx for idx,name in enumerate(header)}

# Group by pair (f1,f2) using header indices
from collections import defaultdict
pairs = defaultdict(list)
for r in rows:
    try:
        f1 = int(r[hmap['f1']]); f2 = int(r[hmap['f2']])
    except Exception:
        continue
    pairs[(f1,f2)].append(r)

# helper to get field by name using header map, returns None if not present or empty
def rget(row,name):
    if name not in hmap: return None
    idx = hmap[name]
    if idx >= len(row): return None
    v = row[idx]
    if v == '': return None
    return v

report = []

# Analyze each pair
for pair, rlist in pairs.items():
    f1,f2 = pair
    issues = []
    # Get planes row
    planes = [r for r in rlist if r[hmap.get('stage',5)]=='planes']
    if not planes:
        issues.append('no_planes_row')
        report.append((pair, issues))
        continue
    p = planes[0]
    # parse plane coefficients robustly
    # preferred: use named columns if present
    def pfloat(x):
        try: return float(x)
        except Exception: return None
    Ap_c = pfloat(rget(p,'Ap'))
    if Ap_c is None:
        # fallback: extract numeric tokens after the 'stage' column
        # find index of stage in header
        stage_idx = hmap.get('stage', 5)
        tokens = [t for t in p[stage_idx+1:] if t != '']
        nums = []
        for t in tokens:
            try: nums.append(float(t))
            except: pass
        if len(nums) >= 8:
            Ap_c,Bp_c,Cp_c,Dp_c,Aq_c,Bq_c,Cq_c,Dq_c = nums[:8]
        else:
            Ap_c = Bp_c = Cp_c = Dp_c = Aq_c = Bq_c = Cq_c = Dq_c = None
    else:
        Bp_c = pfloat(rget(p,'Bp')); Cp_c = pfloat(rget(p,'Cp')); Dp_c = pfloat(rget(p,'Dp'))
        Aq_c = pfloat(rget(p,'Aq')); Bq_c = pfloat(rget(p,'Bq')); Cq_c = pfloat(rget(p,'Cq')); Dq_c = pfloat(rget(p,'Dq'))
    Ap_e,Bp_e,Cp_e,Dp_e = plane_map.get(f1,(None,None,None,None))
    Aq_e,Bq_e,Cq_e,Dq_e = plane_map.get(f2,(None,None,None,None))
    # Compare planes
    for name, cval, evalv in [('Ap',Ap_c,Ap_e),('Bp',Bp_c,Bp_e),('Cp',Cp_c,Cp_e),('Dp',Dp_c,Dp_e)]:
        if cval is None or evalv is None: continue
        if abs(cval-evalv) > EPS_PLANES:
            issues.append(f'plane_mismatch_{name}:C={cval:.9f},PY={evalv:.9f}')
    # Now per-vertex checks
    v_pq = [r for r in rlist if r[hmap.get('stage',5)]=='vertex_pq']
    v_qp = [r for r in rlist if r[hmap.get('stage',5)]=='vertex_qp']
    def getf(row,name):
        v = rget(row,name)
        if v is None: return None
        try: return float(v)
        except: return None
    def check_vertex_rows(vrows, plane_for_eval):
        for r in vrows:
            # parse by header indices
            try:
                vid = int(r[hmap['idx']])
            except Exception:
                # fallback: first token
                vid = int(r[0]) if r and r[0].isdigit() else None
            xo_v = getf(r,'xo'); yo_v = getf(r,'yo'); zo_v = getf(r,'zo')
            wx = getf(r,'wx'); wy = getf(r,'wy'); wz = getf(r,'wz')
            ax_obs_c = getf(r,'ax_obs'); by_obs_c = getf(r,'by_obs'); cz_obs_c = getf(r,'cz_obs'); raw_obs_c = getf(r,'raw_obs'); tempo_obs_c = getf(r,'tempo_obs')
            ax_w_c = getf(r,'ax_w'); by_w_c = getf(r,'by_w'); cz_w_c = getf(r,'cz_w'); raw_w_c = getf(r,'raw_w'); tempo_w_c = getf(r,'tempo_world')
            if None in (xo_v,yo_v,zo_v,wx,wy,wz):
                # attempt positional fallback: tokens after stage
                stage_idx = hmap.get('stage',5)
                tokens = [t for t in r[stage_idx+1:] if t!='']
                nums = []
                for t in tokens:
                    try: nums.append(float(t))
                    except: pass
                if len(nums) >= 11:
                    # mapping: idx,xo,yo,zo,wx,wy,wz,ax_obs,by_obs,cz_obs,raw_obs,tempo_obs,ax_w,by_w,cz_w,raw_w,tempo_world
                    # but some rows may omit some fields; we map first 11
                    idxp = int(nums[0])
                    if vid is None: vid = idxp
                    xo_v,yo_v,zo_v,wx,wy,wz = nums[1:7]
                    ax_obs_c,by_obs_c,cz_obs_c,raw_obs_c,tempo_obs_c = (nums[7:12] + [None]*5)[:5]
                    # next fields as world comps
                    if len(nums) >= 17:
                        ax_w_c,by_w_c,cz_w_c,raw_w_c,tempo_w_c = nums[12:17]
            # compute expected using python planes
            Ap,Bp,Cp,Dp = plane_map[plane_for_eval]
            # If any of xo_v,wx... missing, skip
            if None in (xo_v,yo_v,zo_v,wx,wy,wz):
                issues.append(f'vid{vid}_missing_coords')
                continue
            ax_obs_e = Ap * xo_v; by_obs_e = Bp * yo_v; cz_obs_e = Cp * zo_v; raw_obs_e = ax_obs_e + by_obs_e + cz_obs_e + Dp; tempo_obs_e = arr(raw_obs_e)
            ax_w_e = Ap * wx; by_w_e = Bp * wy; cz_w_e = Cp * wz; raw_w_e = ax_w_e + by_w_e + cz_w_e + Dp; tempo_w_e = arr(raw_w_e)
            if ax_obs_c is not None and abs(ax_obs_c - ax_obs_e) > EPS_COMPONENT:
                issues.append(f'vid{vid}_ax_obs_mismatch_C={ax_obs_c:.6f}_PY={ax_obs_e:.6f}')
            if raw_obs_c is not None and abs(raw_obs_c - raw_obs_e) > EPS_COMPONENT:
                issues.append(f'vid{vid}_raw_obs_mismatch_C={raw_obs_c:.6f}_PY={raw_obs_e:.6f}')
            if tempo_obs_c is not None and abs(tempo_obs_c - tempo_obs_e) > EPS_COMPONENT:
                issues.append(f'vid{vid}_tempo_obs_mismatch_C={tempo_obs_c:.6f}_PY={tempo_obs_e:.6f}')
            if ax_w_c is not None and abs(ax_w_c - ax_w_e) > EPS_COMPONENT:
                issues.append(f'vid{vid}_ax_w_mismatch_C={ax_w_c:.6f}_PY={ax_w_e:.6f}')
            if tempo_w_c is not None and abs(tempo_w_c - tempo_w_e) > EPS_COMPONENT:
                issues.append(f'vid{vid}_tempo_w_mismatch_C={tempo_w_c:.6f}_PY={tempo_w_e:.6f}')
    check_vertex_rows(v_pq, f1)
    check_vertex_rows(v_qp, f2)
    # valobs compare
    valrows = [r for r in rlist if r[hmap.get('stage',5)]=='valobs']
    for vr in valrows:
        stage_idx = hmap.get('stage',5)
        tokens = [t for t in vr[stage_idx+1:] if t!='']
        nums = []
        for t in tokens:
            try: nums.append(float(t))
            except: pass
        # expect at least 10 numbers: ax_p,by_p,cz_p,raw_p,val_p,ax_q,by_q,cz_q,raw_q,val_q
        if len(nums) >= 10:
            ax_p_c,by_p_c,cz_p_c,raw_p_c,val_p_c,ax_q_c,by_q_c,cz_q_c,raw_q_c,val_q_c = nums[:10]
        else:
            ax_p_c=by_p_c=cz_p_c=raw_p_c=val_p_c=ax_q_c=by_q_c=cz_q_c=raw_q_c=val_q_c=None
        # expected
        xobs = D*math.cos(math.radians(H))*math.cos(math.radians(V))
        yobs = D*math.sin(math.radians(H))*math.cos(math.radians(V))
        zobs = D*math.sin(math.radians(V))
        xobs_a = 0.0 if abs(xobs)<0.01 else xobs
        yobs_a = 0.0 if abs(yobs)<0.01 else yobs
        zobs_a = 0.0 if abs(zobs)<0.01 else zobs
        Ap, Bp, Cp, Dp = plane_map[f1]
        Aq, Bq, Cq, Dq = plane_map[f2]
        raw_p_e = Ap * xobs_a + Bp * yobs_a + Cp * zobs_a + Dp
        val_p_e = arr(raw_p_e)
        raw_q_e = Aq * xobs_a + Bq * yobs_a + Cq * zobs_a + Dq
        val_q_e = arr(raw_q_e)
        if val_p_c is not None and abs(val_p_c - val_p_e) > EPS_COMPONENT:
            issues.append(f'valobs_p_mismatch_C={val_p_c:.6f}_PY={val_p_e:.6f}')
        if val_q_c is not None and abs(val_q_c - val_q_e) > EPS_COMPONENT:
            issues.append(f'valobs_q_mismatch_C={val_q_c:.6f}_PY={val_q_e:.6f}')
    report.append((pair, list(set(issues))))

# Print summary
for pair, issues in report:
    f1,f2 = pair
    if not issues:
        print(f'Pair ({f1},{f2}): OK')
    else:
        print(f'Pair ({f1},{f2}):')
        for it in issues:
            print('  -', it)

# helper
def _looks_like_number(s):
    try:
        float(s); return True
    except Exception:
        return False
