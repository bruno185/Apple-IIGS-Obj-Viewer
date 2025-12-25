#!/usr/bin/env python3
"""Pair debug: detailed per-pair diagnostics reproducing tests 4..7 in Fixed32
Usage: python tools/pair_debug.py MODEL.obj H V W D faceA faceB
"""
import sys
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float

EPSILON_FLOAT = 0.01


def per_pair_report(path, H, V, W, D, fA, fB):
    verts, faces = parse_obj(path)
    xo, yo, zo, x2d, y2d = transform_and_project(verts, H, V, W, D)
    face_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)
    if fA < 0 or fA >= len(face_info) or fB < 0 or fB >= len(face_info):
        print('Invalid face indices')
        return
    A = face_info[fA]
    B = face_info[fB]
    print(f'PAIR {fA} vs {fB}')
    def print_face(tag, F):
        print(f" {tag}: n={len(F['idxs'])} display={F['display']} zmin={to_float(F['zmin']):.6f} zmean={to_float(F['zmean']):.6f} zmax={to_float(F['zmax']):.6f}")
        print(f"   a={to_float(F['a']):.6f}, b={to_float(F['b']):.6f}, c={to_float(F['c']):.6f}, d={to_float(F['d']):.6f}")
        print(f"   vertices (zero-based): {F['idxs']}")
    print_face('FaceA', A)
    print_face('FaceB', B)

    # Test 1 depth overlap
    print('\nTest 1: depth overlap')
    print(f" FaceA z_max={to_float(A['zmax']):.6f} FaceB z_min={to_float(B['zmin']):.6f}")
    print(f" FaceB z_max={to_float(B['zmax']):.6f} FaceA z_min={to_float(A['zmin']):.6f}")

    # Test 2/3 bbox overlap
    print('\nTest 2/3: bbox overlap')
    print(f" A minx,maxx,miny,maxy: {A['minx']},{A['maxx']},{A['miny']},{A['maxy']}")
    print(f" B minx,maxx,miny,maxy: {B['minx']},{B['maxx']},{B['miny']},{B['maxy']}")

    a1,b1,c1,d1 = to_float(A['a']), to_float(A['b']), to_float(A['c']), to_float(A['d'])
    a2,b2,c2,d2 = to_float(B['a']), to_float(B['b']), to_float(B['c']), to_float(B['d'])

    # Helper to compute per-vertex tv values for (plane of P) applied to vertices of Q
    def values_for_plane_on_vertices(a,b,c,d, verts_idx):
        vals = []
        for vi in verts_idx:
            tv = a*to_float(xo[vi]) + b*to_float(yo[vi]) + c*to_float(zo[vi]) + d
            vals.append(tv)
        return vals

    # Test 4
    print('\nTest 4: f2 on same side as observer w.r.t plane f1')
    # obs_side1 from d1
    obs_side1 = 0
    if d1 > EPSILON_FLOAT: obs_side1 = 1
    elif d1 < -EPSILON_FLOAT: obs_side1 = -1
    print(f" obs_side1 (from d1) = {obs_side1} (d1={d1:.6f})")
    if obs_side1 != 0:
        vals = values_for_plane_on_vertices(a1,b1,c1,d1, B['idxs'])
        print('  per-vertex test values (plane f1 on f2 vertices):')
        for idx,val in zip(B['idxs'], vals):
            print(f'   v{idx}: {val:.6f}')
        same = all(((v > EPSILON_FLOAT and obs_side1==1) or (v < -EPSILON_FLOAT and obs_side1==-1)) for v in vals)
        print('  all_same_side =', same)
    else:
        print('  obs on plane -> skip Test 4 (d1 approx 0)')

    # Test 5
    print('\nTest 5: f1 on opposite side of observer w.r.t plane f2 (if true -> f1 in front, no swap)')
    obs_side2 = 0
    if d2 > EPSILON_FLOAT: obs_side2 = 1
    elif d2 < -EPSILON_FLOAT: obs_side2 = -1
    print(f" obs_side2 (from d2) = {obs_side2} (d2={d2:.6f})")
    if obs_side2 != 0:
        vals = values_for_plane_on_vertices(a2,b2,c2,d2, A['idxs'])
        print('  per-vertex test values (plane f2 on f1 vertices):')
        for idx,val in zip(A['idxs'], vals):
            print(f'   v{idx}: {val:.6f}')
        all_opposite = all(((val < -EPSILON_FLOAT and obs_side2==1) or (val > EPSILON_FLOAT and obs_side2==-1)) for val in vals)
        print('  all_opposite =', all_opposite)
    else:
        print('  obs on plane -> skip Test 5 (d2 approx 0)')

    # Test 6 / 7 quick checks
    print('\nTest 6: f2 opposite of observer wrt plane f1 -> swap')
    if d1 > EPSILON_FLOAT or d1 < -EPSILON_FLOAT:
        vals = values_for_plane_on_vertices(a1,b1,c1,d1, B['idxs'])
        print('  per-vertex tv:', [f'{v:.6f}' for v in vals])
    else:
        print('  skipped (d1 approx 0)')

    print('\nTest 7: f1 same side of observer wrt plane f2 -> swap')
    if d2 > EPSILON_FLOAT or d2 < -EPSILON_FLOAT:
        vals = values_for_plane_on_vertices(a2,b2,c2,d2, A['idxs'])
        print('  per-vertex tv:', [f'{v:.6f}' for v in vals])
    else:
        print('  skipped (d2 approx 0)')

    # Provide final quick judgement using our single-shot function logic
    # Derive single-shot decision following tests 4/5/6/7 logic
    same_side1 = None
    if abs(d1) <= EPSILON_FLOAT:
        same_side1 = None
    else:
        vals = values_for_plane_on_vertices(a1,b1,c1,d1, B['idxs'])
        same_side1 = all(((v > EPSILON_FLOAT and d1>0) or (v < -EPSILON_FLOAT and d1<0)) for v in vals)
    all_opposite2 = None
    if abs(d2) <= EPSILON_FLOAT:
        all_opposite2 = None
    else:
        vals2 = values_for_plane_on_vertices(a2,b2,c2,d2, A['idxs'])
        all_opposite2 = all(((v < -EPSILON_FLOAT and d2>0) or (v > EPSILON_FLOAT and d2<0)) for v in vals2)
    print('\nDecision heuristics:')
    print(f' same_side1 = {same_side1} (Test4)')
    print(f' all_opposite2 = {all_opposite2} (Test5)')
    # Additional checks: if either test yields definitive result, give it
    if same_side1:
        print(' -> test4 says f2 behind f1 (no swap)')
    elif all_opposite2:
        print(' -> test5 says f1 in front (no swap)')
    else:
        # fallback: ambiguous
        print(' -> inconclusive by tests 4/5; original algorithm must fall back to other heuristics (bbox, zmean, tie-break)')
    # Provide tie-breaker suggestions
    print('\nTie-breaker candidates (z_mean, z_min, z_max, face index):')
    print(f"  A.zmean={to_float(A['zmean']):.6f} B.zmean={to_float(B['zmean']):.6f}")
    print(f"  A.zmin={to_float(A['zmin']):.6f} B.zmin={to_float(B['zmin']):.6f}")
    print(f"  idx order A={fA} B={fB}")

if __name__ == '__main__':
    # Hard-coded run for given problematic pairs and params
    path='cube.obj'
    H=30; V=20; W=0; D=180.0
    pairs = [(4,0),(4,1),(5,0),(5,1),(0,1)]
    for p in pairs:
        print('\n' + '='*60)
        per_pair_report(path,H,V,W,D,p[0],p[1])

if __name__ == '__main__':
    # Hard-coded run for given problematic pairs and params
    path='cube.obj'
    H=30; V=20; W=0; D=180.0
    pairs = [(4,0),(4,1),(5,0),(5,1),(0,1)]
    for p in pairs:
        print('\n' + '='*60)
        per_pair_report(path,H,V,W,D,p[0],p[1])
