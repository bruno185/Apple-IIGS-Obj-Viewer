#!/usr/bin/env python3
"""
pair_diag.py

Diagnostic script that reproduces the relevant parts of the painter (Newell/Sancha)
algorithm for a single pair of faces, comparing plane computation using:
 - 3-vertex formula (current C implementation), and
 - Newell's method (robust for non-planar polygons)

Outputs a GS-friendly CSV (default: PAIRPY.CSV) with detailed per-face and per-vertex
information so you can compare with the in-app PAIRDIAG.CSV.

Usage:
  python scripts/pair_diag.py q1.obj --angle_h 0 --angle_v 0 --angle_w 0 --distance 300

Options:
  --outfile NAME   : output CSV name (default PAIRPY.CSV, must be <=15 chars, no underscore)
  --verbose        : print summary to console

"""
import sys
import math
import argparse
from pathlib import Path


def parse_obj(path):
    verts = []
    faces = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('v '):
                parts = line.split()
                x = float(parts[1]); y = float(parts[2]); z = float(parts[3])
                # Match C code: OBJ often Z-up; they swap (y,z)
                verts.append((x, z, y))
            elif line.startswith('f '):
                parts = line.split()
                idx = []
                for p in parts[1:]:
                    # support formats like v, v//vn, v/vt/vn
                    v = p.split('/')[0]
                    try:
                        vi = int(v) - 1
                    except:
                        vi = 0
                    idx.append(vi)
                faces.append(idx)
    return verts, faces


def deg2rad(d):
    return d * math.pi / 180.0


def transform_vertex(v, angle_h, angle_v, angle_w, distance):
    # Follow C implementation (float version)
    x, y, z = v
    ch = math.cos(deg2rad(angle_h)); sh = math.sin(deg2rad(angle_h))
    cv = math.cos(deg2rad(angle_v)); sv = math.sin(deg2rad(angle_v))
    # products
    ch_cv = ch * cv
    sh_cv = sh * cv
    ch_sv = ch * sv
    sh_sv = sh * sv
    term1 = x * ch_cv
    term2 = y * sh_cv
    term3 = z * sv
    zo = -term1 - term2 - term3 + distance
    if zo > 0:
        xo = - x * sh + y * ch
        yo = - x * ch_sv - y * sh_sv + z * cv
    else:
        xo = 0.0; yo = 0.0
    return xo, yo, zo


def plane_from_first3(v0, v1, v2):
    # compute a,b,c,d using 3-vertex formula (no normalization)
    x1,y1,z1 = v0
    x2,y2,z2 = v1
    x3,y3,z3 = v2
    a = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
    b = -x1*(z2 - z3) + x2*(z1 - z3) - x3*(z1 - z2)
    c = x1*(y2 - y3) - x2*(y1 - y3) + x3*(y1 - y2)
    d = -x1*(y2*z3 - y3*z2) + x2*(y1*z3 - y3*z1) - x3*(y1*z2 - y2*z1)
    return (a,b,c,d)


def plane_newell(poly):
    # poly: list of (x,y,z) observer-space coords
    nx = ny = nz = 0.0
    n = len(poly)
    for i in range(n):
        x1,y1,z1 = poly[i]
        x2,y2,z2 = poly[(i+1)%n]
        nx += (y1 - y2) * (z1 + z2)
        ny += (z1 - z2) * (x1 + x2)
        nz += (x1 - x2) * (y1 + y2)
    # d = -n . p0
    d = -(nx * poly[0][0] + ny * poly[0][1] + nz * poly[0][2])
    return (nx, ny, nz, d)


def bbox_2d(vlist):
    xs = [v[0] for v in vlist]
    ys = [v[1] for v in vlist]
    return int(min(xs)), int(max(xs)), int(min(ys)), int(max(ys))


def cmp_faces_by_zmean(fa, fb):
    if fa['z_mean'] > fb['z_mean']: return -1
    if fa['z_mean'] < fb['z_mean']: return 1
    return 0


def pair_should_swap_py(faces, vtx_obs, f1_idx, f2_idx, use_newell=False):
    f1 = faces[f1_idx]; f2 = faces[f2_idx]
    # Test1: depth
    if f2['z_max'] <= f1['z_min']: return -1
    if f1['z_max'] <= f2['z_min']: return 1
    # bbox quick reject (we'll use 2D bbox on projected x2d/y2d if available)
    if f1['maxx'] <= f2['minx'] or f2['maxx'] <= f1['minx']: return -1
    if f1['maxy'] <= f2['miny'] or f2['maxy'] <= f1['miny']: return -1

    # plane coefficients
    if use_newell:
        a1,b1,c1,d1 = f1['plane_newell']
        a2,b2,c2,d2 = f2['plane_newell']
    else:
        a1,b1,c1,d1 = f1['plane_first3']
        a2,b2,c2,d2 = f2['plane_first3']

    n1 = len(f1['vertices']); n2 = len(f2['vertices'])

    # compute maxabs and epsol
    maxabs1 = abs(d1)
    for v in f2['vertices']:
        tv = a1*v['xo'] + b1*v['yo'] + c1*v['zo'] + d1
        if abs(tv) > maxabs1: maxabs1 = abs(tv)
    eps_rel1 = max(1e-4, maxabs1 * 1e-6)
    obs1 = 0
    if d1 > eps_rel1: obs1 = 1
    elif d1 < -eps_rel1: obs1 = -1
    if obs1 != 0:
        pos = neg = 0
        for v in f2['vertices']:
            tv = a1*v['xo'] + b1*v['yo'] + c1*v['zo'] + d1
            if abs(tv) <= eps_rel1: pass
            elif tv > 0: pos += 1
            else: neg += 1
        thr = (3*n2 + 3)//4
        if (obs1==1 and pos>=thr) or (obs1==-1 and neg>=thr):
            return -1
        if (obs1==1 and neg>=thr) or (obs1==-1 and pos>=thr):
            return 1

    # symmetric for obs2
    maxabs2 = abs(d2)
    for v in f1['vertices']:
        tv = a2*v['xo'] + b2*v['yo'] + c2*v['zo'] + d2
        if abs(tv) > maxabs2: maxabs2 = abs(tv)
    eps_rel2 = max(1e-4, maxabs2 * 1e-6)
    obs2 = 0
    if d2 > eps_rel2: obs2 = 1
    elif d2 < -eps_rel2: obs2 = -1
    if obs2 != 0:
        pos2 = neg2 = 0
        for v in f1['vertices']:
            tv = a2*v['xo'] + b2*v['yo'] + c2*v['zo'] + d2
            if abs(tv) <= eps_rel2: pass
            elif tv > 0: pos2 += 1
            else: neg2 += 1
        thr2 = (3*n1 + 3)//4
        if (obs2==1 and neg2>=thr2) or (obs2==-1 and pos2>=thr2):
            return -1
        if (obs2==1 and pos2>=thr2) or (obs2==-1 and neg2>=thr2):
            return 1

    # Test6 and 7 (when obs1 or obs2 == 0)
    if obs1 == 0:
        pos = neg = 0
        for v in f2['vertices']:
            tv = a1*v['xo'] + b1*v['yo'] + c1*v['zo'] + d1
            if abs(tv) <= eps_rel1: pass
            elif tv > 0: pos += 1
            else: neg += 1
        thr = (3*n2 + 3)//4
        if (obs1==1 and neg>=thr) or (obs1==-1 and pos>=thr): return 1
    if obs2 == 0:
        pos2 = neg2 = 0
        for v in f1['vertices']:
            tv = a2*v['xo'] + b2*v['yo'] + c2*v['zo'] + d2
            if abs(tv) <= eps_rel2: pass
            elif tv > 0: pos2 += 1
            else: neg2 += 1
        thr2 = (3*n1 + 3)//4
        if (obs2==1 and pos2>=thr2) or (obs2==-1 and neg2>=thr2): return 1

    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('objfile', help='OBJ filename')
    parser.add_argument('--angle_h', type=float, default=0.0)
    parser.add_argument('--angle_v', type=float, default=0.0)
    parser.add_argument('--angle_w', type=float, default=0.0)
    parser.add_argument('--distance', type=float, default=300.0)
    parser.add_argument('--obsfile', default=None, help='Optional existing PAIRDIAG/PAIRPY CSV with observer-space coords')
    parser.add_argument('--outfile', default='PAIRPY.CSV')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    verts, faces_idx = parse_obj(args.objfile)
    if len(faces_idx) < 2:
        print('Need at least 2 faces to compare'); return

    # Build observer-space vertices: either by transforming 3D verts, or by importing
    # existing observer-space coords from a PAIRDIAG/PAIRPY file (args.obsfile)
    v_obs = []
    if args.obsfile:
        # Parse lines 'eval,...' expecting format: eval,face_src,vertex_idx,tv,xo,yo,zo,...
        obs_map = {}
        with open(args.obsfile, 'r', encoding='utf-8', errors='ignore') as of:
            for line in of:
                line = line.strip()
                if not line: continue
                if not line.startswith('eval,'): continue
                parts = line.split(',')
                # parts: ['eval', face_src, vertex_idx, tv, xo, yo, zo, ...]
                try:
                    vidx = int(parts[2])
                    xo = float(parts[4]); yo = float(parts[5]); zo = float(parts[6])
                    obs_map[vidx] = {'xo': xo, 'yo': yo, 'zo': zo}
                except Exception:
                    continue
        max_idx = max(max(range(len(verts))), max(obs_map.keys()) if obs_map else -1)
        for i in range(max_idx+1):
            if i in obs_map:
                v_obs.append(obs_map[i])
            else:
                # fall back: transform original 3D vertex
                if i < len(verts):
                    xo, yo, zo = transform_vertex(verts[i], args.angle_h, args.angle_v, args.angle_w, args.distance)
                    v_obs.append({'xo': xo, 'yo': yo, 'zo': zo})
                else:
                    v_obs.append({'xo': 0.0, 'yo': 0.0, 'zo': -99999.0})
    else:
        for v in verts:
            xo, yo, zo = transform_vertex(v, args.angle_h, args.angle_v, args.angle_w, args.distance)
            v_obs.append({'xo': xo, 'yo': yo, 'zo': zo})

    # Build face structures
    faces = []
    for fi, fidx in enumerate(faces_idx):
        face = {'idx': fi, 'indices': fidx}
        verts_obs = []
        xs = []
        ys = []
        zs = []
        xo_arr = []
        yo_arr = []
        zo_arr = []
        for vi in fidx:
            x,y,z = verts[vi]
            xs.append(x); ys.append(y); zs.append(z)
            xo_arr.append(v_obs[vi]['xo']); yo_arr.append(v_obs[vi]['yo']); zo_arr.append(v_obs[vi]['zo'])
            verts_obs.append({'xo': v_obs[vi]['xo'], 'yo': v_obs[vi]['yo'], 'zo': v_obs[vi]['zo'], 'idx': vi})
        face['vertices'] = verts_obs
        face['z_min'] = min(zo_arr)
        face['z_max'] = max(zo_arr)
        face['z_mean'] = sum(zo_arr)/len(zo_arr)
        # 2D bbox on projected coords (approx using xo/zo etc.) - simple projection: x2d = xo/zo, y2d = yo/zo
        x2d = [(xo / zo) if zo>0 else -99999 for xo,zo in zip(xo_arr, zo_arr)]
        y2d = [(yo / zo) if zo>0 else -99999 for yo,zo in zip(yo_arr, zo_arr)]
        face['minx'], face['maxx'], face['miny'], face['maxy'] = min(x2d), max(x2d), min(y2d), max(y2d)
        # planes
        # first-3 method
        idx0 = fidx[0]; idx1 = fidx[1]; idx2 = fidx[2]
        v0 = (v_obs[idx0]['xo'], v_obs[idx0]['yo'], v_obs[idx0]['zo'])
        v1 = (v_obs[idx1]['xo'], v_obs[idx1]['yo'], v_obs[idx1]['zo'])
        v2 = (v_obs[idx2]['xo'], v_obs[idx2]['yo'], v_obs[idx2]['zo'])
        face['plane_first3'] = plane_from_first3(v0, v1, v2)
        # Newell
        poly_obs = [(vv['xo'], vv['yo'], vv['zo']) for vv in verts_obs]
        face['plane_newell'] = plane_newell(poly_obs)
        faces.append(face)

    # Compare pair 0 and 1
    outp = Path(args.outfile)
    with open(outp, 'w', encoding='utf-8') as df:
        df.write('section,face,a3,b3,c3,d3,aN,bN,cN,dN,z_min,z_mean,z_max,vertex_indices\n')
        for fi in range(2):
            f = faces[fi]
            a3,b3,c3,d3 = f['plane_first3']
            aN,bN,cN,dN = f['plane_newell']
            df.write(f'face,{fi},{a3:.6f},{b3:.6f},{c3:.6f},{d3:.6f},{aN:.6f},{bN:.6f},{cN:.6f},{dN:.6f},{f["z_min"]:.6f},{f["z_mean"]:.6f},{f["z_max"]:.6f},"')
            for v in f['indices']:
                df.write(str(v) + ' ')
            df.write('"\n')

        # Simulate triangulation (fan) and evaluate triangle ordering to check if triangulation resolves ambiguity
        tri_faces = []
        orig_to_tri = []
        for fi in range(2):
            idxs = faces[fi]['indices']
            if len(idxs) <= 3:
                tri_faces.append({'indices': idxs, 'orig': fi})
            else:
                for k in range(1, len(idxs)-1):
                    tri_faces.append({'indices': [idxs[0], idxs[k], idxs[k+1]], 'orig': fi})
        # compute plane and stats per triangle
        for ti,tf in enumerate(tri_faces):
            idxs = tf['indices']
            # build temporary face structure
            tmp = {'indices': idxs, 'vertices': [{'xo': v_obs[i]['xo'], 'yo': v_obs[i]['yo'], 'zo': v_obs[i]['zo'], 'idx': i} for i in idxs]}
            xo_arr = [v['xo'] for v in tmp['vertices']]; yo_arr = [v['yo'] for v in tmp['vertices']]; zo_arr = [v['zo'] for v in tmp['vertices']]
            tmp['z_min'] = min(zo_arr); tmp['z_mean'] = sum(zo_arr)/len(zo_arr); tmp['z_max'] = max(zo_arr)
            tmp['plane_first3'] = plane_from_first3((tmp['vertices'][0]['xo'],tmp['vertices'][0]['yo'],tmp['vertices'][0]['zo']),(tmp['vertices'][1]['xo'],tmp['vertices'][1]['yo'],tmp['vertices'][1]['zo']),(tmp['vertices'][2]['xo'],tmp['vertices'][2]['yo'],tmp['vertices'][2]['zo']))
            tmp['plane_newell'] = plane_newell([(v['xo'],v['yo'],v['zo']) for v in tmp['vertices']])
            # compute 2D bbox projection used by tests
            x2d = [(v['xo']/v['zo']) if v['zo']>0 else -99999.0 for v in tmp['vertices']]
            y2d = [(v['yo']/v['zo']) if v['zo']>0 else -99999.0 for v in tmp['vertices']]
            tmp['minx'], tmp['maxx'], tmp['miny'], tmp['maxy'] = min(x2d), max(x2d), min(y2d), max(y2d)
            tri_faces[ti] = tmp
        # Evaluate ordering among triangles
        # We'll compute pair_should_swap for each adjacent pair in tri list (as if they were faces)
        df.write('triangulation_check,triangle_count,{}\n'.format(len(tri_faces)))
        # Detailed triangle-pair comparisons and majority voting
        # We'll compare every triangle of face0 with every triangle of face1
        # and tally votes for keep(-1), swap(1), undetermined(0)
        tri_votes_first3 = {'keep':0, 'swap':0, 'und':0}
        tri_votes_newell = {'keep':0, 'swap':0, 'und':0}
        # map triangle index back to original face (0 or 1)
        tri_orig = [0 if (ti < (len(faces[0]['indices'])-2)) else 1 for ti in range(len(tri_faces))]
        # Build lists of triangles per face
        tris_f0 = [i for i,tf in enumerate(tri_faces) if tf['vertices'][0]['idx'] in faces[0]['indices']]
        tris_f1 = [i for i,tf in enumerate(tri_faces) if tf['vertices'][0]['idx'] in faces[1]['indices']]
        for i in tris_f0:
            for j in tris_f1:
                r3 = pair_should_swap_py(tri_faces, v_obs, i, j, use_newell=False)
                rN = pair_should_swap_py(tri_faces, v_obs, i, j, use_newell=True)
                df.write('tri_pair,{},{},{},{}\n'.format(i, j, r3, rN))
                # tally
                if r3 == -1: tri_votes_first3['keep'] += 1
                elif r3 == 1: tri_votes_first3['swap'] += 1
                else: tri_votes_first3['und'] += 1
                if rN == -1: tri_votes_newell['keep'] += 1
                elif rN == 1: tri_votes_newell['swap'] += 1
                else: tri_votes_newell['und'] += 1
        # majority decision
        def decide(v):
            if v['swap'] > v['keep']: return 1
            if v['keep'] > v['swap']: return -1
            return 0
        dec3 = decide(tri_votes_first3)
        decN = decide(tri_votes_newell)
        df.write('tri_vote,first3,keep={},swap={},und={};decision={}\n'.format(tri_votes_first3['keep'], tri_votes_first3['swap'], tri_votes_first3['und'], dec3))
        df.write('tri_vote,newell,keep={},swap={},und={};decision={}\n'.format(tri_votes_newell['keep'], tri_votes_newell['swap'], tri_votes_newell['und'], decN))
        # Always attempt local occlusion check when projected bboxes overlap for extra robustness
        if tri_votes_first3['keep'] != tri_votes_first3['swap']:
            loc = local_occlusion_global(tri_faces, tris_f0, tris_f1)
            if loc[0] != 0:
                skeep, sswap, sund, svalid, mean0, mean1 = loc[1]
                df.write('TRI_LOCAL_OCCLUSION_FIRST3,keep={},swap={},und={},samples={},mean_d0={:.6f},mean_d1={:.6f},decision={}\n'.format(skeep, sswap, sund, svalid, mean0, mean1, loc[0]))
                dec3 = loc[0]
        if tri_votes_newell['keep'] != tri_votes_newell['swap']:
            loc = local_occlusion_global(tri_faces, tris_f0, tris_f1)
            if loc[0] != 0:
                skeep, sswap, sund, svalid, mean0, mean1 = loc[1]
                df.write('TRI_LOCAL_OCCLUSION_NEWELL,keep={},swap={},und={},samples={},mean_d0={:.6f},mean_d1={:.6f},decision={}\n'.format(skeep, sswap, sund, svalid, mean0, mean1, loc[0]))
                decN = loc[0]

        # If the triangle voting is inconclusive, attempt overlap-based sampling to break ties
        def sample_overlap(tri_faces, tris_f0, tris_f1):
            sample_keep = sample_swap = sample_und = 0
            def point_in_tri_pt(p, a, b, c):
                def orient(a,b,c): return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
                o1 = orient(a,b,p); o2 = orient(b,c,p); o3 = orient(c,a,p)
                s1 = (o1 >= 0); s2 = (o2 >= 0); s3 = (o3 >= 0)
                s1n = (o1 <= 0); s2n = (o2 <= 0); s3n = (o3 <= 0)
                return (s1 and s2 and s3) or (s1n and s2n and s3n)
            for i in tris_f0:
                A = tri_faces[i]
                Ax = [v['xo']/v['zo'] for v in A['vertices']]
                Ay = [v['yo']/v['zo'] for v in A['vertices']]
                Az = [v['zo'] for v in A['vertices']]
                for j in tris_f1:
                    B = tri_faces[j]
                    Bx = [v['xo']/v['zo'] for v in B['vertices']]
                    By = [v['yo']/v['zo'] for v in B['vertices']]
                    Bz = [v['zo'] for v in B['vertices']]
                    # bbox quick-skip
                    if max(Ax) < min(Bx) or max(Bx) < min(Ax) or max(Ay) < min(By) or max(By) < min(Ay):
                        continue
                    # vertex-in-triangle tests
                    found = False
                    sx = sy = None
                    for vi in range(3):
                        if point_in_tri_pt((Ax[vi],Ay[vi]), (Bx[0],By[0]), (Bx[1],By[1]), (Bx[2],By[2])):
                            sx, sy = Ax[vi], Ay[vi]; found = True; break
                    for vj in range(3):
                        if not found and point_in_tri_pt((Bx[vj],By[vj]), (Ax[0],Ay[0]), (Ax[1],Ay[1]), (Ax[2],Ay[2])):
                            sx, sy = Bx[vj], By[vj]; found = True; break
                    # edge intersection fallback
                    if not found:
                        inter = None
                        for e1 in range(3):
                            p1 = (Ax[e1], Ay[e1]); p2 = (Ax[(e1+1)%3], Ay[(e1+1)%3])
                            for e2 in range(3):
                                p3 = (Bx[e2], By[e2]); p4 = (Bx[(e2+1)%3], By[(e2+1)%3])
                                if seg_intersect(p1,p2,p3,p4):
                                    # approximate intersection sample by segment midpoint
                                    inter = ((p1[0]+p2[0])/2.0, (p1[1]+p2[1])/2.0)
                                    break
                            if inter: break
                        if inter:
                            sx, sy = inter[0], inter[1]; found = True
                    if not found:
                        continue
                    # barycentric in 2d for A
                    def bary_weights(px,py, tx,ty):
                        denom = (ty[1]-ty[2])*(tx[0]-tx[2]) + (tx[2]-tx[1])*(ty[0]-ty[2])
                        if abs(denom) < 1e-9: return None
                        w0 = ((ty[1]-ty[2])*(px-tx[2]) + (tx[2]-tx[1])*(py-ty[2])) / denom
                        w1 = ((ty[2]-ty[0])*(px-tx[2]) + (tx[0]-tx[2])*(py-ty[2])) / denom
                        w2 = 1.0 - w0 - w1
                        return (w0,w1,w2)
                    wA = bary_weights(sx,sy, Ax,Ay)
                    wB = bary_weights(sx,sy, Bx,By)
                    if (not wA) or (not wB): sample_und += 1; continue
                    zoA = wA[0]*Az[0] + wA[1]*Az[1] + wA[2]*Az[2]
                    zoB = wB[0]*Bz[0] + wB[1]*Bz[1] + wB[2]*Bz[2]
                    if abs(zoA - zoB) < 1e-6: sample_und += 1
                    elif zoA < zoB: sample_keep += 1
                    else: sample_swap += 1
            return (sample_keep, sample_swap, sample_und)

        # Global local occlusion helper (callable before triangle vote logic)
        def local_occlusion_global(tri_faces, tris_f0, tris_f1):
            # compute projected bbox overlap and do 3x3 grid sampling
            def bbox_from_tris(orig):
                xs = []; ys = []
                for ti,tf in enumerate(tri_faces):
                    if tf.get('orig', None) != orig: continue
                    for vid in tf['indices']:
                        xo = v_obs[vid]['xo']; yo = v_obs[vid]['yo']; zo = v_obs[vid]['zo']
                        if zo <= 0: continue
                        xs.append(xo/zo); ys.append(yo/zo)
                if not xs: return None
                return min(xs), max(xs), min(ys), max(ys)
            if not tris_f0 or not tris_f1: return (0, None)
            b0 = bbox_from_tris(0); b1 = bbox_from_tris(1)
            if not b0 or not b1: return (0, None)
            x0min,x0max,y0min,y0max = b0; x1min,x1max,y1min,y1max = b1
            oxmin,oxmax = max(x0min,x1min), min(x0max,x1max)
            oymin,oymax = max(y0min,y1min), min(y0max,y1max)
            if oxmax <= oxmin or oymax <= oymin: return (0, None)
            cx, cy = (oxmin+oxmax)/2.0, (oymin+oymax)/2.0
            dx = oxmax - oxmin; dy = oymax - oymin
            def point_in_poly_pt(px,py, poly_vertices):
                crossings = 0
                n = len(poly_vertices)
                for i in range(n):
                    a = poly_vertices[i]; b = poly_vertices[(i+1)%n]
                    x0,y0 = a['xo']/a['zo'], a['yo']/a['zo']
                    x1,y1 = b['xo']/b['zo'], b['yo']/b['zo']
                    cond = ((y0 <= py) and (y1 > py)) or ((y0 > py) and (y1 <= py))
                    if cond:
                        xint = x0 + (py - y0) * (x1 - x0) / (y1 - y0)
                        if xint > px: crossings += 1
                return (crossings & 1)
            def depth_at_point(poly_vertices, px, py):
                n = len(poly_vertices)
                for k in range(1, n-1):
                    a = poly_vertices[0]; b = poly_vertices[k]; c = poly_vertices[k+1]
                    Ax,Ay,Az = a['xo']/a['zo'], a['yo']/a['zo'], a['zo']
                    Bx,By,Bz = b['xo']/b['zo'], b['yo']/b['zo'], b['zo']
                    Cx,Cy,Cz = c['xo']/c['zo'], c['yo']/c['zo'], c['zo']
                    o1 = (Bx-Ax)*(py-Ay) - (By-Ay)*(px-Ax)
                    o2 = (Cx-Bx)*(py-By) - (Cy-By)*(px-Bx)
                    o3 = (Ax-Cx)*(py-Cy) - (Ay-Cy)*(px-Cx)
                    if not ((o1>=0 and o2>=0 and o3>=0) or (o1<=0 and o2<=0 and o3<=0)):
                        continue
                    denom = (By - Cy)*(Ax - Cx) + (Cx - Bx)*(Ay - Cy)
                    if abs(denom) < 1e-9: return None
                    wA = ((By - Cy)*(px - Cx) + (Cx - Bx)*(py - Cy)) / denom
                    wB = ((Cy - Ay)*(px - Cx) + (Ax - Cx)*(py - Cy)) / denom
                    wC = 1.0 - wA - wB
                    if wA < -1e-5 or wB < -1e-5 or wC < -1e-5: return None
                    return wA*Az + wB*Bz + wC*Cz
                return None
            # build poly vertices arrays (collect unique vertices from tri list)
            def collect_poly(orig):
                verts = []
                vis = set()
                for ti,tf in enumerate(tri_faces):
                    if tf.get('orig', None) != orig: continue
                    for v in tf['indices']:
                        if v not in vis:
                            verts.append({'xo': v_obs[v]['xo'], 'yo': v_obs[v]['yo'], 'zo': v_obs[v]['zo'], 'idx': v})
                            vis.add(v)
                return verts
            poly0 = collect_poly(0)
            poly1 = collect_poly(1)
            skeep = sswap = sund = svalid = 0
            sumd0 = sumd1 = 0.0
            for ix in (-1,0,1):
                for iy in (-1,0,1):
                    px = cx + ix * (dx * 0.25)
                    py = cy + iy * (dy * 0.25)
                    if not point_in_poly_pt(px,py,poly0):
                        continue
                    if not point_in_poly_pt(px,py,poly1):
                        continue
                    d0 = depth_at_point(poly0, px, py); d1 = depth_at_point(poly1, px, py)
                    if d0 is None or d1 is None: continue
                    svalid += 1; sumd0 += d0; sumd1 += d1
                    epsd = 1e-4
                    if d0 < d1 - epsd: skeep += 1
                    elif d1 < d0 - epsd: sswap += 1
                    else: sund += 1
            if svalid > 0:
                mean0 = sumd0 / svalid; mean1 = sumd1 / svalid
                if sswap > skeep:
                    return (1, (skeep, sswap, sund, svalid, mean0, mean1))
                if skeep >= sswap:
                    return (-1, (skeep, sswap, sund, svalid, mean0, mean1))
            return (0, None)
        def local_occlusion(tri_faces, tris_f0, tris_f1):
            # Grid sampling (3x3) inside projected bbox intersection; cheaper and more robust than single centroid
            def bbox(face):
                xs = [v['xo']/v['zo'] for v in face['vertices']]
                ys = [v['yo']/v['zo'] for v in face['vertices']]
                return min(xs), max(xs), min(ys), max(ys)
            # collect full face polygons from triangles
            def face_vertices_from_tris(orig):
                verts = []
                vis = set()
                for ti,tf in enumerate(tri_faces):
                    if tf.get('orig', None) != orig: continue
                    for v in tf['indices']:
                        if v not in vis:
                            verts.append({'xo': v_obs[v]['xo'], 'yo': v_obs[v]['yo'], 'zo': v_obs[v]['zo'], 'idx': v})
                            vis.add(v)
                return verts
            if not tris_f0 or not tris_f1: return (0, None)
            poly0 = {'vertices': face_vertices_from_tris(0)}
            poly1 = {'vertices': face_vertices_from_tris(1)}
            if not poly0['vertices'] or not poly1['vertices']: return (0, None)
            x0min,x0max,y0min,y0max = bbox(poly0)
            x1min,x1max,y1min,y1max = bbox(poly1)
            oxmin,oxmax = max(x0min,x1min), min(x0max,x1max)
            oymin,oymax = max(y0min,y1min), min(y0max,y1max)
            if oxmax <= oxmin or oymax <= oymin: return (0, None)
            cx, cy = (oxmin+oxmax)/2.0, (oymin+oymax)/2.0
            dx = oxmax - oxmin; dy = oymax - oymin
            def point_in_poly_pt(px,py, poly):
                crossings = 0
                n = len(poly['vertices'])
                for i in range(n):
                    a = poly['vertices'][i]; b = poly['vertices'][(i+1)%n]
                    x0,y0 = a['xo']/a['zo'], a['yo']/a['zo']
                    x1,y1 = b['xo']/b['zo'], b['yo']/b['zo']
                    cond = ((y0 <= py) and (y1 > py)) or ((y0 > py) and (y1 <= py))
                    if cond:
                        xint = x0 + (py - y0) * (x1 - x0) / (y1 - y0)
                        if xint > px: crossings += 1
                return (crossings & 1)
            def depth_at_point(poly, px, py):
                n = len(poly['vertices'])
                for k in range(1, n-1):
                    a = poly['vertices'][0]; b = poly['vertices'][k]; c = poly['vertices'][k+1]
                    Ax,Ay,Az = a['xo']/a['zo'], a['yo']/a['zo'], a['zo']
                    Bx,By,Bz = b['xo']/b['zo'], b['yo']/b['zo'], b['zo']
                    Cx,Cy,Cz = c['xo']/c['zo'], c['yo']/c['zo'], c['zo']
                    # point-in-triangle
                    o1 = (Bx-Ax)*(py-Ay) - (By-Ay)*(px-Ax)
                    o2 = (Cx-Bx)*(py-By) - (Cy-By)*(px-Bx)
                    o3 = (Ax-Cx)*(py-Cy) - (Ay-Cy)*(px-Cx)
                    if not ((o1>=0 and o2>=0 and o3>=0) or (o1<=0 and o2<=0 and o3<=0)):
                        continue
                    denom = (By - Cy)*(Ax - Cx) + (Cx - Bx)*(Ay - Cy)
                    if abs(denom) < 1e-9: return None
                    wA = ((By - Cy)*(px - Cx) + (Cx - Bx)*(py - Cy)) / denom
                    wB = ((Cy - Ay)*(px - Cx) + (Ax - Cx)*(py - Cy)) / denom
                    wC = 1.0 - wA - wB
                    if wA < -1e-5 or wB < -1e-5 or wC < -1e-5: return None
                    return wA*Az + wB*Bz + wC*Cz
                return None
            skeep = sswap = sund = svalid = 0
            sumd0 = sumd1 = 0.0
            for ix in (-1,0,1):
                for iy in (-1,0,1):
                    px = cx + ix * (dx * 0.25)
                    py = cy + iy * (dy * 0.25)
                    if not point_in_poly_pt(px,py,poly0):
                        continue
                    if not point_in_poly_pt(px,py,poly1):
                        continue
                    d0 = depth_at_point(poly0, px, py); d1 = depth_at_point(poly1, px, py)
                    if d0 is None or d1 is None: continue
                    svalid += 1; sumd0 += d0; sumd1 += d1
                    epsd = 1e-4
                    if d0 < d1 - epsd: skeep += 1
                    elif d1 < d0 - epsd: sswap += 1
                    else: sund += 1
            if svalid > 0:
                mean0 = sumd0 / svalid; mean1 = sumd1 / svalid
                if sswap > skeep:
                    return (1, (skeep, sswap, sund, svalid, mean0, mean1))
                if skeep >= sswap:
                    return (-1, (skeep, sswap, sund, svalid, mean0, mean1))
            return (0, None)

        if dec3 == 0 or tri_votes_first3['und'] > 0:
            loc = local_occlusion(tri_faces, tris_f0, tris_f1)
            if loc[0] != 0:
                skeep, sswap, sund, svalid, mean0, mean1 = loc[1]
                df.write('TRI_LOCAL_OCCLUSION_FIRST3,keep={},swap={},und={},samples={},mean_d0={:.6f},mean_d1={:.6f},decision={}\n'.format(skeep, sswap, sund, svalid, mean0, mean1, loc[0]))
                dec3 = loc[0]
            else:
                sk, ss, su = sample_overlap(tri_faces, tris_f0, tris_f1)
                df.write('TRI_OVERLAP_VOTE_FIRST3,keep={},swap={},und={}\n'.format(sk, ss, su))
                if ss > sk: dec3 = 1
                elif sk >= ss: dec3 = -1
        if decN == 0 or tri_votes_newell['und'] > 0:
            loc = local_occlusion(tri_faces, tris_f0, tris_f1)
            if loc[0] != 0:
                skeep, sswap, sund, svalid, mean0, mean1 = loc[1]
                df.write('TRI_LOCAL_OCCLUSION_NEWELL,keep={},swap={},und={},samples={},mean_d0={:.6f},mean_d1={:.6f},decision={}\n'.format(skeep, sswap, sund, svalid, mean0, mean1, loc[0]))
                decN = loc[0]
            else:
                sk, ss, su = sample_overlap(tri_faces, tris_f0, tris_f1)
                df.write('TRI_OVERLAP_VOTE_NEWELL,keep={},swap={},und={}\n'.format(sk, ss, su))
                if ss > sk: decN = 1
                elif sk >= ss: decN = -1

        df.write('eval,face_src,vertex_idx,tv3,xo,yo,zo,tvN,behind_camera\n')
        # face0 plane on face1 vertices
        a3,b3,c3,d3 = faces[0]['plane_first3']; aN,bN,cN,dN = faces[0]['plane_newell']
        for v in faces[1]['vertices']:
            tv3 = a3*v['xo'] + b3*v['yo'] + c3*v['zo'] + d3
            tvN = aN*v['xo'] + bN*v['yo'] + cN*v['zo'] + dN
            df.write(f'eval,0,{v["idx"]},{tv3:.6f},{v["xo"]:.6f},{v["yo"]:.6f},{v["zo"]:.6f},{tvN:.6f},{1 if v["zo"]<=0 else 0}\n')
        # face1 plane on face0 vertices
        a3,b3,c3,d3 = faces[1]['plane_first3']; aN,bN,cN,dN = faces[1]['plane_newell']
        for v in faces[0]['vertices']:
            tv3 = a3*v['xo'] + b3*v['yo'] + c3*v['zo'] + d3
            tvN = aN*v['xo'] + bN*v['yo'] + cN*v['zo'] + dN
            df.write(f'eval,1,{v["idx"]},{tv3:.6f},{v["xo"]:.6f},{v["yo"]:.6f},{v["zo"]:.6f},{tvN:.6f},{1 if v["zo"]<=0 else 0}\n')

        # normals from Newell (unscaled) and comparison dot
        nx0,ny0,nz0 = faces[0]['plane_newell'][0:3]
        nx1,ny1,nz1 = faces[1]['plane_newell'][0:3]
        # compute cross normals from first 3 vertices as additional check
        def tri_normal(face):
            i0,i1,i2 = face['indices'][:3]
            v0 = v_obs[i0]; v1 = v_obs[i1]; v2 = v_obs[i2]
            ux,uy,uz = v1['xo']-v0['xo'], v1['yo']-v0['yo'], v1['zo']-v0['zo']
            vx,vy,vz = v2['xo']-v0['xo'], v2['yo']-v0['yo'], v2['zo']-v0['zo']
            nx = uy*vz - uz*vy
            ny = uz*vx - ux*vz
            nz = ux*vy - uy*vx
            return nx,ny,nz
        tn0 = tri_normal(faces[0]); tn1 = tri_normal(faces[1])
        tri_dot0 = nx0*tn0[0] + ny0*tn0[1] + nz0*tn0[2]
        tri_dot1 = nx1*tn1[0] + ny1*tn1[1] + nz1*tn1[2]
        df.write('normal,0,{:.6f},{:.6f},{:.6f},tri_dot={:.6f}\n'.format(nx0,ny0,nz0,tri_dot0))
        df.write('normal,1,{:.6f},{:.6f},{:.6f},tri_dot={:.6f}\n'.format(nx1,ny1,nz1,tri_dot1))

        # planarity checks
        def planarity(face, plane):
            a,b,c,d = plane
            denom = math.sqrt(a*a + b*b + c*c)
            if denom == 0: return {'max': 0.0, 'avg': 0.0}
            vals = []
            for v in face['vertices']:
                tv = (a*v['xo'] + b*v['yo'] + c*v['zo'] + d) / denom
                vals.append(abs(tv))
            return {'max': max(vals), 'avg': sum(vals)/len(vals)}
        def seg_intersect(p1,p2,p3,p4):
            # 2D segment intersection test (proper or touching)
            def orient(a,b,c):
                return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
            def on_seg(a,b,c):
                return (min(a[0],b[0]) <= c[0] <= max(a[0],b[0]) and min(a[1],b[1]) <= c[1] <= max(a[1],b[1]))
            o1 = orient(p1,p2,p3); o2 = orient(p1,p2,p4); o3 = orient(p3,p4,p1); o4 = orient(p3,p4,p2)
            if o1==0 and on_seg(p1,p2,p3): return True
            if o2==0 and on_seg(p1,p2,p4): return True
            if o3==0 and on_seg(p3,p4,p1): return True
            if o4==0 and on_seg(p3,p4,p2): return True
            return (o1>0) != (o2>0) and (o3>0) != (o4>0)

        for fi in range(2):
            p3 = planarity(faces[fi], faces[fi]['plane_first3'])
            pN = planarity(faces[fi], faces[fi]['plane_newell'])
            # model-space planarity using original 3D coordinates (swapped Y/Z to match in-app convention)
            poly_model = []
            for vid in faces[fi]['indices']:
                x,y,z = verts[vid]
                poly_model.append((x,y,z))
            pm = planarity({'vertices': [{'xo':p[0],'yo':p[1],'zo':p[2]} for p in poly_model]}, plane_newell(poly_model))
            # check 2D self-intersection: use model XY (x,y) and observer projection (xo/zo -> x2d, yo/zo -> y2d)
            # model-space simple check
            pts_model = [(p[0], p[1]) for p in poly_model]
            simple_model = 1
            n = len(pts_model)
            for a in range(n):
                for b in range(a+1,n):
                    # skip adjacent edges
                    if abs(a-b) <= 1 or (a==0 and b==n-1) or (b==0 and a==n-1):
                        continue
                    if seg_intersect(pts_model[a], pts_model[(a+1)%n], pts_model[b], pts_model[(b+1)%n]):
                        simple_model = 0
            # observer projected simple check
            pts_obs = []
            for vid in faces[fi]['indices']:
                v = v_obs[vid]
                x2 = v['xo']/v['zo'] if v['zo']>0 else -99999.0
                y2 = v['yo']/v['zo'] if v['zo']>0 else -99999.0
                pts_obs.append((x2,y2))
            simple_obs = 1
            for a in range(n):
                for b in range(a+1,n):
                    if abs(a-b) <= 1 or (a==0 and b==n-1) or (b==0 and a==n-1):
                        continue
                    if seg_intersect(pts_obs[a], pts_obs[(a+1)%n], pts_obs[b], pts_obs[(b+1)%n]):
                        simple_obs = 0
            df.write('planarity,{},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f},simple_model={},simple_obs={}\n'.format(fi, p3['max'], p3['avg'], pN['max'], pN['avg'], pm['max'], pm['avg'], simple_model, simple_obs))

        # run pair_should_swap using both methods in OBSERVER SPACE
        res3 = pair_should_swap_py(faces, v_obs, 0, 1, use_newell=False)
        resN = pair_should_swap_py(faces, v_obs, 0, 1, use_newell=True)
        df.write(f'pair_should_swap_obs_first3,{res3}\n')
        df.write(f'pair_should_swap_obs_newell,{resN}\n')

        # Now compute planes and tests in MODEL SPACE and evaluate at observer position in model coordinates
        # Observer position in model space: p_obs = distance * (cos_h*cos_v, sin_h*cos_v, sin_v)
        ch = math.cos(deg2rad(args.angle_h)); sh = math.sin(deg2rad(args.angle_h))
        cv = math.cos(deg2rad(args.angle_v)); sv = math.sin(deg2rad(args.angle_v))
        dist = args.distance
        obs_model = (dist * ch * cv, dist * sh * cv, dist * sv)
        df.write('observer_model,{:.6f},{:.6f},{:.6f}\n'.format(obs_model[0], obs_model[1], obs_model[2]))

        # Build face structures in model space
        faces_model = []
        for fi, fidx in enumerate(faces_idx):
            face = {'idx': fi, 'indices': fidx}
            verts_model = []
            for vi in fidx:
                # using parse_obj coords (already matching model-space orientation used by C)
                x,y,z = verts[vi]
                verts_model.append({'x': x, 'y': y, 'z': z, 'idx': vi})
            face['vertices'] = verts_model
            # planes in model space
            # first-3 method
            idx0 = fidx[0]; idx1 = fidx[1]; idx2 = fidx[2]
            v0 = (verts[idx0][0], verts[idx0][1], verts[idx0][2])
            v1 = (verts[idx1][0], verts[idx1][1], verts[idx1][2])
            v2 = (verts[idx2][0], verts[idx2][1], verts[idx2][2])
            face['plane_first3'] = plane_from_first3(v0, v1, v2)
            face['plane_newell'] = plane_newell([(vv['x'], vv['y'], vv['z']) for vv in verts_model])
            faces_model.append(face)

        # helper: evaluate tv in model-space plane
        def eval_tv_model(plane, pt):
            a,b,c,d = plane
            x,y,z = pt
            return a*x + b*y + c*z + d

        # write model-space per-vertex evals and observer evals
        df.write('model_eval,face_src,vertex_idx,tv3_model,tvN_model,x,y,z\n')
        for v in faces_model[1]['vertices']:
            pt = (v['x'], v['y'], v['z'])
            tv3 = eval_tv_model(faces_model[0]['plane_first3'], pt)
            tvN = eval_tv_model(faces_model[0]['plane_newell'], pt)
            df.write('model_eval,0,{},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}\n'.format(v['idx'], tv3, tvN, pt[0], pt[1], pt[2]))
        for v in faces_model[0]['vertices']:
            pt = (v['x'], v['y'], v['z'])
            tv3 = eval_tv_model(faces_model[1]['plane_first3'], pt)
            tvN = eval_tv_model(faces_model[1]['plane_newell'], pt)
            df.write('model_eval,1,{},{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}\n'.format(v['idx'], tv3, tvN, pt[0], pt[1], pt[2]))

        obs_tv_f0_3 = eval_tv_model(faces_model[0]['plane_first3'], obs_model)
        obs_tv_f0_N = eval_tv_model(faces_model[0]['plane_newell'], obs_model)
        obs_tv_f1_3 = eval_tv_model(faces_model[1]['plane_first3'], obs_model)
        obs_tv_f1_N = eval_tv_model(faces_model[1]['plane_newell'], obs_model)
        df.write('observer_eval,face0,tv3_model={:.6f},tvN_model={:.6f}\n'.format(obs_tv_f0_3, obs_tv_f0_N))
        df.write('observer_eval,face1,tv3_model={:.6f},tvN_model={:.6f}\n'.format(obs_tv_f1_3, obs_tv_f1_N))

        # run pair tests in model space
        # reuse pair_should_swap_py but need a wrapper that accepts face structures in model space
        def pair_model_wrapper(faces_m, obs_pt, f1i, f2i, use_newell=False):
            # build faces-like structure with vertices containing xo/yo/zo replaced by x,y,z
            tmp_faces = []
            for fm in faces_m:
                tmp = {'z_min':None, 'z_max':None, 'z_mean':None, 'minx':None, 'maxx':None, 'miny':None, 'maxy':None, 'vertex_count': len(fm['vertices']), 'vertex_indices_buffer': [v['idx'] for v in fm['vertices']], 'vertex_indices_ptr':[0], 'plane_a':0, 'plane_b':0, 'plane_c':0, 'plane_d':0}
                # store vertices as dicts with xo/yo/zo fields
                tmp_v = []
                for v in fm['vertices']:
                    tmp_v.append({'xo': v['x'], 'yo': v['y'], 'zo': v['z'], 'idx': v['idx']})
                tmp['vertices'] = tmp_v
                # attach plane used
                if use_newell:
                    tmp['plane_a'], tmp['plane_b'], tmp['plane_c'], tmp['plane_d'] = fm['plane_newell']
                else:
                    tmp['plane_a'], tmp['plane_b'], tmp['plane_c'], tmp['plane_d'] = fm['plane_first3']
                # compute bbox on projected 2D using x,y coords (model space projection isn't meaningful, but we only need bbox relative)
                xs = [v['xo'] for v in tmp_v]
                ys = [v['yo'] for v in tmp_v]
                tmp['minx'], tmp['maxx'], tmp['miny'], tmp['maxy'] = min(xs), max(xs), min(ys), max(ys)
                tmp_faces.append(tmp)
            # observer point as vertex with 'xo','yo','zo' fields
            # reuse pair_should_swap_py but need faces as list and vtx mapping; we pass vtx as tmp_faces via face indices
            # hack: create a fake vtx list mapping indices to xo/yo/zo from faces
            # Instead, we'll implement a simplified test here using same logic as pair_should_swap_py but adapted
            return pair_should_swap_py_adapted(tmp_faces, obs_pt, f1i, f2i)

        # simplified adapted tests for model wrapper
        def pair_should_swap_py_adapted(faces_tmp, obs_pt, f1i, f2i):
            # similar to pair_should_swap_py but using faces_tmp where vertices are in .vertices with xo/yo/zo
            f1 = faces_tmp[f1i]; f2 = faces_tmp[f2i]
            # depth test using z values (model-space z) via vertices' zo
            zvals_f1 = [v['zo'] for v in f1['vertices']]
            zvals_f2 = [v['zo'] for v in f2['vertices']]
            if max(zvals_f2) <= min(zvals_f1): return -1
            if max(zvals_f1) <= min(zvals_f2): return 1
            # bbox quick reject
            if f1['maxx'] <= f2['minx'] or f2['maxx'] <= f1['minx']: return -1
            if f1['maxy'] <= f2['miny'] or f2['maxy'] <= f1['miny']: return -1
            # load planes
            a1,b1,c1,d1 = f1['plane_a'], f1['plane_b'], f1['plane_c'], f1['plane_d']
            a2,b2,c2,d2 = f2['plane_a'], f2['plane_b'], f2['plane_c'], f2['plane_d']
            n1 = f1['vertex_count']; n2 = f2['vertex_count']
            # Test 4 mimic
            maxabs1 = abs(d1)
            for v in f2['vertices']:
                tv = a1*v['xo'] + b1*v['yo'] + c1*v['zo'] + d1
                if abs(tv) > maxabs1: maxabs1 = abs(tv)
            eps_rel1 = max(1e-4, maxabs1 * 1e-6)
            obs1 = 0
            if d1 > eps_rel1: obs1 = 1
            elif d1 < -eps_rel1: obs1 = -1
            if obs1 != 0:
                pos=neg=0
                for v in f2['vertices']:
                    tv = a1*v['xo'] + b1*v['yo'] + c1*v['zo'] + d1
                    if abs(tv) <= eps_rel1: pass
                    elif tv > 0: pos+=1
                    else: neg+=1
                thr = (3*n2 + 3)//4
                if (obs1==1 and pos>=thr) or (obs1==-1 and neg>=thr): return -1
                if (obs1==1 and neg>=thr) or (obs1==-1 and pos>=thr): return 1
            # symmetric
            maxabs2 = abs(d2)
            for v in f1['vertices']:
                tv = a2*v['xo'] + b2*v['yo'] + c2*v['zo'] + d2
                if abs(tv) > maxabs2: maxabs2 = abs(tv)
            eps_rel2 = max(1e-4, maxabs2 * 1e-6)
            obs2 = 0
            if d2 > eps_rel2: obs2 = 1
            elif d2 < -eps_rel2: obs2 = -1
            if obs2 != 0:
                pos2=neg2=0
                for v in f1['vertices']:
                    tv = a2*v['xo'] + b2*v['yo'] + c2*v['zo'] + d2
                    if abs(tv) <= eps_rel2: pass
                    elif tv > 0: pos2+=1
                    else: neg2+=1
                thr2 = (3*n1 + 3)//4
                if (obs2==1 and neg2>=thr2) or (obs2==-1 and pos2>=thr2): return -1
                if (obs2==1 and pos2>=thr2) or (obs2==-1 and neg2>=thr2): return 1
            return 0

        # compute model-space pair results
        pm_res3 = pair_model_wrapper(faces_model, obs_model, 0, 1, use_newell=False)
        pm_resN = pair_model_wrapper(faces_model, obs_model, 0, 1, use_newell=True)
        df.write(f'pair_should_swap_model_first3,{pm_res3}\n')
        df.write(f'pair_should_swap_model_newell,{pm_resN}\n')


    if args.verbose:
        print(f'Wrote {outp} with comparison (first3 vs Newell).')


if __name__ == '__main__':
    main()
