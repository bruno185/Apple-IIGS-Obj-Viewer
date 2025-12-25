#!/usr/bin/env python3
"""Analyze PAIRDIAG.CSV and compare C decisions with Python replay of pair_should_swap."""
import math
import csv
import sys

def deg2rad(d): return d*math.pi/180.0

def transform_vertex(v, angle_h, angle_v, distance):
    x,y,z = v
    ch = math.cos(deg2rad(angle_h)); sh = math.sin(deg2rad(angle_h))
    cv = math.cos(deg2rad(angle_v)); sv = math.sin(deg2rad(angle_v))
    term1 = x * ch * cv
    term2 = y * sh * cv
    term3 = z * sv
    zo = -term1 - term2 - term3 + distance
    if zo > 0:
        xo = - x * sh + y * ch
        yo = - x * ch * sv - y * sh * sv + z * cv
    else:
        xo = 0.0; yo = 0.0
    return xo, yo, zo

# minimal pair_should_swap_py (use Newell planes or first3 based on use_newell param)
def pair_should_swap_py_face_structs(f1, f2, use_newell=True):
    # f1,f2: dict with 'vertices' list of model-space vertex tuples (x,y,z)
    # compute observer-space vertices first using provided observer in outer scope
    def plane_from_first3(v0,v1,v2):
        x1,y1,z1=v0; x2,y2,z2=v1; x3,y3,z3=v2
        a = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
        b = -x1*(z2 - z3) + x2*(z1 - z3) - x3*(z1 - z2)
        c = x1*(y2 - y3) - x2*(y1 - y3) + x3*(y1 - y2)
        d = -x1*(y2*z3 - y3*z2) + x2*(y1*z3 - y3*z1) - x3*(y1*z2 - y2*z1)
        return a,b,c,d
    def plane_newell(poly):
        nx=ny=nz=0.0
        n=len(poly)
        for i in range(n):
            x1,y1,z1 = poly[i]
            x2,y2,z2 = poly[(i+1)%n]
            nx += (y1 - y2) * (z1 + z2)
            ny += (z1 - z2) * (x1 + x2)
            nz += (x1 - x2) * (y1 + y2)
        d = -(nx * poly[0][0] + ny * poly[0][1] + nz * poly[0][2])
        return nx,ny,nz,d

    # build observer-space vertices
    vobs1 = [transform_vertex(v, global_angle_h, global_angle_v, global_distance) for v in f1['vertices']]
    vobs2 = [transform_vertex(v, global_angle_h, global_angle_v, global_distance) for v in f2['vertices']]
    if use_newell:
        a1,b1,c1,d1 = plane_newell(vobs1)
        a2,b2,c2,d2 = plane_newell(vobs2)
    else:
        a1,b1,c1,d1 = plane_from_first3(vobs1[0], vobs1[1], vobs1[2])
        a2,b2,c2,d2 = plane_from_first3(vobs2[0], vobs2[1], vobs2[2])
    # run simplified tests similar to C version (observer-space only for speed)
    # Test1: depth (use zo values)
    zmin1 = min([v[2] for v in vobs1]); zmax1 = max([v[2] for v in vobs1])
    zmin2 = min([v[2] for v in vobs2]); zmax2 = max([v[2] for v in vobs2])
    if zmax2 <= zmin1: return -1
    if zmax1 <= zmin2: return 1
    # bbox 2D
    xs1 = [v[0]/v[2] if v[2]>0 else -99999 for v in vobs1]
    ys1 = [v[1]/v[2] if v[2]>0 else -99999 for v in vobs1]
    xs2 = [v[0]/v[2] if v[2]>0 else -99999 for v in vobs2]
    ys2 = [v[1]/v[2] if v[2]>0 else -99999 for v in vobs2]
    if max(xs1) <= min(xs2) or max(xs2) <= min(xs1): return -1
    if max(ys1) <= min(ys2) or max(ys2) <= min(ys1): return -1
    # plane tests (observer-space)
    maxabs1 = abs(d1)
    for v in vobs2:
        tv = a1*v[0] + b1*v[1] + c1*v[2] + d1
        if abs(tv) > maxabs1: maxabs1 = abs(tv)
    eps1 = max(1e-4, maxabs1 * 1e-6)
    obs1 = 1 if d1 > eps1 else (-1 if d1 < -eps1 else 0)
    if obs1 != 0:
        pos=neg=0
        for v in vobs2:
            tv = a1*v[0] + b1*v[1] + c1*v[2] + d1
            if abs(tv) <= eps1: pass
            elif tv > 0: pos +=1
            else: neg +=1
        thr = (3*len(vobs2)+3)//4
        if (obs1==1 and pos>=thr) or (obs1==-1 and neg>=thr): return -1
        if (obs1==1 and neg>=thr) or (obs1==-1 and pos>=thr): return 1
    # symmetric
    maxabs2 = abs(d2)
    for v in vobs1:
        tv = a2*v[0] + b2*v[1] + c2*v[2] + d2
        if abs(tv) > maxabs2: maxabs2 = abs(tv)
    eps2 = max(1e-4, maxabs2 * 1e-6)
    obs2 = 1 if d2 > eps2 else (-1 if d2 < -eps2 else 0)
    if obs2 != 0:
        pos2=neg2=0
        for v in vobs1:
            tv = a2*v[0] + b2*v[1] + c2*v[2] + d2
            if abs(tv) <= eps2: pass
            elif tv > 0: pos2 +=1
            else: neg2 +=1
        thr2 = (3*len(vobs1)+3)//4
        if (obs2==1 and neg2>=thr2) or (obs2==-1 and pos2>=thr2): return -1
        if (obs2==1 and pos2>=thr2) or (obs2==-1 and neg2>=thr2): return 1
    # triangulation fallback (observer-space) - simplified votes
    tri_keep = tri_swap = tri_und = 0
    if len(vobs1)>=3 and len(vobs2)>=3:
        t1_count = len(vobs1)-2; t2_count = len(vobs2)-2
        for a in range(t1_count):
            i0,i1,i2 = 0,a+1,a+2
            for b in range(t2_count):
                j0,j1,j2 = 0,b+1,b+2
                # compute triangle plane using observer-space coords and do simplified test
                # implement simple heuristic: if majority of tri2 verts are behind tri1 plane -> keep
                # This is enough for reproduction and comparison
                # compute normal by cross product (observer-space)
                ux = vobs1[i1][0]-vobs1[i0][0]; uy=vobs1[i1][1]-vobs1[i0][1]; uz=vobs1[i1][2]-vobs1[i0][2]
                vx = vobs1[i2][0]-vobs1[i0][0]; vy=vobs1[i2][1]-vobs1[i0][1]; vz=vobs1[i2][2]-vobs1[i0][2]
                ta = uy*vz - uz*vy
                tb = uz*vx - ux*vz
                tc = ux*vy - uy*vx
                td = -(ta*vobs1[i0][0] + tb*vobs1[i0][1] + tc*vobs1[i0][2])
                pos=neg=0
                for jj in (j0,j1,j2):
                    tv = ta*vobs2[jj][0] + tb*vobs2[jj][1] + tc*vobs2[jj][2] + td
                    if tv > 0: pos+=1
                    elif tv < 0: neg+=1
                if pos>=2: tri_keep+=1
                elif neg>=2: tri_swap+=1
                else: tri_und+=1
    if tri_swap > tri_keep: return 1
    if tri_keep >= tri_swap: return -1
    # fallback zmean
    zm1 = sum([v[2] for v in vobs1])/len(vobs1)
    zm2 = sum([v[2] for v in vobs2])/len(vobs2)
    if zm1 > zm2: return -1
    if zm1 < zm2: return 1
    return -1

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: analyze_pair_logs.py PAIRDIAG.CSV')
        sys.exit(1)
    infile = sys.argv[1]
    # parse
    pairs = []
    verts = { }  # model-space verts per face
    global_angle_h = None; global_angle_v = None; global_distance = None
    with open(infile, 'r') as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith('NEWELL_V,'):
                # NEWELL_V,face,vi,x,y,z
                parts=line.split(',')
                fi=int(parts[1]); vi=int(parts[2]); x=float(parts[3]); y=float(parts[4]); z=float(parts[5])
                verts.setdefault(fi,[]).append((x,y,z))
            if line.startswith('PAIR_DETAIL,'):
                parts=line.split(',')
                # fields: PAIR_DETAIL,angle_h,angle_v,distance,f1,f2,a1,b1,c1,d1,a2,b2,c2,d2,zmean1,zmean2,pos1,neg1,pos2,neg2,tri_keep,tri_swap,tri_und,decision
                try:
                    angle_h=float(parts[1]); angle_v=float(parts[2]); distance=float(parts[3])
                    f1=int(parts[4]); f2=int(parts[5])
                    a1=float(parts[6]); b1=float(parts[7]); c1=float(parts[8]); d1=float(parts[9])
                    a2=float(parts[10]); b2=float(parts[11]); c2=float(parts[12]); d2=float(parts[13])
                    zmean1=float(parts[14]); zmean2=float(parts[15])
                    tri_keep=int(parts[20]); tri_swap=int(parts[21]); tri_und=int(parts[22]); decision=int(parts[23])
                except Exception as e:
                    print('Parse error', e, line)
                    continue
                pairs.append({'angle_h':angle_h,'angle_v':angle_v,'distance':distance,'f1':f1,'f2':f2,'a1':a1,'b1':b1,'c1':c1,'d1':d1,'a2':a2,'b2':b2,'c2':c2,'d2':d2,'z1':zmean1,'z2':zmean2,'tri_keep':tri_keep,'tri_swap':tri_swap,'tri_und':tri_und,'decision':decision})
            if line.startswith('observer_model,'):
                parts=line.split(',')
                global_distance = float(parts[1])
    # analyze
    mismatches = []
    for p in pairs:
        global_angle_h = p['angle_h']; global_angle_v = p['angle_v']; global_distance = p['distance']
        # build face dicts from verts
        f1 = {'vertices': verts.get(p['f1'], [])}
        f2 = {'vertices': verts.get(p['f2'], [])}
        py_dec = pair_should_swap_py_face_structs(f1,f2,use_newell=True)
        if py_dec != p['decision']:
            mismatches.append((p,py_dec))
    print('Total pairs:', len(pairs), 'mismatches:', len(mismatches))
    for pm,py in mismatches[:50]:
        print('Angle H,V,D:', pm['angle_h'], pm['angle_v'], pm['distance'], 'faces', pm['f1'], pm['f2'], 'C_dec', pm['decision'], 'PY_dec', py, 'tri_keep', pm['tri_keep'], 'tri_swap', pm['tri_swap'])
    if mismatches:
        sys.exit(2)
    else:
        print('No mismatches found')
