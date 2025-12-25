#!/usr/bin/env python3
"""Pascal simulation diagnostic script
Reproduces Fixed32 16.16 arithmetic and tests 1..7 used by painter_newell_sancha.
Usage: python tools/pascal_sim.py cube.obj H V W distance
"""
import sys
import math

FIXED_SHIFT = 16
FIXED_SCALE = 1 << FIXED_SHIFT

# Fixed32 helpers
def to_fixed(f):
    return int(round(f * FIXED_SCALE))

def to_float(x):
    return float(x) / FIXED_SCALE

def fixed_mul_64(a, b):
    return int((int(a) * int(b)) >> FIXED_SHIFT)

def fixed_add(a, b):
    return int(a + b)

def fixed_sub(a, b):
    return int(a - b)

def fixed_neg(a):
    return -a

# Parse a simple OBJ

def parse_obj(path):
    verts = []
    faces = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            if parts[0] == 'v':
                x = float(parts[1]); y = float(parts[2]); z = float(parts[3])
                verts.append((x,y,z))
            elif parts[0] == 'f':
                # faces may be like "f 1 2 3 4" or with slashes
                idxs = []
                for p in parts[1:]:
                    if '/' in p:
                        p = p.split('/')[0]
                    idxs.append(int(p)-1)
                faces.append(idxs)
    return verts, faces

# Trig in fixed (16.16)
def fixed_sin_deg(deg):
    r = math.radians(deg)
    return to_fixed(math.sin(r))

def fixed_cos_deg(deg):
    r = math.radians(deg)
    return to_fixed(math.cos(r))

# Reproduce processModelFast transform and projection

def transform_and_project(verts, H, V, W, distance):
    n = len(verts)
    xo = [0]*n; yo=[0]*n; zo=[0]*n; x2d=[0]*n; y2d=[0]*n
    cos_h = fixed_cos_deg(H); sin_h = fixed_sin_deg(H)
    cos_v = fixed_cos_deg(V); sin_v = fixed_sin_deg(V)
    cos_w = fixed_cos_deg(W); sin_w = fixed_sin_deg(W)

    cos_h_cos_v = fixed_mul_64(cos_h, cos_v)
    sin_h_cos_v = fixed_mul_64(sin_h, cos_v)
    cos_h_sin_v = fixed_mul_64(cos_h, sin_v)
    sin_h_sin_v = fixed_mul_64(sin_h, sin_v)

    scale = to_fixed(100)  # projection scale
    centre_x_f = to_fixed(160)
    centre_y_f = to_fixed(100)
    # follow code: distance = params.distance ; distance = distance * 4
    d = to_fixed(distance)
    d = fixed_mul_64(d, to_fixed(4))

    for i,(xf,yf,zf) in enumerate(verts):
        x = to_fixed(xf); y = to_fixed(yf); z = to_fixed(zf)
        term1 = fixed_mul_64(x, cos_h_cos_v)
        term2 = fixed_mul_64(y, sin_h_cos_v)
        term3 = fixed_mul_64(z, sin_v)
        zo_i = fixed_add(fixed_sub(fixed_sub(fixed_neg(term1), term2), term3), d)
        zo[i] = zo_i
        if zo_i != 0:
            xo_i = fixed_add(fixed_neg(fixed_mul_64(x, sin_h)), fixed_mul_64(y, cos_h))
            yo_i = fixed_add(fixed_sub(fixed_neg(fixed_mul_64(x, cos_h_sin_v)), fixed_mul_64(y, sin_h_sin_v)), fixed_mul_64(z, cos_v))
            xo[i] = xo_i
            yo[i] = yo_i
            inv_zo = int((int(scale) << FIXED_SHIFT) // zo_i)  # scale / zo, in fixed
            x2d_temp = fixed_add(fixed_mul_64(xo_i, inv_zo), centre_x_f)
            y2d_temp = fixed_sub(centre_y_f, fixed_mul_64(yo_i, inv_zo))
            # rotation W (zero often)
            x2d_i = fixed_add(fixed_sub(fixed_mul_64(cos_w, fixed_sub(x2d_temp, centre_x_f)), fixed_mul_64(sin_w, fixed_sub(centre_y_f, y2d_temp))), centre_x_f)
            y2d_i = fixed_sub(centre_y_f, fixed_add(fixed_mul_64(sin_w, fixed_sub(x2d_temp, centre_x_f)), fixed_mul_64(cos_w, fixed_sub(centre_y_f, y2d_temp))))
            x2d[i] = int(round(to_float(x2d_i)))
            y2d[i] = int(round(to_float(y2d_i)))
        else:
            xo[i]=yo[i]=0
            x2d[i]=y2d[i]=0
    return xo, yo, zo, x2d, y2d

# Compute face planes and depths

def compute_faces(verts, faces, xo, yo, zo, x2d, y2d):
    face_info = []
    for fi, idxs in enumerate(faces):
        n = len(idxs)
        display_flag = 1
        zmin = None; zmax=None; sumz=0
        minx=miny=99999; maxx=maxy=-99999
        for idx in idxs:
            zz = zo[idx]
            if zz < 0: display_flag = 0
            if zmin is None or zz < zmin: zmin = zz
            if zmax is None or zz > zmax: zmax = zz
            sumz += zz
            xx2 = x2d[idx]; yy2 = y2d[idx]
            if xx2 < minx: minx = xx2
            if xx2 > maxx: maxx = xx2
            if yy2 < miny: miny = yy2
            if yy2 > maxy: maxy = yy2
        if not display_flag or n < 3:
            a=b=c=d=0
        else:
            idx0=idxs[0]; idx1=idxs[1]; idx2=idxs[2]
            x1,y1,z1 = xo[idx0], yo[idx0], zo[idx0]
            x2,y2,z2 = xo[idx1], yo[idx1], zo[idx1]
            x3,y3,z3 = xo[idx2], yo[idx2], zo[idx2]
            a = fixed_add(fixed_add(fixed_mul_64(y1, fixed_sub(z2,z3)), fixed_mul_64(y2, fixed_sub(z3,z1))), fixed_mul_64(y3, fixed_sub(z1,z2)))
            b = fixed_sub(fixed_add(fixed_neg(fixed_mul_64(x1, fixed_sub(z2,z3))), fixed_mul_64(x2, fixed_sub(z1,z3))), fixed_mul_64(x3, fixed_sub(z1,z2)))
            c = fixed_add(fixed_sub(fixed_mul_64(x1, fixed_sub(y2,y3)), fixed_mul_64(x2, fixed_sub(y1,y3))), fixed_mul_64(x3, fixed_sub(y1,y2)))
            t1 = fixed_sub(fixed_mul_64(y2,z3), fixed_mul_64(y3,z2))
            t2 = fixed_sub(fixed_mul_64(y1,z3), fixed_mul_64(y3,z1))
            t3 = fixed_sub(fixed_mul_64(y1,z2), fixed_mul_64(y2,z1))
            d = fixed_add(fixed_add(fixed_neg(fixed_mul_64(x1,t1)), fixed_mul_64(x2,t2)), fixed_neg(fixed_mul_64(x3,t3)))
        zmean = int(sumz // n) if n>0 else 0
        face_info.append({'a':a,'b':b,'c':c,'d':d,'zmin':zmin,'zmax':zmax,'zmean':zmean,'display':display_flag,'minx':minx,'maxx':maxx,'miny':miny,'maxy':maxy,'idxs':idxs})
    return face_info

# Tests 1..7 per painter code

def test_pair(faces, f1i, f2i):
    f1=faces[f1i]; f2=faces[f2i]
    out = {'pair':(f1i,f2i)}
    # Test1 depth overlap
    if f2['zmax'] <= f1['zmin']:
        out['result']='ok_order_depth1'
        return out
    if f1['zmax'] <= f2['zmin']:
        out['result']='swap_due_depth'
        return out
    # bbox quick reject
    if f1['maxx'] <= f2['minx'] or f2['maxx'] <= f1['minx']:
        out['result']='ok_no_overlap_x'
        return out
    if f1['maxy'] <= f2['miny'] or f2['maxy'] <= f1['miny']:
        out['result']='ok_no_overlap_y'
        return out
    # prepare plane coeffs
    a1,b1,c1,d1 = f1['a'], f1['b'], f1['c'], f1['d']
    a2,b2,c2,d2 = f2['a'], f2['b'], f2['c'], f2['d']
    epsilon = to_fixed(0.01)
    # Test4
    obs_side1 = 0
    if d1 > epsilon: obs_side1 = 1
    elif d1 < -epsilon: obs_side1 = -1
    else:
        out['test4']='skip_obs_on_plane'
    if obs_side1 != 0:
        all_same_side = True
        for idx in f2['idxs']:
            tv = a1*xo[idx] + b1*yo[idx] + c1*zo[idx] + d1
            if tv > epsilon: side=1
            elif tv < -epsilon: side=-1
            else: side=0
            if side != obs_side1:
                all_same_side = False; break
        if all_same_side:
            out['result']='ok_test4'
            return out
    # Test5
    obs_side2 = 0
    if d2 > epsilon: obs_side2 = 1
    elif d2 < -epsilon: obs_side2 = -1
    else:
        out['test5']='skip_obs_on_plane'
    if obs_side2 != 0:
        all_opposite = True
        for idx in f1['idxs']:
            tv = a2*xo[idx] + b2*yo[idx] + c2*zo[idx] + d2
            if tv > epsilon: side=1
            elif tv < -epsilon: side=-1
            else: side=0
            if side == obs_side2:
                all_opposite = False; break
        if all_opposite:
            out['result']='ok_test5'
            return out
    # Tests 6 and 7 similarly
    # Test6
    if d1 > epsilon: obs_side1 = 1
    elif d1 < -epsilon: obs_side1 = -1
    else: out['test6']='skip_obs_on_plane'
    if obs_side1 != 0:
        all_opposite = True
        for idx in f2['idxs']:
            tv = a1*xo[idx] + b1*yo[idx] + c1*zo[idx] + d1
            if tv > epsilon: side=1
            else: side=-1
            if side == obs_side1:
                all_opposite = False; break
        if all_opposite:
            out['result']='swap_test6'
            return out
    # Test7
    if d2 > epsilon: obs_side2 = 1
    elif d2 < -epsilon: obs_side2 = -1
    else: out['test7']='skip_obs_on_plane'
    if obs_side2 != 0:
        all_same = True
        for idx in f1['idxs']:
            tv = a2*xo[idx] + b2*yo[idx] + c2*zo[idx] + d2
            if tv > epsilon: side=1
            else: side=-1
            if side != obs_side2:
                all_same = False; break
        if all_same:
            out['result']='swap_test7'
            return out
    out['result']='undetermined'
    return out


def test_pair_rel(faces, f1i, f2i, xo, yo, zo, rel_factor=1e-6, min_eps=0.0001):
    """Relative-epsilon + majority vote variant of the test_pair logic
    xo,yo,zo: per-vertex transformed coordinates (fixed) required to evaluate planes
    rel_factor: fraction of max absolute plane test value used as eps
    min_eps: minimal epsilon in float units
    """
    f1=faces[f1i]; f2=faces[f2i]
    out = {'pair':(f1i,f2i)}
    # Test1 depth overlap
    if f2['zmax'] <= f1['zmin']:
        out['result']='ok_order_depth1'
        return out
    if f1['zmax'] <= f2['zmin']:
        out['result']='swap_due_depth'
        return out
    # bbox quick reject
    if f1['maxx'] <= f2['minx'] or f2['maxx'] <= f1['minx']:
        out['result']='ok_no_overlap_x'
        return out
    if f1['maxy'] <= f2['miny'] or f2['maxy'] <= f1['miny']:
        out['result']='ok_no_overlap_y'
        return out
    # prepare plane coeffs
    a1,b1,c1,d1 = f1['a'], f1['b'], f1['c'], f1['d']
    a2,b2,c2,d2 = f2['a'], f2['b'], f2['c'], f2['d']

    # Test4 relative eps + majority
    # build maxabs from plane f1 values on vertices of f2 and d1
    tvals_f1_on_f2 = [a1*xo[idx] + b1*yo[idx] + c1*zo[idx] + d1 for idx in f2['idxs']]
    maxabs1 = max(abs(x) for x in tvals_f1_on_f2 + [d1])
    eps1 = max(min_eps, abs(maxabs1) * rel_factor)
    # observer side
    obs_side1 = 0
    if d1 > eps1: obs_side1 = 1
    elif d1 < -eps1: obs_side1 = -1
    else:
        out['test4']='skip_obs_near_plane'
    if obs_side1 != 0:
        pos = sum(1 for v in tvals_f1_on_f2 if v > eps1)
        neg = sum(1 for v in tvals_f1_on_f2 if v < -eps1)
        n = len(tvals_f1_on_f2)
        thr = (3*n + 3)//4
        if (obs_side1 == 1 and pos >= thr) or (obs_side1 == -1 and neg >= thr):
            out['result']='ok_test4'
            return out

    # Test5 relative eps + majority
    tvals_f2_on_f1 = [a2*xo[idx] + b2*yo[idx] + c2*zo[idx] + d2 for idx in f1['idxs']]
    maxabs2 = max(abs(x) for x in tvals_f2_on_f1 + [d2])
    eps2 = max(min_eps, abs(maxabs2) * rel_factor)
    obs_side2 = 0
    if d2 > eps2: obs_side2 = 1
    elif d2 < -eps2: obs_side2 = -1
    else:
        out['test5']='skip_obs_near_plane'
    if obs_side2 != 0:
        pos2 = sum(1 for v in tvals_f2_on_f1 if v > eps2)
        neg2 = sum(1 for v in tvals_f2_on_f1 if v < -eps2)
        n2 = len(tvals_f2_on_f1)
        thr2 = (3*n2 + 3)//4
        if (obs_side2 == 1 and neg2 >= thr2) or (obs_side2 == -1 and pos2 >= thr2):
            out['result']='ok_test5'
            return out

    # Test6: f2 opposite of observer wrt plane f1 -> swap
    if obs_side1 == 0:
        # recompute if not already computed
        pos=neg=0
        pos = sum(1 for v in tvals_f1_on_f2 if v > eps1)
        neg = sum(1 for v in tvals_f1_on_f2 if v < -eps1)
        n = len(tvals_f1_on_f2)
        thr = (3*n + 3)//4
    if (obs_side1 == 1 and neg >= thr) or (obs_side1 == -1 and pos >= thr):
        out['result']='swap_test6'
        return out

    # Test7: f1 same side of observer wrt plane f2 -> swap
    if obs_side2 == 0:
        pos2 = sum(1 for v in tvals_f2_on_f1 if v > eps2)
        neg2 = sum(1 for v in tvals_f2_on_f1 if v < -eps2)
        n2 = len(tvals_f2_on_f1)
        thr2 = (3*n2 + 3)//4
    if (obs_side2 == 1 and pos2 >= thr2) or (obs_side2 == -1 and neg2 >= thr2):
        out['result']='swap_test7'
        return out

    out['result']='undetermined'
    return out

# Main
if __name__ == '__main__':
    if len(sys.argv) < 6:
        print('Usage: pascal_sim.py MODEL.obj H V W distance')
        sys.exit(1)
    path = sys.argv[1]
    H_arg = sys.argv[2]
    if H_arg.lower() == 'sweep':
        # sweep H 0..359 and report problematic angles
        V = int(sys.argv[3]); W = int(sys.argv[4]); D = float(sys.argv[5])
        verts, faces = parse_obj(path)
        problematic = []
        for H in range(0,360):
            xo, yo, zo, x2d, y2d = transform_and_project(verts, H, V, W, D)
            face_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)
            # run painter sim
            def simulate_painter_for_report(face_info):
                n = len(face_info)
                order = list(range(n))
                order.sort(key=lambda ii: (-face_info[ii]['zmean'], ii))
                ordered_pairs = set()
                checksum_history = []
                max_passes = 50
                for pass_no in range(1, max_passes+1):
                    swapped = 0; swap_count=0; undetermined=0
                    for i in range(n-1):
                        f1 = order[i]; f2 = order[i+1]
                        if (f1,f2) in ordered_pairs: continue
                        res = test_pair(face_info, f1, f2)
                        r = res['result']
                        if r == 'swap_due_depth' or r.startswith('swap'):
                            order[i], order[i+1] = order[i+1], order[i]
                            swapped = 1; swap_count += 1; ordered_pairs.add((f2,f1))
                        elif r.startswith('ok'):
                            ordered_pairs.add((f1,f2))
                        else:
                            undetermined += 1; ordered_pairs.add((f2,f1))
                    csum = tuple(order)
                    if checksum_history.count(csum) > 1:
                        return ('oscillation', pass_no, swap_count, undetermined)
                    checksum_history.append(csum)
                    if swapped == 0:
                        return ('stabilized', pass_no, swap_count, undetermined)
                return ('maxed', max_passes, swap_count, undetermined)
            rep = simulate_painter_for_report(face_info)
            if rep[0] != 'stabilized' or rep[3] != 0:
                problematic.append((H, rep))
        if not problematic:
            print('Sweep result: no problematic angles found (no undetermined/oscillation)')
        else:
            print('Sweep found problematic angles:')
            for H,rep in problematic:
                print(f' H={H} -> {rep}')
        sys.exit(0)
    else:
        H = int(H_arg); V = int(sys.argv[3]); W = int(sys.argv[4]); D = float(sys.argv[5])

        verts, faces = parse_obj(path)
        xo, yo, zo, x2d, y2d = transform_and_project(verts, H, V, W, D)
        face_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)
        # Print per-face summary
        print(f'Model: {path} H={H} V={V} W={W} D={D}\nFaces: {len(faces)} vertices: {len(verts)}')
        for i,f in enumerate(face_info):
            print(f'Face {i}: n={len(f["idxs"])}, display={f["display"]}, zmin={to_float(f["zmin"]) if f["zmin"] is not None else None}, zmax={to_float(f["zmax"]) if f["zmax"] is not None else None}, a={to_float(f["a"])}, b={to_float(f["b"])}, c={to_float(f["c"])}, d={to_float(f["d"])})')
        print('\nPairwise tests (single-shot diagnostics):')
        for i in range(len(faces)):
            for j in range(i+1, len(faces)):
                res = test_pair(face_info, i, j)
                if res['result'] != 'ok_no_overlap_x' and res['result'] != 'ok_no_overlap_y' and res['result'] != 'ok_test4' and res['result'] != 'ok_test5':
                    print(f'Pair ({i},{j}): {res}')
                else:
                    print(f'Pair ({i},{j}): {res["result"]}')

        # Now simulate the painter correction loop (adjacent swaps)
        def simulate_painter(face_info):
            n = len(face_info)
            # initial order: sort by zmean descending, tie-breaker index
            order = list(range(n))
            order.sort(key=lambda ii: (-face_info[ii]['zmean'], ii))
            print('\nInitial order (by zmean desc):', order)
            ordered_pairs = set()  # pairs definitively ordered (face1,face2)
            checksum_history = []
            max_passes = 50
            for pass_no in range(1, max_passes+1):
                swapped = 0
                swap_count = 0
                undetermined = 0
                for i in range(n-1):
                    f1 = order[i]; f2 = order[i+1]
                    if (f1,f2) in ordered_pairs:
                        continue
                    res = test_pair(face_info, f1, f2)
                    r = res['result']
                    if r == 'swap_due_depth' or r.startswith('swap'):
                        # swap
                        order[i], order[i+1] = order[i+1], order[i]
                        swapped = 1; swap_count += 1
                        ordered_pairs.add((f2,f1))
                    elif r.startswith('ok'):
                        # ordered correctly
                        ordered_pairs.add((f1,f2))
                    else:
                        undetermined += 1
                        ordered_pairs.add((f2,f1))
                # compute checksum
                csum = tuple(order)
                checksum_history.append(csum)
                if checksum_history.count(csum) > 1:
                    print(f'[WARN] Oscillation detected at pass {pass_no}; order repeats: {order}')
                    return order, pass_no, swap_count, undetermined
                print(f'PASS {pass_no}: swaps={swap_count} undetermined={undetermined} order={order}')
                if swapped == 0:
                    print('[INFO] Painter stabilized')
                    return order, pass_no, swap_count, undetermined
            print('[WARN] Reached max passes without stabilization')
            return order, max_passes, swap_count, undetermined

    def compute_pascal_order(face_info):
        # Full pairwise insertion: repeatedly ensure for every pair (i<j) that order respects tests
        n = len(face_info)
        order = list(range(n))
        changed = True
        passes = 0
        while changed and passes < 500:
            changed = False
            passes += 1
            for i in range(n):
                for j in range(i+1, n):
                    f_i = order[i]; f_j = order[j]
                    # test pair as in painter but deciding if we need to invert
                    res = test_pair(face_info, f_i, f_j)
                    r = res['result']
                    # if res indicates f_j should be before f_i, move it
                    if r == 'swap_due_depth' or r.startswith('swap'):
                        # move f_j before f_i
                        order.pop(j)
                        order.insert(i, f_j)
                        changed = True
                        break
                if changed: break
        return order, passes

    final_order, passes, swaps, und = simulate_painter(face_info)
    pascal_order, pascal_passes = compute_pascal_order(face_info)
    print('\nPainter simulation result:')
    print(f' final_order={final_order} passes={passes} swaps={swaps} undetermined={und}')
    print('Pascal-sim order:', pascal_order, 'passes:', pascal_passes)
