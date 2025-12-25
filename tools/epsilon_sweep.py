#!/usr/bin/env python3
"""Try different EPS thresholds to see effect on painter vs full-pair ordering
Usage: python tools/epsilon_sweep.py MODEL.obj H V W D
"""
import sys
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float, to_fixed


def test_pair_eps(face_info, xo, yo, zo, f1i, f2i, eps_fixed):
    f1 = face_info[f1i]; f2 = face_info[f2i]
    # Test1 depth
    if f2['zmax'] <= f1['zmin']: return 'ok_depth1'
    if f1['zmax'] <= f2['zmin']: return 'swap_depth'
    # bbox
    if f1['maxx'] <= f2['minx'] or f2['maxx'] <= f1['minx']: return 'ok_bbox_x'
    if f1['maxy'] <= f2['miny'] or f2['maxy'] <= f1['miny']: return 'ok_bbox_y'
    a1,b1,c1,d1 = f1['a'], f1['b'], f1['c'], f1['d']
    a2,b2,c2,d2 = f2['a'], f2['b'], f2['c'], f2['d']
    obs1 = 0
    if d1 > eps_fixed: obs1=1
    elif d1 < -eps_fixed: obs1=-1
    else: obs1=0
    if obs1 != 0:
        all_same = True
        for idx in f2['idxs']:
            tv = a1*xo[idx] + b1*yo[idx] + c1*zo[idx] + d1
            if tv > eps_fixed: side=1
            elif tv < -eps_fixed: side=-1
            else: side=0
            if side != obs1:
                all_same=False; break
        if all_same: return 'ok_test4'
    obs2=0
    if d2 > eps_fixed: obs2=1
    elif d2 < -eps_fixed: obs2=-1
    else: obs2=0
    if obs2 != 0:
        all_opposite=True
        for idx in f1['idxs']:
            tv = a2*xo[idx] + b2*yo[idx] + c2*zo[idx] + d2
            if tv > eps_fixed: side=1
            elif tv < -eps_fixed: side=-1
            else: side=0
            if side == obs2:
                all_opposite=False; break
        if all_opposite: return 'ok_test5'
    # Test 6
    if obs1 != 0:
        all_opposite=True
        for idx in f2['idxs']:
            tv = a1*xo[idx] + b1*yo[idx] + c1*zo[idx] + d1
            if tv > eps_fixed: side=1
            else: side=-1
            if side == obs1:
                all_opposite=False; break
        if all_opposite: return 'swap_test6'
    if obs2 != 0:
        all_same=True
        for idx in f1['idxs']:
            tv = a2*xo[idx] + b2*yo[idx] + c2*zo[idx] + d2
            if tv > eps_fixed: side=1
            else: side=-1
            if side != obs2:
                all_same=False; break
        if all_same: return 'swap_test7'
    return 'undetermined'


def simulate_painter_insert(face_info, xo, yo, zo, eps_fixed):
    """Insertion-based: when swap/undetermined occurs, insert the right element earlier by scanning left
    This often mimics full-pairwise insertion but with lower cost normally."""
    n = len(face_info)
    order = list(range(n))
    order.sort(key=lambda ii: (-face_info[ii]['zmean'], ii))
    # one pass inserting elements where needed, repeat until stable
    for pass_n in range(50):
        changed = False
        i = 0
        while i < len(order)-1:
            f1, f2 = order[i], order[i+1]
            res = test_pair_eps(face_info, xo, yo, zo, f1, f2, eps_fixed)
            if res.startswith('swap') or res=='swap_depth' or res=='undetermined':
                # remove f2 and insert it at the earliest position where it should be
                cand = order.pop(i+1)
                j = i
                # move left while cand should be before order[j-1]
                while j > 0:
                    prev = order[j-1]
                    # if prev should be after cand, i.e., test_pair(prev,cand) startswith 'swap', then we should put cand before prev
                    r = test_pair_eps(face_info, xo, yo, zo, prev, cand, eps_fixed)
                    if r.startswith('swap') or r=='swap_depth' or r=='undetermined':
                        j -= 1
                    else:
                        break
                order.insert(j, cand)
                changed = True
                # continue scanning after inserted element
                i = j + 1
            else:
                i += 1
        if not changed:
            return order
    return order

def simulate_painter_insert(face_info, xo, yo, zo, eps_fixed):
    """Insertion-based: when swap/undetermined occurs, insert the right element earlier by scanning left
    This often mimics full-pairwise insertion but with lower cost normally."""
    n = len(face_info)
    order = list(range(n))
    order.sort(key=lambda ii: (-face_info[ii]['zmean'], ii))
    # one pass inserting elements where needed, repeat until stable
    for pass_n in range(50):
        changed = False
        i = 0
        while i < len(order)-1:
            f1, f2 = order[i], order[i+1]
            res = test_pair_eps(face_info, xo, yo, zo, f1, f2, eps_fixed)
            if res.startswith('swap') or res=='swap_depth' or res=='undetermined':
                # remove f2 and insert it at the earliest position where it should be
                cand = order.pop(i+1)
                j = i
                # move left while cand should be before order[j-1]
                while j > 0:
                    prev = order[j-1]
                    # if prev should be after cand, i.e., test_pair(prev,cand) startswith 'swap', then we should put cand before prev
                    r = test_pair_eps(face_info, xo, yo, zo, prev, cand, eps_fixed)
                    if r.startswith('swap') or r=='swap_depth' or r=='undetermined':
                        j -= 1
                    else:
                        break
                order.insert(j, cand)
                changed = True
                # continue scanning after inserted element
                i = j + 1
            else:
                i += 1
        if not changed:
            return order
    return order


def full_pascal_order(face_info, xo, yo, zo, eps_fixed):
    # a simple O(n^2) comparator + insertion sort using test_pair_eps
    n = len(face_info)
    order = [0]
    for i in range(1,n):
        placed=False
        for j in range(len(order)):
            res = test_pair_eps(face_info, xo, yo, zo, order[j], i, eps_fixed)
            if res.startswith('swap'): # i should be before order[j]
                order.insert(j, i)
                placed=True
                break
        if not placed:
            order.append(i)
    return order

if __name__=='__main__':
    if len(sys.argv) < 6:
        print('Usage: epsilon_sweep.py MODEL.obj H V W D')
        sys.exit(1)
    path=sys.argv[1]; H=int(sys.argv[2]); V=int(sys.argv[3]); W=int(sys.argv[4]); D=float(sys.argv[5])
    verts, faces = parse_obj(path)
    xo, yo, zo, x2d, y2d = transform_and_project(verts, H, V, W, D)
    face_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)
    eps_list = [0.01, 0.02, 0.05, 0.1, 0.5, 1.0]
    for eps in eps_list:
        eps_fixed = to_fixed(eps)
        painter_order = simulate_painter(face_info, xo, yo, zo, eps_fixed)
        painter_insert_order = simulate_painter_insert(face_info, xo, yo, zo, eps_fixed)
        pascal_order = full_pascal_order(face_info, xo, yo, zo, eps_fixed)
        print(f'eps {eps:.3f}: painter_adj {painter_order} painter_insert {painter_insert_order} pascal {pascal_order}')
        # print decisions for problematic pairs
        test_pairs = [(4,0),(4,1),(5,0),(5,1),(0,1)]
        for p in test_pairs:
            res = test_pair_eps(face_info, xo, yo, zo, p[0], p[1], eps_fixed)
            print(f'   pair {p}: {res}')

    print('\nDone sweep.')
