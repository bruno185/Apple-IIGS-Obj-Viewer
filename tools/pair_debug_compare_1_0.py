#!/usr/bin/env python3
import csv
from pascal_sim import parse_obj, transform_and_project, compute_faces, to_float
import math

PVF = 'Windows\\build_vs\\bin\\Release\\pair_values_v4.csv'
OBJ = 'Windows\\3d objects\\q2.obj'
H=30;V=20;W=0;D=39.214455

def arr(r): return 0.0 if abs(r) < 0.01 else r

verts, faces = parse_obj(OBJ)
xo,yo,zo,x2d,y2d = transform_and_project(verts,H,V,W,D)
faces_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)
plane_map = {fi: (arr(to_float(f['a'])), arr(to_float(f['b'])), arr(to_float(f['c'])), arr(to_float(f['d']))) for fi,f in enumerate(faces_info)}

# Read CSV with header
with open(PVF,newline='') as f:
    rdr = csv.reader(f)
    header = next(rdr)
    rows = list(rdr)

hmap = {name:idx for idx,name in enumerate(header)}

pair = (1,0)
f1,f2 = pair
rows_pair = [r for r in rows if len(r)> max(hmap['f1'],hmap['f2']) and int(r[hmap['f1']])==f1 and int(r[hmap['f2']])==f2]

# find planes row
for r in rows_pair:
    if r[hmap['stage']]=='planes':
        # find numeric tokens after stage
        stage_idx = hmap['stage']
        tokens = [t for t in r[stage_idx+1:] if t!='']
        nums=[]
        for t in tokens:
            try: nums.append(float(t))
            except: pass
        print('C planes tokens:',nums[:8])
        break

Ap_py,Bp_py,Cp_py,Dp_py = plane_map[f1]
Aq_py,Bq_py,Cq_py,Dq_py = plane_map[f2]
print('PY plane f1:',(Ap_py,Bp_py,Cp_py,Dp_py))
print('PY plane f2:',(Aq_py,Bq_py,Cq_py,Dq_py))

# print per-vertex comparison
print('\nPer-vertex comparison for vertex_pq (f2 verts on plane f1):')
for r in rows_pair:
    if r[hmap['stage']]=='vertex_pq':
        # tokens after stage
        stage_idx = hmap['stage']
        tokens = [t for t in r[stage_idx+1:] if t!='']
        nums=[]
        for t in tokens:
            try: nums.append(float(t))
            except: nums.append(None)
        # Try to map robustly
        # Some rows are: idx,xo,yo,zo,wx,wy,wz,ax_obs,by_obs,cz_obs,raw_obs,tempo_obs,ax_w,by_w,cz_w,raw_w,tempo_world
        if len(nums) >= 17:
            vid = int(nums[0]); xo_c,yo_c,zo_c = nums[1],nums[2],nums[3]; wx_c,wy_c,wz_c = nums[4],nums[5],nums[6]
            ax_obs_c,by_obs_c,cz_obs_c,raw_obs_c,tempo_obs_c = nums[7],nums[8],nums[9],nums[10],nums[11]
            ax_w_c,by_w_c,cz_w_c,raw_w_c,tempo_w_c = nums[12],nums[13],nums[14],nums[15],nums[16]
        else:
            # fallback using header indices
            def g(n):
                try: return float(r[hmap[n]])
                except: return None
            vid = int(r[hmap['idx']])
            xo_c,yo_c,zo_c = g('xo'), g('yo'), g('zo')
            wx_c,wy_c,wz_c = g('wx'), g('wy'), g('wz')
            ax_obs_c,by_obs_c,cz_obs_c = g('ax_obs'), g('by_obs'), g('cz_obs')
            raw_obs_c,tempo_obs_c = g('raw_obs'), g('tempo_obs')
            ax_w_c,by_w_c,cz_w_c = g('ax_w'), g('by_w'), g('cz_w')
            raw_w_c,tempo_w_c = g('raw_w'), g('tempo_world')
        # Python expected values
        xo_py,yo_py,zo_py = to_float(xo[vid]), to_float(yo[vid]), to_float(zo[vid])
        wx_py,wy_py,wz_py = verts[vid]
        ax_obs_py = Ap_py * xo_py if xo_py is not None else None
        raw_obs_py = None
        if None not in (ax_obs_py, yo_c, zo_c):
            by_obs_py = Bp_py * yo_c
            cz_obs_py = Cp_py * zo_c
            raw_obs_py = ax_obs_py + by_obs_py + cz_obs_py + Dp_py
            tempo_obs_py = arr(raw_obs_py)
        else:
            by_obs_py = Bp_py * yo_c if yo_c is not None else None
            cz_obs_py = Cp_py * zo_c if zo_c is not None else None
            tempo_obs_py = None
        ax_w_py = Ap_py * wx_py
        raw_w_py = ax_w_py + Bp_py * wy_py + Cp_py * wz_py + Dp_py
        tempo_w_py = arr(raw_w_py)
        print(f'VID {vid}: C xo={xo_c:.6f} {yo_c:.6f} {zo_c:.6f} | PY xo={xo_py:.6f} {yo_py:.6f} {zo_py:.6f}')
        print(f'       C ax_obs={ax_obs_c}, raw_obs={raw_obs_c}, tempo_obs={tempo_obs_c} | PY ax_obs={ax_obs_py:.6f}, raw_obs={raw_obs_py:.6f}, tempo_obs={tempo_obs_py}')
        print(f'       C ax_w={ax_w_c}, raw_w={raw_w_c}, tempo_w={tempo_w_c} | PY ax_w={ax_w_py:.6f}, raw_w={raw_w_py:.6f}, tempo_w={tempo_w_py}')

# same for vertex_qp
print('\nPer-vertex comparison for vertex_qp (f1 verts on plane f2):')
for r in rows_pair:
    if r[hmap['stage']]=='vertex_qp':
        stage_idx = hmap['stage']
        tokens = [t for t in r[stage_idx+1:] if t!='']
        nums=[]
        for t in tokens:
            try: nums.append(float(t))
            except: nums.append(None)
        if len(nums) >= 17:
            vid = int(nums[0]); xo_c,yo_c,zo_c = nums[1],nums[2],nums[3]; wx_c,wy_c,wz_c = nums[4],nums[5],nums[6]
            ax_obs_c,by_obs_c,cz_obs_c,raw_obs_c,tempo_obs_c = nums[7],nums[8],nums[9],nums[10],nums[11]
            ax_w_c,by_w_c,cz_w_c,raw_w_c,tempo_w_c = nums[12],nums[13],nums[14],nums[15],nums[16]
        else:
            def g(n):
                try: return float(r[hmap[n]])
                except: return None
            vid = int(r[hmap['idx']])
            xo_c,yo_c,zo_c = g('xo'), g('yo'), g('zo')
            wx_c,wy_c,wz_c = g('wx'), g('wy'), g('wz')
            ax_obs_c,by_obs_c,cz_obs_c = g('ax_obs'), g('by_obs'), g('cz_obs')
            raw_obs_c,tempo_obs_c = g('raw_obs'), g('tempo_obs')
            ax_w_c,by_w_c,cz_w_c = g('ax_w'), g('by_w'), g('cz_w')
            raw_w_c,tempo_w_c = g('raw_w'), g('tempo_world')
        # Python expected
        xo_py,yo_py,zo_py = to_float(xo[vid]), to_float(yo[vid]), to_float(zo[vid])
        ax_obs_py = Aq_py * xo_py
        raw_obs_py = ax_obs_py + Bq_py * yo_c + Cq_py * zo_c + Dq_py
        tempo_obs_py = arr(raw_obs_py)
        wx_py,wy_py,wz_py = verts[vid]
        ax_w_py = Aq_py * wx_py
        raw_w_py = ax_w_py + Bq_py * wy_py + Cq_py * wz_py + Dq_py
        tempo_w_py = arr(raw_w_py)
        print(f'VID {vid}: C xo={xo_c:.6f} {yo_c:.6f} {zo_c:.6f} | PY xo={xo_py:.6f} {yo_py:.6f} {zo_py:.6f}')
        print(f'       C ax_obs={ax_obs_c}, raw_obs={raw_obs_c}, tempo_obs={tempo_obs_c} | PY ax_obs={ax_obs_py:.6f}, raw_obs={raw_obs_py:.6f}, tempo_obs={tempo_obs_py}')
        print(f'       C ax_w={ax_w_c}, raw_w={raw_w_c}, tempo_w={tempo_w_c} | PY ax_w={ax_w_py:.6f}, raw_w={raw_w_py:.6f}, tempo_w={tempo_w_py}')

print('\nDone')
