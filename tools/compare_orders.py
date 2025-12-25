from pascal_sim import parse_obj, transform_and_project, compute_faces, test_pair

verts, faces = parse_obj('cube.obj')
H=30; V=20; W=0; D=180.0
xo,yo,zo,x2d,y2d = transform_and_project(verts,H,V,W,D)
face_info = compute_faces(verts,faces,xo,yo,zo,x2d,y2d)

# compute painter one-pass then stabilize
order = list(range(len(face_info)))
order.sort(key=lambda ii: (-face_info[ii]['zmean'], ii))
# simulate painter passes until stable
changed=True
while changed:
    changed=False
    # define local test_pair using xo,yo,zo captured above
    def local_test_pair(faces, f1i, f2i):
        f1=faces[f1i]; f2=faces[f2i]
        # Depth quick test
        if f2['zmax'] <= f1['zmin']: return {'result':'ok_order_depth'}
        if f1['zmax'] <= f2['zmin']: return {'result':'swap_due_depth'}
        # bbox
        if f1['maxx'] <= f2['minx'] or f2['maxx'] <= f1['minx']: return {'result':'ok_no_overlap_x'}
        if f1['maxy'] <= f2['miny'] or f2['maxy'] <= f1['miny']: return {'result':'ok_no_overlap_y'}
        a1,b1,c1,d1 = f1['a'], f1['b'], f1['c'], f1['d']
        a2,b2,c2,d2 = f2['a'], f2['b'], f2['c'], f2['d']
        epsilon = 1 << 10  # approx 0.015625 in fixed scale but we operate floats here
        # compute obs_side1
        obs_side1 = 0
        if d1 > 0.01: obs_side1 = 1
        elif d1 < -0.01: obs_side1 = -1
        if obs_side1 != 0:
            all_same = True
            for idx in f2['idxs']:
                # idxs are 1-based in CSV or 0-based? face_info uses 0-based vertex indices
                vid = idx
                tv = a1*xo[vid] + b1*yo[vid] + c1*zo[vid] + d1
                if tv > 0.01: side = 1
                elif tv < -0.01: side = -1
                else: side = 0
                if side != obs_side1:
                    all_same = False; break
            if all_same: return {'result':'ok_test4'}
        # Test5
        obs_side2 = 0
        if d2 > 0.01: obs_side2 = 1
        elif d2 < -0.01: obs_side2 = -1
        if obs_side2 != 0:
            all_opposite = True
            for idx in f1['idxs']:
                vid = idx
                tv = a2*xo[vid] + b2*yo[vid] + c2*zo[vid] + d2
                if tv > 0.01: side = 1
                elif tv < -0.01: side = -1
                else: side = 0
                if side == obs_side2:
                    all_opposite = False; break
            if all_opposite: return {'result':'ok_test5'}
        # tests 6 & 7 omitted for brevity; treat as undetermined
        return {'result':'undetermined'}

    for i in range(len(order)-1):
        res = local_test_pair(face_info, order[i], order[i+1])
        if res['result']=='swap_due_depth' or res['result'].startswith('swap'):
            order[i], order[i+1] = order[i+1], order[i]
            changed=True

# compute pascal full pairwise insertion (iterative)
pascal = list(range(len(face_info)))
changed=True
while changed:
    changed=False
    for i in range(len(pascal)):
        for j in range(i+1, len(pascal)):
            r = local_test_pair(face_info, pascal[i], pascal[j])
            if r['result']=='swap_due_depth' or r['result'].startswith('swap'):
                val = pascal.pop(j)
                pascal.insert(i, val)
                changed=True
                break
        if changed: break

print('painter order', order)
print('pascal order', pascal)

# list pairwise disagreements
pos_painter={f:i for i,f in enumerate(order)}
pos_pascal={f:i for f,i in pos_painter.items()}
# Note: pos_pascal mapping
pos_pascal={f:i for i,f in enumerate(pascal)}
print('\nPairwise disagreements:')
for i in range(len(face_info)):
    for j in range(i+1,len(face_info)):
        p_pos = pos_painter[i] < pos_painter[j]
        pas_pos = pos_pascal[i] < pos_pascal[j]
        if p_pos != pas_pos:
            print(f' Faces {i} and {j}: painter says {"i before j" if p_pos else "j before i"}, pascal says {"i before j" if pas_pos else "j before i"}')

# Print more detailed reason for first few disagreements
count=0
for i in range(len(face_info)):
    for j in range(i+1,len(face_info)):
        p_pos = pos_painter[i] < pos_painter[j]
        pas_pos = pos_pascal[i] < pos_pascal[j]
        if p_pos != pas_pos:
            print('\nDetailed pair:', i, j)
            print(' face', i, 'a,b,c,d =', face_info[i]['a'], face_info[i]['b'], face_info[i]['c'], face_info[i]['d'])
            print(' face', j, 'a,b,c,d =', face_info[j]['a'], face_info[j]['b'], face_info[j]['c'], face_info[j]['d'])
            print(' tests single-shot:', local_test_pair(face_info, i, j))
            count+=1
            if count>=5: break
    if count>=5: break
