#!/usr/bin/env python3
from pascal_sim import parse_obj, transform_and_project, compute_faces, test_pair_rel

path='cube.obj'
H=30; V=20; W=0; D=180.0
verts, faces = parse_obj(path)
xo, yo, zo, x2d, y2d = transform_and_project(verts, H, V, W, D)
face_info = compute_faces(verts, faces, xo, yo, zo, x2d, y2d)

# painter adjacency using test_pair_rel
n = len(face_info)
order = list(range(n))
order.sort(key=lambda ii: (-face_info[ii]['zmean'], ii))
changed = True
passes=0
while changed and passes < 50:
    changed=False
    passes+=1
    for i in range(n-1):
        f1=order[i]; f2=order[i+1]
        res = test_pair_rel(face_info, f1, f2, xo, yo, zo)
        if res['result'].startswith('swap') or res['result']=='undetermined':
            order[i], order[i+1] = order[i+1], order[i]
            changed=True
print('painter_rel order', order)

# full pascal using test_pair_rel as comparator (simple insertion)
pascal_order=[0]
for i in range(1,n):
    placed=False
    for j in range(len(pascal_order)):
        res = test_pair_rel(face_info, pascal_order[j], i, xo, yo, zo)
        if res['result'].startswith('swap'):
            pascal_order.insert(j, i)
            placed=True
            break
    if not placed:
        pascal_order.append(i)
print('pascal_rel order', pascal_order)

# print pair results for problematic pairs
pairs=[(4,0),(4,1),(5,0),(5,1),(0,1)]
for p in pairs:
    r = test_pair_rel(face_info, p[0], p[1], xo, yo, zo)
    print('pair',p, r)