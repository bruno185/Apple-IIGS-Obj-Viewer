from pascal_sim import parse_obj, transform_and_project, compute_faces

verts, faces = parse_obj('Windows\\3d objects\\q2.obj')
H=30; V=20; W=0; D=300.0
xo,yo,zo,x2d,y2d = transform_and_project(verts,H,V,W,D)
face_info = compute_faces(verts,faces,xo,yo,zo,x2d,y2d)

# compute painter (simple one-pass stabilization as in compare_orders)
order = list(range(len(face_info)))
order.sort(key=lambda ii: (-face_info[ii]['zmean'], ii))
changed=True
while changed:
    changed=False
    for i in range(len(order)-1):
        f1=face_info[order[i]]; f2=face_info[order[i+1]]
        if f2['zmax'] <= f1['zmin']:
            continue
        if f1['zmax'] <= f2['zmin']:
            order[i],order[i+1]=order[i+1],order[i]; changed=True; continue

# pascal insertion approximation (depth-only fallback)
pascal = list(range(len(face_info)))
changed=True
while changed:
    changed=False
    for i in range(len(pascal)):
        for j in range(i+1,len(pascal)):
            f1=face_info[pascal[i]]; f2=face_info[pascal[j]]
            if f1['zmax'] <= f2['zmin']:
                continue
            if f2['zmax'] <= f1['zmin']:
                val = pascal.pop(j); pascal.insert(i,val); changed=True; break
        if changed: break

print('painter order', order)
print('pascal order', pascal)

print('\nPairwise disagreements:')
for i in range(len(face_info)):
    for j in range(i+1,len(face_info)):
        p_pos = order.index(i) < order.index(j)
        pas_pos = pascal.index(i) < pascal.index(j)
        if p_pos != pas_pos:
            print(' Faces', i, 'and', j, ': painter says', 'i before j' if p_pos else 'j before i', ', pascal says', 'i before j' if pas_pos else 'j before i')
