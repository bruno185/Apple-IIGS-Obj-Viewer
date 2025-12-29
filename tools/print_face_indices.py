from pascal_sim import parse_obj
verts,faces = parse_obj('Windows\\3d objects\\q2.obj')
print('num verts', len(verts), 'num faces', len(faces))
print('face 1 indices', faces[1])
print('face 1 indices 1-based', [i+1 for i in faces[1]])
