import math
FIXED_SCALE = 1<<16
sin_vals = [str(int(round(math.sin(math.radians(d))*FIXED_SCALE))) for d in range(360)]
cos_vals = [str(int(round(math.cos(math.radians(d))*FIXED_SCALE))) for d in range(360)]

def format_arr(vals, name):
    lines = []
    line = ""
    for i,v in enumerate(vals):
        if i%10==0:
            if line:
                lines.append(line.rstrip()+',')
            line = '    '+v
        else:
            line += ', '+v
    lines.append(line.rstrip()+',')
    return 'static const Fixed32 %s[360] = {\n%s\n};\n' % (name, '\n'.join(lines))

print(format_arr(sin_vals, 'sin_table'))
print(format_arr(cos_vals, 'cos_table'))
