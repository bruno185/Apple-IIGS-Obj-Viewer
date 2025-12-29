#!/usr/bin/env python3
import csv
# read C values from pair_debug_v4.csv
cvals = {}
with open('Windows/build_vs/bin/Release/pair_debug_v4.csv') as f:
    for line in f:
        if line.startswith('#'):
            parts=line.strip().split(',')
            if parts[0]=='#DETAIL' and parts[3]=='1' and parts[4].startswith('vertex_pq'):
                # not robust; we'll extract Ap,Bp,... from parent pair line
                pass
# simpler: parse pair_debug to get plane row first
planes = {}
with open('Windows/build_vs/bin/Release/pair_debug_v4.csv') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts=line.strip().split(',')
        if len(parts) > 12 and parts[0].isdigit():
            f1=int(parts[0]); f2=int(parts[1]); Ap=float(parts[4]); Bp=float(parts[5]); Cp=float(parts[6]); Dp=float(parts[7])
            planes[(f1,f2)]=(Ap,Bp,Cp,Dp)
# get g_obsv from pascal file (pair_values_pascal.csv has xo,yo,zo per vertex)
obs = {}
with open('Windows/build_vs/bin/Release/pair_values_pascal.csv') as f:
    rdr=csv.DictReader(f)
    for r in rdr:
        if r['stage'].startswith('vertex_'):
            vid=int(r['idx']); xo=float(r['xo']); yo=float(r['yo']); zo=float(r['zo'])
            obs[vid]=(xo,yo,zo)
# pick pair (0,1) vid 1
Ap,Bp,Cp,Dp = planes[(0,1)]
xo,yo,zo = obs[1]
val = Ap*xo + Bp*yo + Cp*zo + Dp
print('Using C plane coeffs and observer coords:')
print('Ap,Bp,Cp,Dp =',Ap,Bp,Cp,Dp)
print('xo,yo,zo =',xo,yo,zo)
print('sum =',val)
# now Pascal plane coeffs from pair_values_pascal
pplanes = {}
with open('Windows/build_vs/bin/Release/pair_values_pascal.csv') as f:
    rdr=csv.DictReader(f)
    for r in rdr:
        if r['stage']=='planes':
            f1=int(r['i']); f2=int(r['j'])
            pplanes[(f1,f2)]=(float(r['Ap']),float(r['Bp']),float(r['Cp']),float(r['Dp']))
Ap_p,Bp_p,Cp_p,Dp_p = pplanes[(0,1)]
valp = Ap_p*xo + Bp_p*yo + Cp_p*zo + Dp_p
print('Using Pascal plane coeffs:',Ap_p,Bp_p,Cp_p,Dp_p)
print('sum_pascal =',valp)
