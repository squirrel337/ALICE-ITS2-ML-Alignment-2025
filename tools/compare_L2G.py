#!/usr/bin/env python3
"""Compare O2's getMatrixL2G against the ROOT-only cache, chip by chip.

    python3 compare_L2G.py o2_L2G.txt cache_L2G.txt

Reports the rotation difference and, separately, the translation difference
projected on the sensor normal (local y) and within the plane -- the normal
component is the one that decides whether the cache sits on the right plane.
"""
import sys, math

def load(p):
    d={}
    for ln in open(p):
        if ln.startswith('#') or not ln.strip(): continue
        f=ln.split(); d[int(f[0])]=[float(x) for x in f[1:13]]
    return d

a,b = load(sys.argv[1]), load(sys.argv[2])
common = sorted(set(a) & set(b))
print(f"chips: O2={len(a)}  cache={len(b)}  common={len(common)}")
if not common: sys.exit("no overlap")

LAY=[0,108,252,432,3120,6480,14712,24120]
def layer(i):
    for l in range(7):
        if LAY[l] <= i < LAY[l+1]: return l
    return -1

stat={l:{'n':0,'rot':0.0,'nrm':[],'inp':[]} for l in range(7)}
worst=(0,None)
for i in common:
    A,B=a[i],b[i]
    rot=max(abs(A[k]-B[k]) for k in range(9))
    # A's rows are (Rxx,Rxy,Rxz),(Ryx,..),(Rzx,..); local-y axis in global = column 1
    ny=(A[1],A[4],A[7])
    dT=[B[9+k]-A[9+k] for k in range(3)]
    n=sum(dT[k]*ny[k] for k in range(3))                 # along the sensor normal
    ip=math.sqrt(max(sum(x*x for x in dT)-n*n,0.0))      # within the plane
    L=layer(i); s=stat[L]
    s['n']+=1; s['rot']=max(s['rot'],rot); s['nrm'].append(n); s['inp'].append(ip)
    if rot>worst[0]: worst=(rot,i)

def ms(v):
    m=sum(v)/len(v); return m, math.sqrt(sum((x-m)**2 for x in v)/len(v))

print(f"\n{'lay':<5}{'chips':>7}{'max|dR|':>12}{'dT_normal mean':>17}{'rms':>9}{'dT_inplane mean':>18}")
for l in range(7):
    s=stat[l]
    if not s['n']: continue
    nm,nr=ms(s['nrm']); im,_=ms(s['inp'])
    print(f"L{l:<4}{s['n']:>7}{s['rot']:>12.2e}{1e4*nm:>14.3f} um{1e4*nr:>7.3f}{1e4*im:>15.3f} um")

alln=[x for l in stat for x in stat[l]['nrm']]
m,r=ms(alln)
print(f"\nmax rotation difference : {worst[0]:.2e} (chip {worst[1]})")
print(f"normal offset overall   : {1e4*m:+.3f} um  rms {1e4*r:.3f} um")
print("\ninterpretation")
print("  |dR| ~ 1e-14 and dT_normal ~ 0      -> cache reproduces O2, correction was right")
print("  dT_normal ~ -1 um uniformly         -> drop kEffLayer in export_geometry_cache.C")
print("  dT_normal ~ +4 um (IB) / -21 um (OB)-> O2 uses the chip frame; revert the whole fix")
