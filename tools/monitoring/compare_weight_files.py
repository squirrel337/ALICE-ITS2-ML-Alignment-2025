#!/usr/bin/env python3
"""Compare two weights.txt files sensor by sensor and say WHERE they differ.

    compare_weight_files.py A.txt B.txt [--top N]

A weights file is three normalisation/header sections followed by one line per
sensor: the sensor id, then the network parameters. This walks both files in
lockstep and reports the difference broken down by parameter column and by
detector layer, because the shape of the disagreement is what identifies its
cause -- a constant offset everywhere points at a global constant, a spread that
tracks |value| points at rounding, and a difference confined to some sensors or
some layers points at a different set of tracks reaching them.
"""
import sys

CHIP_BOUNDARY = [0, 108, 252, 432, 3120, 6480, 14712, 24120]

def layer_of(sid):
    for l in range(7):
        if CHIP_BOUNDARY[l] <= sid < CHIP_BOUNDARY[l + 1]:
            return l
    return -1

def read(path):
    head, rows = [], {}
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if line.startswith('#') or len(t) < 4:
                head.append(line.rstrip('\n')); continue
            try:
                sid = int(t[0]); vals = [float(x) for x in t[1:]]
            except ValueError:
                head.append(line.rstrip('\n')); continue
            rows[sid] = vals
    return head, rows

def main():
    a_path, b_path = sys.argv[1], sys.argv[2]
    top = 10
    if '--top' in sys.argv:
        top = int(sys.argv[sys.argv.index('--top') + 1])

    ha, A = read(a_path); hb, B = read(b_path)
    print("A = %s   (%d sensors)" % (a_path, len(A)))
    print("B = %s   (%d sensors)" % (b_path, len(B)))
    if ha != hb:
        print("\n!! header sections differ")
        for x, y in zip(ha, hb):
            if x != y: print("   A: %s\n   B: %s" % (x, y))
    only = set(A) ^ set(B)
    if only:
        print("\n!! %d sensor ids present in only one file" % len(only))

    ids = sorted(set(A) & set(B))
    if not ids:
        print("no common sensors"); return
    ncol = min(len(A[ids[0]]), len(B[ids[0]]))

    # ---- per column -------------------------------------------------------
    print("\ncol :   n_diff        max|dA|        max rel      RMS rel   at sensor")
    worst = []
    for c in range(ncol):
        nd = 0; mx = 0.0; mrel = 0.0; s2 = 0.0; where = -1
        for sid in ids:
            x, y = A[sid][c], B[sid][c]
            d = abs(x - y)
            if d == 0.0: continue
            nd += 1
            r = d / max(abs(x), abs(y)) if max(abs(x), abs(y)) > 0 else float('inf')
            s2 += r * r
            if d > mx: mx, where = d, sid
            if r > mrel: mrel = r
        if nd:
            print("%3d : %8d  %13.4e  %13.4e  %11.4e  %7d"
                  % (c, nd, mx, mrel, (s2 / nd) ** 0.5, where))
            worst.append((mx, c))
        else:
            print("%3d : %8d  %13s" % (c, 0, "identical"))

    # ---- per layer --------------------------------------------------------
    print("\nlayer : n_sensors  n_with_diff     max|dA|      mean|dA|")
    for l in range(7):
        sel = [s for s in ids if layer_of(s) == l]
        nd = 0; mx = 0.0; tot = 0.0; cnt = 0
        for sid in sel:
            dm = max(abs(A[sid][c] - B[sid][c]) for c in range(ncol))
            if dm > 0: nd += 1
            mx = max(mx, dm); tot += dm; cnt += 1
        print("%5d : %9d  %11d  %11.4e  %12.4e"
              % (l, len(sel), nd, mx, tot / cnt if cnt else 0))

    # ---- worst sensors ----------------------------------------------------
    rank = sorted(((max(abs(A[s][c] - B[s][c]) for c in range(ncol)), s) for s in ids),
                  reverse=True)[:top]
    print("\nlargest %d sensor differences:" % top)
    print("  sensor  layer      max|dA|   column   A                 B")
    for dm, sid in rank:
        c = max(range(ncol), key=lambda k: abs(A[sid][k] - B[sid][k]))
        print("  %6d  %5d  %11.4e  %6d   %-17.10g %-17.10g"
              % (sid, layer_of(sid), dm, c, A[sid][c], B[sid][c]))

    ndiff = sum(1 for s in ids if any(A[s][c] != B[s][c] for c in range(ncol)))
    print("\n%d of %d sensors differ (%.2f%%)" % (ndiff, len(ids), 100.0 * ndiff / len(ids)))

    # Per-sensor, does ONE scale factor explain all of that sensor's columns? This is
    # the question that separates causes. If it does, the sensor's gradient points the
    # same way in both runs and only its length moved -- the same tracks contributed,
    # and something scaled the result. If it does not, the per-track quantities differ,
    # so tracks entered or left the sum.
    cols = [c for c in range(ncol) if any(A[s][c] for s in ids)]
    print("\nlayer   n   median scale   5%..95%%      median residual after scaling")
    for l in range(7):
        ks, rr = [], []
        for sid in ids:
            if layer_of(sid) != l: continue
            num = sum(A[sid][c] * B[sid][c] for c in cols)
            den = sum(A[sid][c] ** 2 for c in cols)
            d0 = sum((B[sid][c] - A[sid][c]) ** 2 for c in cols)
            if den <= 0 or d0 <= 0: continue
            k = num / den
            d1 = sum((B[sid][c] - k * A[sid][c]) ** 2 for c in cols)
            ks.append(k); rr.append((d1 / d0) ** 0.5)
        if not ks: continue
        ks.sort(); rr.sort(); n = len(ks)
        print("%5d %5d %13.6f  %.4f..%-7.4f %15.4f"
              % (l, n, ks[n // 2], ks[n // 20], ks[-max(1, n // 20)], rr[n // 2]))
    print("   residual near 0 -> same direction, different magnitude")
    print("   residual near 1 -> the difference is orthogonal to the update")

main()
