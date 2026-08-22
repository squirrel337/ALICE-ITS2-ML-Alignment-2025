#!/usr/bin/env python3
"""Per-parameter run-to-run spread of the epoch-0 weights.

    weight_spread.py OUT.png LABEL=costs.tsv [LABEL=costs.tsv ...]

Takes the costs.tsv written by run_to_run.sh, whose fifth column is the path to
that repetition's weights_Epoch_At_0.txt, loads every one of them, and measures
how much each of the 24120 x 17 sensor parameters moves between repetitions of
the *same* configuration.

This is the quantity that matters: the cost is one scalar summary, but the
alignment is the parameter set, and a comparison between two versions of the
module is only meaningful outside the band that repeating one version already
produces. Spread is reported per network-parameter column, since the columns are
on wildly different scales and a single pooled number hides which ones move.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NPAR = 17
# The eleven network parameters that map to the six rigid-body degrees of
# freedom, plus the bookkeeping columns the file carries alongside them.
SKIP_HEADER = "#"


def load_one(path):
    rows = []
    started = False
    for line in open(path):
        if line.startswith(SKIP_HEADER):
            started = "synapses weights" in line
            continue
        if not started:
            continue
        f = line.split()
        if len(f) == NPAR + 1:
            rows.append([float(v) for v in f[1:]])
    return np.asarray(rows)


def load_arm(tsv):
    paths = []
    for line in open(tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5 and f[3] == "ok" and f[4] != "NA":
            paths.append(f[4])
    mats = []
    for p in paths:
        try:
            m = load_one(p)
        except OSError:
            continue
        if m.size:
            mats.append(m)
    if not mats:
        return None
    shape = mats[0].shape
    mats = [m for m in mats if m.shape == shape]
    return np.stack(mats)          # (nrun, nsensor, NPAR)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out = sys.argv[1]
    arms = []
    for spec in sys.argv[2:]:
        label, _, tsv = spec.partition("=")
        a = load_arm(tsv)
        if a is None:
            print(f"  {label}: no weight files listed in {tsv}, skipped")
            continue
        arms.append((label, a))
        print(f"  {label}: {a.shape[0]} runs x {a.shape[1]} sensors x {a.shape[2]} parameters")
    if not arms:
        sys.exit("nothing to analyse")

    print()
    for label, a in arms:
        sd = a.std(axis=0, ddof=1)          # (nsensor, NPAR)
        mean = a.mean(axis=0)
        moved = (sd > 0).sum()
        tot = sd.size
        print(f"=== {label} : {a.shape[0]} repetitions of the identical configuration")
        print(f"    parameters that differ between repetitions : "
              f"{moved}/{tot}  ({100*moved/tot:.1f}%)")
        # relative spread where the parameter is not identically zero
        nz = np.abs(mean) > 0
        rel = np.zeros_like(sd)
        rel[nz] = sd[nz] / np.abs(mean[nz])
        print(f"    absolute sd : median {np.median(sd[sd>0]) if moved else 0:.3e}   "
              f"p95 {np.percentile(sd[sd>0], 95) if moved else 0:.3e}   "
              f"max {sd.max():.3e}")
        print(f"    relative sd : median {np.median(rel[rel>0]) if (rel>0).any() else 0:.3e}   "
              f"p95 {np.percentile(rel[rel>0], 95) if (rel>0).any() else 0:.3e}")
        print("    per column (sd across repetitions, median over sensors):")
        for c in range(a.shape[2]):
            s = sd[:, c]
            if (s > 0).any():
                print(f"      col {c:2d}  median {np.median(s[s>0]):.3e}  max {s.max():.3e}  "
                      f"({100*(s>0).mean():.0f}% of sensors move)")
            else:
                print(f"      col {c:2d}  identical in every repetition")
        print()

    # ---------------------------------------------------------------- plot
    colors = ["#2b6cb0", "#c05621", "#2f855a"]
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6))

    for i, (label, a) in enumerate(arms):
        sd = a.std(axis=0, ddof=1)
        s = sd[sd > 0]
        if s.size:
            ax[0].hist(np.log10(s), bins=60, alpha=.55, color=colors[i],
                       label=f"{label} (n={a.shape[0]})")
    ax[0].set_xlabel(r"$\log_{10}$ (per-parameter sd across repetitions)")
    ax[0].set_ylabel("parameters")
    ax[0].set_title("같은 설정을 반복했을 때 파라미터가 흔들리는 크기")
    ax[0].legend(frameon=False, fontsize=9)
    ax[0].grid(alpha=.25)

    width = .8 / max(len(arms), 1)
    for i, (label, a) in enumerate(arms):
        sd = a.std(axis=0, ddof=1)
        med = [np.median(sd[:, c][sd[:, c] > 0]) if (sd[:, c] > 0).any() else np.nan
               for c in range(a.shape[2])]
        x = np.arange(a.shape[2]) + i * width - .4 + width / 2
        ax[1].bar(x, med, width=width, color=colors[i], label=label, alpha=.85)
    ax[1].set_yscale("log")
    ax[1].set_xlabel("network parameter column")
    ax[1].set_ylabel("median sd across repetitions")
    ax[1].set_title("컬럼별 run-to-run 폭")
    ax[1].set_xticks(range(arms[0][1].shape[2]))
    ax[1].legend(frameon=False, fontsize=9)
    ax[1].grid(alpha=.25, axis="y")

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
