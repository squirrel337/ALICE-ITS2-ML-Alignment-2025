#!/usr/bin/env python3
"""Run-to-run spread of the epoch-0 weights, for one or two backends.

    weight_spread.py OUT.png LABEL=costs.tsv [LABEL=costs.tsv ...]

Takes the costs.tsv written by run_to_run.sh, whose fifth column is the path to
that repetition's weights_Epoch_At_0.txt, loads every one of them, and measures
how far the 24120 x 17 sensor parameters move between repetitions of the *same*
configuration.

Why the weights and not the cost: the cost is one scalar summary, while the
alignment *is* the parameter set, so a comparison between two versions of the
module is a comparison of these numbers. Use plot_run_to_run.py alongside this if
the cost is wanted too -- both are recorded in costs.tsv.

Three numbers come out, in increasing order of usefulness:

  scalar    one RMS per run over all parameters. Cheap to plot, and its spread
            is the run-to-run band as a single figure.
  pairwise  the RMS difference between two runs of the same configuration,
            over all pairs. This is the band as it is actually used: "two runs
            of one version land this far apart".
  per-parameter  the sd of each parameter across repetitions, kept per column
            because the columns sit on very different scales.

With two arms the widths are compared, not the centres -- two backends computing
the same arithmetic should be equally unstable -- so the tests are on the
variances (F-test, and Brown-Forsythe which assumes no normality), with the mean
reported separately as a systematic offset.
"""
import sys
import math
import itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from _plotfont import setup_font, log_ticks

NPAR = 17


def _labels():
    """Korean titles where a Korean font exists, English where it does not.

    matplotlib does not fall back per glyph: with no Hangul in the active font
    every Korean character renders as an empty box and the figure ships broken.
    So the label set is chosen from what is actually installed rather than
    assumed, and the tool stays portable to machines without a Korean font.
    """
    if setup_font():
        return dict(
            scatter="run별 파라미터 RMS (epoch 0)",
            dev="arm 평균 대비 편차  [ppm]",
            rep="반복 회차",
            hist="파라미터별 흔들림 크기",
            histx=r"$\log_{10}$ (반복 간 파라미터 표준편차)",
            histy="파라미터 수",
            col="컬럼별 run-to-run 폭",
            colx="네트워크 파라미터 컬럼",
            coly="반복 간 표준편차의 중앙값",
        )
    return dict(
        scatter="Per-run parameter RMS (epoch 0)",
        dev="deviation from arm mean  [ppm]",
        rep="repetition",
        hist="Per-parameter instability",
        histx=r"$\log_{10}$ (parameter sd across repetitions)",
        histy="parameters",
        col="Run-to-run width by column",
        colx="network parameter column",
        coly="median sd across repetitions",
    )


def load_one(path):
    rows, started = [], False
    for line in open(path):
        if line.startswith("#"):
            started = "synapses weights" in line
            continue
        if not started:
            continue
        f = line.split()
        if len(f) == NPAR + 1:
            rows.append([float(v) for v in f[1:]])
    return np.asarray(rows, dtype=np.float64)


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
    return np.stack([m for m in mats if m.shape == shape])   # (nrun, nsensor, NPAR)


def sd_ci(x, conf=.95):
    n = len(x)
    if n < 2:
        return (math.nan, math.nan)
    s2 = x.var(ddof=1)
    a = (1 - conf) / 2
    return (math.sqrt((n - 1) * s2 / stats.chi2.ppf(1 - a, n - 1)),
            math.sqrt((n - 1) * s2 / stats.chi2.ppf(a, n - 1)))


def pairwise_rms(a):
    """RMS difference between every pair of runs."""
    n = a.shape[0]
    flat = a.reshape(n, -1)
    return np.array([np.sqrt(((flat[i] - flat[j]) ** 2).mean())
                     for i, j in itertools.combinations(range(n), 2)])


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

    scalars, pw = {}, {}
    for label, a in arms:
        scalars[label] = np.sqrt((a.reshape(a.shape[0], -1) ** 2).mean(axis=1))
        pw[label] = pairwise_rms(a)

    print("\n" + "=" * 78)
    print("RUN-TO-RUN SPREAD  (identical configuration, repeated)")
    print("=" * 78)
    for label, a in arms:
        s = scalars[label]
        lo, hi = sd_ci(s)
        print(f"\n--- {label}  n={a.shape[0]}")
        print(f"  per-run scalar (RMS over all parameters)")
        print(f"    mean {s.mean():.6e}   sd {s.std(ddof=1):.3e}   "
              f"relative sd {s.std(ddof=1)/s.mean():.3e}")
        print(f"    range {s.max()-s.min():.3e}   relative range "
              f"{(s.max()-s.min())/s.mean():.3e}")
        print(f"    sd 95% CI [{lo:.3e}, {hi:.3e}]")
        d = pw[label]
        print(f"  pairwise RMS difference between two runs  ({len(d)} pairs)")
        print(f"    median {np.median(d):.3e}   min {d.min():.3e}   max {d.max():.3e}")
        sd = a.std(axis=0, ddof=1)
        moved, tot = int((sd > 0).sum()), sd.size
        print(f"  per parameter")
        print(f"    differ between repetitions: {moved}/{tot} ({100*moved/tot:.1f}%)")
        if moved:
            print(f"    sd  median {np.median(sd[sd>0]):.3e}  "
                  f"p95 {np.percentile(sd[sd>0],95):.3e}  max {sd.max():.3e}")

    if len(arms) == 2:
        (la, a), (lb, b) = arms
        sa, sb = scalars[la], scalars[lb]
        Fa, Fb = sa.std(ddof=1), sb.std(ddof=1)
        F = Fa**2 / Fb**2
        dfa, dfb = len(sa) - 1, len(sb) - 1
        p_f = 2 * min(stats.f.cdf(F, dfa, dfb), 1 - stats.f.cdf(F, dfa, dfb))
        _, p_bf = stats.levene(sa, sb, center="median")
        _, p_t = stats.ttest_ind(sa, sb, equal_var=False)
        lo = math.sqrt(F / stats.f.ppf(.975, dfa, dfb))
        hi = math.sqrt(F / stats.f.ppf(.025, dfa, dfb))
        print("\n" + "-" * 78)
        print(f"WIDTH   sd({la}) / sd({lb}) = {Fa/Fb:.3f}   95% CI [{lo:.3f}, {hi:.3f}]")
        print(f"        F-test p = {p_f:.3g}   Brown-Forsythe p = {p_bf:.3g}")
        print(f"        (CI containing 1 means the two backends are equally unstable)")
        print(f"CENTRE  mean({la}) - mean({lb}) = {sa.mean()-sb.mean():+.3e}   "
              f"Welch p = {p_t:.3g}")
        print(f"CROSS   median pairwise RMS  {la} {np.median(pw[la]):.3e}   "
              f"{lb} {np.median(pw[lb]):.3e}")
    print("=" * 78)

    # ---------------------------------------------------------------- plot
    L = _labels()
    colors = ["#2b6cb0", "#c05621", "#2f855a"]
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))

    for i, (label, a) in enumerate(arms):
        s = scalars[label]
        d = (s - s.mean()) / s.mean() * 1e6
        ax[0].plot(np.arange(1, len(s) + 1), d, "o", ms=5, color=colors[i],
                   alpha=.85, label=f"{label} (n={len(s)})")
        ax[0].axhspan(-d.std(ddof=1), d.std(ddof=1), color=colors[i], alpha=.10)
    ax[0].axhline(0, color="#444", lw=.8)
    ax[0].set_xlabel(L["rep"])
    ax[0].set_ylabel(L["dev"])
    ax[0].set_title(L["scatter"])
    ax[0].legend(frameon=False, fontsize=9)
    ax[0].grid(alpha=.25)

    for i, (label, a) in enumerate(arms):
        sd = a.std(axis=0, ddof=1)
        s = sd[sd > 0]
        if s.size:
            ax[1].hist(np.log10(s), bins=60, alpha=.55, color=colors[i], label=label)
    ax[1].set_xlabel(L["histx"])
    ax[1].set_ylabel(L["histy"])
    ax[1].set_title(L["hist"])
    ax[1].legend(frameon=False, fontsize=9)
    ax[1].grid(alpha=.25)

    width = .8 / max(len(arms), 1)
    for i, (label, a) in enumerate(arms):
        sd = a.std(axis=0, ddof=1)
        med = [np.median(sd[:, c][sd[:, c] > 0]) if (sd[:, c] > 0).any() else np.nan
               for c in range(a.shape[2])]
        x = np.arange(a.shape[2]) + i * width - .4 + width / 2
        ax[2].bar(x, med, width=width, color=colors[i], label=label, alpha=.85)
    ax[2].set_yscale("log")
    log_ticks(ax[2].yaxis)
    ax[2].set_xlabel(L["colx"])
    ax[2].set_ylabel(L["coly"])
    ax[2].set_title(L["col"])
    ax[2].set_xticks(range(arms[0][1].shape[2]))
    ax[2].tick_params(axis="x", labelsize=7)
    ax[2].legend(frameon=False, fontsize=9)
    ax[2].grid(alpha=.25, axis="y")

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
