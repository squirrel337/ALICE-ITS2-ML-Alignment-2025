#!/usr/bin/env python3
"""Compare the run-to-run spread of two (or one) sets of repeated runs.

    plot_run_to_run.py OUT.png LABEL=costs.tsv [LABEL=costs.tsv ...]

Each costs.tsv comes from tools/monitoring/run_to_run.sh and holds one line per
repetition: index, fit cost, CHSYM, status.

What is being measured: the job is run repeatedly from pristine copies with
nEPOCH=1, every repetition seeing exactly the same events from the same seed, so
any difference between repetitions is non-determinism in the module itself. The
quantity plotted is the epoch 0 training cost.

nEPOCH=0 would measure nothing at all -- the epoch loop at
YMultiLayerPerceptron.cxx:1748 is "for (iepoch = 0; iepoch < nEpoch; ...)", so at
zero it is never entered and the job exits without computing anything.

For the parameters themselves rather than this one scalar, use weight_spread.py.

The comparison of interest is the *width* of each distribution, not its centre:
two backends that compute the same thing should be equally unstable. So the
headline test is on the variances (F-test, plus Brown-Forsythe which does not
assume normality), and the mean is reported separately as a systematic offset.
"""
import sys
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from _plotfont import setup_font


def read(path):
    v = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and f[3] == "ok":
            v.append(float(f[1]))
    return np.array(v)


def sd_ci(x, conf=0.95):
    """Two-sided CI on the standard deviation, chi-square, normal assumption."""
    n = len(x)
    if n < 2:
        return (math.nan, math.nan)
    s2 = x.var(ddof=1)
    a = (1 - conf) / 2
    lo = (n - 1) * s2 / stats.chi2.ppf(1 - a, n - 1)
    hi = (n - 1) * s2 / stats.chi2.ppf(a, n - 1)
    return (math.sqrt(lo), math.sqrt(hi))


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out = sys.argv[1]
    arms = []
    for spec in sys.argv[2:]:
        label, _, path = spec.partition("=")
        x = read(path)
        if len(x) == 0:
            print(f"  {label}: no completed runs in {path}, skipped")
            continue
        arms.append((label, x))
    if not arms:
        sys.exit("nothing to plot")

    print("=" * 74)
    print(f"{'arm':10} {'n':>3} {'mean':>12} {'sd':>11} {'rel sd':>10} "
          f"{'range':>11} {'rel range':>10}")
    print("-" * 74)
    for label, x in arms:
        n, m, s = len(x), x.mean(), x.std(ddof=1)
        rng = x.max() - x.min()
        print(f"{label:10} {n:3d} {m:12.6f} {s:11.3e} {s/m:10.3e} "
              f"{rng:11.3e} {rng/m:10.3e}")
    print("-" * 74)
    for label, x in arms:
        lo, hi = sd_ci(x)
        print(f"{label:10} sd 95% CI = [{lo:.3e}, {hi:.3e}]   "
              f"(n={len(x)}, so the width is known to about "
              f"{100*(hi-lo)/(2*x.std(ddof=1)):.0f}%)")

    verdict = []
    if len(arms) == 2:
        (la, a), (lb, b) = arms
        sa, sb = a.std(ddof=1), b.std(ddof=1)
        # F-test on the variance ratio (normal assumption)
        F = sa**2 / sb**2
        dfa, dfb = len(a) - 1, len(b) - 1
        p_f = 2 * min(stats.f.cdf(F, dfa, dfb), 1 - stats.f.cdf(F, dfa, dfb))
        # Brown-Forsythe: Levene centred on the median, no normality assumed
        bf, p_bf = stats.levene(a, b, center="median")
        # Welch on the means: a systematic offset, separate from the width
        t, p_t = stats.ttest_ind(a, b, equal_var=False)
        # Two-sided CI on the variance ratio -> sd ratio
        lo = F / stats.f.ppf(0.975, dfa, dfb)
        hi = F / stats.f.ppf(0.025, dfa, dfb)
        print("-" * 74)
        print(f"width      sd({la})/sd({lb}) = {sa/sb:.3f}  "
              f"95% CI [{math.sqrt(lo):.3f}, {math.sqrt(hi):.3f}]")
        print(f"           F-test p = {p_f:.3g}    Brown-Forsythe p = {p_bf:.3g}")
        print(f"centre     mean({la}) - mean({lb}) = {a.mean()-b.mean():+.3e}   "
              f"Welch p = {p_t:.3g}")
        verdict = [
            f"sd ratio {sa/sb:.3f}, 95% CI [{math.sqrt(lo):.3f}, {math.sqrt(hi):.3f}]",
            f"F p={p_f:.3g}, Brown-Forsythe p={p_bf:.3g}",
            f"mean offset {a.mean()-b.mean():+.2e} (Welch p={p_t:.3g})",
        ]
    print("=" * 74)

    # ---------------------------------------------------------------- plot
    kr = setup_font()
    T = (lambda k, e: {"scatter":"run별 cost (epoch 0)","dist":"분포","width":"run-to-run 폭 (95% CI)"}[k] if kr else e)
    colors = ["#2b6cb0", "#c05621", "#2f855a"]
    ncol = 3 if len(arms) == 2 else 2
    fig, ax = plt.subplots(1, ncol, figsize=(5.2 * ncol, 4.4))
    if ncol == 2:
        ax = list(ax)

    # 1. every repetition, as deviation from that arm's own mean in ppm
    for i, (label, x) in enumerate(arms):
        d = (x - x.mean()) / x.mean() * 1e6
        ax[0].plot(np.arange(1, len(x) + 1), d, "o", ms=5, color=colors[i],
                   label=f"{label} (n={len(x)})", alpha=.85)
        s = d.std(ddof=1)
        ax[0].axhspan(-s, s, color=colors[i], alpha=.10)
    ax[0].axhline(0, color="#444", lw=.8)
    ax[0].set_xlabel("repetition")
    ax[0].set_ylabel("deviation from arm mean  [ppm]")
    ax[0].set_title(T("scatter","Per-run cost (epoch 0)"))
    ax[0].legend(frameon=False, fontsize=9)
    ax[0].grid(alpha=.25)

    # 2. the distributions themselves
    allc = np.concatenate([(x - x.mean()) / x.mean() * 1e6 for _, x in arms])
    bins = np.linspace(allc.min(), allc.max(), 13)
    for i, (label, x) in enumerate(arms):
        d = (x - x.mean()) / x.mean() * 1e6
        ax[1].hist(d, bins=bins, alpha=.55, color=colors[i], label=label)
    ax[1].set_xlabel("deviation from arm mean  [ppm]")
    ax[1].set_ylabel("runs")
    ax[1].set_title(T("dist","Distribution"))
    ax[1].legend(frameon=False, fontsize=9)
    ax[1].grid(alpha=.25)

    # 3. the width, with the uncertainty on the width
    if len(arms) == 2:
        for i, (label, x) in enumerate(arms):
            s = x.std(ddof=1) / x.mean() * 1e6
            lo, hi = sd_ci(x)
            lo, hi = lo / x.mean() * 1e6, hi / x.mean() * 1e6
            ax[2].errorbar([i], [s], yerr=[[s - lo], [hi - s]], fmt="o",
                           ms=9, capsize=7, color=colors[i], lw=2)
        ax[2].set_xticks(range(len(arms)))
        ax[2].set_xticklabels([l for l, _ in arms])
        ax[2].set_xlim(-.6, len(arms) - .4)
        ax[2].set_ylabel("relative sd  [ppm]")
        ax[2].set_title(T("width","Run-to-run width (95% CI)"))
        ax[2].grid(alpha=.25, axis="y")

    fig.tight_layout()
    if verdict:
        # A figure-level caption, not text inside the axes: anchored in the panel
        # it lands on top of whichever marker sits lowest.
        fig.subplots_adjust(bottom=.26)
        fig.text(.5, .04, "   |   ".join(verdict), ha="center", va="bottom",
                 fontsize=9, color="#333")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
