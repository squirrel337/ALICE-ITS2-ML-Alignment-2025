"""Font setup shared by the run-to-run plotting scripts.

matplotlib does not fall back glyph by glyph. If the active font has no Hangul,
every Korean character in a title renders as an empty box and the figure ships
broken, so a script that wants Korean labels has to ask whether a Korean face is
actually installed rather than assume one is.

Two things have to be right at once, and they pull in opposite directions:

  - the Korean face has to be reachable for the text, and
  - it must NOT be reachable for mathtext. Log-scale ticks come out of the
    formatter as "$\\mathdefault{10^{-15}}$", which resolves through the text
    font; NanumGothic has no U+2212, so every exponent renders as 10¤15.

So the Korean face goes into the sans-serif list with DejaVu behind it rather
than into font.family directly, and any log axis should additionally use
log_ticks() below, which formats in plain ASCII and avoids mathtext entirely.
"""
import matplotlib
import matplotlib.font_manager as fm
from matplotlib.ticker import FuncFormatter

_CANDIDATES = ("NanumGothic", "NanumBarunGothic", "Noto Sans CJK KR",
               "NanumSquare", "Malgun Gothic", "AppleGothic")


def setup_font():
    """Configure a Korean-capable font if one exists. Returns True if it did."""
    for name in _CANDIDATES:
        try:
            if not fm.findfont(fm.FontProperties(family=name),
                               fallback_to_default=False):
                continue
        except Exception:
            continue
        matplotlib.rcParams["font.family"] = "sans-serif"
        matplotlib.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
        matplotlib.rcParams["mathtext.fontset"] = "dejavusans"
        matplotlib.rcParams["axes.unicode_minus"] = False
        return True
    return False


def log_ticks(axis):
    """Label a log axis as 1e-15 rather than through mathtext."""
    import numpy as np
    axis.set_major_formatter(
        FuncFormatter(lambda v, _: f"1e{int(round(np.log10(v)))}" if v > 0 else ""))
