"""Stage 11 (figures): Fig. 1, the primary cross-phase figure.

One model, fit twice; the figure shows the result and its two caveats in one view.

(a) Phase 1 (IFP interaction types) and (b) phase 2 (BRICS fragment chemistry): posterior
    total-variation distance between the profiles of all 15 subpocket pairs, against each pair's
    same-profile null. Per row: grey span = the null's 97.5% point (TV not resolved above
    same-profile noise at that pair's effective sample size N/R); thin line = the 95% CrI widened
    about its mean by sqrt(R); thick line = the unwidened 95% CrI; dot = posterior mean, hollow
    where the widened lower bound does not clear the null. Rows share the phase 1 ordering so the
    two panels align, but the x-axes are NOT comparable: the vocabularies differ (7 interaction
    types vs 4 chemical categories), which puts TV on different scales.
(c) The FP-SE pair resolved: per-category posterior differences pi_FP - pi_SE (phase 2), mean and
    95% CrI, with each category's share of the pair's TV.

  figures/fig1_main.{pdf,png}
"""
import argparse
import importlib

import numpy as np
import pandas as pd

_f7 = importlib.import_module("07_figures")
plt, GridSpec = _f7.plt, _f7.GridSpec
INK, INK_2, GRID, SURFACE = _f7.INK, _f7.INK_2, _f7.GRID, _f7.SURFACE
SERIES, REFERENCE = _f7.SERIES, _f7.REFERENCE
style, save = _f7.style, _f7.save

from common import ALPHAS, FRAGMENT_LABELS, PHASES, RESULTS, log  # noqa: E402

HILITE = "#f2e6c9"   # row highlight for the marginal pair (neutral, not a series colour)
PAIR = ("FP", "SE")
PHASE_TITLES = {"ifp": "(a) Phase 1: IFP interaction types",
                "fragment": "(b) Phase 2: BRICS fragment chemistry"}


def load_null(phase, analysis, prior):
    d = pd.read_csv(RESULTS / f"pairwise_tv_null{PHASES[phase]['suffix']}.csv")
    d = d[(d.analysis == analysis) & (d.prior == prior) & d.both_subpockets].copy()
    d["pair"] = d.region_a + "-" + d.region_b
    return d.set_index("pair")


def panel_pairs(ax, d, order, title, show_labels):
    """One phase: 15 pair rows, TV intervals against each pair's null."""
    xmax = max(0.05, float(np.ceil(d["cri_hi_inflated"].max() * 20) / 20))
    for i, pair in enumerate(order):
        r = d.loc[pair]
        if pair == "-".join(PAIR):
            ax.axhspan(i - 0.46, i + 0.46, color=HILITE, lw=0, zorder=0)
        # The null: TV at or below this is not resolved above same-profile noise at N/R.
        ax.add_patch(plt.Rectangle((0, i - 0.3), r.null_eff_q975, 0.6,
                                   color=REFERENCE, alpha=0.3, lw=0, zorder=1))
        ax.hlines(i, r.cri_lo_inflated, r.cri_hi_inflated, color=SERIES, lw=1, alpha=0.45, zorder=2)
        ax.hlines(i, r.cri_lo, r.cri_hi, color=SERIES, lw=2.4, zorder=3)
        ax.plot(r.tv_mean, i, "o", ms=4.2, zorder=4, color=SERIES if r.distinguishable else SURFACE,
                mec=SERIES, mew=0.9)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlim(0, xmax)
    ax.set_yticks(range(len(order)))
    if show_labels:
        ax.set_yticklabels(order, fontsize=6.3, color=INK)
        for t in ax.get_yticklabels():
            if t.get_text() == "-".join(PAIR):
                t.set_fontweight("bold")
    else:
        ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Total-variation distance between profiles", fontsize=6.5, color=INK_2)
    n_dist = int(d["distinguishable"].sum())
    ax.set_title(f"{title}\n{n_dist} of {len(order)} pairs resolved above the null, $R$ = {d.R.mean():.2f}",
                 loc="left", fontsize=7.5, color=INK, pad=4, linespacing=1.45)
    ax.title.set_color(INK)
    ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)


def legend(fig, y):
    """Figure-level key, in the gap between the pair panels and panel (c)."""
    ax = fig.add_axes([0.085, y, 0.9, 0.055]); ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.add_patch(plt.Rectangle((0.0, 0.3), 0.055, 0.42, color=REFERENCE, alpha=0.3, lw=0))
    ax.text(0.068, 0.5, "same-profile null, 97.5% point, simulated at $N/R$", va="center",
            fontsize=6.2, color=INK_2)
    ax.hlines(0.5, 0.40, 0.475, color=SERIES, lw=1, alpha=0.45)
    ax.hlines(0.5, 0.418, 0.457, color=SERIES, lw=2.4)
    ax.plot(0.4375, 0.5, "o", ms=4.2, color=SERIES)
    ax.text(0.488, 0.5, "posterior mean, 95% CrI, and the CrI widened by $\\sqrt{R}$",
            va="center", fontsize=6.2, color=INK_2)
    ax.plot(0.862, 0.5, "o", ms=4.2, color=SURFACE, mec=SERIES, mew=0.9)
    ax.text(0.878, 0.5, "not resolved", va="center", fontsize=6.2, color=INK_2)


def panel_pair_detail(ax, prior, tv_frag):
    """Panel (c): what resolves FP-SE, per category (phase 2)."""
    d = pd.read_csv(RESULTS / "pairwise_diffs_fragments.csv")
    d = d[(d.analysis == "main") & (d.prior == prior)
          & (d.region_a == PAIR[0]) & (d.region_b == PAIR[1])].copy()
    d["contrib"] = 0.5 * d.diff_mean.abs()
    d["share"] = d.contrib / d.contrib.sum()
    d = d.sort_values("contrib")                       # largest contributor at the top
    ax.axvline(0, color=INK_2, lw=0.7, zorder=2)
    for i, (_, r) in enumerate(d.iterrows()):
        ax.hlines(i, r.cri_lo, r.cri_hi, color=SERIES, lw=2.4, zorder=3)
        ax.plot(r.diff_mean, i, "o", ms=4.2, color=SERIES, zorder=4)
        ax.text(0.985, i / max(len(d) - 1, 1) * 0 + i, f"{r.share:.0%} of TV",
                transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=5.8, color=INK_2)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([FRAGMENT_LABELS[v] for v in d.interaction], fontsize=6.3, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, len(d) - 0.4)
    m = float(np.abs(d[["cri_lo", "cri_hi"]].to_numpy()).max()) * 1.45
    ax.set_xlim(-m, m)
    ax.set_xlabel(f"$\\leftarrow$ more in {PAIR[1]}          posterior difference  "
                  f"$\\pi_{{{PAIR[0]}}} - \\pi_{{{PAIR[1]}}}$  (fragment share)"
                  f"          more in {PAIR[0]} $\\rightarrow$", fontsize=6.5, color=INK_2)
    tv = tv_frag.loc["-".join(PAIR)]
    ax.set_title(f"(c) What resolves {PAIR[0]}-{PAIR[1]} (phase 2): TV = {tv.tv_mean:.3f}, "
                 f"no single category dominates", loc="left", fontsize=7.5, color=INK, pad=4)
    ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    return d


def fig1(analysis, prior):
    ifp, frag = load_null("ifp", analysis, prior), load_null("fragment", analysis, prior)
    order = list(ifp.sort_values("tv_mean").index)   # smallest TV (FP-SE) at the bottom
    missing = set(order) ^ set(frag.index)
    if missing:
        raise SystemExit(f"pair sets differ between phases: {sorted(missing)}")
    style()
    fig = plt.figure(figsize=(7.0, 4.9), facecolor=SURFACE)
    top = GridSpec(1, 2, figure=fig, left=0.085, right=0.985, top=0.90, bottom=0.44, wspace=0.075)
    panel_pairs(fig.add_subplot(top[0, 0]), ifp, order, PHASE_TITLES["ifp"], True)
    panel_pairs(fig.add_subplot(top[0, 1]), frag, order, PHASE_TITLES["fragment"], False)
    legend(fig, 0.295)
    bot = GridSpec(1, 1, figure=fig, left=0.145, right=0.985, top=0.235, bottom=0.10)
    det = panel_pair_detail(fig.add_subplot(bot[0, 0]), prior, frag)
    save(fig, "fig1_main" if (analysis, prior) == ("main", "uniform")
         else f"fig1_main_{analysis.replace(':', '_')}_{prior}")

    log.info("row order (phase 1 TV, descending): %s", order)
    for name, d in (("phase 1", ifp), ("phase 2", frag)):
        log.info("%s: %d of %d resolved; R = %.2f; unresolved: %s", name,
                 int(d.distinguishable.sum()), len(d), d.R.mean(),
                 list(d.index[~d.distinguishable]))
    log.info("panel (c) shares of TV:\n%s",
             det[["interaction", "diff_mean", "cri_lo", "cri_hi", "share"]].round(4).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--analysis", default="main")
    ap.add_argument("--prior", default="uniform", choices=list(ALPHAS))
    a = ap.parse_args()
    fig1(a.analysis, a.prior)
