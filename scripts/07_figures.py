"""Stage 7 (figures): primary figure, posterior interaction-type proportions by subpocket.

One panel per interaction type; x = KinFragLib subpocket; point = posterior mean, thick bar =
50% CrI, thin bar = 95% CrI (main analysis, alpha = 1). The grey band is the 'other' reference
distribution's 95% CrI, so each subpocket is read against the contrast class. Hollow points mark
cells with n < 10 (prior-dominated). Each panel has its own y-scale starting at 0.

  figures/fig1_posteriors.pdf / .png
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import (ALPHAS, FIGS, INTERACTION_LABELS, INTERACTION_TYPES, OTHER, PRIOR_DOMINATED_N, PROC,  # noqa: E402
                    RESULTS, SUBPOCKETS, log)

INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SERIES = "#2a78d6"      # categorical slot 1
REFERENCE = "#898781"   # neutral: 'other' is a reference band, not a series
N_DRAWS, SEED = 40_000, 7


def draws_for(counts, alpha, rng):
    return rng.dirichlet(alpha + counts.astype(float), size=N_DRAWS)


def main(analysis, prior):
    alpha = ALPHAS[prior]
    long = pd.read_csv(PROC / "counts_long.csv")
    g = long[long.analysis == analysis].pivot(index="region", columns="interaction", values="n")
    g = g.reindex(columns=INTERACTION_TYPES).fillna(0).astype(int)
    rng = np.random.default_rng(SEED)
    q = {}
    for r in SUBPOCKETS + [OTHER]:
        d = draws_for(g.loc[r].to_numpy(), alpha, rng)
        q[r] = np.quantile(d, [0.025, 0.25, 0.75, 0.975], axis=0)
    post = pd.read_csv(RESULTS / "posterior_summary.csv")
    post = post[(post.analysis == analysis) & (post.prior == prior)].set_index(["region", "interaction"])

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.edgecolor": INK_2,
                         "axes.labelcolor": INK_2, "xtick.color": INK_2, "ytick.color": INK_2,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
    fig, axes = plt.subplots(2, 4, figsize=(7.0, 3.5), facecolor=SURFACE)
    x = np.arange(len(SUBPOCKETS))
    for k, (ax, it) in enumerate(zip(axes.flat, INTERACTION_TYPES)):
        ax.set_facecolor(SURFACE)
        lo95, lo50, hi50, hi95 = (q[OTHER][j][k] for j in range(4))
        ax.axhspan(lo95, hi95, color=REFERENCE, alpha=0.22, lw=0, zorder=0)
        ax.axhline(post.loc[(OTHER, it), "post_mean"], color=REFERENCE, lw=1, zorder=1)
        for i, r in enumerate(SUBPOCKETS):
            ax.vlines(i, q[r][0][k], q[r][3][k], color=SERIES, lw=1, zorder=2)
            ax.vlines(i, q[r][1][k], q[r][2][k], color=SERIES, lw=3, zorder=2)
            small = post.loc[(r, it), "n"] < PRIOR_DOMINATED_N
            ax.plot(i, post.loc[(r, it), "post_mean"], "o", ms=4.5, zorder=3, mec=SERIES, mew=1,
                    mfc=SURFACE if small else SERIES)
        ax.set_title(INTERACTION_LABELS[it], fontsize=7, color=INK, loc="left", pad=3)
        ax.set_xticks(x, SUBPOCKETS)
        ax.set_xlim(-0.6, len(SUBPOCKETS) - 0.4)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
        ax.tick_params(length=2, pad=1.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        if k % 4 == 0:
            ax.set_ylabel("Posterior proportion")

    key = axes.flat[-1]
    key.axis("off")
    key.vlines(0.08, 0.62, 0.92, color=SERIES, lw=1)
    key.vlines(0.08, 0.70, 0.84, color=SERIES, lw=3)
    key.plot(0.08, 0.77, "o", ms=4.5, color=SERIES)
    key.text(0.2, 0.77, "Subpocket: mean,\n50% and 95% CrI", va="center", color=INK, fontsize=6.5)
    key.plot(0.08, 0.47, "o", ms=4.5, mec=SERIES, mfc=SURFACE)
    key.text(0.2, 0.47, f"n < {PRIOR_DOMINATED_N}: prior-dominated", va="center", color=INK, fontsize=6.5)
    key.add_patch(plt.Rectangle((0.02, 0.17), 0.12, 0.12, color=REFERENCE, alpha=0.22, lw=0))
    key.plot([0.02, 0.14], [0.23, 0.23], color=REFERENCE, lw=1)
    key.text(0.2, 0.23, "'other' positions:\nreference, 95% CrI", va="center", color=INK, fontsize=6.5)
    key.set_xlim(0, 1)
    key.set_ylim(0, 1)

    n_struct = pd.read_parquet(PROC / "structures_main.parquet").shape[0]
    fig.tight_layout(w_pad=0.8, h_pad=1.0)
    stem = FIGS / ("fig1_posteriors" if (analysis, prior) == ("main", "uniform")
                   else f"fig_posteriors_{analysis.replace(':', '_')}_{prior}")
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=300, facecolor=SURFACE)
    log.info("wrote %s.{pdf,png} (analysis=%s, prior=%s, alpha=%s; main set n=%d before assignment drops)",
             stem.relative_to(FIGS.parent), analysis, prior, alpha, n_struct)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--analysis", default="main")
    ap.add_argument("--prior", default="uniform", choices=list(ALPHAS))
    a = ap.parse_args()
    main(a.analysis, a.prior)
