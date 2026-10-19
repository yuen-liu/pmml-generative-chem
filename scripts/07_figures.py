"""Stage 7 (figures).

fig1_posteriors (primary, two panels; main analysis, alpha = 1)
  (a) Posterior total-variation distance between region interaction-type profiles, mean and 95%
      CrI. Dagger: a region in the pair has prior share 7*alpha/(7*alpha + N_r) > 1%.
  (b) Interaction-type profiles (posterior mean; 50% and 95% CrI) for the four subpockets with
      the most main-set observations (rule fixed before the full run), on a shared log scale,
      against the 'other' reference (grey band = its 95% CrI, tick = its mean). Hollow points:
      n < 10 (prior-dominated).
figS1_posteriors_by_type (supplementary): all six subpockets, one panel per interaction type.

  figures/fig1_posteriors.{pdf,png}, figures/figS1_posteriors_by_type.{pdf,png}
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

from common import (ALPHAS, FIGS, OTHER, PHASES, PRIOR_DOMINATED_N, PROC,  # noqa: E402
                    REGIONS, RESULTS, SUBPOCKETS, log)

INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SERIES = "#2a78d6"      # categorical slot 1
REFERENCE = "#898781"   # neutral: 'other' is a reference band, not a series
# Sequential blue ramp, steps 100-700 (reference palette).
SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
       "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
N_DRAWS, SEED = 40_000, 7
N_PROFILE_REGIONS = 4
PRIOR_SHARE_FLAG = 0.01
LOG_FLOOR = 1e-5  # values below are drawn at the floor (zero-count cells reach ~3e-5)
INTERACTION_TYPES, INTERACTION_LABELS, SFX, COUNTS_FILE = [], {}, "", ""
SHORT_IFP = {
    "apolar": "Apolar",
    "aromatic_face_to_face": "Arom. face-face",
    "aromatic_edge_to_face": "Arom. edge-face",
    "hbond_donor_protein": "H-bond, prot. donor",
    "hbond_acceptor_protein": "H-bond, prot. acceptor",
    "cation_protein_anion_ligand": "Prot. cation : lig. anion",
    "anion_protein_cation_ligand": "Prot. anion : lig. cation",
}


def set_phase(phase):
    """Bind the vocabulary, labels and output suffix for this phase (module-level, read by panels)."""
    global INTERACTION_TYPES, INTERACTION_LABELS, SHORT, SFX, COUNTS_FILE
    cfg = PHASES[phase]
    INTERACTION_TYPES, INTERACTION_LABELS, SFX = cfg["vocab"], cfg["labels"], cfg["suffix"]
    COUNTS_FILE = cfg["counts"]
    SHORT = SHORT_IFP if phase == "ifp" else dict(cfg["labels"])


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.edgecolor": INK_2,
                         "axes.labelcolor": INK_2, "xtick.color": INK_2, "ytick.color": INK_2,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                         "xtick.minor.width": 0.4})


def load(analysis, prior):
    alpha = ALPHAS[prior]
    long = pd.read_csv(PROC / COUNTS_FILE)
    g = long[long.analysis == analysis].pivot(index="region", columns="interaction", values="n")
    g = g.reindex(index=REGIONS, columns=INTERACTION_TYPES).fillna(0).astype(int)
    rng = np.random.default_rng(SEED)
    q = {r: np.quantile(rng.dirichlet(alpha + g.loc[r].to_numpy(float), size=N_DRAWS),
                        [0.025, 0.25, 0.75, 0.975], axis=0) for r in REGIONS}
    mean = {r: (alpha + g.loc[r].to_numpy()) / (len(INTERACTION_TYPES) * alpha + g.loc[r].sum()) for r in REGIONS}
    return alpha, g, q, mean


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(FIGS / f"{stem}.{ext}", dpi=300, facecolor=SURFACE)
    plt.close(fig)
    log.info("wrote figures/%s.{pdf,png}", stem)


def panel_tv(ax, analysis, prior, g, alpha):
    tv = pd.read_csv(RESULTS / f"pairwise_tv{SFX}.csv")
    tv = tv[(tv.analysis == analysis) & (tv.prior == prior)]
    V = len(INTERACTION_TYPES)
    share = {r: V * alpha / (V * alpha + g.loc[r].sum()) for r in REGIONS}
    n = len(REGIONS)
    vmax = max(0.05, float(np.ceil(tv["tv_mean"].max() * 20) / 20))
    cmap = LinearSegmentedColormap.from_list("seq", SEQ)
    ax.set_facecolor(SURFACE)
    for _, row in tv.iterrows():
        i, j = REGIONS.index(row.region_b), REGIONS.index(row.region_a)  # row i > column j
        if i < j:
            i, j = j, i
        frac = row.tv_mean / vmax
        ax.add_patch(plt.Rectangle((j - 0.5 + 0.03, i - 0.5 + 0.03), 0.94, 0.94, color=cmap(frac), lw=0))
        dagger = "†" if max(share[row.region_a], share[row.region_b]) > PRIOR_SHARE_FLAG else ""
        txt_color = SURFACE if frac > 0.55 else INK
        ax.text(j, i - 0.12, f"{row.tv_mean:.2f}{dagger}", ha="center", va="center", fontsize=6.3,
                color=txt_color, fontweight="bold")
        ax.text(j, i + 0.22, f"{row.cri_lo:.2f}–{row.cri_hi:.2f}", ha="center", va="center",
                fontsize=4.6, color=txt_color)
    labels = [f"{r}\n{g.loc[r].sum():,}" for r in REGIONS]
    ax.set_xticks(range(n - 1), labels[:-1], fontsize=6)
    ax.set_yticks(range(1, n), labels[1:], fontsize=6)
    ax.set_xlim(-0.5, n - 1.5)
    ax.set_ylim(n - 0.5, 0.5)
    ax.tick_params(length=0, pad=2)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("(a) Posterior TV distance between profiles", loc="left", fontsize=7.5, color=INK, pad=4)
    ax.text(n - 1.55, 0.6, f"cell: mean, 95% CrI\naxis: region, N obs.\n† prior share > {PRIOR_SHARE_FLAG:.0%}",
            ha="right", va="top", fontsize=5.5, color=INK_2)


def panel_profiles(axes, g, q, mean):
    top = g.loc[SUBPOCKETS].sum(axis=1).sort_values(ascending=False).index[:N_PROFILE_REGIONS].tolist()
    y = np.arange(len(INTERACTION_TYPES))
    for c, (ax, r) in enumerate(zip(axes, top)):
        ax.set_facecolor(SURFACE)
        for k, it in enumerate(INTERACTION_TYPES):
            lo = max(q[OTHER][0][k], LOG_FLOOR)
            ax.fill_betweenx([k - 0.32, k + 0.32], lo, q[OTHER][3][k], color=REFERENCE, alpha=0.22, lw=0, zorder=0)
            ax.vlines(max(mean[OTHER][k], LOG_FLOOR), k - 0.32, k + 0.32, color=REFERENCE, lw=1, zorder=1)
            ax.hlines(k, max(q[r][0][k], LOG_FLOOR), q[r][3][k], color=SERIES, lw=1, zorder=2)
            ax.hlines(k, max(q[r][1][k], LOG_FLOOR), max(q[r][2][k], LOG_FLOOR), color=SERIES, lw=3, zorder=2)
            small = g.loc[r, it] < PRIOR_DOMINATED_N
            ax.plot(max(mean[r][k], LOG_FLOOR), k, "o", ms=4, mec=SERIES, mew=1,
                    mfc=SURFACE if small else SERIES, zorder=3)
        ax.set_xscale("log")
        ax.set_xlim(LOG_FLOOR / 2, 2.0)
        ax.set_xticks([1e-4, 1e-1], ["10⁻⁴", "10⁻¹"], fontsize=6)
        ax.set_xticks([1e-5, 1e-3, 1e-2, 1], minor=True)
        ax.tick_params(axis="x", which="minor", length=1.5, labelbottom=False)
        ax.set_ylim(len(INTERACTION_TYPES) - 0.5, -0.5)
        ax.set_yticks(y, [SHORT[it] for it in INTERACTION_TYPES] if c == 0 else [])
        ax.grid(axis="x", which="major", color=GRID, lw=0.5, zorder=0)
        ax.tick_params(axis="y", length=0, pad=2)
        ax.tick_params(axis="x", length=2, pad=1.5)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.set_title(f"{r}\nN = {g.loc[r].sum():,}", fontsize=6.5, color=INK, pad=3, linespacing=1.1)
    axes[0].text(0, 1.13, f"(b) Profiles, {N_PROFILE_REGIONS} best-populated subpockets",
                 transform=axes[0].transAxes, fontsize=7.5, color=INK, ha="left", va="bottom")
    fig = axes[0].figure
    fig.text(0.985, 0.02, "posterior proportion (log; values < 10⁻⁵ drawn at 10⁻⁵)   "
             "● mean, 50% / 95% CrI   ○ n < 10   ▮ 'other' reference, 95% CrI",
             ha="right", va="bottom", fontsize=5.3, color=INK_2)
    return top


def fig1(analysis, prior):
    alpha, g, q, mean = load(analysis, prior)
    style()
    fig = plt.figure(figsize=(7.0, 3.3), facecolor=SURFACE)
    gs = GridSpec(1, 1 + N_PROFILE_REGIONS, figure=fig, width_ratios=[3.0] + [1.0] * N_PROFILE_REGIONS,
                  left=0.06, right=0.985, top=0.86, bottom=0.12, wspace=0.12)
    panel_tv(fig.add_subplot(gs[0, 0]), analysis, prior, g, alpha)
    sub = GridSpec(1, N_PROFILE_REGIONS, figure=fig, left=0.59, right=0.985, top=0.84, bottom=0.12, wspace=0.18)
    top = panel_profiles([fig.add_subplot(sub[0, c]) for c in range(N_PROFILE_REGIONS)], g, q, mean)
    stem = f"fig1_posteriors{SFX}" if (analysis, prior) == ("main", "uniform") else \
        f"fig1_posteriors{SFX}_{analysis.replace(':', '_')}_{prior}"
    save(fig, stem)
    log.info("panel (b) subpockets (most observations): %s", top)


def figS1(analysis, prior):
    """Supplementary: all six subpockets, one panel per interaction type (independent y-scales)."""
    alpha, g, q, mean = load(analysis, prior)
    style()
    fig, axes = plt.subplots(2, 4, figsize=(7.0, 3.5), facecolor=SURFACE)
    x = np.arange(len(SUBPOCKETS))
    for k, (ax, it) in enumerate(zip(axes.flat, INTERACTION_TYPES)):
        ax.set_facecolor(SURFACE)
        ax.axhspan(q[OTHER][0][k], q[OTHER][3][k], color=REFERENCE, alpha=0.22, lw=0, zorder=0)
        ax.axhline(mean[OTHER][k], color=REFERENCE, lw=1, zorder=1)
        for i, r in enumerate(SUBPOCKETS):
            ax.vlines(i, q[r][0][k], q[r][3][k], color=SERIES, lw=1, zorder=2)
            ax.vlines(i, q[r][1][k], q[r][2][k], color=SERIES, lw=3, zorder=2)
            small = g.loc[r, it] < PRIOR_DOMINATED_N
            ax.plot(i, mean[r][k], "o", ms=4.5, zorder=3, mec=SERIES, mew=1, mfc=SURFACE if small else SERIES)
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
    fig.tight_layout(w_pad=0.8, h_pad=1.0)
    save(fig, f"figS1_posteriors_by_type{SFX}" if (analysis, prior) == ("main", "uniform") else
         f"figS1_posteriors_by_type{SFX}_{analysis.replace(':', '_')}_{prior}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--analysis", default="main")
    ap.add_argument("--prior", default="uniform", choices=list(ALPHAS))
    ap.add_argument("--phase", default="ifp", choices=list(PHASES))
    a = ap.parse_args()
    set_phase(a.phase)
    fig1(a.analysis, a.prior)
    figS1(a.analysis, a.prior)
