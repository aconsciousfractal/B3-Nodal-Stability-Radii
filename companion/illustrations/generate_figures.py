"""Generate the two illustrative vector figures for the nodal-face stability article.

Run `python illustrations/generate_figures.py` from the companion directory.
The support locations come from the atomic certificate. The fibre curves
sample the exact field formulas for plotting.
"""

from math import atan, cos, pi, sin, sqrt
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt


OUTPUT = Path(__file__).resolve().parent
ETA = atan(1 / 100)
OMEGA = atan(1 / 2)
V_ETA = (sin(ETA) * cos(ETA)) ** 2
DELTA = atan(sqrt(V_ETA) / 4)
SIN_DELTA = sin(DELTA)

INK = "#26333e"
MUTED = "#526371"
EDGE = "#9aaab7"
GUIDE = "#ccd5dc"
TUBE = "#f2f5f7"
WINDOW = "#f3d8a6"
AMBER = "#95600f"
BLUE = "#205b87"
RUST = "#a34225"


def save_figure(fig, name):
    fig.savefig(
        OUTPUT / name,
        metadata={"Creator": Path(__file__).name, "CreationDate": None},
    )
    plt.close(fig)


def tube_support():
    """Plot the tube chart and the actual 20 atoms of one supported jet."""
    alpha = 3 * DELTA
    jet = 1 / 10
    nodes = np.cos(np.arange(10) * pi / 9)
    ct = sqrt(1 - (SIN_DELTA * jet) ** 2)
    w = np.exp(1j * alpha) - 1
    b = (1 - SIN_DELTA**2 * jet * nodes) / ct
    q = (SIN_DELTA * (nodes - jet) / (ct * abs(w))) ** 2
    discriminant = np.sqrt(b**2 - 4 * q)
    roots = np.column_stack((b - discriminant, b + discriminant)) / 2
    offsets = np.angle(b[:, None] + roots * w) / DELTA

    fig, (overview, detail) = plt.subplots(
        2, 1, figsize=(5.906, 4.55), gridspec_kw={"height_ratios": [1, 2.05]}
    )
    fig.subplots_adjust(left=0.126, right=0.982, bottom=0.194, top=0.938, hspace=0.67)

    delta_mrad = 1000 * DELTA
    overview.fill_between([ETA, OMEGA], -delta_mrad, delta_mrad, color=TUBE, lw=0)
    overview.fill_between(
        [ETA, ETA + alpha], -delta_mrad, delta_mrad, color=WINDOW, lw=0
    )
    overview.hlines([-delta_mrad, delta_mrad], ETA, OMEGA, color=INK, lw=1.05)
    overview.vlines([ETA, OMEGA], -delta_mrad, delta_mrad, color=INK, lw=1.25)
    overview.hlines(0, ETA, OMEGA, color=BLUE, lw=1.05)
    overview.annotate(
        "Detail below",
        xy=(ETA + alpha / 2, 0.75),
        xytext=(0.062, 1.32),
        color=AMBER,
        va="center",
        arrowprops={"arrowstyle": "-", "color": AMBER, "lw": 0.85},
    )
    overview.text(0.447, 0.26, r"$r=0$", ha="right", color=BLUE)
    overview.text(0.447, 2.08, r"$r=+\delta$", ha="right", va="top", color=INK)
    overview.text(0.447, -2.08, r"$r=-\delta$", ha="right", va="bottom", color=INK)
    overview.set(
        xlim=(ETA - 0.012, OMEGA + 0.012),
        ylim=(-3.1, 3.1),
        xticks=[ETA, 0.1, 0.2, 0.3, 0.4, OMEGA],
        xticklabels=[r"$\eta$", "0.1", "0.2", "0.3", "0.4", r"$\omega$"],
        yticks=[-delta_mrad, 0, delta_mrad],
        yticklabels=[r"$-2.50$", "0", "2.50"],
        xlabel=r"Angular coordinate $\theta$ (rad)",
        ylabel=r"Normal $r$ (mrad)",
    )
    overview.set_title("(a)  Closed tube", loc="left")

    detail.axhspan(-1, 1, color=TUBE, lw=0)
    detail.axhline(0, color=BLUE, lw=0.7, alpha=0.55)
    detail.axhline(-1, color=INK, lw=0.9)
    detail.axhline(1, color=INK, lw=0.9)
    detail.axvline(0, color=INK, lw=1.05)
    detail.axvline(alpha / DELTA, color=AMBER, lw=0.9, ls="--")
    detail.hlines(nodes, offsets[:, 0], offsets[:, 1], color=GUIDE, lw=0.77)
    detail.scatter(
        offsets[::2, :].ravel(), np.repeat(nodes[::2], 2), s=19, zorder=4,
        color=BLUE, marker="o", label=r"Contact sign $+1$",
    )
    detail.scatter(
        offsets[1::2, :].ravel(), np.repeat(nodes[1::2], 2), s=21, zorder=4,
        facecolors="white", edgecolors=RUST, linewidths=0.9, marker="D",
        label=r"Contact sign $-1$",
    )
    detail.scatter(
        [0], [jet], s=22, zorder=5, color=INK, marker="+", linewidths=1.15,
        label="Evaluation point",
    )
    detail.set(
        xlim=(-0.10, 3.10),
        ylim=(-1.12, 1.12),
        xticks=[0, 0.5, 1, 1.5, 2, 2.5, 3],
        yticks=[-1, -0.5, 0, 0.5, 1],
        xlabel=r"Angular offset $(\theta-\eta)/\delta$",
        ylabel=r"Normal coordinate $s=\sin r/a$",
    )
    detail.set_title(r"(b)  Atomic support near $\eta$", loc="left")
    fig.legend(
        *detail.get_legend_handles_labels(), loc="lower center", ncols=3,
        bbox_to_anchor=(0.5, 0.02), frameon=False,
    )
    save_figure(fig, "tube_support.pdf")


def fibre_threshold():
    """Plot the target and its critical and amplified perturbations at eta."""
    t = np.linspace(-1, 1, 1601)
    q = (SIN_DELTA * t) ** 2
    # Exact target F(eta,r)/(a beta) from the factored R_eta(q).
    target = t * (1 - q) ** 2 * ((1 - q) ** 2 - q * (1 - q) / V_ETA + q**2 / V_ETA)
    # T_9(t)/9 in odd polynomial form, preserving its exact zero at t=0.
    u = t**2
    perturbation = -t * (9 + u * (-120 + u * (432 + u * (-576 + u * 256)))) / 9
    t_zero = sin(pi / 9)

    fig, (full, close) = plt.subplots(1, 2, figsize=(5.906, 2.72))
    fig.subplots_adjust(left=0.119, right=0.98, bottom=0.285, top=0.895, wspace=0.30)
    full.set_title("(a)  Full normal fibre", loc="left")
    close.set_title("(b)  Central detail", loc="left")
    full.set(xlim=(-1.025, 1.025), ylim=(-1.24, 1.24), xticks=[-1, -0.5, 0, 0.5, 1])
    close.set(
        xlim=(-0.375, 0.375), ylim=(-0.46, 0.46),
        xticks=[-t_zero, -0.2, 0, 0.2, t_zero],
        xticklabels=[r"$-t_0$", r"$-0.2$", "0", "0.2", r"$t_0$"],
        yticks=[-0.4, -0.2, 0, 0.2, 0.4],
    )
    close.axvline(-t_zero, color=EDGE, lw=0.8, ls=(0, (2, 2)), zorder=1)
    close.axvline(t_zero, color=EDGE, lw=0.8, ls=(0, (2, 2)), zorder=1)
    full.set_ylabel(r"Field divided by $a\beta$")
    for axis in (full, close):
        axis.axhline(0, color=EDGE, lw=0.8, zorder=1)
        axis.plot(t, target, color=MUTED, lw=1.35, ls=(0, (2, 1.6)), label=r"Target $F$")
        axis.plot(t, target + perturbation, color=BLUE, lw=1.65, label=r"Critical $F+E_C$")
        axis.plot(
            t, target + 2 * perturbation, color=RUST, lw=1.45,
            ls=(0, (5, 2.35)), label=r"Amplified $F+2E_C$",
        )
    close.scatter([0], [0], s=14, color=INK, zorder=5)
    fig.supxlabel(r"Normal coordinate $t=\sin r/a$", y=0.137, size=10)
    fig.legend(
        *full.get_legend_handles_labels(), loc="lower center", ncols=3,
        bbox_to_anchor=(0.5, 0.015), frameon=False,
    )
    save_figure(fig, "fibre_threshold.pdf")


def main():
    with plt.rc_context({
        "font.family": "STIXGeneral",
        "font.size": 10,
        "mathtext.fontset": "stix",
        "axes.titlesize": 10.5,
        "axes.labelsize": 10,
        "axes.labelcolor": INK,
        "axes.edgecolor": EDGE,
        "axes.linewidth": 0.65,
        "text.color": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "xtick.major.size": 2.7,
        "ytick.major.size": 2.7,
        "legend.fontsize": 10,
        "pdf.fonttype": 42,
    }):
        tube_support()
        fibre_threshold()


if __name__ == "__main__":
    main()
