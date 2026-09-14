#!/usr/bin/env python3
# FIGURE GENERATION FOR THE CHEBYSHEV/SMOLYAK METHODS WRITE-UP.
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

sys.path.insert(0, "/Users/Huawei/Quantitative_Model/code/bocola2016")
from smolyak import SmolyakGrid, _level_points, chebyshev_basis_1d  # faithful grid

OUT = "/private/tmp/claude-501/-Users-Huawei-Quantitative-Model/d4031c6e-4ebf-4ea1-9802-95747d474540/scratchpad/"

plt.rcParams.update({
    "font.size": 11, "axes.grid": True, "grid.alpha": 0.25,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.bbox": "tight",
})
BLUE, RED, GREEN, GREY = "#2c5aa0", "#c0392b", "#218c5a", "#7f8c8d"


def runge():
    # RUNGE PHENOMENON: EQUISPACED vs CHEBYSHEV INTERPOLATION.
    f = lambda x: 1.0 / (1.0 + 25.0 * x ** 2)
    xx = np.linspace(-1, 1, 1000)
    n = 14
    # equispaced nodes
    xe = np.linspace(-1, 1, n + 1)
    ce = np.polyfit(xe, f(xe), n)
    ye = np.polyval(ce, xx)
    # chebyshev-extrema nodes (same count)
    xc = -np.cos(np.pi * np.arange(n + 1) / n)
    cc = np.polyfit(xc, f(xc), n)
    yc = np.polyval(cc, xx)

    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.6), sharey=True)
    for a in ax:
        a.plot(xx, f(xx), color=GREY, lw=2.2, label=r"$f(x)=1/(1+25x^2)$")
    ax[0].plot(xx, ye, color=RED, lw=1.6, label=f"degree-{n} interpolant")
    ax[0].plot(xe, f(xe), "o", color=RED, ms=5, mfc="white")
    ax[0].set_title("Equispaced nodes  —  Runge blow-up")
    ax[0].set_ylim(-0.6, 1.4)
    ax[1].plot(xx, yc, color=BLUE, lw=1.6, label=f"degree-{n} interpolant")
    ax[1].plot(xc, f(xc), "o", color=BLUE, ms=5, mfc="white")
    ax[1].set_title("Chebyshev–Lobatto nodes  —  stable")
    for a in ax:
        a.set_xlabel("$x$"); a.legend(loc="upper center", fontsize=8.5, framealpha=0.9)
    ax[0].set_ylabel("$f,\\; p_n$")
    fig.tight_layout()
    fig.savefig(OUT + "runge.pdf")
    plt.close(fig)


def nodes():
    # CHEBYSHEV NODES AS PROJECTED EQUISPACED POINTS ON THE SEMICIRCLE.
    n = 12
    th = np.pi * np.arange(n + 1) / n
    x = -np.cos(th)
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    tt = np.linspace(0, np.pi, 300)
    ax.plot(-np.cos(tt), np.sin(tt), color=GREY, lw=1.4)
    for xi, ti in zip(x, th):
        ax.plot([xi, xi], [0, np.sin(ti)], color=BLUE, lw=0.8, alpha=0.6)
        ax.plot(xi, np.sin(ti), "o", color=BLUE, ms=5)
        ax.plot(xi, 0, "o", color=RED, ms=6)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_title("Chebyshev–Lobatto nodes $x_k=-\\cos(\\pi k/n)$ "
                 "= equispaced angles projected down")
    ax.set_xlabel("$x$"); ax.set_yticks([])
    ax.set_ylim(-0.12, 1.15); ax.set_aspect("equal")
    ax.grid(False)
    fig.tight_layout(); fig.savefig(OUT + "nodes.pdf"); plt.close(fig)


def smolyak_grid():
    # 2D: FULL TENSOR PRODUCT vs SMOLYAK SPARSE GRIDS (FAITHFUL SmolyakGrid).
    g2 = SmolyakGrid([-1, -1], [1, 1], mu=2)
    g3 = SmolyakGrid([-1, -1], [1, 1], mu=3)
    # tensor product of level-4 1D extrema (9 points) -> 81 nodes, comparable degree
    p = _level_points(4)
    TX, TY = np.meshgrid(p, p)
    fig, ax = plt.subplots(1, 3, figsize=(10.2, 3.5))
    ax[0].plot(TX.ravel(), TY.ravel(), "o", color=GREY, ms=4)
    ax[0].set_title(f"Full tensor grid\n$9\\times 9={TX.size}$ points")
    ax[1].plot(g2.points_unit[:, 0], g2.points_unit[:, 1], "o", color=BLUE, ms=5)
    ax[1].set_title(f"Smolyak $\\mu=2$\n{g2.n} points")
    ax[2].plot(g3.points_unit[:, 0], g3.points_unit[:, 1], "o", color=RED, ms=5)
    ax[2].set_title(f"Smolyak $\\mu=3$\n{g3.n} points")
    for a in ax:
        a.set_xlim(-1.15, 1.15); a.set_ylim(-1.15, 1.15); a.set_aspect("equal")
        a.set_xlabel("$x_1$"); a.grid(alpha=0.2)
    ax[0].set_ylabel("$x_2$")
    fig.tight_layout(); fig.savefig(OUT + "smolyak2d.pdf"); plt.close(fig)


def growth():
    # NODE-COUNT GROWTH: TENSOR vs SMOLYAK ACROSS DIMENSION.
    dims = np.arange(1, 11)
    tensor = 5.0 ** dims                     # 5 pts / dim (crude, illustrative)
    smol2 = [SmolyakGrid([-1]*d, [1]*d, mu=2).n for d in dims]
    smol3 = [SmolyakGrid([-1]*d, [1]*d, mu=3).n for d in dims]
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.semilogy(dims, tensor, "o-", color=GREY, label=r"tensor product ($5^d$)")
    ax.semilogy(dims, smol3, "s-", color=RED, label=r"Smolyak $\mu=3$")
    ax.semilogy(dims, smol2, "^-", color=BLUE, label=r"Smolyak $\mu=2$")
    ax.axvline(6, color="k", ls=":", lw=1)
    ax.text(6.1, 3e4, "model: $d=6$", fontsize=9)
    ax.set_xlabel("state dimension $d$"); ax.set_ylabel("collocation nodes")
    ax.set_title("Sparse grids defeat the tensor-product curse of dimensionality")
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(OUT + "growth.pdf"); plt.close(fig)


def convergence():
    # SPECTRAL CONVERGENCE ON CHEBYSHEV NODES (ANALYTIC vs KINKED TARGET).
    xx = np.linspace(-1, 1, 2000)
    analytic = lambda x: np.exp(np.sin(3 * x))          # entire -> geometric
    kinked = lambda x: np.abs(x - 0.2)                   # kink -> algebraic
    degs = np.arange(2, 41, 2)
    ea, ek = [], []
    for n in degs:
        xc = -np.cos(np.pi * np.arange(n + 1) / n)
        for f, store in ((analytic, ea), (kinked, ek)):
            c = np.linalg.solve(chebyshev_basis_1d(xc, n), f(xc))
            approx = chebyshev_basis_1d(xx, n) @ c
            store.append(np.max(np.abs(approx - f(xx))))
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.semilogy(degs, ea, "o-", color=BLUE, label="analytic $e^{\\sin 3x}$ (spectral)")
    ax.semilogy(degs, ek, "s-", color=RED, label="kinked $|x-0.2|$ (algebraic)")
    ax.set_xlabel("polynomial degree $n$")
    ax.set_ylabel(r"$\max_x |f-p_n|$")
    ax.set_title("Smoothness governs the rate: geometric vs algebraic")
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(OUT + "convergence.pdf"); plt.close(fig)


def branch_tree():
    # SCHEMATIC OF THE TWO-BRANCH QUADRATURE + LOGISTIC DEFAULT PROBABILITY.
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.8),
                           gridspec_kw={"width_ratios": [1.5, 1]})
    a = ax[0]; a.axis("off"); a.set_xlim(0, 10); a.set_ylim(0, 10)

    def box(x, y, w, h, text, fc):
        a.add_patch(plt.Rectangle((x, y), w, h, fc=fc, ec="k", lw=1.2, alpha=0.9))
        a.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9)

    def arrow(x0, y0, x1, y1, txt="", col="k"):
        a.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                    mutation_scale=12, lw=1.2, color=col))
        if txt:
            a.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.25, txt, fontsize=8.5,
                   ha="center", color=col)

    box(0.2, 4.2, 2.2, 1.6, "state $S_t$\n$(K,B,P,\\Delta z,g,s)$", "#dce6f2")
    box(3.4, 4.2, 2.3, 1.6, "GH nodes\n$3^3=27$\n$(\\epsilon_z,\\epsilon_g,\\epsilon_s)$", "#eae5f2")
    arrow(2.4, 5.0, 3.4, 5.0, "AR(1)")
    box(7.0, 6.7, 2.7, 1.5, "no default $d'{=}0$\nfull payoff", "#d7ede0")
    box(7.0, 1.9, 2.7, 1.5, "default $d'{=}1$\nhaircut $1-D$", "#f2dcdc")
    arrow(5.7, 5.4, 7.0, 7.4, "$1-p^d(s)$", GREEN)
    arrow(5.7, 4.6, 7.0, 2.6, "$p^d(s)$", RED)
    a.text(5.0, 0.6, r"$\mathbb{E}_t[\cdot]=\sum_j w_j\,[(1-p^d)\,(\cdot)^{d'=0}"
                     r"+p^d\,(\cdot)^{d'=1}]$", fontsize=9, ha="center")
    a.set_title("Two-branch quadrature (expectations_full.py)", fontsize=11)

    b = ax[1]
    s = np.linspace(-11, -3, 400)
    b.plot(s, 1.0 / (1.0 + np.exp(-s)), color=RED, lw=2)
    for sv, lab in [(-7.06, "$s^\\star$"), (-3.66, "stress")]:
        b.axvline(sv, color=GREY, ls=":", lw=1)
        b.text(sv + 0.1, 0.5, lab, fontsize=8, rotation=90, va="center")
    b.set_xlabel("risk factor $s$")
    b.set_ylabel("$p^d(s)=1/(1+e^{-s})$")
    b.set_title("Priced default probability")
    fig.subplots_adjust(left=0.02, right=0.95, bottom=0.14, top=0.90, wspace=0.25)
    # the axis("off") schematic confuses the "tight" bbox estimator; save the
    # figure at its true 10x3.8in extent instead of cropping to artist bounds
    old = matplotlib.rcParams["savefig.bbox"]
    matplotlib.rcParams["savefig.bbox"] = "standard"
    fig.savefig(OUT + "branch.pdf")
    matplotlib.rcParams["savefig.bbox"] = old
    plt.close(fig)


if __name__ == "__main__":
    runge(); nodes(); smolyak_grid(); growth(); convergence(); branch_tree()
    print("figures written to", OUT)
