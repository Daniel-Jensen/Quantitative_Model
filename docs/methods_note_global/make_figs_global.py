#!/usr/bin/env python3
# FIGURES FOR THE code/global RECURSIVE-SOLVER METHODS NOTE.
# Grids are the ACTUAL solver_recursive/state_grid.SmolyakGrid so every point
# count and node layout is faithful to the running code.
import importlib.util as _u
import itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

# import the real grid module by path (no package side effects)
_spec = _u.spec_from_file_location(
    "sg", "/Users/Huawei/Quantitative_Model/code/global/solver_recursive/state_grid.py")
sg = _u.module_from_spec(_spec); _spec.loader.exec_module(sg)

OUT = "/private/tmp/claude-501/-Users-Huawei-Quantitative-Model/d4031c6e-4ebf-4ea1-9802-95747d474540/scratchpad/"

plt.rcParams.update({
    "font.size": 11, "axes.grid": True, "grid.alpha": 0.25,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 160, "savefig.bbox": "tight",
})
BLUE, RED, GREEN, GREY, PURP, ORAN = ("#2c5aa0", "#c0392b", "#218c5a",
                                      "#7f8c8d", "#7b4ea3", "#d98c1f")


def runge():
    f = lambda x: 1.0 / (1.0 + 25.0 * x ** 2)
    xx = np.linspace(-1, 1, 1000)
    n = 14
    xe = np.linspace(-1, 1, n + 1)
    ye = np.polyval(np.polyfit(xe, f(xe), n), xx)
    xc = -np.cos(np.pi * np.arange(n + 1) / n)
    yc = np.polyval(np.polyfit(xc, f(xc), n), xx)
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.6), sharey=True)
    for a in ax:
        a.plot(xx, f(xx), color=GREY, lw=2.2, label=r"true $f(x)=1/(1+25x^2)$")
    ax[0].plot(xx, ye, color=RED, lw=1.6, label=f"degree-{n} fit")
    ax[0].plot(xe, f(xe), "o", color=RED, ms=5, mfc="white")
    ax[0].set_title("Evenly spaced nodes  —  the fit blows up")
    ax[0].set_ylim(-0.6, 1.4)
    ax[1].plot(xx, yc, color=BLUE, lw=1.6, label=f"degree-{n} fit")
    ax[1].plot(xc, f(xc), "o", color=BLUE, ms=5, mfc="white")
    ax[1].set_title("Chebyshev nodes  —  the fit is stable")
    for a in ax:
        a.set_xlabel("$x$"); a.legend(loc="upper center", fontsize=8.5, framealpha=0.9)
    ax[0].set_ylabel("$f,\\; p_n$")
    fig.tight_layout(); fig.savefig(OUT + "runge.pdf"); plt.close(fig)


def nodes():
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
    ax.set_title("Chebyshev nodes = evenly spaced angles dropped onto the line")
    ax.set_xlabel("$x_k=-\\cos(\\pi k/n)$"); ax.set_yticks([])
    ax.set_ylim(-0.12, 1.15); ax.set_aspect("equal"); ax.grid(False)
    fig.tight_layout(); fig.savefig(OUT + "nodes.pdf"); plt.close(fig)


def convergence():
    xx = np.linspace(-1, 1, 2000)
    analytic = lambda x: np.exp(np.sin(3 * x))
    kinked = lambda x: np.abs(x - 0.2)
    degs = np.arange(2, 41, 2)
    ea, ek = [], []
    for n in degs:
        xc = -np.cos(np.pi * np.arange(n + 1) / n)
        for f, store in ((analytic, ea), (kinked, ek)):
            c = np.linalg.solve(sg.chebyshev_basis_1d(xc, n), f(xc))
            store.append(np.max(np.abs(sg.chebyshev_basis_1d(xx, n) @ c - f(xx))))
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.semilogy(degs, ea, "o-", color=BLUE, label="smooth $e^{\\sin 3x}$ (fast)")
    ax.semilogy(degs, ek, "s-", color=RED, label="kinked $|x-0.2|$ (slow)")
    ax.set_xlabel("polynomial degree $n$")
    ax.set_ylabel(r"largest error $\max_x|f-p_n|$")
    ax.set_title("How fast the error falls depends on smoothness")
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(OUT + "convergence.pdf"); plt.close(fig)


def growth():
    dims = np.arange(1, 11)
    tensor = 5.0 ** dims
    smol2 = [sg.SmolyakGrid([-1]*d, [1]*d, mu=2).n for d in dims]
    smol3 = [sg.SmolyakGrid([-1]*d, [1]*d, mu=3).n for d in dims]
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.semilogy(dims, tensor, "o-", color=GREY, label=r"full grid ($5^d$)")
    ax.semilogy(dims, smol3, "s-", color=RED, label=r"Smolyak $\mu=3$")
    ax.semilogy(dims, smol2, "^-", color=BLUE, label=r"Smolyak $\mu=2$")
    ax.axvline(6, color="k", ls=":", lw=1)
    ax.text(5.0, 4e4, "our model:\n$d=6$ states", fontsize=9)
    ax.set_xlabel("number of state variables $d$")
    ax.set_ylabel("grid points to solve at")
    ax.set_title("Sparse grids dodge the explosion")
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(OUT + "growth.pdf"); plt.close(fig)


def smolyak2d():
    g2 = sg.SmolyakGrid([-1, -1], [1, 1], mu=2)
    g3 = sg.SmolyakGrid([-1, -1], [1, 1], mu=3)
    p = sg._level_points(4)
    TX, TY = np.meshgrid(p, p)
    fig, ax = plt.subplots(1, 3, figsize=(10.2, 3.5))
    ax[0].plot(TX.ravel(), TY.ravel(), "o", color=GREY, ms=4)
    ax[0].set_title(f"full grid\n$9\\times 9={TX.size}$ points")
    ax[1].plot(g2.points_unit[:, 0], g2.points_unit[:, 1], "o", color=BLUE, ms=5)
    ax[1].set_title(f"Smolyak $\\mu=2$\n{g2.n} points")
    ax[2].plot(g3.points_unit[:, 0], g3.points_unit[:, 1], "o", color=RED, ms=5)
    ax[2].set_title(f"Smolyak $\\mu=3$\n{g3.n} points")
    for a in ax:
        a.set_xlim(-1.15, 1.15); a.set_ylim(-1.15, 1.15); a.set_aspect("equal")
        a.set_xlabel("$x_1$"); a.grid(alpha=0.2)
    ax[0].set_ylabel("$x_2$")
    fig.tight_layout(); fig.savefig(OUT + "smolyak2d.pdf"); plt.close(fig)


def smolyak_blocks():
    # d=2, mu=2 -> 13 points, coloured by which multi-index block they come from.
    idx = sg._multi_indices(2, 2, np.array([2, 2]))
    cols = [BLUE, GREEN, ORAN, PURP, RED, "#c0398c"]
    mk = ["o", "s", "^", "D", "P", "X"]
    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    for b, iv in enumerate(idx):
        axes = [sg._new_points(i) for i in iv]
        pts = np.array(list(itertools.product(*axes)))
        lab = f"$i=({iv[0]},{iv[1]})$: {len(pts)} pt" + ("s" if len(pts) > 1 else "")
        ax.plot(pts[:, 0], pts[:, 1], mk[b], color=cols[b], ms=12, mfc=cols[b],
                mec="white", label=lab, alpha=0.9)
    ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.25, 1.4); ax.set_aspect("equal")
    ax.set_xlabel("$x_1$"); ax.set_ylabel("$x_2$")
    ax.set_title("Building the $\\mu=2$ grid in 2-D: 6 blocks, 13 points")
    ax.legend(fontsize=8.5, ncol=2, loc="upper center", framealpha=0.95)
    fig.tight_layout(); fig.savefig(OUT + "smolyak_blocks.pdf"); plt.close(fig)


def branch_global():
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.9),
                           gridspec_kw={"width_ratios": [1.55, 1]})
    a = ax[0]; a.axis("off"); a.set_xlim(0, 10); a.set_ylim(0, 10)

    def box(x, y, w, h, text, fc):
        a.add_patch(plt.Rectangle((x, y), w, h, fc=fc, ec="k", lw=1.2, alpha=0.92))
        a.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8.6)

    def arr(x0, y0, x1, y1, txt="", col="k"):
        a.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                    mutation_scale=12, lw=1.2, color=col))
        if txt:
            a.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.28, txt, fontsize=8.5,
                   ha="center", color=col)

    box(0.1, 4.1, 2.5, 1.8, "state today\n$S=(K_D,K_F,$\n$P_D,P_F,B_D,s)$", "#dce6f2")
    box(3.5, 4.3, 2.2, 1.5, "7 Gauss-\nHermite draws\nof $\\epsilon_s$", "#eae5f2")
    arr(2.6, 5.0, 3.5, 5.0, "AR(1) for $s$")
    box(7.1, 6.9, 2.7, 1.6, "$d'{=}0$: no default\nrules$(0,S')$,  SDF $\\beta$", "#d7ede0")
    box(7.1, 1.7, 2.7, 1.6, "$d'{=}1$: default\nrules$(1,S')$, haircut,\nincome SDF", "#f2dcdc")
    arr(5.7, 5.5, 7.1, 7.6, "$1-p^d(s)$", GREEN)
    arr(5.7, 4.6, 7.1, 2.5, "$p^d(s)$", RED)
    a.text(5.0, 0.5, r"$\mathbb{E}[\cdot]=\sum_{j}w_j[(1-p^d)\,(\cdot)^{d'=0}"
                     r"+p^d\,(\cdot)^{d'=1}]$", fontsize=9.5, ha="center")
    a.set_title("Two-branch quadrature — point\\_map.py $\\;(\\_cont,\\ \\_E)$",
                fontsize=10.5)

    b = ax[1]
    sp = sg.s_process_params({})
    s = np.linspace(-9.0, -3.5, 400)
    b.plot(s, 1.0 / (1.0 + np.exp(-s)), color=RED, lw=2)
    for sv, lab in [(sp["s_star"], "$s^\\star$ (0.1%)"),
                    (sp["s_star"] + 2 * sp["sigma_s"], "$+2\\sigma$ (2%)")]:
        b.axvline(sv, color=GREY, ls=":", lw=1)
        b.plot(sv, 1/(1+np.exp(-sv)), "o", color=BLUE, ms=6)
        b.annotate(lab, (sv, 1/(1+np.exp(-sv))), textcoords="offset points",
                   xytext=(6, -2), fontsize=8)
    b.set_xlabel("risk factor $s$")
    b.set_ylabel("$p^d(s)=1/(1+e^{-s})$")
    b.set_title("Priced default probability")
    old = matplotlib.rcParams["savefig.bbox"]
    matplotlib.rcParams["savefig.bbox"] = "standard"
    fig.subplots_adjust(left=0.02, right=0.96, bottom=0.15, top=0.89, wspace=0.28)
    fig.savefig(OUT + "branch_global.pdf")
    matplotlib.rcParams["savefig.bbox"] = old
    plt.close(fig)


def loop():
    # THE TIME-ITERATION LOOP, HIGHLIGHTING WHERE CHEBYSHEV/SMOLYAK DO THE WORK.
    fig, a = plt.subplots(figsize=(8.6, 5.4))
    a.axis("off"); a.set_xlim(0, 10); a.set_ylim(0, 8)

    def box(cx, cy, text, fc, w=3.0, h=1.5):
        a.add_patch(plt.Rectangle((cx - w/2, cy - h/2), w, h, fc=fc, ec="k",
                                  lw=1.3, alpha=0.95, zorder=2))
        a.text(cx, cy, text, ha="center", va="center", fontsize=9, zorder=3)

    def arr(p0, p1, txt="", col="k", dx=0.0, dy=0.35):
        a.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=15,
                    lw=1.5, color=col, shrinkA=42, shrinkB=42, zorder=1))
        a.text((p0[0]+p1[0])/2 + dx, (p0[1]+p1[1])/2 + dy, txt, fontsize=8.2,
               ha="center", color=col)

    A = (5.0, 6.9); B = (8.4, 4.0); C = (5.0, 1.1); D = (1.6, 4.0)
    box(*A, "decision rules:\n30 coefficient vectors\n(15 quantities $\\times$ 2 regimes)\non the Smolyak grid", "#dce6f2")
    box(*B, "solve the 7 equilibrium\nequations at each of the\n85 grid points; the\nexpectation evaluates the\nrules at off-grid $S'$", "#eae5f2", w=3.2, h=2.0)
    box(*C, "refit: one linear solve\n$\\Phi\\,c=f$ turns point\nvalues back into coefficients", "#d7ede0", w=3.4)
    box(*D, "damp (blend new\n& old) and check\nthe residual", "#f2ead7")

    arr(A, B, "freeze as next\nperiod's behaviour", dy=0.55)
    arr(B, C, "new values at\nall 85 points", dx=0.9, dy=0.15)
    arr(C, D, "converged?\nno $\\rightarrow$ loop", dy=-0.6)
    arr(D, A, "updated rules", dx=-0.6)

    # callouts: the two places the Chebyshev/Smolyak machinery is load-bearing
    a.annotate("Chebyshev eval + clip\n(interpolate off-grid)", (B[0], B[1]-1.05),
               (B[0]+0.1, B[1]-2.5), fontsize=8, color=RED, ha="center",
               arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
    a.annotate("Smolyak fit\n(the $\\Phi c=f$ step)", (C[0]-1.6, C[1]),
               (1.6, 0.55), fontsize=8, color=GREEN, ha="center",
               arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.2))
    a.set_title("The solve loop, and where Chebyshev / Smolyak are used",
                fontsize=11)
    old = matplotlib.rcParams["savefig.bbox"]
    matplotlib.rcParams["savefig.bbox"] = "standard"
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.03, top=0.92)
    fig.savefig(OUT + "loop.pdf")
    matplotlib.rcParams["savefig.bbox"] = old
    plt.close(fig)


if __name__ == "__main__":
    runge(); nodes(); convergence(); growth(); smolyak2d()
    smolyak_blocks(); branch_global(); loop()
    print("global figures written to", OUT)
