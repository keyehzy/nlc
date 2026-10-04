#!/usr/bin/env python3
"""Figures 1-5 of Rigol, Bryant and Singh, PRL 97, 187202 (2006), from the
CSV files written by `nlc run` (see tools/run_all.sh).

    python3 tools/plot_figures.py [--data data] [--out figures]

Thin black lines are the lower order, thick coloured lines the higher order.
As in the paper, a bare-sum curve is drawn only down to the temperature where
the two orders stop agreeing (--tol, relative to the curve's maximum).
"""
import argparse
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def fname(data, prefix, jp=0, hx=0, hz=0):
    g = lambda x: f"{x:g}"
    return os.path.join(data, f"{prefix}_jp{g(jp)}_hx{g(hx)}_hz{g(hz)}.csv")


def converged_mask(lo, hi, tol):
    """True from high T down to the first point where |hi - lo| > tol * max|hi|."""
    bad = np.abs(hi - lo) > tol * np.max(np.abs(hi))
    mask = np.ones_like(hi, dtype=bool)
    # Temperatures increase along the arrays; cut everything below the
    # highest temperature where the orders disagree.
    idx = np.nonzero(bad)[0]
    if idx.size:
        mask[: idx[-1] + 1] = False
    return mask


STYLES = [
    dict(color="tab:blue", ls="-"),
    dict(color="tab:green", ls="--"),
    dict(color="c", ls=":"),
    dict(color="m", ls="-."),
    dict(color="r", ls=(0, (6, 2, 1, 2, 1, 2))),
]


def two_orders(ax, d, prop, lo, hi, style, label, tol, extend):
    """Lower order thin black, higher order thick coloured; both cut where
    they part (plus `extend` points of the divergence for visibility)."""
    T = d["T"]
    a, b = d[f"{prop}_{lo}"], d[f"{prop}_{hi}"]
    m = converged_mask(a, b, tol)
    first = np.argmax(m)
    m[max(first - extend, 0) :] = True
    ax.plot(T[m], a[m], color="k", lw=0.7)
    ax.plot(T[m], b[m], lw=2.0, label=label, **style)


def fig_two_panel(data, out, prefix, key, values, lo, hi, tmin, tmax, ylims, labels, tol, name):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.2, 5.2), sharex=True, gridspec_kw=dict(hspace=0.05))
    for v, st, lab in zip(values, STYLES, labels):
        d = load(fname(data, prefix, **{key: v}))
        two_orders(ax1, d, "S", lo, hi, st, lab, tol, 6)
        two_orders(ax2, d, "C", lo, hi, st, lab, tol, 6)
    for ax, yl, lab in ((ax1, ylims[0], "$S$"), (ax2, ylims[1], "$C_v$")):
        ax.set_xscale("log")
        ax.set_xlim(tmin, tmax)
        ax.set_ylim(*yl)
        ax.set_ylabel(lab)
    ax1.text(0.93, 0.35, "(a)", transform=ax1.transAxes, ha="right")
    ax2.text(0.93, 0.85, "(b)", transform=ax2.transAxes, ha="right")
    ax2.set_xlabel("$T$")
    h, l = ax1.get_legend_handles_labels()
    ax1.legend(h[:2], l[:2], loc="upper left", frameon=False, fontsize=8)
    ax2.legend(h[2:], l[2:], loc="upper left", frameon=False, fontsize=8)
    fig.savefig(os.path.join(out, name), bbox_inches="tight", dpi=150)
    plt.close(fig)


def fig3(data, out):
    d = load(fname(data, "fig2_kagome_xxz", jp=1))
    T = d["T"]
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    # The paper shows the two orders down to T = 0.28, just past where they part.
    m = T >= 0.28
    ax.plot(T[m], d["C_7"][m], color="tab:green", ls="-.", lw=2, label="NLC 7")
    ax.plot(T[m], d["C_8"][m], color="r", lw=2, label="NLC 8")
    ax.set_xscale("log")
    ax.set_xlim(0.1, 10)
    ax.set_ylim(0, 0.25)
    ax.set_xlabel("$T$")
    ax.set_ylabel("$C_v$")
    ax.legend(frameon=False)
    fig.savefig(os.path.join(out, "fig3_kagome_heisenberg_cv.png"), bbox_inches="tight", dpi=150)
    plt.close(fig)


def fig4(data, out):
    d = load(fname(data, "fig4_kagome_tlfim", hx=0.5, hz=0.25))
    T = d["T"]
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    m = T >= 0.1
    for n, ls in zip((3, 4, 5, 6), (":", "--", "-.", "-")):
        ax.plot(T[m], d[f"C_{n}"][m], color="k", lw=0.8, ls=ls, label=f"bare ({n})")
    ax.plot(T, d["C_wynn2"], color="r", lw=2, label="Wynn$_2$")
    ax.plot(T, d["C_euler"], color="tab:green", ls="--", lw=2, label="Euler")
    ax.set_xscale("log")
    ax.set_xlim(0.01, 100)
    ax.set_ylim(0, 0.4)
    ax.set_xlabel("$T$")
    ax.set_ylabel("$C_v$")
    ax.legend(frameon=False, loc="center right")
    fig.savefig(os.path.join(out, "fig4_kagome_tlfim_cv.png"), bbox_inches="tight", dpi=150)
    plt.close(fig)


def fig5(data, out):
    tr = load(fname(data, "fig5_triangular", jp=1))
    sq = load(fname(data, "fig5_square", jp=1))
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    ax.plot(tr["T"], tr["S_wynn3"], color="purple", ls=":", lw=2, label="TR Wynn$_3$")
    ax.plot(tr["T"], tr["S_wynn4"], color="b", ls="--", lw=2, label="TR Wynn$_4$")
    ax.plot(sq["T"], sq["S_brez2"], color="lime", lw=2, label="SQ Brez$_2$")
    ax.plot(sq["T"], sq["S_wynn2"], color="r", ls="-.", lw=2, label="SQ Wynn$_2$")
    ax.set_xscale("log")
    ax.set_xlim(0.1, 10)
    ax.set_ylim(0, 0.8)
    ax.set_xlabel("$T$")
    ax.set_ylabel("$S$")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.savefig(os.path.join(out, "fig5_triangular_square_entropy.png"), bbox_inches="tight", dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data", default="data")
    p.add_argument("--out", default="figures")
    p.add_argument("--tol", type=float, default=0.02)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    jobs = {
        "fig1": lambda: fig_two_panel(
            a.data, a.out, "fig1_kagome_tfim", "hx", [0.0, 0.01, 0.25, 0.5, 1.0], 5, 6, 1e-3, 10,
            [(0, 0.8), (0, 0.4)], [f"$h_x$={v}" for v in ("0.0", "0.01", "0.25", "0.5", "1.0")], a.tol,
            "fig1_kagome_tfim.png"),
        "fig2": lambda: fig_two_panel(
            a.data, a.out, "fig2_kagome_xxz", "jp", [0.0, 0.01, 0.1, 0.25, 1.0], 7, 8, 1e-3, 10,
            [(0, 0.8), (0, 0.2)], [rf"$J_\perp$={v}" for v in ("0.0", "0.01", "0.1", "0.25", "1.0")], a.tol,
            "fig2_kagome_xxz.png"),
        "fig3": lambda: fig3(a.data, a.out),
        "fig4": lambda: fig4(a.data, a.out),
        "fig5": lambda: fig5(a.data, a.out),
    }
    for name, job in jobs.items():
        try:
            job()
            print("wrote", name)
        except FileNotFoundError as e:
            print("skipped", name, "(missing", os.path.basename(e.filename) + ")")


if __name__ == "__main__":
    main()
