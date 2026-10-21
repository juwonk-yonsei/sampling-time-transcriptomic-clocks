#!/usr/bin/env python3
"""Manuscript figures for TOD instructions 11 and 12. Stored tables only. No new tests."""

import math
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, FancyBboxPatch, Rectangle
from matplotlib.text import Annotation, Text

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
OUT = ROOT / "figures" / "manuscript"
SUB = ROOT / "figures" / "submission"
SRC = OUT / "source_data"
MM = 1 / 25.4
WIDTH_MM = 174

for face in (
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-BoldItalic.ttf",
):
    font_manager.fontManager.addfont(face)

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Liberation Sans",
    "mathtext.it": "Liberation Sans:italic",
    "mathtext.bf": "Liberation Sans:bold",
    "mathtext.cal": "Liberation Sans:italic",
    "mathtext.sf": "Liberation Sans",
    "mathtext.tt": "Liberation Sans",
    "lines.linewidth": 0.8,
})

PRIMARY = ["TAms_chrono", "TAms_mort", "RA_default"]
SECONDARY = ["RA_Peters", "RA_GTEx"]
CLOCK_NAME = {
    "TAms_chrono": "tAge-chronological",
    "TAms_mort": "tAge-mortality",
    "RA_default": "RNAAgeCalc-DESeq2",
    "RA_Peters": "RNAAgeCalc-Peters",
    "RA_GTEx": "RNAAgeCalc-GTEx",
}
DATA_COLOR = {
    "GSE39445": "#56B4E9",
    "GSE56931": "#0072B2",
    "GSE113883": "#009E73",
    "GSE48113": "#D55E00",
    "GSE295758": "#E69F00",
    "GSE122725": "#CC79A7",
}
DATA_MARKER = {
    "GSE39445": "o",
    "GSE56931": "o",
    "GSE113883": "o",
    "GSE48113": "s",
    "GSE295758": "s",
    "GSE122725": "D",
}
GEO = {"L1": "GSE48113", "L2": "GSE122725", "L3": "GSE295758"}
MODULE_SHORT = {
    "TM_turquoise_mortality": "TM-turquoise",
    "CMP_TAms_chrono_turquoise": "chrono-turquoise",
    "CMP_TAms_chrono_yellow": "chrono-yellow",
    "CMP_TAms_chrono_orange": "chrono-orange",
    "CMP_TAms_chrono_brown": "chrono-brown",
    "CMP_TAms_chrono_red": "chrono-red",
    "CMP_TAms_mort_turquoise": "mort-turquoise",
    "CMP_TAms_mort_darkred": "mort-darkred",
}
MODULES = list(MODULE_SHORT)
POPULATIONS = [
    "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells",
    "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils",
]
SCENARIO_LS = {"lab": "-", "free-living": "--", "free-living upper": ":"}
PHI = r"$\bar{\mathrm{\varphi}}$"
SCENARIO_LW = {"lab": 1.0, "free-living": 1.0, "free-living upper": 0.6}
SCENARIO_NAME = {
    "lab": "Laboratory",
    "free-living": "Free-living",
    "free-living upper": "Free-living, upper bound",
}
DEV_ORDER = ["GSE39445", "GSE56931", "GSE113883", "L1"]
AUDIT = []


def load(path):
    return pd.read_csv(ROOT / path, sep="\t", dtype=str, keep_default_na=False)


def write_src(frame, name):
    SRC.mkdir(parents=True, exist_ok=True)
    frame.to_csv(SRC / name, sep="\t", index=False)


def num(text):
    return float(text)


def geo_name(dataset):
    return GEO.get(dataset, dataset)


def one(frame, **eq):
    matched = frame
    for key, value in eq.items():
        matched = matched.loc[matched[key] == value]
    if len(matched) != 1:
        raise SystemExit(f"{len(matched)} rows for {eq}")
    return matched.iloc[0]


def fmt_p(text):
    p = float(text)
    if p < 0.001:
        exp = math.floor(math.log10(p))
        mant = float(f"{p / 10 ** exp:.1f}")
        if mant >= 10:
            mant = 1.0
            exp += 1
        return rf"p = ${mant:.1f} \times 10^{{{exp}}}$"
    return f"p = {p:.2g}"


def fmt_peak(text):
    return f"{float(text):.1f} h"


def yes_no(text):
    if text == "TRUE":
        return "Yes"
    if text == "FALSE":
        return "No"
    raise SystemExit(f"yes/no value {text!r}")


def sign_char(text):
    if text == "+":
        return "+"
    if text == "-":
        return "\u2212"
    raise SystemExit(f"sign {text!r}")


def use_legend(leg):
    frame = leg.get_frame()
    frame.set_alpha(None)
    frame.set_facecolor("none")
    frame.set_edgecolor("none")
    frame.set_linewidth(0)
    return leg


def panel(ax, letter):
    ax.text(0.0, 1.01, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left", clip_on=False)


def new_fig(height_mm):
    if height_mm > 234:
        raise SystemExit(f"height {height_mm} mm exceeds 234")
    return plt.figure(figsize=(WIDTH_MM * MM, height_mm * MM))


def style_open(open_mark, color):
    return {"markerfacecolor": ("white" if open_mark else color), "markeredgecolor": color, "markeredgewidth": 0.6}


def draw_interval(ax, y, x, lo, hi, xmin, xmax, color, marker, open_mark=False):
    show_lo = min(max(lo, xmin), xmax)
    show_hi = max(min(hi, xmax), xmin)
    if show_hi >= show_lo:
        ax.plot([show_lo, show_hi], [y, y], color=color, lw=0.8, zorder=2, solid_capstyle="butt")
    span = xmax - xmin
    clipped = lo < xmin - 1e-9 or hi > xmax + 1e-9
    if clipped:
        ax.figure._tod_clips = getattr(ax.figure, "_tod_clips", 0) + 1
    if lo < xmin - 1e-9:
        ax.annotate("", xy=(xmin, y), xytext=(xmin + 0.04 * span, y), arrowprops={"arrowstyle": "->", "color": color, "lw": 0.6})
    if hi > xmax + 1e-9:
        ax.annotate("", xy=(xmax, y), xytext=(xmax - 0.04 * span, y), arrowprops={"arrowstyle": "->", "color": color, "lw": 0.6})
    xx = min(max(x, xmin), xmax)
    ax.plot([xx], [y], linestyle="None", marker=marker, ms=4.5, zorder=3, **style_open(open_mark, color))


def axis_from_18(hour):
    hour = float(hour) % 24
    return hour if hour >= 18 else hour + 24


def shape_legend():
    return [
        Line2D([0], [0], marker="o", color="black", linestyle="None", ms=5, label="Development"),
        Line2D([0], [0], marker="s", color="0.35", linestyle="None", ms=5, label="Confirmation"),
        Line2D([0], [0], marker="D", color="0.35", linestyle="None", ms=4.5, label="Post hoc (GSE122725)"),
        Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor="black", linestyle="None", ms=5, label="Open symbols: LRT p \u2265 0.05"),
    ]


def group_labels(ax, groups):
    for clock, ypos in groups:
        ax.text(1.02, ypos, CLOCK_NAME[clock], transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=8, clip_on=False)


def forest_ratio(ax, clocks):
    effects = load("results/stage8/effect_summary.tsv")
    effects = effects.loc[effects["kind"] == "cosinor"].copy()
    order = ["GSE39445", "GSE56931", "GSE113883", "L1", "L2"]
    rows = []
    y = 0
    yticks, ylabels, groups = [], [], []
    for clock in clocks:
        y0 = y
        for dataset in order:
            row = one(effects, dataset=dataset, clock=clock)
            geo = geo_name(dataset)
            open_mark = num(row.p) >= 0.05
            ax.plot([num(row.low), num(row.high)], [y, y], color=DATA_COLOR[geo], lw=0.8, zorder=2)
            ax.plot([num(row.estimate)], [y], linestyle="None", marker=DATA_MARKER[geo], ms=5, zorder=3, **style_open(open_mark, DATA_COLOR[geo]))
            yticks.append(y)
            ylabels.append(geo)
            rows.append({"dataset": geo, "internal_dataset": dataset, "clock": clock, "estimate": row.estimate, "low": row.low, "high": row.high, "lrt_p": row.p, "open_marker": open_mark, "role": row.role, "y": y})
            y += 1
        groups.append((clock, (y0 + y - 1) / 2))
        y += 0.7
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_ylim(y - 0.35, -0.8)
    ax.set_xlim(left=0)
    ax.set_xlabel("2A/\u03c3")
    group_labels(ax, groups)
    return pd.DataFrame(rows)


def paired_panel(ax, clocks):
    effects = load("results/stage8/effect_summary.tsv")
    paired = effects.loc[(effects["dataset"] == "L3") & (effects["kind"] == "paired")]
    f2 = load("prereg/frozen/F2.tsv")
    rows = []
    labels = []
    for i, clock in enumerate(clocks):
        obs = one(paired, clock=clock)
        pred = one(f2, clock=clock)
        dev = [num(pred[c]) for c in ("d_GSE39445", "d_GSE56931", "d_GSE113883")]
        ax.plot([min(dev), max(dev)], [i - 0.18, i - 0.18], color="0.45", lw=1.2, solid_capstyle="butt", zorder=1)
        ax.plot([num(pred.d_exp_median)], [i - 0.18], linestyle="None", marker="D", ms=5, markerfacecolor="white", markeredgecolor="black", markeredgewidth=0.6, zorder=3)
        ax.plot([num(obs.low), num(obs.high)], [i, i], color=DATA_COLOR["GSE295758"], lw=0.8, zorder=2)
        ax.plot([num(obs.estimate)], [i], linestyle="None", marker="s", ms=5, color=DATA_COLOR["GSE295758"], zorder=4)
        name = CLOCK_NAME[clock].replace("RNAAgeCalc-", "RNAAgeCalc-\n")
        labels.append(name)
        rows.append({"clock": clock, "d": obs.estimate, "low": obs.low, "high": obs.high, "d_exp_median": pred.d_exp_median, "d_GSE39445": pred.d_GSE39445, "d_GSE56931": pred.d_GSE56931, "d_GSE113883": pred.d_GSE113883})
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(range(len(clocks)))
    ax.set_yticklabels(labels)
    ax.set_xlim(-0.2, 0.8)
    ax.set_ylim(len(clocks) - 0.45, -0.55)
    ax.set_xlabel("Effect size d (morning \u2212 evening)")
    use_legend(ax.legend(handles=[
        Line2D([0], [0], marker="s", color=DATA_COLOR["GSE295758"], linestyle="None", ms=5, label="Observed, GSE295758 (95% CI)"),
        Line2D([0], [0], marker="D", markerfacecolor="white", markeredgecolor="black", linestyle="None", ms=5, label="Predicted (median of development)"),
        Line2D([0], [0], color="0.45", lw=1.2, label="Range of development predictions"),
    ], frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8))
    return pd.DataFrame(rows)


def peak_panel(ax):
    effects = load("results/stage8/effect_summary.tsv")
    point = effects.loc[effects["kind"] == "cosinor"]
    f1 = load("prereg/frozen/F1.tsv")
    order = ["GSE39445", "GSE56931", "GSE113883", "L1"]
    rows, yticks, ylabels, groups = [], [], [], []
    y = 0
    for clock in PRIMARY:
        y0 = y
        phi = one(f1.loc[f1["dataset"] == "GSE39445"], clock=clock).phi_bar
        ys = []
        for dataset in order:
            row = one(point, dataset=dataset, clock=clock)
            geo = geo_name(dataset)
            center = num(row.peak_center) % 24
            half = num(row.peak_halfwidth_95)
            open_mark = num(row.p) >= 0.05
            draw_interval(ax, y, num(row.peak_h) % 24, center - half, center + half, 0, 24, DATA_COLOR[geo], DATA_MARKER[geo], open_mark)
            yticks.append(y)
            ylabels.append(geo)
            ys.append(y)
            rows.append({"clock": clock, "dataset": geo, "peak_h": row.peak_h, "peak_center": row.peak_center, "peak_halfwidth_95": row.peak_halfwidth_95, "lrt_p": row.p, "phi_bar": phi, "open_marker": open_mark, "y": y})
            y += 1
        ax.fill_betweenx([ys[0] - 0.45, ys[-1] + 0.45], num(phi) - 3, num(phi) + 3, color="#E6E6E6", zorder=0)
        ax.plot([num(phi), num(phi)], [ys[0] - 0.45, ys[-1] + 0.45], color="black", lw=0.6, ls=":", zorder=1)
        groups.append((clock, (y0 + y - 1) / 2))
        y += 0.8
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_ylim(y - 0.4, -0.7)
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.set_xlabel("Clock time (h)")
    group_labels(ax, groups)
    return pd.DataFrame(rows)


def fig1():
    src = load("figures/manuscript/source_data/Fig1_A.tsv")
    counts = {}
    for _, row in src.loc[src["item"] == "sample_size"].iterrows():
        counts[row.dataset] = (row.n_subjects, row.n_samples)
    sha = one(src, item="manifest_sha256").value
    when = one(src, item="freeze_commit_time").value[:10]
    nmod = one(src, item="f3_rows").n_subjects
    gates = {}
    for geo in ("GSE48113", "GSE122725", "GSE295758"):
        genes = src.loc[(src["item"] == "gate") & (src["dataset"] == geo), "gene"].tolist()
        flag = src.loc[(src["item"] == "gate") & (src["dataset"] == geo), "dataset_pass"].iloc[0]
        gates[geo] = (genes, flag)
    fig = new_fig(120)
    ax = fig.add_axes([0.01, 0.02, 0.98, 0.96])
    ax.set_xlim(0, 174)
    ax.set_ylim(0, 116)
    ax.axis("off")

    def box(x, y, w, h, text, fc, ec="0.2", ls="-"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.2", facecolor=fc, edgecolor=ec, linestyle=ls, linewidth=0.8))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8)

    def arrow(x1, y1, x2, y2, ls="-"):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops={"arrowstyle": "->", "lw": 0.7, "color": "0.2", "linestyle": ls})

    def lock(x, y):
        ax.add_patch(Rectangle((x, y - 3.2), 3.4, 2.3, facecolor="white", edgecolor="0.15", linewidth=0.6, zorder=4))
        ax.add_patch(Arc((x + 1.7, y - 1.1), 2.2, 2.2, theta1=0, theta2=180, edgecolor="0.15", linewidth=0.6, zorder=4))

    c394, c569, c113 = counts["GSE39445"], counts["GSE56931"], counts["GSE113883"]
    ax.text(2, 112, "Development", fontsize=8, fontweight="bold", color="#009E73", va="top")
    box(2, 86, 52, 22, f"(participants / samples)\nGSE39445 {c394[0]}/{c394[1]}\nGSE56931 {c569[0]}/{c569[1]}\nGSE113883 {c113[0]}/{c113[1]}", "#D9F0E3")
    arrow(54, 97, 58, 97)
    box(58, 88, 30, 18, "5 clocks\n(3 primary,\n2 secondary)", "#D9F0E3")
    arrow(88, 97, 92, 97)
    box(92, 88, 28, 18, "Mixed cosinor\nmodel", "#D9F0E3")
    arrow(120, 97, 124, 97)
    box(124, 84, 48, 26, "Frozen predictions:\npeak " + PHI + f", amplitude,\ncomposition mediation,\nmorning\u2212evening d,\n{nmod} module outcomes", "#D9F0E3")
    arrow(148, 84, 148, 83)
    arrow(148, 83, 36, 83)
    arrow(36, 83, 36, 72)

    ax.text(50, 77.6, "Preregistration and locked confirmation", fontsize=8, fontweight="bold", color="#D55E00", va="center")
    box(2, 46, 40, 26, f"Freeze\n{when}\nmanifest SHA-256\n{sha[:8]}\u2026", "#F6D6C9")
    arrow(42, 59, 48, 59)
    c481, c122, c295 = counts["GSE48113"], counts["GSE122725"], counts["GSE295758"]
    box(48, 44, 64, 30, f"Locked until freeze\nGSE48113 {c481[0]}/{c481[1]}\nGSE122725 {c122[0]}/{c122[1]}\nGSE295758 {c295[0]} pairs", "#F6D6C9", ec="#D55E00", ls="--")
    lock(49.2, 72.6)
    arrow(112, 59, 118, 59)
    gtxt = []
    for geo, word in (("GSE48113", "pass"), ("GSE122725", "fail"), ("GSE295758", "pass")):
        genes, flag = gates[geo]
        expect = "TRUE" if word == "pass" else "FALSE"
        if flag != expect:
            raise SystemExit(f"{geo} dataset_pass {flag}")
        gtxt.append(f"{geo} {word} ({', '.join(genes)})")
    box(118, 42, 54, 34, "Gate\n" + "\n".join(gtxt) + "\nHolm primary family\nReplication: \u22652 datasets", "#F6D6C9")

    arrow(87, 44, 87, 36, ls="--")
    ax.text(90, 40, "after confirmation", fontsize=8, va="center", ha="left")
    ax.text(2, 34, "Post hoc", fontsize=8, fontweight="bold", color="#CC79A7", va="center")
    box(2, 4, 40, 24, "Time-axis\ncorrection\n(GSE122725)", "#F7D6EA")
    box(46, 4, 40, 24, "Missing-value\nsensitivity", "#F7D6EA")
    box(90, 4, 40, 24, "Bias\nsimulation", "#F7D6EA")
    box(134, 4, 38, 24, "Mitigation", "#F7D6EA")
    return fig


def fig2():
    letters = ["A", "B", "C", "D", "E", "F"]
    fig = new_fig(142)
    axes = fig.subplots(2, 3)
    fig.subplots_adjust(left=0.16, right=0.97, top=0.892, bottom=0.088, hspace=0.38, wspace=0.52)
    summary = []
    for i, letter in enumerate(letters):
        ax = axes.flat[i]
        src = load(f"figures/manuscript/source_data/Fig2_{letter}.tsv")
        dataset = src["dataset"].iloc[0]
        clock = src["clock"].iloc[0]
        points = src.loc[src["record"] == "point"].copy()
        segs = src.loc[src["record"] == "segment"]
        for subject, person in points.groupby("subject_id"):
            person = person.copy()
            person["_order"] = person["order_h"].map(num)
            person = person.sort_values(["_order", "gsm"])
            for _, seg in segs.loc[segs["subject_id"] == subject].iterrows():
                sl = person.iloc[int(float(seg.i_start)):int(float(seg.i_end)) + 1]
                ax.plot(sl["clock_time_h"].map(num), sl["centered"].map(num), color="#BFBFBF", lw=0.5, marker="o", ms=1.8)
        bins = src.loc[src["record"] == "bin"]
        yerr = [0 if row["sem"] == "" else num(row["sem"]) for _, row in bins.iterrows()]
        ax.errorbar(bins["bin_center"].map(num), bins["mean"].map(num), yerr=yerr, color=DATA_COLOR[dataset], lw=1.0, marker=DATA_MARKER[dataset], ms=3.5, capsize=1.5, zorder=3)
        curve = src.loc[src["record"] == "curve"].copy()
        curve["_t"] = curve["curve_t"].map(num)
        curve = curve.sort_values("_t")
        ax.plot(curve["_t"], curve["curve_y"].map(num), color="black", lw=1.0, zorder=4)
        ann = one(src, record="annotation")
        ax.set_title(f"{CLOCK_NAME[clock]}\n{fmt_p(ann.p)}; peak {fmt_peak(ann.peak_h)}", fontsize=8, pad=3)
        ax.text(0.03, 0.97, letter.lower(), transform=ax.transAxes, fontsize=10, fontweight="bold", va="top", ha="left", zorder=6)
        ax.set_xlim(0, 24)
        ax.set_xticks([0, 6, 12, 18, 24])
        if i == 0:
            ax.set_yticks([-0.075, -0.050, -0.025, 0.000, 0.025, 0.050, 0.075, 0.100])
        ax.set_ylabel("Centered estimate (years)" if clock == "RA_default" else "Centered score")
        if i >= 3:
            ax.set_xlabel("Clock time (h)")
        def whole(value):
            text = str(value)
            if text.endswith(".0"):
                return text[:-2]
            return text
        summary.append({"panel": letter, "dataset": dataset, "clock": clock, "n_subjects": whole(ann.n_subjects), "n_samples": whole(ann.n_samples), "n_bins": whole(ann.n_bins), "p": ann.p, "peak_h": ann.peak_h, "beta_source": ann.beta_source})
    for row, text in ((0, "GSE113883 (development)"), (1, "GSE48113 Day 1 (confirmation)")):
        pos = axes[row, 0].get_position()
        fig.text(pos.x0, pos.y1 + 0.058, text, fontsize=9, fontweight="bold", ha="left", va="bottom")
    write_src(pd.DataFrame(summary), "Fig2_summary.tsv")
    return fig


def fig3():
    fig = new_fig(228)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.35, 0.72, 1.25], hspace=0.38, left=0.16, right=0.75, top=0.97, bottom=0.10)
    axa = fig.add_subplot(gs[0])
    panel(axa, "a")
    write_src(forest_ratio(axa, PRIMARY), "Fig3_A.tsv")
    axb = fig.add_subplot(gs[1])
    panel(axb, "b")
    pos = axb.get_position()
    axb.set_position([pos.x0, pos.y0, 0.56 - pos.x0, pos.height])
    write_src(paired_panel(axb, PRIMARY), "Fig3_B.tsv")
    axc = fig.add_subplot(gs[2])
    panel(axc, "c")
    write_src(peak_panel(axc), "Fig3_C.tsv")
    use_legend(fig.legend(handles=shape_legend(), loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.47, 0.005)))
    return fig


def figs1():
    fig = new_fig(188)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.45, 0.7, 0.42], hspace=0.42, left=0.18, right=0.75, top=0.93, bottom=0.12)
    axa = fig.add_subplot(gs[0])
    panel(axa, "a")
    write_src(forest_ratio(axa, SECONDARY), "FigS1_A.tsv")
    axb = fig.add_subplot(gs[1])
    panel(axb, "b")
    pos = axb.get_position()
    axb.set_position([pos.x0, pos.y0, 0.56 - pos.x0, pos.height])
    write_src(paired_panel(axb, SECONDARY), "FigS1_B.tsv")
    axc = fig.add_subplot(gs[2])
    panel(axc, "c")
    dev = load("results/stage3/h2/projection_ci.tsv")
    dev = dev.loc[dev["method"] == "lin"]
    row = one(dev, dataset="GSE56931", clock="RA_Peters")
    if not num(row.lrt_p_unadj) < 0.05:
        raise SystemExit("S1c row is not an original rhythm")
    open_mark = num(row.lrt_p_unadj) >= 0.05
    draw_interval(axc, 0, num(row.p_med), num(row.p_med_p025), num(row.p_med_p975), -0.5, 1.5, DATA_COLOR["GSE56931"], "o", open_mark)
    axc.axvline(0, color="black", lw=0.5)
    axc.axvline(1, color="0.4", lw=0.5, ls="--")
    axc.set_xlim(-0.5, 1.5)
    axc.set_yticks([0])
    axc.set_yticklabels(["GSE56931"])
    axc.set_ylim(0.8, -0.8)
    axc.set_xlabel("Proportion mediated by composition")
    axc.text(1.02, 0, CLOCK_NAME["RA_Peters"], transform=axc.get_yaxis_transform(), va="center", ha="left", fontsize=8, clip_on=False)
    write_src(pd.DataFrame([{
        "dataset": "GSE56931", "clock": "RA_Peters", "p_med": row.p_med, "low": row.p_med_p025, "high": row.p_med_p975,
        "lrt_p": row.lrt_p_unadj, "source": "results/stage3/h2/projection_ci.tsv", "y": 0,
    }]), "FigS1_C.tsv")
    use_legend(fig.legend(handles=shape_legend(), loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.46, 0.005)))
    return fig


def mediation_rows():
    pieces = []
    for path, method in (
        ("results/stage5/projection_ci.tsv", None),
        ("results/stage3/h2/projection_ci.tsv", "lin"),
        ("results/stage4/rnaseq/projection_ci.tsv", "lin"),
    ):
        tab = load(path)
        if method:
            tab = tab.loc[tab["method"] == method]
        tab = tab.copy()
        tab["source"] = path
        pieces.append(tab)
    dev = pd.concat(pieces, ignore_index=True)
    dev = dev.loc[dev["lrt_p_unadj"].map(num) < 0.05]
    rows = []
    for clock in PRIMARY:
        for dataset in ("GSE39445", "GSE56931", "GSE113883"):
            hit = dev.loc[(dev["dataset"] == dataset) & (dev["clock"] == clock)]
            if hit.empty:
                continue
            if len(hit) != 1:
                raise SystemExit(f"projection rows {dataset} {clock}")
            row = hit.iloc[0]
            rows.append({"dataset": dataset, "clock": clock, "p_med": row.p_med, "low": row.p_med_p025, "high": row.p_med_p975, "lrt_p": row.lrt_p_unadj, "source": row.source})
    h2 = load("results/confirm/h2.tsv")
    for _, row in h2.iterrows():
        if row.clock not in PRIMARY:
            raise SystemExit(f"h2 clock {row.clock}")
        rows.append({"dataset": "GSE48113", "clock": row.clock, "p_med": row.p_med, "low": row.p_med_p025, "high": row.p_med_p975, "lrt_p": row.lrt_p, "source": "results/confirm/h2.tsv"})
    return rows


def fig4():
    fig = new_fig(214)
    ax = fig.add_axes([0.16, 0.58, 0.58, 0.36])
    panel(ax, "a")
    rows = mediation_rows()
    y = 0
    yticks, ylabels, groups = [], [], []
    out = []
    for clock in PRIMARY:
        block = [r for r in rows if r["clock"] == clock]
        if not block:
            continue
        y0 = y
        for rec in block:
            open_mark = num(rec["lrt_p"]) >= 0.05
            draw_interval(ax, y, num(rec["p_med"]), num(rec["low"]), num(rec["high"]), -0.5, 1.5, DATA_COLOR[rec["dataset"]], DATA_MARKER[rec["dataset"]], open_mark)
            yticks.append(y)
            ylabels.append(rec["dataset"])
            out.append({**rec, "y": y, "open_marker": open_mark})
            y += 1
        groups.append((clock, (y0 + y - 1) / 2))
        y += 0.7
    ax.axvline(0, color="black", lw=0.5)
    ax.axvline(1, color="0.4", lw=0.5, ls="--")
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_ylim(y - 0.35, -0.7)
    ax.set_xlim(-0.5, 1.5)
    ax.set_xlabel("Proportion mediated by composition")
    group_labels(ax, groups)
    write_src(pd.DataFrame(out), "Fig4_A.tsv")

    axb = fig.add_axes([0.30, 0.08, 0.36, 0.40])
    panel(axb, "b")
    peaks = load("results/stage10/fig4b_composition_peaks.tsv")
    datasets = ["GSE39445", "GSE56931", "GSE113883", "GSE48113"]
    brows = []
    for i, population in enumerate(POPULATIONS):
        for j, dataset in enumerate(datasets):
            row = one(peaks, dataset=dataset, population=population)
            yy = i + (j - 1.5) * 0.16
            center = axis_from_18(row.peak_center)
            half = num(row.peak_halfwidth_95)
            open_mark = num(row.lrt_p) >= 0.05
            draw_interval(axb, yy, axis_from_18(row.peak_h), center - half, center + half, 18, 42, DATA_COLOR[dataset], DATA_MARKER[dataset], open_mark)
            brows.append({"dataset": dataset, "population": population, "y": yy, "open_marker": open_mark, "axis_peak": axis_from_18(row.peak_h), "axis_low": center - half, "axis_high": center + half, "axis_rule": "hour mod 24; if < 18 then +24; interval is axis(peak_center) ± peak_halfwidth_95", **{c: row[c] for c in peaks.columns}})
    axb.set_yticks(range(len(POPULATIONS)))
    axb.set_yticklabels(POPULATIONS)
    axb.set_xlim(18, 42)
    axb.set_xticks([18, 24, 30, 36, 42])
    axb.set_xticklabels(["18", "0", "6", "12", "18"])
    axb.set_xlabel("Clock time (h)")
    axb.set_ylim(len(POPULATIONS) - 0.35, -0.7)
    handles = [Line2D([0], [0], marker=DATA_MARKER[d], color=DATA_COLOR[d], linestyle="None", ms=5, label=("GSE48113 Day 1 (post hoc)" if d == "GSE48113" else d)) for d in datasets]
    handles.append(Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor="0.2", linestyle="None", ms=5, label="Open symbols: LRT p \u2265 0.05"))
    use_legend(axb.legend(handles=handles, frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1)))
    write_src(pd.DataFrame(brows), "Fig4_B.tsv")
    return fig


def marker_area(p):
    return 14 * (-math.log10(p))


def ms_from_area(area):
    return 2 * math.sqrt(area / math.pi)


def fig5():
    tab = load("figures/manuscript/source_data/Fig5_A.tsv")
    fig = new_fig(214)
    ax = fig.add_axes([0.08, 0.40, 0.46, 0.54])
    panel(ax, "a")
    xs = tab["phi_bar"].map(num)
    ys = tab["y"].map(num)
    frozen_y = [num(v) if v != "" else None for v in tab["frozen_y"]]
    lo = min(list(xs) + list(ys) + [v for v in frozen_y if v is not None]) - 3
    hi = max(list(xs) + list(ys) + [v for v in frozen_y if v is not None]) + 4
    grid = np.linspace(lo, hi, 80)
    ax.fill_between(grid, grid - 3, grid + 3, color="#E6E6E6", zorder=0)
    ax.plot(grid, grid, color="black", lw=0.6, zorder=1)
    label_rows = []
    # Hand positions in data coordinates. Leaders are drawn to the plotted point.
    name_xy = {
        "TM-turquoise": (6.2, 14.6),
        "chrono-turquoise": (17.2, 8.4),
        "chrono-yellow": (16.6, 11.6),
        "chrono-orange": (6.2, 1.0),
        "chrono-brown": (8.6, -1.6),
        "chrono-red": (6.2, 5.4),
        "mort-turquoise": (17.2, 14.2),
        "mort-darkred": (18.4, 16.8),
    }
    n_xy = {
        "TM-turquoise": (6.2, 21.8),
        "chrono-turquoise": (17.4, 20.6),
        "chrono-yellow": (16.8, -0.6),
        "chrono-orange": (6.2, 8.2),
        "mort-turquoise": (17.4, 24.6),
        "mort-darkred": (18.4, 4.4),
    }
    for _, row in tab.iterrows():
        phi = num(row.phi_bar)
        y = num(row.y)
        filled = row.filled == "True"
        ax.plot([phi], [y], linestyle="None", marker="o", ms=ms_from_area(num(row.marker_area)), zorder=3, **style_open(not filled, "black"))
        lx, ly = name_xy[row.short]
        ax.annotate(row.short, xy=(phi, y), xytext=(lx, ly), fontsize=8, ha="left", va="center", arrowprops={"arrowstyle": "-", "color": "0.45", "lw": 0.5})
        label_rows.append({"clock": row.clock, "kind": "name", "x": phi, "y": y, "label_x": lx, "label_y": ly, "text": row.short})
        if row.frozen_y != "":
            fy = num(row.frozen_y)
            ax.plot([phi, phi], [y, fy], color="0.65", lw=0.5, ls=":", zorder=1)
            ax.plot([phi], [fy], linestyle="None", marker="x", color="0.35", ms=ms_from_area(marker_area(num(row.frozen_p))), zorder=2)
            nx, ny = n_xy[row.short]
            ax.annotate(f"n = {row.frozen_n_entered}", xy=(phi, fy), xytext=(nx, ny), fontsize=8, color="0.25", ha="left", va="center", arrowprops={"arrowstyle": "-", "color": "0.55", "lw": 0.5})
            label_rows.append({"clock": row.clock, "kind": "n", "x": phi, "y": fy, "label_x": nx, "label_y": ny, "text": f"n = {row.frozen_n_entered}"})
    ax.set_xlim(lo, hi + 2)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Predicted peak " + PHI + " (h)")
    ax.set_ylabel("Observed peak, GSE48113 (h)")
    ax.set_aspect("equal", adjustable="box", anchor="W")
    size_handles = [Line2D([0], [0], marker="o", color="black", linestyle="None", ms=ms_from_area(marker_area(p)), label=lab) for p, lab in ((0.05, r"p = 0.05"), (1e-4, r"p = $10^{-4}$"), (1e-8, r"p = $10^{-8}$"))]
    use_legend(fig.legend(handles=[
        Line2D([0], [0], marker="o", color="black", linestyle="None", ms=5, label="Participant-mean imputation (post hoc)"),
        Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor="black", linestyle="None", ms=5, label="Open symbols: LRT p \u2265 0.05"),
        Line2D([0], [0], marker="x", color="0.35", linestyle="None", ms=5, label="Preregistered scoring"),
        *size_handles,
    ], loc="upper left", bbox_to_anchor=(0.58, 0.92), frameon=False))
    write_src(tab, "Fig5_A.tsv")
    write_src(pd.DataFrame(label_rows), "Fig5_A_labels.tsv")

    b = load("figures/manuscript/source_data/Fig5_B.tsv")
    axb = fig.add_axes([0.04, 0.03, 0.94, 0.34])
    panel(axb, "b")
    axb.set_xlim(0, 164)
    axb.set_ylim(-0.8, len(MODULES) + 5.4)
    axb.axis("off")
    group_headers = [
        (54, "Development\n\u2212log10 q"),
        (88, "Preregistered\nscoring"),
        (112, "Participant-mean\nimputation\n(post hoc)"),
        (144, "GSE295758"),
    ]
    for x, text in group_headers:
        axb.text(x, len(MODULES) + 4.5, text, ha="center", va="center", fontsize=8)
    for x, text in ((40, "GSE\n39445"), (54, "GSE\n56931"), (68, "GSE\n113883")):
        axb.text(x, len(MODULES) + 1.6, text, ha="center", va="center", fontsize=8)
    for x, text in ((82, "Features"), (94, "Confirmed"), (106, "Features"), (118, "Confirmed"), (132, "Predicted"), (144, "Observed"), (156, "Agree")):
        axb.text(x, len(MODULES) + 1.9, text, ha="center", va="top", rotation=90, fontsize=8, rotation_mode="anchor")
    qvals = []
    for _, row in b.iterrows():
        for col in ("q_GSE39445", "q_GSE56931", "q_GSE113883"):
            qvals.append(-math.log10(num(row[col])))
    vmin, vmax = min(qvals), max(qvals)
    cmap = plt.get_cmap("Blues")
    brows = []
    for r, row in b.iterrows():
        y = len(MODULES) - r
        axb.text(1, y, row.short, ha="left", va="center", fontsize=8)
        for i, col in enumerate(("q_GSE39445", "q_GSE56931", "q_GSE113883")):
            lq = -math.log10(num(row[col]))
            color = cmap((lq - vmin) / (vmax - vmin if vmax > vmin else 1))
            cx = (40, 54, 68)[i]
            axb.add_patch(Rectangle((cx - 4.2, y - 0.38), 8.4, 0.76, facecolor=color, edgecolor="white", linewidth=0.4))
            rgb = color[:3]
            tone = "white" if sum(rgb) / 3 < 0.55 else "black"
            axb.text(cx, y, f"{lq:.1f}", ha="center", va="center", fontsize=8, color=tone)
        frozen_feat = "\u2014" if row.frozen_n == "0" else row.frozen_n
        d11_feat = "\u2014" if row.d11_n == "0" else row.d11_n
        agree = yes_no(row.sign_agrees) if row.sign_agrees in ("TRUE", "FALSE") else "\u2014"
        for x, text in ((82, frozen_feat), (94, yes_no(row.frozen_confirmed)), (106, d11_feat), (118, yes_no(row.d11_confirmed)), (132, sign_char(row.predicted_sign)), (144, sign_char(row.observed_sign)), (156, agree)):
            axb.text(x, y, text, ha="center", va="center", fontsize=8)
        brows.append({**row.to_dict(), "neg_log10_q_GSE39445": f"{-math.log10(num(row.q_GSE39445)):.1f}", "neg_log10_q_GSE56931": f"{-math.log10(num(row.q_GSE56931)):.1f}", "neg_log10_q_GSE113883": f"{-math.log10(num(row.q_GSE113883)):.1f}", "agree_drawn": agree, "white_text_if_mean_rgb_below": 0.55})
    write_src(pd.DataFrame(brows), "Fig5_B.tsv")
    return fig


def fig6():
    a = load("figures/manuscript/source_data/Fig6_A.tsv")
    fig = new_fig(220)
    ax = fig.add_axes([0.10, 0.62, 0.52, 0.32])
    panel(ax, "a")
    xs = a["x_gse39445"].map(num)
    ys = list(a["y_geo"].map(num)) + list(a["y_corrected"].map(num))
    lo = min(list(xs) + ys) - 3
    hi = max(list(xs) + ys) + 4
    grid = np.linspace(lo, hi, 60)
    ax.fill_between(grid, grid - 2, grid + 2, color="#E6E6E6", zorder=0)
    ax.plot(grid, grid, color="black", lw=0.6, zorder=1)
    # Right-hand names follow the target y order. PER2 sits between PER1 and PER3.
    # Its leader goes up, then right, so the horizontal leg clears the PER3 diamond.
    target_y = {row.gene: num(row.y_corrected) for _, row in a.iterrows()}
    gap = 2.3
    per2_label_y = (target_y["PER3"] + target_y["PER1"]) / 2
    label_y = {"PER2": per2_label_y, "PER1": per2_label_y + gap, "CRY2": target_y["CRY2"]}
    y = per2_label_y
    for gene in ("PER3", "NR1D2", "NR1D1", "CRY1"):
        y = min(target_y[gene], y - gap)
        label_y[gene] = y
    name_xy = {gene: (28.0, label_y[gene]) for gene in label_y}
    name_ha = {gene: "left" for gene in label_y}
    name_xy["TEF"] = (10.5, 16.6)
    name_xy["ARNTL"] = (10.5, 13.4)
    name_xy["DBP"] = (8.5, -4.2)
    name_ha["TEF"] = name_ha["ARNTL"] = name_ha["DBP"] = "left"
    for _, row in a.iterrows():
        x = num(row.x_gse39445)
        y0, y1 = num(row.y_geo), num(row.y_corrected)
        ax.plot([x, x], [y0, y1], color="0.65", lw=0.5, zorder=1)
        ax.plot([x], [y0], linestyle="None", marker="o", ms=5, markerfacecolor="white", markeredgecolor="0.45", markeredgewidth=0.6, zorder=2)
        filled = row.filled == "True"
        ax.plot([x], [y1], linestyle="None", marker="D", ms=4.5, zorder=3, **style_open(not filled, "#CC79A7"))
        lx, ly = name_xy[row.gene]
        arrowprops = {"arrowstyle": "-", "color": "0.45", "lw": 0.5}
        if row.gene == "PER2":
            arrowprops["connectionstyle"] = "angle,angleA=180,angleB=90,rad=0"
        ax.annotate(row.gene, xy=(x, y1), xytext=(lx, ly), fontsize=8, ha=name_ha[row.gene], va="center", arrowprops=arrowprops)
    ax.set_xlim(min(lo, 0), max(hi, 34))
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Peak, GSE39445 (h since wake)")
    ax.set_ylabel("Peak, GSE122725 (h since wake)")
    use_legend(ax.legend(handles=[
        Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor="0.45", linestyle="None", ms=5, label="GEO metadata definition"),
        Line2D([0], [0], marker="D", color="#CC79A7", linestyle="None", ms=4.5, label="Publication definition (post hoc)"),
        Line2D([0], [0], marker="D", markerfacecolor="white", markeredgecolor="#CC79A7", linestyle="None", ms=4.5, label="Open symbols: LRT p \u2265 0.05"),
    ], frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1)))
    write_src(a, "Fig6_A.tsv")
    write_src(pd.DataFrame([{"gene": g, "label_x": xy[0], "label_y": xy[1]} for g, xy in name_xy.items()]), "Fig6_A_labels.tsv")

    b = load("figures/manuscript/source_data/Fig6_B.tsv")
    axb = fig.add_axes([0.16, 0.28, 0.28, 0.24])
    panel(axb, "b")
    axb.set_title("GSE48113", fontsize=8, pad=8)
    style = {
        "frozen": ("o", "black", "black"),
        "D2_knn": ("s", "#808080", "#808080"),
        "D11_participant_mean": ("D", "white", "black"),
    }
    for _, row in b.iterrows():
        marker, face, edge = style[row.version]
        if row.n_entered == "0":
            axb.annotate("0", xy=(0.4, num(row.y)), xytext=(2, 0), textcoords="offset points", fontsize=8, va="center", ha="left")
        else:
            axb.plot([num(row.n_entered)], [num(row.y)], linestyle="None", marker=marker, ms=5, markerfacecolor=face, markeredgecolor=edge, markeredgewidth=0.6)
    weights = {"TAms_chrono": "1839", "TAms_mort": "10487"}
    for clock, weight in weights.items():
        yy = num(one(b, clock=clock, version="frozen").y)
        axb.plot([num(weight), num(weight)], [yy - 0.32, yy + 0.32], color="0.35", lw=1.0)
    labels, ys = [], []
    for clock in ["TAms_chrono", "TAms_mort"] + MODULES:
        labels.append(CLOCK_NAME[clock] if clock in CLOCK_NAME else MODULE_SHORT[clock])
        ys.append(num(one(b, clock=clock, version="frozen").y))
    axb.set_yticks(ys)
    axb.set_yticklabels(labels)
    axb.set_xscale("log")
    axb.set_xlim(0.4, 30000)
    axb.set_xlabel("Features used (log scale)")
    use_legend(axb.legend(handles=[
        Line2D([0], [0], marker="o", color="black", linestyle="None", ms=5, label="Preregistered scoring"),
        Line2D([0], [0], marker="s", color="#808080", linestyle="None", ms=5, label="kNN imputation"),
        Line2D([0], [0], marker="D", markerfacecolor="white", markeredgecolor="black", linestyle="None", ms=5, label="Participant-mean imputation"),
        Line2D([0], [0], color="0.35", lw=1.0, linestyle="None", marker="|", ms=10, markeredgewidth=1.0, label="Nonzero-weight features of the clock"),
        Line2D([0], [0], color="none", linestyle="None", label="0 = no feature available"),
    ], frameon=False, loc="lower left", bbox_to_anchor=(0.02, 0.012), bbox_transform=fig.transFigure, ncol=1))
    write_src(b, "Fig6_B.tsv")

    c = load("figures/manuscript/source_data/Fig6_C.tsv")
    axc = fig.add_axes([0.56, 0.30, 0.40, 0.20])
    panel(axc, "c")
    axc.set_title("GSE39445, 20% of values masked", fontsize=8, pad=8)
    unmasked = c.loc[(c["treatment"] == "unmasked") & (c["mask"] == "random")]
    x = 2.0
    xpos = [0]
    ticklabels = ["Unmasked"]
    groups = []
    drawn = []
    for _, row in unmasked.iterrows():
        drawn.append({**row.to_dict(), "x_slot": 0})
    for mask, label in (("random", "Random masking"), ("low_expression", "Low-expression masking")):
        start = x
        for treat, short in (("kNN", "kNN"), ("participant mean", "Participant\nmean")):
            block = c.loc[(c["mask"] == mask) & (c["treatment"] == treat)]
            xpos.append(x)
            ticklabels.append(short)
            for _, row in block.iterrows():
                drawn.append({**row.to_dict(), "x_slot": x})
            x += 1.8
        groups.append((start, x - 1.8, label))
        x += 0.8
    markers = {"RA_default": "o", "RA_Peters": "^", "RA_GTEx": "v"}
    shift = {"RA_default": -0.12, "RA_Peters": 0.0, "RA_GTEx": 0.12}
    for rec in drawn:
        rec["x_drawn"] = rec["x_slot"] + shift[rec["clock"]]
        axc.plot([rec["x_drawn"]], [num(rec["neg_log10_p"])], linestyle="None", marker=markers[rec["clock"]], ms=5, color="#404040", zorder=3)
    y05 = -math.log10(0.05)
    axc.axhline(y05, color="0.2", lw=0.6, ls="--", zorder=1)
    axc.text(max(xpos), y05, "p = 0.05", fontsize=8, va="bottom", ha="right")
    axc.set_xticks(xpos)
    axc.set_xticklabels(ticklabels, fontsize=8)
    axc.set_xlim(-0.9, max(xpos) + 0.9)
    axc.set_ylabel("\u2212log10 LRT p")
    axc.set_ylim(bottom=0)
    trans = axc.get_xaxis_transform()
    for x0, x1, label in groups:
        yb = -0.42
        axc.plot([x0, x0, x1, x1], [yb + 0.08, yb, yb, yb + 0.08], transform=trans, color="0.35", lw=0.5, clip_on=False, solid_capstyle="butt")
        axc.text((x0 + x1) / 2, yb - 0.02, label, transform=trans, ha="center", va="top", fontsize=8, clip_on=False)
    use_legend(axc.legend(handles=[Line2D([0], [0], marker=markers[ck], color="#404040", linestyle="None", ms=5, label=CLOCK_NAME[ck]) for ck in ("RA_default", "RA_Peters", "RA_GTEx")], frameon=False, loc="lower left", bbox_to_anchor=(0.62, 0.012), bbox_transform=fig.transFigure, ncol=1))
    write_src(pd.DataFrame(drawn), "Fig6_C.tsv")
    return fig


def reference_lines(ax):
    ax.axhline(0.05, color="#BFBFBF", lw=0.6, ls="-", zorder=0)
    ax.axhline(0.10, color="#808080", lw=0.6, ls="-.", zorder=0)


def fig7(prefix):
    a = load(f"figures/manuscript/source_data/{prefix}_A.tsv")
    b = load(f"figures/manuscript/source_data/{prefix}_B.tsv")
    c = load(f"figures/manuscript/source_data/{prefix}_C.tsv")
    d = load(f"figures/manuscript/source_data/{prefix}_D.tsv")
    e = load(f"figures/manuscript/source_data/{prefix}_E.tsv")
    clocks = list(dict.fromkeys(a["clock"]))
    fig = new_fig(220)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1, 1], hspace=0.52, wspace=0.34, left=0.10, right=0.90, top=0.93, bottom=0.20)
    panels = {
        "a": fig.add_subplot(gs[0, :]),
        "b": fig.add_subplot(gs[1, 0]),
        "c": fig.add_subplot(gs[1, 1]),
        "d": fig.add_subplot(gs[2, 0]),
        "e": fig.add_subplot(gs[2, 1]),
    }
    for clock in clocks:
        for scenario, ls in SCENARIO_LS.items():
            block = a.loc[(a["clock"] == clock) & (a["scenario"] == scenario)].copy()
            if "record" in block.columns:
                block = block.loc[block["record"] != "annotation"]
            block = block.loc[block["b_over_total_sd"] != ""].copy()
            block["_s"] = block["s"].map(num)
            block = block.sort_values("_s")
            panels["a"].plot(block["_s"], block["b_over_total_sd"].map(num), color=DATA_COLOR.get(clock, None) if False else {
                "TAms_chrono": "#0072B2", "TAms_mort": "#E69F00", "RA_default": "#009E73", "RA_Peters": "#CC79A7", "RA_GTEx": "#56B4E9",
            }[clock], ls=ls, lw=SCENARIO_LW[scenario], zorder=3)
            if prefix == "Fig7" and clock == "RA_default" and scenario in ("lab", "free-living"):
                hit = one(block, s="3")
                panels["a"].plot([3], [num(hit.b_over_total_sd)], linestyle="None", marker="o", ms=4.5, markerfacecolor="white", markeredgecolor={
                    "TAms_chrono": "#0072B2", "TAms_mort": "#E69F00", "RA_default": "#009E73",
                }[clock], markeredgewidth=0.6, zorder=4)
            sub = b.loc[(b["clock"] == clock) & (b["scenario"] == scenario)].copy()
            sub["_s"] = sub["s"].map(num)
            sub = sub.sort_values("_s")
            color = {"TAms_chrono": "#0072B2", "TAms_mort": "#E69F00", "RA_default": "#009E73", "RA_Peters": "#CC79A7", "RA_GTEx": "#56B4E9"}[clock]
            panels["b"].plot(sub["_s"], sub["fpr"].map(num), color=color, ls=ls, lw=SCENARIO_LW[scenario], zorder=3)
            sub = c.loc[(c["clock"] == clock) & (c["scenario"] == scenario)].copy()
            sub["_n"] = sub["n"].map(num)
            sub = sub.sort_values("_n")
            panels["c"].plot(sub["_n"], sub["fpr"].map(num), color=color, ls=ls, lw=SCENARIO_LW[scenario], zorder=3)
            sub = d.loc[(d["clock"] == clock) & (d["scenario"] == scenario)].copy()
            sub["_n"] = sub["n_per_group"].map(num)
            sub = sub.sort_values("_n")
            panels["d"].plot(sub["_n"], sub["fpr"].map(num), color=color, ls=ls, lw=SCENARIO_LW[scenario], zorder=3)
            if scenario == "lab":
                sub = e.loc[e["clock"] == clock].copy()
                sub["_w"] = sub["W"].map(num)
                sub = sub.sort_values("_w")
                ys = sub["power"].map(num).to_numpy()
                lo = sub["wilson_low"].map(num).to_numpy()
                hiw = sub["wilson_high"].map(num).to_numpy()
                panels["e"].errorbar(sub["_w"], ys, yerr=[ys - lo, hiw - ys], color=color, lw=1.0, marker="o", ms=3.5, capsize=2)
    color = None
    panels["a"].axhline(0, color="#BFBFBF", lw=0.5, zorder=1)
    panels["a"].set_xlabel("s (h)")
    panels["a"].set_ylabel("Bias (fraction of total SD)")
    panels["a"].set_xlim(0, 12)
    ys = []
    for line in panels["a"].get_lines():
        ys.extend(float(v) for v in line.get_ydata() if np.isfinite(v))
    lo, hi = min(ys + [0.0]), max(ys + [0.0])
    pad = 0.05 * (hi - lo if hi > lo else 1.0)
    panels["a"].set_ylim(lo - pad, hi + pad)
    if prefix == "Fig7":
        panels["a"].text(0.03, 0.04, "RNAAgeCalc-DESeq2, s = 3 h:\n1.1 years (laboratory), 0.6 years (free-living)", transform=panels["a"].transAxes, fontsize=8, va="bottom", ha="left")
    for key, text in (("b", "Design 1, n = 100"), ("c", "Design 1, s = 3 h"), ("d", "Design 2, s = 4 h"), ("e", "Design 3, laboratory, n = 50")):
        panels[key].text(0.02, 0.98, text, transform=panels[key].transAxes, fontsize=8, va="top", ha="left")
    for key in ("b", "c", "d"):
        reference_lines(panels[key])
        panels[key].set_ylim(bottom=0)
    panels["b"].set_xlabel("s (h)")
    panels["b"].set_ylabel("False-positive rate")
    panels["c"].set_xlabel("n")
    panels["c"].set_ylabel("False-positive rate")
    panels["d"].set_xlabel("n per group")
    panels["d"].set_ylabel("False-positive rate")
    panels["e"].set_xlabel("W (h)")
    panels["e"].set_ylabel("Power")
    panels["e"].set_ylim(0, 1)
    for letter, axp in panels.items():
        panel(axp, letter)
    clock_color = {"TAms_chrono": "#0072B2", "TAms_mort": "#E69F00", "RA_default": "#009E73", "RA_Peters": "#CC79A7", "RA_GTEx": "#56B4E9"}
    handles = [Line2D([0], [0], color=clock_color[c], lw=1.2, label=CLOCK_NAME[c]) for c in clocks]
    handles += [Line2D([0], [0], color="black", lw=SCENARIO_LW[s], ls=SCENARIO_LS[s], label=SCENARIO_NAME[s]) for s in SCENARIO_LS]
    handles.append(Line2D([0], [0], color="#BFBFBF", lw=0.6, label="Nominal α = 0.05"))
    handles.append(Line2D([0], [0], color="#808080", lw=0.6, ls="-.", label="10% false-positive rate"))
    use_legend(fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.005)))
    for frame, name in ((a, "A"), (b, "B"), (c, "C"), (d, "D"), (e, "E")):
        write_src(frame, f"{prefix}_{name}.tsv")
    return fig


def place_ratio_labels(ax, items):
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    ymin, ymax = ax.get_ylim()
    height = ax.get_window_extent(renderer).height
    data_per_pt = (ymax - ymin) / (height * 72 / ax.figure.dpi) if height else 0
    # Clear gap of one 8 pt line between label boxes.
    gap = 16 * data_per_pt
    order = sorted(range(len(items)), key=lambda i: items[i][1])
    placed = [items[i][1] for i in order]
    for i in range(1, len(placed)):
        if placed[i] < placed[i - 1] + gap:
            placed[i] = placed[i - 1] + gap
    top = ymax - 0.02 * (ymax - ymin)
    if placed and placed[-1] > top:
        shift = placed[-1] - top
        placed = [v - shift for v in placed]
    out = [None] * len(items)
    for pos, idx in zip(placed, order):
        out[idx] = pos
    for (x, y, text, color), yy in zip(items, out):
        if abs(yy - y) > 0.15 * gap:
            ax.annotate("", xy=(x, y), xytext=(1.28, yy), arrowprops={"arrowstyle": "-", "color": "0.55", "lw": 0.5})
        ax.text(1.32, yy, text, fontsize=8, va="center", ha="left", color=color, clip_on=False)


def fig8(clocks, prefix):
    cos = load("results/stage9/composition_cosinor.tsv")
    l3 = load("results/stage9/composition_l3.tsv")
    fig = new_fig(168)
    gs = fig.add_gridspec(2, len(clocks), height_ratios=[1.15, 1], hspace=0.62, wspace=0.55, left=0.12, right=0.96, top=0.94, bottom=0.18)
    arows = []
    datasets = ["GSE39445", "GSE56931", "GSE113883", "L1"]
    for i, clock in enumerate(clocks):
        ax = fig.add_subplot(gs[0, i])
        items = []
        marks = []
        offsets = (-0.06, -0.02, 0.02, 0.06)
        for j, dataset in enumerate(datasets):
            row = one(cos, dataset=dataset, clock=clock)
            geo = geo_name(dataset)
            dx = offsets[j]
            before_open = num(row.lrt_p_before) >= 0.05
            after_open = num(row.lrt_p_after) >= 0.05
            y0, y1 = num(row.two_a_before), num(row.two_a_after)
            ax.plot([dx, 1 + dx], [y0, y1], color=DATA_COLOR[geo], lw=0.8, zorder=2)
            marks.append((dx, y0, DATA_MARKER[geo], DATA_COLOR[geo], before_open))
            marks.append((1 + dx, y1, DATA_MARKER[geo], DATA_COLOR[geo], after_open))
            items.append((1 + dx, y1, f"{num(row.two_a_ratio):.2f}", DATA_COLOR[geo]))
            arows.append({"clock": clock, "dataset": geo, "x_offset": dx, "two_a_before": row.two_a_before, "two_a_after": row.two_a_after, "two_a_ratio": row.two_a_ratio, "lrt_p_before": row.lrt_p_before, "lrt_p_after": row.lrt_p_after, "open_before": before_open, "open_after": after_open, "source": "results/stage9/composition_cosinor.tsv"})
        for x, y, marker, color, open_mark in marks:
            if open_mark:
                ax.plot([x], [y], linestyle="None", marker=marker, ms=5, zorder=3, **style_open(True, color))
        for x, y, marker, color, open_mark in marks:
            if not open_mark:
                ax.plot([x], [y], linestyle="None", marker=marker, ms=5, zorder=4, **style_open(False, color))
        ax.set_xlim(-0.45, 1.90)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Unadjusted", "Composition-adjusted"], rotation=30, ha="right")
        ax.set_ylim(bottom=0)
        ax.set_title(CLOCK_NAME[clock], fontsize=8)
        ax.set_ylabel("2A (years)" if clock.startswith("RA_") else "2A (score)")
        if i == 0:
            panel(ax, "a")
        place_ratio_labels(ax, items)
    brows = []
    for i, clock in enumerate(clocks):
        ax = fig.add_subplot(gs[1, i])
        row = one(l3, clock=clock)
        series = [
            ("Unadjusted", row.mean_delta_before, row.ci_low_before, row.ci_high_before),
            ("Composition-\nadjusted", row.mean_delta_composition, row.ci_low_composition, row.ci_high_composition),
        ]
        if row.mean_delta_time != "":
            series.append(("Time-\ncorrected", row.mean_delta_time, row.ci_low_time, row.ci_high_time))
        for j, (lab, est, lo, hi) in enumerate(series):
            ax.errorbar(j, num(est), yerr=[[num(est) - num(lo)], [num(hi) - num(est)]], color=DATA_COLOR["GSE295758"], marker="s", ms=5, capsize=2, lw=0.8, markeredgewidth=0.6)
            brows.append({"clock": clock, "adjustment": lab.replace("\n", " "), "estimate": est, "low": lo, "high": hi, "time_source": row.time_source, "marker": "s", "color": DATA_COLOR["GSE295758"]})
        ax.axhline(0, color="0.45", lw=0.5)
        ax.set_xticks(range(len(series)))
        ax.set_xticklabels([lab for lab, *_ in series], rotation=30, ha="right")
        ax.set_xlim(-0.65, len(series) - 0.35)
        ax.set_title(CLOCK_NAME[clock], fontsize=8)
        ax.set_ylabel("Mean \u0394, morning \u2212 evening (years)" if clock.startswith("RA_") else "Mean \u0394, morning \u2212 evening (score)")
        if i == 0:
            panel(ax, "b")
    handles = [Line2D([0], [0], marker=DATA_MARKER[geo_name(d)], color=DATA_COLOR[geo_name(d)], lw=0.8, label=geo_name(d)) for d in datasets]
    handles.append(Line2D([0], [0], marker="s", color=DATA_COLOR["GSE295758"], lw=0.8, label="GSE295758"))
    handles.append(Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor="0.2", linestyle="None", ms=5, label="Open symbols: LRT p \u2265 0.05"))
    use_legend(fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0)))
    write_src(pd.DataFrame(arows), f"{prefix}_A.tsv")
    write_src(pd.DataFrame(brows), f"{prefix}_B.tsv")
    return fig


def text_items(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    items = []
    for text in fig.findobj(Text):
        if not text.get_visible() or not str(text.get_text()).strip():
            continue
        items.append((text, text.get_window_extent(renderer).expanded(0, 0)))
    return items, renderer


def same_axis_ticks(fig, left, right):
    for ax in fig.axes:
        ticks = list(ax.get_xticklabels()) + list(ax.get_yticklabels())
        ids = {id(t) for t in ticks}
        if id(left) in ids and id(right) in ids:
            return True
    return False


def marker_points(fig, renderer):
    legend_ids = set()
    legends = list(fig.legends)
    for ax in fig.axes:
        if ax.get_legend() is not None:
            legends.append(ax.get_legend())
    for leg in legends:
        legend_ids.update(id(obj) for obj in leg.findobj())
    points = []
    for line in fig.findobj(Line2D):
        if line.axes is None or id(line) in legend_ids or line.get_marker() in (None, "None", "none"):
            continue
        xy = line.get_xydata()
        if len(xy) == 0:
            continue
        disp = line.get_transform().transform(xy)
        radius = line.get_markersize() * fig.dpi / 72 / 2
        for pt in disp:
            points.append((pt[0], pt[1], radius, line))
    return points, legends


def boxes_overlap(a, b):
    return a.x0 < b.x1 - 0.6 and b.x0 < a.x1 - 0.6 and a.y0 < b.y1 - 0.6 and b.y0 < a.y1 - 0.6


def segment_hits_box(p0, p1, box):
    steps = int(max(abs(p1[0] - p0[0]), abs(p1[1] - p0[1]))) + 1
    for i in range(steps + 1):
        t = i / steps
        x = p0[0] + t * (p1[0] - p0[0])
        y = p0[1] + t * (p1[1] - p0[1])
        if box.x0 + 0.4 < x < box.x1 - 0.4 and box.y0 + 0.4 < y < box.y1 - 0.4:
            return True
    return False


def data_outside(fig):
    hits = []
    for ax in fig.axes:
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        if x0 > x1:
            x0, x1 = x1, x0
        if y0 > y1:
            y0, y1 = y1, y0
        ex = 1e-6 * (abs(x1 - x0) + 1)
        ey = 1e-6 * (abs(y1 - y0) + 1)
        for line in ax.get_lines():
            if line.axes is None or line.get_visible() is False:
                continue
            if line.get_transform() is not ax.transData:
                continue
            xy = line.get_xydata()
            if len(xy) == 0:
                continue
            for x, y in xy:
                if not (np.isfinite(x) and np.isfinite(y)):
                    continue
                if x < x0 - ex or x > x1 + ex or y < y0 - ey or y > y1 + ey:
                    hits.append(f"{ax.get_label() or id(ax)}:{x:.4g},{y:.4g}")
    return hits


def leader_hits(fig, items, points):
    text_hits = []
    point_hits = []
    for ann in fig.findobj(Annotation):
        patch = ann.arrow_patch
        if patch is None or not str(ann.get_text()).strip():
            continue
        verts = patch.get_path().transformed(patch.get_transform()).vertices
        if len(verts) < 2:
            continue
        own = ann
        for i in range(len(verts) - 1):
            p0, p1 = verts[i], verts[i + 1]
            for text, ext in items:
                if text is own:
                    continue
                if segment_hits_box(p0, p1, ext):
                    text_hits.append(str(text.get_text()).replace("\n", " ")[:40])
            for x, y, radius, line in points:
                if min((verts[0][0] - x) ** 2 + (verts[0][1] - y) ** 2, (verts[-1][0] - x) ** 2 + (verts[-1][1] - y) ** 2) <= (radius + 1.5) ** 2:
                    continue
                steps = 12
                for k in range(steps + 1):
                    t = k / steps
                    px = p0[0] + t * (p1[0] - p0[0])
                    py = p0[1] + t * (p1[1] - p0[1])
                    if (px - x) ** 2 + (py - y) ** 2 <= (radius + 0.6) ** 2:
                        point_hits.append(f"{x:.0f},{y:.0f}")
                        break
    return text_hits, point_hits


def proper_cross(a, b, c, d):
    def side(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    return side(a, b, c) * side(a, b, d) < 0 and side(c, d, a) * side(c, d, b) < 0


def leader_geometry(fig, points):
    arrows = []
    for ann in fig.findobj(Annotation):
        patch = ann.arrow_patch
        if patch is None or not str(ann.get_text()).strip():
            continue
        verts = patch.get_path().transformed(patch.get_transform()).vertices
        if len(verts) < 2:
            continue
        arrows.append((ann, verts))
    cross = []
    for i in range(len(arrows)):
        _, vi = arrows[i]
        for j in range(i + 1, len(arrows)):
            _, vj = arrows[j]
            hit = False
            for a in range(len(vi) - 1):
                for b in range(len(vj) - 1):
                    if proper_cross(vi[a], vi[a + 1], vj[b], vj[b + 1]):
                        hit = True
            if hit:
                cross.append(str(arrows[i][0].get_text()).strip() + " x " + str(arrows[j][0].get_text()).strip())
    foreign = []
    for ann, verts in arrows:
        seen = set()
        for i in range(len(verts) - 1):
            p0, p1 = verts[i], verts[i + 1]
            for x, y, radius, line in points:
                if id(line) in seen:
                    continue
                if min((verts[0][0] - x) ** 2 + (verts[0][1] - y) ** 2, (verts[-1][0] - x) ** 2 + (verts[-1][1] - y) ** 2) <= (radius + 1.5) ** 2:
                    continue
                steps = 16
                for k in range(steps + 1):
                    t = k / steps
                    px = p0[0] + t * (p1[0] - p0[0])
                    py = p0[1] + t * (p1[1] - p0[1])
                    if (px - x) ** 2 + (py - y) ** 2 <= radius ** 2:
                        seen.add(id(line))
                        foreign.append(str(ann.get_text()).strip())
                        break
    return cross, foreign


def other_axes_text(fig, items, renderer):
    hits = []
    for text, ext in items:
        for ax in fig.axes:
            if ax is text.axes or not ax.axison:
                continue
            if boxes_overlap(ext, ax.get_window_extent(renderer)):
                hits.append(str(text.get_text()).replace("\n", " ")[:40])
                break
    return hits


def ylabel_neighbor_mm(fig, renderer):
    axes = [ax for ax in fig.axes if ax.axison and ax.get_visible()]
    best = None
    for ax in axes:
        box = ax.get_window_extent(renderer)
        spines = []
        for other in axes:
            if other is ax:
                continue
            ob = other.get_window_extent(renderer)
            if ob.x1 <= box.x0 + 1 and ob.y0 < box.y1 and ob.y1 > box.y0:
                spines.append(ob.x1)
        if not spines:
            continue
        spine_x = max(spines)
        texts = [ax.yaxis.label, *ax.get_yticklabels()]
        for text in texts:
            if not text.get_visible() or not str(text.get_text()).strip():
                continue
            dist_mm = (text.get_window_extent(renderer).x0 - spine_x) / fig.dpi * 25.4
            if best is None or dist_mm < best:
                best = dist_mm
    return best


def text_on_line_or_marker(fig, items, points):
    """Text boxes that meet a connector, leader, diagonal, or marker. Own leader excluded."""
    own_arrow = {}
    for ann in fig.findobj(Annotation):
        patch = ann.arrow_patch
        if patch is None:
            continue
        own_arrow[id(ann)] = patch.get_path().transformed(patch.get_transform()).vertices
    lines = []
    for line in fig.findobj(Line2D):
        if line.axes is None or not line.get_visible() or line.get_transform() is not line.axes.transData:
            continue
        if line.get_linestyle() in ("None", "none", None):
            continue
        xy = line.get_transform().transform(line.get_xydata())
        if len(xy) >= 2:
            lines.append(xy)
    hits = []
    for text, ext in items:
        if ext.width < 1 and ext.height < 1:
            continue
        hit = False
        for xy in lines:
            for i in range(len(xy) - 1):
                if segment_hits_box(xy[i], xy[i + 1], ext):
                    hit = True
                    break
            if hit:
                break
        if not hit:
            for ann_id, other in own_arrow.items():
                if ann_id == id(text) or len(other) < 2:
                    continue
                for i in range(len(other) - 1):
                    if segment_hits_box(other[i], other[i + 1], ext):
                        hit = True
                        break
                if hit:
                    break
        if not hit:
            for x, y, radius, line in points:
                box = ext.expanded(radius, radius)
                if box.x0 <= x <= box.x1 and box.y0 <= y <= box.y1:
                    hit = True
                    break
        if hit:
            hits.append(str(text.get_text()).replace("\n", " ")[:40])
    return hits


def audit_fig(fig, stem):
    items, renderer = text_items(fig)
    text_overlap = []
    for i in range(len(items)):
        ti, ei = items[i]
        for j in range(i + 1, len(items)):
            tj, ej = items[j]
            if boxes_overlap(ei, ej):
                text_overlap.append((str(ti.get_text()).replace("\n", " "), str(tj.get_text()).replace("\n", " ")))
    points, legends = marker_points(fig, renderer)
    text_marker = []
    for text, ext in items:
        for x, y, radius, line in points:
            box = ext.expanded(radius, radius)
            if box.x0 <= x <= box.x1 and box.y0 <= y <= box.y1:
                text_marker.append(str(text.get_text()).replace("\n", " "))
                break
    legend_hits = []
    for leg in legends:
        box = leg.get_window_extent(renderer)
        for x, y, radius, line in points:
            if box.x0 <= x <= box.x1 and box.y0 <= y <= box.y1:
                legend_hits.append(leg.get_texts()[0].get_text() if leg.get_texts() else "legend")
                break
    outside = []
    fb = fig.bbox
    for text, ext in items:
        if ext.x0 < fb.x0 or ext.x1 > fb.x1 or ext.y0 < fb.y0 or ext.y1 > fb.y1:
            outside.append(str(text.get_text()).replace("\n", " "))
    outside_data = data_outside(fig)
    leader_text, leader_point = leader_hits(fig, items, points)
    leader_cross, leader_foreign = leader_geometry(fig, points)
    axes_text = other_axes_text(fig, items, renderer)
    neighbor_mm = ylabel_neighbor_mm(fig, renderer)
    geom_hits = text_on_line_or_marker(fig, items, points)
    line_text = []
    seen_line_text = set()
    for ax in fig.axes:
        for line in ax.get_lines():
            if line.get_transform() is not ax.transData or not line.get_visible():
                continue
            xy = line.get_transform().transform(line.get_xydata())
            if len(xy) < 2:
                continue
            for i in range(len(xy) - 1):
                for text, ext in items:
                    if id(text) in seen_line_text:
                        continue
                    if segment_hits_box(xy[i], xy[i + 1], ext):
                        seen_line_text.add(id(text))
                        line_text.append(str(text.get_text()).replace("\n", " ")[:50])
    sizes = [t.get_fontsize() for t, _ in items]
    alpha_n = 0
    for artist in fig.findobj():
        alpha = artist.get_alpha()
        if alpha is not None and abs(alpha - 1) > 1e-8:
            alpha_n += 1
    AUDIT.append({
        "stem": stem,
        "n_text": len(items),
        "text_overlap": len(text_overlap),
        "outside": len(outside),
        "leader_text": len(leader_text),
        "leader_point": len(leader_point),
        "leader_cross": len(leader_cross),
        "leader_foreign": len(leader_foreign),
        "axes_text": len(axes_text),
        "text_geom": len(geom_hits),
        "text_geom_detail": " | ".join(geom_hits[:20]),
        "ylabel_gap_mm": "" if neighbor_mm is None else f"{neighbor_mm:.2f}",
        "leader_cross_detail": " | ".join(leader_cross[:12]),
        "leader_foreign_detail": " | ".join(leader_foreign[:12]),
        "axes_text_detail": " | ".join(axes_text[:12]),
        "line_text": len(line_text),
        "line_text_detail": " | ".join(line_text[:15]),
        "data_outside": len(outside_data),
        "intentional_clip": int(getattr(fig, "_tod_clips", 0)),
        "text_marker": len(text_marker),
        "legend_data": len(legend_hits),
        "min_fontsize": min(sizes) if sizes else "",
        "alpha_not_1": alpha_n,
        "text_overlap_detail": " | ".join(f"{a} || {b}" for a, b in text_overlap[:20]),
        "outside_detail": " | ".join(outside[:20]),
        "leader_text_detail": " | ".join(leader_text[:20]),
        "leader_point_detail": " | ".join(leader_point[:20]),
        "data_outside_detail": " | ".join(outside_data[:20]),
        "text_marker_detail": " | ".join(text_marker[:12]),
    })


def pdf_edge_words(path):
    raw = subprocess.check_output(["pdftotext", "-bbox", str(path), "-"], text=True, errors="replace")
    inset = 72 / 25.4
    hits = []
    for page in re.finditer(r'<page width="([\d.]+)" height="([\d.]+)"[^>]*>(.*?)</page>', raw, re.S):
        width, height = float(page.group(1)), float(page.group(2))
        body = page.group(3)
        for word in re.finditer(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)"[^>]*>([^<]*)</word>', body):
            x0, y0, x1, y1 = (float(word.group(i)) for i in range(1, 5))
            if x0 < inset or y0 < inset or x1 > width - inset or y1 > height - inset:
                hits.append(word.group(5))
    return hits


def save_fig(fig, stem):
    audit_fig(fig, stem)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.eps")
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    fig.savefig(OUT / f"{stem}.tiff", dpi=1200, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    hits = pdf_edge_words(OUT / f"{stem}.pdf")
    AUDIT[-1]["pdf_clip"] = len(hits)
    AUDIT[-1]["pdf_clip_detail"] = " | ".join(hits[:30])


def submission_copies():
    SUB.mkdir(parents=True, exist_ok=True)
    main = {
        "Fig1_design": "Fig1", "Fig2_trajectories": "Fig2", "Fig3_effects": "Fig3", "Fig4_composition": "Fig4",
        "Fig5_modules": "Fig5", "Fig6_public": "Fig6", "Fig7_bias": "Fig7", "Fig8_mitigation": "Fig8",
    }
    for src, dest in main.items():
        shutil.copyfile(OUT / f"{src}.eps", SUB / f"{dest}.eps")
        shutil.copyfile(OUT / f"{src}.tiff", SUB / f"{dest}.tif")
    for src, dest in (("FigS1_effects_secondary", "FigS1"), ("FigS2_bias_secondary", "FigS2"), ("FigS3_mitigation_secondary", "FigS3")):
        shutil.copyfile(OUT / f"{src}.pdf", SUB / f"{dest}.pdf")


def merge_review():
    stems = [
        "Fig1_design", "Fig2_trajectories", "Fig3_effects", "Fig4_composition", "Fig5_modules",
        "Fig6_public", "Fig7_bias", "Fig8_mitigation", "FigS1_effects_secondary", "FigS2_bias_secondary", "FigS3_mitigation_secondary",
    ]
    cmd = ["gs", "-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", "-sOutputFile=" + str(OUT / "all_figures_review.pdf")]
    cmd += [str(OUT / f"{stem}.pdf") for stem in stems]
    subprocess.check_call(cmd)


def main():
    save_fig(fig1(), "Fig1_design")
    save_fig(fig2(), "Fig2_trajectories")
    save_fig(fig3(), "Fig3_effects")
    save_fig(fig4(), "Fig4_composition")
    save_fig(fig5(), "Fig5_modules")
    save_fig(fig6(), "Fig6_public")
    save_fig(fig7("Fig7"), "Fig7_bias")
    save_fig(fig8(PRIMARY, "Fig8"), "Fig8_mitigation")
    save_fig(figs1(), "FigS1_effects_secondary")
    save_fig(fig7("FigS2"), "FigS2_bias_secondary")
    save_fig(fig8(SECONDARY, "FigS3"), "FigS3_mitigation_secondary")
    merge_review()
    submission_copies()
    write_src(pd.DataFrame(AUDIT), "audit.tsv")
    cols = ["stem", "pdf_clip", "outside", "text_overlap", "leader_text", "data_outside", "leader_cross", "leader_foreign", "axes_text", "text_geom", "ylabel_gap_mm", "min_fontsize", "alpha_not_1"]
    print(pd.DataFrame(AUDIT)[cols].to_string(index=False))


if __name__ == "__main__":
    main()
