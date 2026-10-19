#!/usr/bin/env python3
"""Stage 9 figures. Same clock names and default font as stage 8."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
FIG = ROOT / "figures"
STAGE = ROOT / "results" / "stage9"
NAMES = {
    "TAms_chrono": "tAge chronological",
    "TAms_mort": "tAge mortality",
    "RA_default": "RNAAgeCalc blood",
}
CLOCKS = ["TAms_chrono", "TAms_mort", "RA_default", "RA_Peters", "RA_GTEx"]
RA = ["RA_default", "RA_Peters", "RA_GTEx"]
SCENARIO_COLOR = {
    "lab": "#1b9e77",
    "free-living": "#d95f02",
    "free-living upper": "#7570b3",
}
CLOCK_COLOR = {
    "TAms_chrono": "#1b9e77",
    "TAms_mort": "#d95f02",
    "RA_default": "#7570b3",
    "RA_Peters": "#e7298a",
    "RA_GTEx": "#66a61e",
}


def clock_label(clock: str) -> str:
    return NAMES.get(clock, clock)


def save(fig, stem: str) -> None:
    fig.savefig(FIG / f"{stem}.pdf")
    fig.savefig(FIG / f"{stem}.png", dpi=300)
    plt.close(fig)


def fig6() -> list[dict]:
    bias = pd.read_csv(STAGE / "sim_design1_bias.tsv", sep="\t")
    fpr1 = pd.read_csv(STAGE / "sim_design1_fpr.tsv", sep="\t")
    fpr2 = pd.read_csv(STAGE / "sim_design2_fpr.tsv", sep="\t")
    power = pd.read_csv(STAGE / "sim_design3_power.tsv", sep="\t")
    scenarios = ["lab", "free-living", "free-living upper"]
    styles = {"lab": "-", "free-living": "--", "free-living upper": ":"}
    fig = plt.figure(figsize=(12, 14))
    gs = fig.add_gridspec(4, 4, height_ratios=[1, 1, 1, 0.85], hspace=0.45, wspace=0.35)
    for i, scenario in enumerate(scenarios):
        ax = fig.add_subplot(gs[0, i])
        block = bias.loc[bias["scenario"] == scenario]
        for clock in CLOCKS:
            one = block.loc[block["clock"] == clock].sort_values("s")
            ax.plot(one["s"], one["b_over_total_sd"], marker="o", ms=3, color=CLOCK_COLOR[clock])
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_title(f"A  {scenario}", fontsize=9)
        ax.set_xlabel("s (h)")
        if i == 0:
            ax.set_ylabel("b(s) / total SD")
        ax = fig.add_subplot(gs[1, i])
        block = bias.loc[(bias["scenario"] == scenario) & (bias["clock"].isin(RA))]
        for clock in RA:
            one = block.loc[block["clock"] == clock].sort_values("s")
            ax.plot(one["s"], one["b"], marker="o", ms=3, color=CLOCK_COLOR[clock])
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_xlabel("s (h)")
        if i == 0:
            ax.set_ylabel("b(s), years")
    fpr_axes = [
        (fig.add_subplot(gs[2, 0]), fpr1, "n", 2, "B  design 1, s = 2"),
        (fig.add_subplot(gs[2, 1]), fpr1, "n", 4, "B  design 1, s = 4"),
        (fig.add_subplot(gs[2, 2]), fpr2, "n_per_group", 2, "B  design 2, s = 2"),
        (fig.add_subplot(gs[2, 3]), fpr2, "n_per_group", 4, "B  design 2, s = 4"),
    ]
    for ax, tab, ncol, s_value, title_text in fpr_axes:
        block = tab.loc[tab["s"] == s_value]
        for scenario in scenarios:
            for clock in CLOCKS:
                one = block.loc[(block["scenario"] == scenario) & (block["clock"] == clock)].sort_values(ncol)
                ax.plot(one[ncol], one["fpr"], color=CLOCK_COLOR[clock], linestyle=styles[scenario], marker="o", ms=3)
        ax.axhline(0.05, color="grey", lw=0.5)
        ax.set_title(title_text, fontsize=8)
        ax.set_xlabel("n" if ncol == "n" else "n per group")
        ax.set_ylabel("False positive rate")
    ax = fig.add_subplot(gs[3, 1:3])
    for clock in CLOCKS:
        one = power.loc[power["clock"] == clock].sort_values("W")
        ax.plot(one["W"], one["power"], marker="o", color=CLOCK_COLOR[clock], label=clock_label(clock))
    ax.set_xlabel("W (h)")
    ax.set_ylabel("Power")
    ax.set_title("C  design 3", fontsize=9)
    handles = [plt.Line2D([0], [0], color=CLOCK_COLOR[c], label=clock_label(c)) for c in CLOCKS]
    handles += [plt.Line2D([0], [0], color="black", linestyle=ls, label=name) for name, ls in styles.items()]
    fig.legend(handles=handles, frameon=False, fontsize=8, loc="lower center", ncol=4)
    fig.suptitle("Bias simulation (post hoc)", y=0.995)
    save(fig, "Fig6_simulation")
    return [{
        "figure": "Fig6",
        "panel": "A b(s); B false positive rate at s=2 and 4; C design 3 power",
        "n_subjects": "design 1 n 20-200; design 2 n per group 25-200; design 3 n 50",
        "n_samples": "design 3 Monte Carlo 2000",
        "role": "post hoc",
        "source": "results/stage9/sim_design1_bias.tsv; results/stage9/sim_design1_fpr.tsv; results/stage9/sim_design2_fpr.tsv; results/stage9/sim_design3_power.tsv",
    }]


def fig7() -> list[dict]:
    cos = pd.read_csv(STAGE / "composition_cosinor.tsv", sep="\t")
    l3 = pd.read_csv(STAGE / "composition_l3.tsv", sep="\t")
    datasets = ["GSE39445", "GSE56931", "GSE113883", "L1"]
    dataset_color = {"GSE39445": "#1b9e77", "GSE56931": "#d95f02", "GSE113883": "#7570b3", "L1": "#e7298a"}
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 5, hspace=0.45, wspace=0.4)
    for i, clock in enumerate(CLOCKS):
        ax = fig.add_subplot(gs[0, i])
        block = cos.loc[cos["clock"] == clock]
        for dataset in datasets:
            one = block.loc[block["dataset"] == dataset]
            if one.empty:
                continue
            ax.plot([0, 1], [float(one["two_a_before"].iloc[0]), float(one["two_a_after"].iloc[0])], marker="o", color=dataset_color[dataset], label=dataset)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["before", "after"])
        ax.set_title(clock_label(clock), fontsize=8)
        ax.set_xlim(-0.2, 1.2)
        if i == 0:
            ax.set_ylabel("A  2A")
            ax.legend(frameon=False, fontsize=6)
    kinds = [
        ("before", "mean_delta_before", "ci_low_before", "ci_high_before"),
        ("composition", "mean_delta_composition", "ci_low_composition", "ci_high_composition"),
        ("time", "mean_delta_time", "ci_low_time", "ci_high_time"),
    ]
    for i, clock in enumerate(CLOCKS):
        ax = fig.add_subplot(gs[1, i])
        one = l3.loc[l3["clock"] == clock].iloc[0]
        for k, (label, est, lo, hi) in enumerate(kinds):
            if pd.isna(one[est]):
                continue
            y = float(one[est])
            low = float(one[lo]) if pd.notna(one[lo]) else y
            high = float(one[hi]) if pd.notna(one[hi]) else y
            ax.errorbar(k, y, yerr=[[max(0.0, y - low)], [max(0.0, high - y)]], fmt="o", color="#d95f02")
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["before", "composition", "time"], rotation=30, ha="right", fontsize=7)
        if i == 0:
            ax.set_ylabel("B  mean delta")
    fig.suptitle("Composition and time adjustment (post hoc)")
    save(fig, "Fig7_mitigation")
    n_bits = []
    for dataset, sub in cos.groupby("dataset", sort=False):
        n_bits.append(f"{dataset} {int(sub['n_subjects'].iloc[0])}/{int(sub['n_samples'].iloc[0])}")
    return [{
        "figure": "Fig7",
        "panel": "A cosinor 2A before and after; B L3 mean delta",
        "n_subjects": "; ".join(n_bits) + "; L3 58",
        "n_samples": "L3 116; L1 drops non-finite composition",
        "role": "post hoc",
        "source": "results/stage9/composition_cosinor.tsv; results/stage9/composition_l3.tsv; results/confirm/h5.tsv; results/confirm/family_P.tsv; results/confirm/family_S1.tsv",
    }]


def main() -> None:
    sources = fig6() + fig7()
    path = FIG / "figure_sources.tsv"
    old = pd.read_csv(path, sep="\t", dtype=str) if path.exists() else pd.DataFrame()
    if not old.empty:
        old = old.loc[~old["figure"].isin(["Fig6", "Fig6A", "Fig6B", "Fig6C", "Fig7", "Fig7A", "Fig7B"])]
    out = pd.concat([old, pd.DataFrame(sources)], ignore_index=True)
    out = out.drop_duplicates(subset=["figure", "panel"], keep="last")
    out.to_csv(path, sep="\t", index=False)


if __name__ == "__main__":
    main()
