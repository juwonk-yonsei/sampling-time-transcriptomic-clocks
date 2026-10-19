#!/usr/bin/env python3
"""Manuscript figures for stage 8. Post hoc / exploratory."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

NAMES = {
    "TAms_chrono": "tAge chronological",
    "TAms_mort": "tAge mortality",
    "RA_default": "RNAAgeCalc blood",
}
PRIMARY = ["TAms_chrono", "TAms_mort", "RA_default"]
ROLE_COLOR = {"development": "#1b9e77", "confirmation": "#d95f02", "post hoc": "#7570b3"}
ROLE_MARKER = {"development": "o", "confirmation": "s", "post hoc": "D"}
BLOOD = [
    "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells",
    "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils",
]


def save(fig, stem: str) -> None:
    fig.savefig(FIG / f"{stem}.pdf")
    fig.savefig(FIG / f"{stem}.png", dpi=300)
    plt.close(fig)


def clock_label(clock: str) -> str:
    return NAMES.get(clock, clock)


def finite_xerr(est: float, lo: float, hi: float):
    return [[max(0.0, est - lo)], [max(0.0, hi - est)]]


def fig2() -> list[dict]:
    sources = []
    l1_samples = pd.read_csv(ROOT / "results/confirm/L1/samples_scoring.tsv", sep="\t", dtype=str)
    day1 = l1_samples.loc[l1_samples["period"] == "D1"].copy()
    day1["clock_time_h"] = day1["clock_time_h"].astype(float)
    tage = pd.read_csv(ROOT / "results/confirm/L1/tage_scores.tsv", sep="\t")
    ra = pd.read_csv(ROOT / "results/confirm/L1/ra_scores.tsv", sep="\t")
    l1_scores = pd.concat([tage, ra], ignore_index=True)
    l1_family = pd.read_csv(ROOT / "results/confirm/family_P.tsv", sep="\t", dtype=str)
    s113 = pd.read_csv(ROOT / "results/stage3/GSE113883_samples.tsv", sep="\t")
    s113 = s113.loc[s113["condition"] == "control"].copy()
    t113 = pd.read_csv(ROOT / "results/stage5/clock_scores.tsv", sep="\t")
    t113 = t113.loc[(t113["dataset"] == "GSE113883") & (t113["clock"].isin(PRIMARY)), ["gsm", "clock", "value"]]
    r113 = pd.read_csv(ROOT / "results/stage4/rnaseq/clock_scores.tsv", sep="\t")
    r113 = r113.loc[(r113["method"] == "lin") & (r113["clock"] == "RA_default"), ["gsm", "clock", "value"]]
    dev_scores = pd.concat([t113, r113], ignore_index=True)
    dev_fit = pd.read_csv(ROOT / "results/stage8/cosinor_effects_boot.tsv", sep="\t")
    panels = [
        ("GSE113883", "development", s113, dev_scores, dev_fit.loc[dev_fit["dataset"] == "GSE113883"]),
        ("L1 Day 1", "confirmation", day1, l1_scores, l1_family.loc[l1_family["dataset"] == "L1"]),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11, 7), sharex=True)
    hours = np.linspace(0, 24, 200)
    ang = 2 * math.pi * hours / 24
    for row, (dataset, role, samples, scores, fits) in enumerate(panels):
        samples = samples.copy()
        samples["gsm"] = samples["gsm"].astype(str)
        for col, clock in enumerate(PRIMARY):
            ax = axes[row, col]
            block = scores.loc[scores["clock"] == clock, ["gsm", "value"]].copy()
            block["gsm"] = block["gsm"].astype(str)
            merged = samples.merge(block, on="gsm")
            fit = fits.loc[fits["clock"] == clock].iloc[0]
            beta_cos = float(fit["beta_cos"])
            beta_sin = float(fit["beta_sin"])
            t = merged["clock_time_h"].astype(float).to_numpy()
            y = merged["value"].astype(float).to_numpy()
            wave = beta_cos * np.cos(2 * math.pi * t / 24) + beta_sin * np.sin(2 * math.pi * t / 24)
            offset = float(np.mean(y - wave))
            for _, one in merged.groupby("subject_id"):
                one = one.sort_values("clock_time_h")
                ax.plot(one["clock_time_h"], one["value"], color=ROLE_COLOR[role], alpha=0.35, lw=0.8)
            ax.plot(hours, offset + beta_cos * np.cos(ang) + beta_sin * np.sin(ang), color="black", lw=1.5)
            if row == 0:
                ax.set_title(clock_label(clock))
            if col == 0:
                ax.set_ylabel(dataset)
            ax.set_xlim(0, 24)
            sources.append({
                "figure": "Fig2",
                "panel": f"{dataset} {clock}",
                "n_subjects": merged["subject_id"].nunique(),
                "n_samples": len(merged),
                "role": role,
                "source": "results/stage5/clock_scores.tsv; results/stage4/rnaseq/clock_scores.tsv; results/stage3/GSE113883_samples.tsv"
                if dataset.startswith("GSE") else "results/confirm/L1/tage_scores.tsv; results/confirm/L1/ra_scores.tsv; results/confirm/L1/samples_scoring.tsv; results/confirm/family_P.tsv",
            })
    for ax in axes[1]:
        ax.set_xlabel("Wall-clock time (h)")
    handles = [
        plt.Line2D([0], [0], color=ROLE_COLOR["development"], label="development"),
        plt.Line2D([0], [0], color=ROLE_COLOR["confirmation"], label="confirmation"),
        plt.Line2D([0], [0], color="black", label="group cosinor"),
    ]
    fig.legend(handles=handles, frameon=False, fontsize=8, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Participant trajectories and group cosinor (post hoc, exploratory)")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "Fig2_trajectories")
    return sources


def fig3(effects: pd.DataFrame) -> list[dict]:
    sources = []
    fig, axes = plt.subplots(1, 2, figsize=(11, 7))
    cos = effects.loc[effects["kind"] == "cosinor"].copy()
    clocks = PRIMARY + ["RA_Peters", "RA_GTEx"]
    datasets = ["GSE39445", "GSE56931", "GSE113883", "L1", "L2"]
    ax = axes[0]
    y = 0
    yticks = []
    ylabels = []
    for clock in clocks:
        for dataset in datasets:
            row = cos.loc[(cos["clock"] == clock) & (cos["dataset"] == dataset)]
            if row.empty or not np.isfinite(float(row["estimate"].iloc[0])):
                continue
            role = row["role"].iloc[0]
            est = float(row["estimate"].iloc[0])
            lo = float(row["low"].iloc[0])
            hi = float(row["high"].iloc[0])
            ax.errorbar(est, y, xerr=finite_xerr(est, lo, hi), fmt=ROLE_MARKER[role], color=ROLE_COLOR[role], ms=5)
            yticks.append(y)
            ylabels.append(f"{clock_label(clock)} | {dataset}")
            y += 1
        y += 0.6
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("2A/σ")
    ax.set_title("Cosinor datasets")
    ax.axvline(0, color="grey", lw=0.5)
    ax.invert_yaxis()
    paired = effects.loc[effects["kind"] == "paired"].copy()
    ax = axes[1]
    y = 0
    yticks = []
    ylabels = []
    for clock in clocks:
        row = paired.loc[paired["clock"] == clock]
        if row.empty:
            continue
        est = float(row["estimate"].iloc[0])
        lo = float(row["low"].iloc[0])
        hi = float(row["high"].iloc[0])
        ax.errorbar(est, y, xerr=finite_xerr(est, lo, hi), fmt="s", color=ROLE_COLOR["confirmation"], ms=5)
        yticks.append(y)
        ylabels.append(clock_label(clock))
        y += 1
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_xlabel("d")
    ax.set_title("L3")
    ax.axvline(0, color="grey", lw=0.5)
    ax.invert_yaxis()
    handles = [
        plt.Line2D([0], [0], marker=ROLE_MARKER[k], color=ROLE_COLOR[k], linestyle="", label=k)
        for k in ["development", "confirmation", "post hoc"]
    ]
    axes[0].legend(handles=handles, frameon=False, fontsize=8)
    fig.suptitle("Effect summary (post hoc, exploratory)")
    fig.tight_layout()
    save(fig, "Fig3A_effect_forest")
    counts = []
    for dataset, sub in effects.groupby("dataset", sort=False):
        counts.append(f"{dataset}:{sub['n_subjects'].iloc[0]}/{sub['n_samples'].iloc[0]}")
    sources.append({
        "figure": "Fig3A", "panel": "2A/sigma and d",
        "n_subjects": "; ".join(counts), "n_samples": "subjects/samples in n_subjects",
        "role": "mixed", "source": "results/stage8/effect_summary.tsv",
    })

    fig, axes = plt.subplots(1, 3, figsize=(10, 4), subplot_kw={"projection": "polar"})
    f1 = pd.read_csv(ROOT / "prereg/frozen/F1.tsv", sep="\t")
    for ax, clock in zip(axes, PRIMARY):
        phi = float(f1.loc[f1["clock"] == clock, "phi_bar"].iloc[0])
        ax.plot([phi * math.pi / 12], [1], marker="*", color="black", ms=12, linestyle="")
        block = cos.loc[(cos["clock"] == clock) & (cos["dataset"] != "L2")]
        for _, row in block.iterrows():
            if not np.isfinite(float(row["peak_h"])):
                continue
            ang = float(row["peak_h"]) * math.pi / 12
            ax.plot([ang], [1], marker=ROLE_MARKER[row["role"]], color=ROLE_COLOR[row["role"]], ms=7, linestyle="")
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_xticks([h * math.pi / 12 for h in range(0, 24, 6)])
        ax.set_xticklabels(["0", "6", "12", "18"])
        ax.set_yticks([])
        ax.set_title(clock_label(clock))
    handles = [
        plt.Line2D([0], [0], marker="*", color="black", linestyle="", ms=10, label="development phi_bar"),
        plt.Line2D([0], [0], marker=ROLE_MARKER["development"], color=ROLE_COLOR["development"], linestyle="", label="development"),
        plt.Line2D([0], [0], marker=ROLE_MARKER["confirmation"], color=ROLE_COLOR["confirmation"], linestyle="", label="confirmation"),
    ]
    axes[2].legend(handles=handles, frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.05, 1))
    fig.suptitle("Peak hour and development phi_bar (post hoc, exploratory)")
    fig.tight_layout()
    save(fig, "Fig3B_peak_clock")
    counts_b = []
    for dataset, sub in effects.loc[effects["dataset"] != "L2"].groupby("dataset", sort=False):
        counts_b.append(f"{dataset}:{sub['n_subjects'].iloc[0]}/{sub['n_samples'].iloc[0]}")
    sources.append({
        "figure": "Fig3B", "panel": "primary peaks, L2 omitted",
        "n_subjects": "; ".join(counts_b), "n_samples": "subjects/samples in n_subjects",
        "role": "mixed", "source": "prereg/frozen/F1.tsv; results/stage8/effect_summary.tsv",
    })
    return sources


def fig4() -> list[dict]:
    rows = []
    for path in [
        ROOT / "results/stage5/projection_ci.tsv",
        ROOT / "results/stage3/h2/projection_ci.tsv",
        ROOT / "results/stage4/rnaseq/projection_ci.tsv",
    ]:
        tab = pd.read_csv(path, sep="\t", dtype=str)
        tab["source_file"] = str(path.relative_to(ROOT))
        rows.append(tab)
    dev = pd.concat(rows, ignore_index=True)
    dev = dev.loc[dev["clock"].isin(PRIMARY)]
    l1 = pd.read_csv(ROOT / "results/confirm/h2.tsv", sep="\t", dtype=str)
    fig = plt.figure(figsize=(11, 6))
    ax = fig.add_subplot(1, 2, 1)
    y = 0
    yticks, ylabels = [], []
    for clock in PRIMARY:
        for dataset in ["GSE39445", "GSE56931", "GSE113883"]:
            row = dev.loc[(dev["dataset"] == dataset) & (dev["clock"] == clock)]
            if row.empty:
                continue
            est, lo, hi = map(float, [row["p_med"].iloc[0], row["p_med_p025"].iloc[0], row["p_med_p975"].iloc[0]])
            ax.errorbar(est, y, xerr=finite_xerr(est, lo, hi), fmt="o", color=ROLE_COLOR["development"], ms=5)
            yticks.append(y)
            ylabels.append(f"{clock_label(clock)} | {dataset}")
            y += 1
        row = l1.loc[l1["clock"] == clock]
        if not row.empty:
            est, lo, hi = map(float, [row["p_med"].iloc[0], row["p_med_p025"].iloc[0], row["p_med_p975"].iloc[0]])
            ax.errorbar(est, y, xerr=finite_xerr(est, lo, hi), fmt="s", color=ROLE_COLOR["confirmation"], ms=5)
            yticks.append(y)
            ylabels.append(f"{clock_label(clock)} | L1")
            y += 1
        y += 0.4
    ax.axvline(0, color="grey", lw=0.5)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("p_med")
    ax.set_title("Composition mediation")
    ax.invert_yaxis()
    med_handles = [
        plt.Line2D([0], [0], marker="o", color=ROLE_COLOR["development"], linestyle="", label="development"),
        plt.Line2D([0], [0], marker="s", color=ROLE_COLOR["confirmation"], linestyle="", label="confirmation"),
    ]
    ax.legend(handles=med_handles, frameon=False, fontsize=7)
    peaks = pd.concat([
        pd.read_csv(ROOT / "results/stage2/h2_population_cosinor.tsv", sep="\t"),
        pd.read_csv(ROOT / "results/stage4/rnaseq/h2_population_cosinor.tsv", sep="\t"),
    ], ignore_index=True)
    peaks = peaks.loc[peaks["population"].isin(BLOOD)]
    ax = fig.add_subplot(1, 2, 2)
    datasets = ["GSE39445", "GSE56931", "GSE113883"]
    for i, dataset in enumerate(datasets):
        block = peaks.loc[peaks["dataset"] == dataset]
        xs, ys = [], []
        for j, pop in enumerate(BLOOD):
            one = block.loc[block["population"] == pop]
            if one.empty:
                continue
            xs.append(float(one["peak_h"].iloc[0]))
            ys.append(j + (i - 1) * 0.18)
        ax.plot(xs, ys, marker="o", linestyle="", color=ROLE_COLOR["development"], alpha=0.4 + 0.25 * i, ms=6, label=dataset)
    ax.set_yticks(range(len(BLOOD)))
    ax.set_yticklabels(BLOOD, fontsize=7)
    ax.set_xlim(-0.4, 24.4)
    ax.set_xlabel("Peak hour")
    ax.set_title("MCPcounter peak hour")
    ax.legend(frameon=False, fontsize=7)
    fig.suptitle("Composition mediation (post hoc, exploratory)")
    fig.tight_layout()
    save(fig, "Fig4_composition")
    n_bits = []
    for dataset, sub in peaks.groupby("dataset", sort=False):
        n_bits.append(f"{dataset} {int(sub['n_subjects'].iloc[0])}/{int(sub['n_samples'].iloc[0])}")
    if not l1.empty:
        n_bits.append(f"L1 H2 {l1['n_samples'].iloc[0]} samples")
    return [{
        "figure": "Fig4", "panel": "p_med and cell peaks",
        "n_subjects": "; ".join(n_bits), "n_samples": "subjects/samples in n_subjects",
        "role": "development and confirmation",
        "source": "results/stage5/projection_ci.tsv; results/stage3/h2/projection_ci.tsv; results/stage4/rnaseq/projection_ci.tsv; results/confirm/h2.tsv; results/stage2/h2_population_cosinor.tsv; results/stage4/rnaseq/h2_population_cosinor.tsv",
    }]


def fig5() -> list[dict]:
    f3 = pd.read_csv(ROOT / "prereg/frozen/F3.tsv", sep="\t")
    s4 = pd.read_csv(ROOT / "results/confirm/s4.tsv", sep="\t", dtype=str)
    d11 = pd.read_csv(ROOT / "results/stage8/l1_f3_d11_model_m.tsv", sep="\t", dtype=str)
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    for i, (_, row) in enumerate(f3.iterrows()):
        hour = float(row["phi_bar"]) % 24
        ang = hour * math.pi / 12
        frozen = s4.loc[(s4["dataset"] == "L1") & (s4["clock"] == row["clock"]), "confirmed"]
        frozen_yes = bool(len(frozen) and frozen.iloc[0] == "TRUE")
        d11_yes = d11.loc[d11["clock"] == row["clock"], "posthoc_confirmed"].iloc[0] == "TRUE"
        l3 = s4.loc[(s4["dataset"] == "L3") & (s4["clock"] == row["clock"])]
        sign_yes = bool(len(l3) and l3["sign_agrees"].iloc[0] == "TRUE")
        ax.plot([ang], [1], marker="o", color="#1b9e77" if frozen_yes else "#bbbbbb", ms=9)
        ax.plot([ang], [0.78], marker="s", color="#d95f02" if d11_yes else "#bbbbbb", ms=7)
        ax.plot([ang], [0.56], marker="D", color="#7570b3" if sign_yes else "#bbbbbb", ms=6)
        screen = math.radians(90 - hour * 15)
        reach = 28 + 10 * (i % 4)
        dx = reach * math.cos(screen)
        dy = reach * math.sin(screen)
        label = row["clock"].replace("CMP_", "").replace("TAms_", "")
        ax.annotate(
            label, xy=(ang, 1), xytext=(dx, dy), textcoords="offset points",
            fontsize=6, ha="center", va="center",
            arrowprops={"arrowstyle": "-", "lw": 0.4, "color": "#888888"},
        )
    ax.set_xticks([h * math.pi / 12 for h in range(0, 24, 3)])
    ax.set_xticklabels([str(h) for h in range(0, 24, 3)])
    ax.set_yticks([])
    ax.plot([], [], marker="o", color="#1b9e77", linestyle="", label="L1 frozen confirmed, outer")
    ax.plot([], [], marker="s", color="#d95f02", linestyle="", label="L1 D11 post-hoc confirmed, middle")
    ax.plot([], [], marker="D", color="#7570b3", linestyle="", label="L3 sign agrees, inner")
    ax.plot([], [], marker="o", color="#bbbbbb", linestyle="", label="status false")
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.08, 1.02))
    ax.set_title("F3 module peaks (post hoc, exploratory)")
    fig.tight_layout(rect=(0, 0, 0.72, 1))
    save(fig, "Fig5_modules")
    return [{
        "figure": "Fig5", "panel": "F3", "n_subjects": "L1 22; L3 58", "n_samples": "L1 147; L3 116",
        "role": "confirmation and post hoc",
        "source": "prereg/frozen/F3.tsv; results/confirm/s4.tsv; results/stage8/l1_f3_d11_model_m.tsv",
    }]


def effect_table() -> pd.DataFrame:
    boot = pd.read_csv(ROOT / "results/stage8/cosinor_effects_boot.tsv", sep="\t", dtype=str)
    rows = []
    for _, row in boot.iterrows():
        rows.append({
            "dataset": row["dataset"], "clock": row["clock"], "role": row["family"], "kind": "cosinor",
            "n_subjects": row["n_subjects"], "n_samples": row["n_samples"],
            "estimate": row["two_a_over_sigma"], "low": row["ratio_p025"], "high": row["ratio_p975"],
            "r2_within": row["r2_within"], "r2_low": row["r2_p025"], "r2_high": row["r2_p975"],
            "peak_h": row["peak_h"], "peak_center": row["peak_center"], "peak_halfwidth_95": row["peak_halfwidth_95"],
            "p": row["p"], "title": "post hoc, exploratory",
        })
    fam = pd.concat([
        pd.read_csv(ROOT / "results/confirm/family_P.tsv", sep="\t", dtype=str),
        pd.read_csv(ROOT / "results/confirm/family_S1.tsv", sep="\t", dtype=str),
    ], ignore_index=True)
    for _, row in fam.loc[fam["dataset"] == "L1"].iterrows():
        rows.append({
            "dataset": "L1", "clock": row["clock"], "role": "confirmation", "kind": "cosinor",
            "n_subjects": row["n_subjects"], "n_samples": row["n_samples"],
            "estimate": row["two_a_over_sigma"], "low": row["ratio_p025"], "high": row["ratio_p975"],
            "r2_within": row["r2_within"], "r2_low": row["r2_p025"], "r2_high": row["r2_p975"],
            "peak_h": row["peak_h"], "peak_center": row["peak_center"], "peak_halfwidth_95": row["peak_halfwidth_95"],
            "p": row["p"], "title": "post hoc, exploratory",
        })
    dtab = pd.read_csv(ROOT / "results/stage8/l3_d_vs_prediction.tsv", sep="\t", dtype=str)
    for _, row in fam.loc[fam["dataset"] == "L3"].iterrows():
        drow = dtab.loc[dtab["clock"] == row["clock"]].iloc[0]
        rows.append({
            "dataset": "L3", "clock": row["clock"], "role": "confirmation", "kind": "paired",
            "n_subjects": row["n_subjects"], "n_samples": row["n_samples"],
            "estimate": drow["d"], "low": drow["d_p025"], "high": drow["d_p975"],
            "mean_delta": row["mean_delta"], "mean_low": row["ci_low"], "mean_high": row["ci_high"],
            "p": row["p"], "peak_h": "", "title": "post hoc, exploratory",
        })
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results/stage8/effect_summary.tsv", sep="\t", index=False)
    return out


def main() -> None:
    effects = effect_table()
    sources = []
    sources.extend(fig2())
    sources.extend(fig3(effects))
    sources.extend(fig4())
    sources.extend(fig5())
    pd.DataFrame(sources).to_csv(FIG / "figure_sources.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
