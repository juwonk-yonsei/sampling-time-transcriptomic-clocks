#!/usr/bin/env python3
"""Manuscript figures for TOD instruction 10. No new hypothesis tests."""

import math
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
OUT = ROOT / "figures" / "manuscript"
SRC = OUT / "source_data"
MM = 1 / 25.4
WIDTH_MM = 174

for face in (
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
):
    if Path(face).exists():
        font_manager.fontManager.addfont(face)
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
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
CLOCK_COLOR = {
    "TAms_chrono": "#0072B2",
    "TAms_mort": "#E69F00",
    "RA_default": "#009E73",
    "RA_Peters": "#CC79A7",
    "RA_GTEx": "#56B4E9",
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
SCENARIO_LW = {"lab": 1.0, "free-living": 1.0, "free-living upper": 0.6}


def load(path):
    return pd.read_csv(ROOT / path, sep="\t", dtype=str, keep_default_na=False)


def write_src(frame, name):
    SRC.mkdir(parents=True, exist_ok=True)
    frame.to_csv(SRC / name, sep="\t", index=False)


def num(text):
    return float(text)


def signed_diff(a, b):
    d = (a - b) % 24
    if d > 12:
        d -= 24
    return d


def panel(ax, letter):
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right", clip_on=False)


def new_fig(height_mm):
    if height_mm > 234:
        raise SystemExit(f"height {height_mm} mm exceeds 234")
    return plt.figure(figsize=(WIDTH_MM * MM, height_mm * MM))


def save_fig(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def geo_name(dataset):
    return GEO.get(dataset, dataset)


def one(frame, **eq):
    matched = frame
    for key, value in eq.items():
        matched = matched.loc[matched[key] == value]
    if len(matched) != 1:
        raise SystemExit(f"{len(matched)} rows for {eq}")
    return matched.iloc[0]


def fig1():
    freeze = load("results/stage7/preunlock_freeze_check.tsv")
    sha = one(freeze, path="prereg/manifest.tsv").status
    commit = "8a203ffdc5a5f1e933cd4f64154c0d5972480e32"
    when = subprocess.check_output(["git", "log", "-1", "--format=%ci", commit], cwd=ROOT, text=True).strip()
    effects = load("results/stage8/effect_summary.tsv")
    gates = {
        "GSE48113": load("results/confirm/L1/gate.tsv"),
        "GSE122725": load("results/confirm/L2/gate.tsv"),
        "GSE295758": load("results/confirm/L3/gate.tsv"),
    }
    f3 = load("prereg/frozen/F3.tsv")
    rows = []
    counts = {}
    for dataset, geo in [("GSE39445", "GSE39445"), ("GSE56931", "GSE56931"), ("GSE113883", "GSE113883"), ("L1", "GSE48113"), ("L2", "GSE122725"), ("L3", "GSE295758")]:
        block = effects.loc[effects["dataset"] == dataset]
        ns = sorted(set(zip(block["n_subjects"], block["n_samples"])))
        if len(ns) != 1:
            raise SystemExit(f"sample size not unique for {dataset}: {ns}")
        counts[geo] = ns[0]
        rows.append({"item": "sample_size", "dataset": geo, "n_subjects": ns[0][0], "n_samples": ns[0][1], "source": "results/stage8/effect_summary.tsv"})
    rows.append({"item": "manifest_sha256", "dataset": "", "n_subjects": "", "n_samples": "", "source": "results/stage7/preunlock_freeze_check.tsv", "value": sha})
    rows.append({"item": "freeze_commit_time", "dataset": "", "n_subjects": "", "n_samples": "", "source": commit, "value": when})
    rows.append({"item": "f3_rows", "dataset": "", "n_subjects": str(len(f3)), "n_samples": "", "source": "prereg/frozen/F3.tsv", "value": str(len(f3))})
    for geo, gate in gates.items():
        for _, grow in gate.iterrows():
            rec = {"item": "gate", "dataset": geo, "source": f"results/confirm gate {geo}"}
            rec.update({c: grow[c] for c in gate.columns})
            rows.append(rec)
    write_src(pd.DataFrame(rows).fillna(""), "Fig1_A.tsv")

    fig = new_fig(168)
    ax = fig.add_axes([0.03, 0.04, 0.94, 0.90])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel(ax, "A")

    def box(x, y, w, h, text, fc, ec="0.2", ls="-"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.008", facecolor=fc, edgecolor=ec, linestyle=ls, linewidth=0.8, mutation_aspect=0.45))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=7)

    def arrow(x1, x2, y):
        ax.annotate("", xy=(x2, y), xytext=(x1, y), arrowprops={"arrowstyle": "->", "lw": 0.7, "color": "0.2"})

    c394, c569, c113 = counts["GSE39445"], counts["GSE56931"], counts["GSE113883"]
    box(0.02, 0.78, 0.22, 0.14, f"Development\nGSE39445 {c394[0]}/{c394[1]}\nGSE56931 {c569[0]}/{c569[1]}\nGSE113883 {c113[0]}/{c113[1]}", "#D9F0E3")
    arrow(0.24, 0.28, 0.85)
    box(0.28, 0.80, 0.16, 0.10, "5 clocks", "#D9F0E3")
    arrow(0.44, 0.48, 0.85)
    box(0.48, 0.80, 0.18, 0.10, "Mixed cosinor\nmodel", "#D9F0E3")
    arrow(0.66, 0.70, 0.85)
    box(0.70, 0.76, 0.28, 0.18, f"Frozen predictions\npeak, amplitude,\nmorning−evening d,\n{len(f3)} module results", "#D9F0E3")
    ax.text(0.02, 0.94, "Development", fontsize=8, fontweight="bold", color="#009E73")

    c481, c122, c295 = counts["GSE48113"], counts["GSE122725"], counts["GSE295758"]
    gtxt = []
    for geo, word in (("GSE48113", "pass"), ("GSE122725", "fail"), ("GSE295758", "pass")):
        flag = gates[geo]["dataset_pass"].iloc[0]
        expect = "TRUE" if word == "pass" else "FALSE"
        if flag != expect:
            raise SystemExit(f"{geo} dataset_pass {flag}")
        genes = ",".join(gates[geo]["gene"].tolist())
        gtxt.append(f"{geo} {word} ({genes})")
    box(0.02, 0.44, 0.28, 0.18, f"Freeze\n{when[:10]}", "#F6D6C9")
    ax.text(0.02, 0.40, sha[:32] + "\n" + sha[32:], fontsize=7, va="top", ha="left", family="DejaVu Sans")
    arrow(0.32, 0.36, 0.54)
    box(0.36, 0.42, 0.28, 0.24, "Opened after freeze\n" + f"GSE48113 {c481[0]}/{c481[1]}\n" + f"GSE122725 {c122[0]}/{c122[1]}\n" + f"GSE295758 {c295[0]} pairs", "#F6D6C9", ec="#D55E00", ls="--")
    arrow(0.64, 0.68, 0.54)
    box(0.68, 0.42, 0.30, 0.24, "Gate\n" + "\n".join(gtxt) + "\nHolm primary family\nReplication: ≥2 datasets", "#F6D6C9")
    ax.text(0.02, 0.68, "Preregistration and locked confirmation", fontsize=8, fontweight="bold", color="#D55E00")

    box(0.02, 0.08, 0.22, 0.16, "Time-axis\ncorrection\n(GSE122725)", "#F7D6EA")
    box(0.27, 0.08, 0.22, 0.16, "Missing-value\nsensitivity", "#F7D6EA")
    box(0.52, 0.08, 0.22, 0.16, "Bias\nsimulation", "#F7D6EA")
    box(0.77, 0.08, 0.21, 0.16, "Mitigation", "#F7D6EA")
    ax.text(0.02, 0.28, "Post hoc", fontsize=8, fontweight="bold", color="#CC79A7")
    save_fig(fig, "Fig1_design")


def score_frame(dataset):
    if dataset == "GSE113883":
        samples = load("results/stage3/GSE113883_samples.tsv")
        samples = samples.loc[samples["condition"] == "control"].copy()
        tage = load("results/stage5/clock_scores.tsv")
        tage = tage.loc[(tage["dataset"] == "GSE113883") & (tage["clock"].isin(["TAms_chrono", "TAms_mort"])), ["gsm", "subject_id", "clock", "value"]]
        ra = load("results/stage4/rnaseq/clock_scores.tsv")
        ra = ra.loc[ra["method"] == "lin", ["gsm", "clock", "value"]]
        scores = pd.concat([tage, ra], ignore_index=True)
        order_col = "elapsed_h"
    elif dataset == "GSE48113":
        samples = load("results/confirm/L1/samples_scoring.tsv")
        samples = samples.loc[(samples["period"] == "D1") & (samples["analysis_set"] == "day1")].copy()
        tage = load("results/confirm/L1/tage_scores.tsv")
        ra = load("results/confirm/L1/ra_scores.tsv")
        scores = pd.concat([
            tage.loc[tage["clock"].isin(PRIMARY), ["gsm", "clock", "value"]],
            ra.loc[ra["clock"].isin(PRIMARY), ["gsm", "clock", "value"]],
        ], ignore_index=True)
        order_col = "elapsed_h"
    else:
        raise SystemExit(dataset)
    samples["gsm"] = samples["gsm"].astype(str)
    scores["gsm"] = scores["gsm"].astype(str)
    out = []
    for clock in PRIMARY:
        sc = scores.loc[scores["clock"] == clock, ["gsm", "value"]]
        merged = samples.merge(sc, on="gsm", how="inner")
        if merged["gsm"].duplicated().any():
            raise SystemExit(f"duplicate gsm {dataset} {clock}")
        merged["clock"] = clock
        merged["dataset"] = dataset
        merged["order_h"] = merged[order_col].map(num)
        merged["clock_time_h"] = merged["clock_time_h"].map(num)
        merged["value"] = merged["value"].map(num)
        merged["centered"] = merged["value"] - merged.groupby("subject_id")["value"].transform("mean")
        out.append(merged)
    return pd.concat(out, ignore_index=True)


def betas_for(dataset, clock):
    if dataset == "GSE113883" and clock in ("TAms_chrono", "TAms_mort"):
        tab = load("results/stage5/h1_point_estimates.tsv")
        row = one(tab, dataset=dataset, clock=clock)
        return row, "results/stage5/h1_point_estimates.tsv"
    if dataset == "GSE113883":
        tab = load("results/stage4/rnaseq/h1_point_estimates.tsv")
        row = one(tab, dataset=dataset, clock=clock, method="lin")
        return row, "results/stage4/rnaseq/h1_point_estimates.tsv"
    tab = load("results/confirm/family_P.tsv")
    row = one(tab, dataset="L1", clock=clock)
    return row, "results/confirm/family_P.tsv"


def fig2():
    frames = {ds: score_frame(ds) for ds in ("GSE113883", "GSE48113")}
    letters = {("GSE113883", clock): letter for letter, clock in zip("ABC", PRIMARY)}
    letters.update({("GSE48113", clock): letter for letter, clock in zip("DEF", PRIMARY)})
    fig = new_fig(150)
    axes = fig.subplots(2, 3)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.93, bottom=0.08, hspace=0.48, wspace=0.32)
    summary = []
    for r, dataset in enumerate(("GSE113883", "GSE48113")):
        block = frames[dataset]
        n_subj = block.groupby("clock")["subject_id"].nunique()
        if dataset == "GSE113883" and set(n_subj) != {11}:
            raise SystemExit(f"GSE113883 subjects {n_subj.to_dict()}")
        if dataset == "GSE48113" and set(n_subj) != {22}:
            raise SystemExit(f"GSE48113 subjects {n_subj.to_dict()}")
        for c, clock in enumerate(PRIMARY):
            ax = axes[r, c]
            sub = block.loc[block["clock"] == clock].copy()
            beta, beta_src = betas_for(dataset, clock)
            records = []
            for subject, person in sub.groupby("subject_id"):
                person = person.sort_values(["order_h", "gsm"])
                times = person["clock_time_h"].to_numpy()
                ys = person["centered"].to_numpy()
                ax.plot(times, ys, color="0.45", alpha=0.25, lw=0.5, marker="o", ms=1.5)
                start = 0
                for i in range(1, len(times) + 1):
                    cut = i == len(times) or times[i] <= times[i - 1]
                    if cut:
                        records.append({"record": "segment", "subject_id": subject, "i_start": start, "i_end": i - 1, "broken_before_next": i != len(times)})
                        start = i
                for _, pt in person.iterrows():
                    records.append({"record": "point", "subject_id": subject, "gsm": pt["gsm"], "clock_time_h": pt["clock_time_h"], "order_h": pt["order_h"], "value": pt["value"], "centered": pt["centered"]})
            # redraw lines with explicit breaks; the first plot connected everything
            ax.cla()
            for subject, person in sub.groupby("subject_id"):
                person = person.sort_values(["order_h", "gsm"])
                times = person["clock_time_h"].to_numpy()
                ys = person["centered"].to_numpy()
                start = 0
                for i in range(1, len(times) + 1):
                    if i == len(times) or times[i] <= times[i - 1]:
                        if i - start >= 2:
                            ax.plot(times[start:i], ys[start:i], color="0.45", alpha=0.25, lw=0.5)
                        ax.plot(times[start:i], ys[start:i], color="0.45", alpha=0.25, lw=0, marker="o", ms=1.5)
                        start = i
            bin_rows = []
            sub = sub.copy()
            sub["bin_start"] = np.floor(np.mod(sub["clock_time_h"], 24) / 2.0) * 2
            for b0, chunk in sub.groupby("bin_start"):
                per = chunk.groupby("subject_id")["centered"].mean()
                mean = float(per.mean())
                n = int(per.shape[0])
                sem = float(per.std(ddof=1) / math.sqrt(n)) if n >= 2 else ""
                bin_rows.append((b0 + 1, mean, sem, n))
                records.append({"record": "bin", "bin_start": b0, "bin_center": b0 + 1, "n_subjects": n, "mean": mean, "sem": sem})
            if bin_rows:
                xs, ys, yerr, ns = zip(*[(a, b, 0 if c == "" else c, d) for a, b, c, d in bin_rows])
                ax.errorbar(xs, ys, yerr=yerr, color=DATA_COLOR[dataset], lw=1.0, marker=DATA_MARKER[dataset], ms=3.5, capsize=1.5, zorder=3)
            t = np.linspace(0, 24, 241)
            curve = num(beta.beta_cos) * np.cos(2 * np.pi * t / 24) + num(beta.beta_sin) * np.sin(2 * np.pi * t / 24)
            curve = curve - curve.mean()
            ax.plot(t, curve, color="black", lw=1.0, zorder=4)
            for tt, yy in zip(t[::8], curve[::8]):
                records.append({"record": "curve", "curve_t": tt, "curve_y": yy, "beta_source": beta_src, "beta_cos": beta.beta_cos, "beta_sin": beta.beta_sin})
            p_text = beta["p"] if "p" in beta.index and beta["p"] != "" else beta["lrt_p"]
            ax.text(0.03, 0.97, letters[(dataset, clock)], transform=ax.transAxes, fontsize=9, fontweight="bold", va="top", ha="left")
            ax.text(0.98, 0.04, f"p = {p_text}\npeak {beta.peak_h} h", transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 0.3})
            records.append({"record": "annotation", "p": p_text, "peak_h": beta.peak_h, "beta_source": beta_src, "n_subjects": int(sub['subject_id'].nunique()), "n_samples": int(len(sub)), "n_bins": len(bin_rows)})
            ylab = "years (centered)" if clock == "RA_default" else "score (centered)"
            ax.set_xlim(0, 24)
            ax.set_xticks([0, 6, 12, 18, 24])
            if r == 1:
                ax.set_xlabel("Clock time (h)")
            if c == 0:
                ax.set_ylabel(ylab)
            ax.set_title(f"{dataset}  {CLOCK_NAME[clock]}", fontsize=8, pad=10)
            src = pd.DataFrame(records).fillna("")
            src.insert(0, "dataset", dataset)
            src.insert(1, "clock", clock)
            write_src(src, f"Fig2_{letters[(dataset, clock)]}.tsv")
            summary.append({"panel": letters[(dataset, clock)], "dataset": dataset, "clock": clock, "n_subjects": int(sub["subject_id"].nunique()), "n_samples": int(len(sub)), "n_bins": len(bin_rows), "p": p_text, "peak_h": beta.peak_h, "beta_source": beta_src})
    save_fig(fig, "Fig2_trajectories")
    write_src(pd.DataFrame(summary), "Fig2_summary.tsv")


def forest_ratio(ax, clocks, include_posthoc):
    effects = load("results/stage8/effect_summary.tsv")
    effects = effects.loc[effects["kind"] == "cosinor"].copy()
    order = ["GSE39445", "GSE56931", "GSE113883", "L1"] + (["L2"] if include_posthoc else [])
    rows = []
    y = 0
    yticks = []
    ylabels = []
    groups = []
    for clock in clocks:
        y0 = y
        for dataset in order:
            row = one(effects, dataset=dataset, clock=clock)
            geo = geo_name(dataset)
            lo, hi, est = num(row.low), num(row.high), num(row.estimate)
            ax.plot([lo, hi], [y, y], color=DATA_COLOR[geo], lw=0.8)
            ax.plot(est, y, marker=DATA_MARKER[geo], color=DATA_COLOR[geo], ms=4)
            yticks.append(y)
            ylabels.append(geo)
            rows.append({"dataset": geo, "internal_dataset": dataset, "clock": clock, "estimate": row.estimate, "low": row.low, "high": row.high, "role": row.role, "y": y})
            y += 1
        groups.append((clock, (y0 + y - 1) / 2))
        y += 0.6
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_ylim(-0.6, y)
    ax.set_xlim(left=0)
    ax.set_xlabel("2A/σ")
    ax.invert_yaxis()
    for clock, ypos in groups:
        ax.text(1.01, ypos, CLOCK_NAME[clock], transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=7, clip_on=False)
    if include_posthoc:
        ax.legend(handles=[Line2D([0], [0], marker="D", color=DATA_COLOR["GSE122725"], lw=0, label="GSE122725 (post hoc)")], frameon=False, loc="lower right")
    return pd.DataFrame(rows)


def fig3(clocks, stem, letters):
    fig = new_fig(220 if clocks == PRIMARY else 170)
    gs = fig.add_gridspec(3 if clocks == PRIMARY else 2, len(clocks), height_ratios=([1.3, 0.7, 1.15] if clocks == PRIMARY else [1.3, 0.8]), hspace=0.45, wspace=0.35, left=0.16, right=0.84, top=0.94, bottom=0.07)
    axa = fig.add_subplot(gs[0, :])
    panel(axa, letters[0])
    src_a = forest_ratio(axa, clocks, include_posthoc=True)
    write_src(src_a, f"{stem.split('_')[0]}_{letters[0]}.tsv" if False else f"Fig{stem}_{letters[0]}.tsv")
    # source names are set by caller via rewrite below
    axa_src_name = f"Fig{ '3' if clocks==PRIMARY else 'S1' }_{letters[0]}.tsv"
    write_src(src_a, axa_src_name)

    effects = load("results/stage8/effect_summary.tsv")
    paired = effects.loc[(effects["dataset"] == "L3") & (effects["kind"] == "paired")]
    f2 = load("prereg/frozen/F2.tsv")
    axb = fig.add_subplot(gs[1, :])
    panel(axb, letters[1])
    brows = []
    for i, clock in enumerate(clocks):
        obs = one(paired, clock=clock)
        pred = one(f2, clock=clock)
        dev = [num(pred[c]) for c in ("d_GSE39445", "d_GSE56931", "d_GSE113883")]
        axb.plot([min(dev), max(dev)], [i - 0.18, i - 0.18], color="0.55", lw=1.2, solid_capstyle="butt", zorder=1)
        axb.plot(num(pred.d_exp_median), i - 0.18, marker="D", markerfacecolor="white", markeredgecolor="black", ms=5, zorder=3)
        axb.plot([num(obs.low), num(obs.high)], [i, i], color=DATA_COLOR["GSE295758"], lw=0.8, zorder=2)
        axb.plot(num(obs.estimate), i, marker="s", color=DATA_COLOR["GSE295758"], ms=4.5, zorder=4)
        brows.append({"clock": clock, "d": obs.estimate, "low": obs.low, "high": obs.high, "d_exp_median": pred.d_exp_median, "d_GSE39445": pred.d_GSE39445, "d_GSE56931": pred.d_GSE56931, "d_GSE113883": pred.d_GSE113883})
    axb.axvline(0, color="black", lw=0.5)
    axb.set_yticks(range(len(clocks)))
    axb.set_yticklabels([CLOCK_NAME[c] for c in clocks])
    axb.set_xlim(-0.2, 0.8)
    axb.set_xlabel("d")
    axb.invert_yaxis()
    write_src(pd.DataFrame(brows), f"Fig{'3' if clocks==PRIMARY else 'S1'}_{letters[1]}.tsv")

    if clocks == PRIMARY:
        f1 = load("prereg/frozen/F1.tsv")
        point = load("results/stage8/effect_summary.tsv")
        point = point.loc[point["kind"] == "cosinor"]
        crows = []
        for i, clock in enumerate(clocks):
            ax = fig.add_subplot(gs[2, i], projection="polar")
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
            ax.set_xticks([k / 24 * 2 * np.pi for k in (0, 6, 12, 18)])
            ax.set_xticklabels(["0", "6", "12", "18"])
            ax.set_ylim(0, 1.35)
            ax.set_yticks([0.8, 1.0])
            ax.set_yticklabels([])
            phi = one(f1.loc[f1["dataset"] == "GSE39445"], clock=clock).phi_bar
            ax.plot([num(phi) / 24 * 2 * np.pi], [1.15], marker="*", color="black", ms=8, zorder=4)
            for dataset, radius in (("GSE39445", 1.0), ("GSE56931", 1.0), ("GSE113883", 1.0), ("L1", 0.8)):
                row = one(point, dataset=dataset, clock=clock)
                geo = geo_name(dataset)
                center = num(row.peak_center)
                half = num(row.peak_halfwidth_95)
                hours = np.linspace(center - half, center + half, 180)
                ax.plot(hours / 24 * 2 * np.pi, np.full(hours.shape, radius), color=DATA_COLOR[geo], lw=1.0)
                ax.plot([num(row.peak_h) / 24 * 2 * np.pi], [radius], marker=DATA_MARKER[geo], color=DATA_COLOR[geo], ms=4, zorder=3)
                crows.append({"clock": clock, "dataset": geo, "radius": radius, "phi_bar": phi, "peak_h": row.peak_h, "peak_center": row.peak_center, "peak_halfwidth_95": row.peak_halfwidth_95})
            ax.set_title(CLOCK_NAME[clock], fontsize=8, pad=10)
            if i == 0:
                ax.text(-0.22, 1.08, "C", transform=ax.transAxes, fontsize=9, fontweight="bold", clip_on=False)
        write_src(pd.DataFrame(crows), "Fig3_C.tsv")
        handles = [
            Line2D([0], [0], marker="o", color="0.3", lw=0, label="development, r = 1.0"),
            Line2D([0], [0], marker="s", color=DATA_COLOR["GSE48113"], lw=0, label="GSE48113, r = 0.8"),
            Line2D([0], [0], marker="*", color="black", lw=0, label="predicted peak, r = 1.15"),
        ]
        fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    save_fig(fig, "Fig3_effects" if clocks == PRIMARY else "FigS1_effects_secondary")


def fig4a():
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
    h2 = load("results/confirm/h2.tsv")
    rows = []
    fig = new_fig(200)
    ax = fig.add_axes([0.28, 0.55, 0.62, 0.38])
    panel(ax, "A")
    y = 0
    for dataset in ("GSE39445", "GSE56931", "GSE113883"):
        for clock in PRIMARY + SECONDARY:
            hit = dev.loc[(dev["dataset"] == dataset) & (dev["clock"] == clock)]
            if hit.empty:
                continue
            if len(hit) != 1:
                raise SystemExit(f"projection rows {dataset} {clock}")
            row = hit.iloc[0]
            lo, hi, est = num(row.p_med_p025), num(row.p_med_p975), num(row.p_med)
            draw_clipped(ax, y, est, lo, hi, -0.5, 1.5, DATA_COLOR[dataset], DATA_MARKER[dataset])
            ax.text(-0.02, y, f"{dataset} {CLOCK_NAME[clock]}", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7)
            rows.append({"panel": "A", "dataset": dataset, "clock": clock, "p_med": row.p_med, "low": row.p_med_p025, "high": row.p_med_p975, "lrt_p": row.lrt_p_unadj, "source": row.source, "y": y})
            y += 1
    for _, row in h2.iterrows():
        lo, hi, est = num(row.p_med_p025), num(row.p_med_p975), num(row.p_med)
        draw_clipped(ax, y, est, lo, hi, -0.5, 1.5, DATA_COLOR["GSE48113"], "s")
        ax.text(-0.02, y, f"GSE48113 {CLOCK_NAME[row.clock]}", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7)
        rows.append({"panel": "A", "dataset": "GSE48113", "clock": row.clock, "p_med": row.p_med, "low": row.p_med_p025, "high": row.p_med_p975, "lrt_p": row.lrt_p, "source": "results/confirm/h2.tsv", "y": y})
        y += 1
    ax.axvline(0, color="black", lw=0.5)
    ax.axvline(1, color="0.4", lw=0.5, ls="--")
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.6, y - 0.4)
    ax.set_yticks([])
    ax.set_xlabel("p_med")
    ax.invert_yaxis()
    return fig, rows


def draw_clipped(ax, y, est, lo, hi, xmin, xmax, color, marker):
    show_lo = min(max(lo, xmin), xmax)
    show_hi = max(min(hi, xmax), xmin)
    ax.plot([show_lo, show_hi], [y, y], color=color, lw=0.8, zorder=2)
    span = xmax - xmin
    if lo < xmin:
        ax.annotate("", xy=(xmin, y), xytext=(xmin + 0.06 * span, y), arrowprops={"arrowstyle": "->", "color": color, "lw": 0.7})
    if hi > xmax:
        ax.annotate("", xy=(xmax, y), xytext=(xmax - 0.06 * span, y), arrowprops={"arrowstyle": "->", "color": color, "lw": 0.7})
    x = min(max(est, xmin), xmax)
    ax.plot(x, y, marker=marker, color=color, ms=4, zorder=3)


def to_clock_axis(hour):
    hour = hour % 24
    return hour if hour >= 12 else hour + 24


def fig4():
    fig, rows = fig4a()
    peaks = load("results/stage10/fig4b_composition_peaks.tsv")
    ax = fig.add_axes([0.28, 0.08, 0.48, 0.38])
    panel(ax, "B")
    datasets = ["GSE39445", "GSE56931", "GSE113883", "GSE48113"]
    brows = []
    for i, population in enumerate(POPULATIONS):
        for j, dataset in enumerate(datasets):
            row = one(peaks, dataset=dataset, population=population)
            y = i + (j - 1.5) * 0.16
            center = num(row.peak_center)
            half = num(row.peak_halfwidth_95)
            cx = to_clock_axis(center)
            lo, hi = cx - half, cx + half
            open_mark = num(row.lrt_p) >= 0.05
            draw_time(ax, y, to_clock_axis(num(row.peak_h)), lo, hi, DATA_COLOR[dataset], DATA_MARKER[dataset], open_mark)
            brows.append({"dataset": dataset, "population": population, "y": y, "open_marker": open_mark, "axis_peak": to_clock_axis(num(row.peak_h)), "axis_low": lo, "axis_high": hi, **{c: row[c] for c in peaks.columns}})
    ax.set_yticks(range(len(POPULATIONS)))
    ax.set_yticklabels(POPULATIONS)
    ax.set_xlim(12, 36)
    ax.set_xticks([12, 18, 24, 30, 36])
    ax.set_xticklabels(["12", "18", "0", "6", "12"])
    ax.set_xlabel("Clock time (h)")
    ax.set_ylim(-0.6, len(POPULATIONS) - 0.4)
    ax.invert_yaxis()
    handles = [Line2D([0], [0], marker=DATA_MARKER[d], color=DATA_COLOR[d], lw=0, label=(d + " (post hoc)" if d == "GSE48113" else d)) for d in datasets]
    handles.append(Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor="0.2", lw=0, label="LRT p ≥ 0.05"))
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1))
    write_src(pd.DataFrame(rows), "Fig4_A.tsv")
    write_src(pd.DataFrame(brows), "Fig4_B.tsv")
    save_fig(fig, "Fig4_composition")


def draw_time(ax, y, peak, lo, hi, color, marker, open_mark):
    show_lo = min(max(lo, 12), 36)
    show_hi = max(min(hi, 36), 12)
    ax.plot([show_lo, show_hi], [y, y], color=color, lw=0.7, zorder=2)
    if lo < 12:
        ax.annotate("", xy=(12, y), xytext=(13.2, y), arrowprops={"arrowstyle": "->", "color": color, "lw": 0.6})
    if hi > 36:
        ax.annotate("", xy=(36, y), xytext=(34.8, y), arrowprops={"arrowstyle": "->", "color": color, "lw": 0.6})
    x = min(max(peak, 12), 36)
    ax.plot(x, y, marker=marker, markerfacecolor=("white" if open_mark else color), markeredgecolor=color, color=color, ms=3.5, zorder=3)


def fig5():
    f3 = load("prereg/frozen/F3.tsv")
    d11 = load("results/stage8/l1_f3_d11_model_m.tsv")
    s4 = load("results/confirm/s4.tsv")
    feat = load("results/stage8/f3_feature_counts.tsv")
    fig = new_fig(210)
    ax = fig.add_axes([0.10, 0.48, 0.82, 0.46])
    panel(ax, "A")
    arows = []
    xs, ys = [], []
    for clock in MODULES:
        phi = num(one(f3, clock=clock).phi_bar)
        post = one(d11, clock=clock)
        diff = signed_diff(num(post.peak_h), phi)
        y = phi + diff
        p = num(post.p)
        filled = p < 0.05
        size = 14 * (-math.log10(p))
        ax.scatter([phi], [y], s=size, marker="o", facecolors=(DATA_COLOR["GSE48113"] if filled else "white"), edgecolors=DATA_COLOR["GSE48113"], zorder=3)
        nudge = {
            "TM-turquoise": (-72, 4), "chrono-turquoise": (6, -12), "chrono-yellow": (6, 8),
            "chrono-orange": (-78, 6), "chrono-red": (6, -10), "chrono-brown": (6, 4),
            "mort-turquoise": (6, 10), "mort-darkred": (6, 4),
        }[MODULE_SHORT[clock]]
        ax.annotate(MODULE_SHORT[clock], (phi, y), textcoords="offset points", xytext=nudge, fontsize=7, color=DATA_COLOR["GSE48113"])
        frozen = s4.loc[(s4["dataset"] == "L1") & (s4["clock"] == clock) & (s4["in_f3"] == "TRUE")]
        if len(frozen) != 1:
            raise SystemExit(clock)
        fr = frozen.iloc[0]
        nfeat = one(feat, version="frozen", clock=clock).n_entered
        if fr.p != "" and fr.peak_h != "":
            fp = num(fr.p)
            fy = phi + signed_diff(num(fr.peak_h), phi)
            ax.scatter([phi], [fy], s=14 * (-math.log10(fp)), marker="x", color="0.35", zorder=2)
            ax.annotate(nfeat, (phi, fy), textcoords="offset points", xytext=(5, -9), fontsize=7, color="0.35")
            frozen_y = fy
        else:
            frozen_y = ""
        xs.append(phi)
        ys.append(y)
        arows.append({"clock": clock, "short": MODULE_SHORT[clock], "phi_bar": one(f3, clock=clock).phi_bar, "d11_peak_h": post.peak_h, "d11_p": post.p, "signed_diff_h": diff, "y": y, "marker_area": size, "filled": filled, "frozen_peak_h": fr.peak_h, "frozen_p": fr.p, "frozen_y": frozen_y, "frozen_n_entered": nfeat, "marker_area_rule": "14 * -log10(p)"})
    plotted_y = [r["y"] for r in arows] + [r["frozen_y"] for r in arows if r["frozen_y"] != ""]
    lo = min(xs + plotted_y) - 2
    hi = max(xs + plotted_y) + 2
    grid = np.linspace(lo, hi, 50)
    ax.fill_between(grid, grid - 3, grid + 3, color="0.85", zorder=0)
    ax.plot(grid, grid, color="black", lw=0.6, zorder=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Predicted peak (h)")
    ax.set_ylabel("Predicted peak + signed circular difference (h)")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(handles=[
        Line2D([0], [0], marker="o", markerfacecolor=DATA_COLOR["GSE48113"], markeredgecolor=DATA_COLOR["GSE48113"], lw=0, label="GSE48113 D11 (post hoc), LRT p < 0.05"),
        Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor=DATA_COLOR["GSE48113"], lw=0, label="GSE48113 D11 (post hoc), LRT p ≥ 0.05"),
        Line2D([0], [0], marker="x", color="0.35", lw=0, label="frozen"),
    ], frameon=False, fontsize=7, loc="upper left")
    write_src(pd.DataFrame(arows), "Fig5_A.tsv")

    axb = fig.add_axes([0.07, 0.06, 0.78, 0.36])
    panel(axb, "B")
    axb.set_xlim(0, 12)
    axb.set_ylim(0, len(MODULES) + 1.6)
    axb.axis("off")
    headers = [(1.0, "module"), (3.2, "−log10 q"), (5.8, "frozen\nn"), (7.7, "D11 (post hoc)\nn"), (10.2, "GSE295758\npred / obs / agree")]
    for x, text in headers:
        axb.text(x, len(MODULES) + 1.2, text, ha="center", va="center", fontsize=7)
    for i, name in enumerate(("39445", "56931", "113883")):
        axb.text(2.35 + i * 0.95, len(MODULES) + 0.55, name, ha="center", va="center", fontsize=7)
    qvals = []
    for clock in MODULES:
        row = one(f3, clock=clock)
        for col in ("q_GSE39445", "q_GSE56931", "q_GSE113883"):
            qvals.append(-math.log10(num(row[col])))
    vmin, vmax = min(qvals), max(qvals)
    cmap = plt.get_cmap("Blues")
    brows = []
    for r, clock in enumerate(MODULES):
        y = len(MODULES) - r
        row = one(f3, clock=clock)
        axb.text(1.0, y, MODULE_SHORT[clock], ha="center", va="center", fontsize=7)
        for i, col in enumerate(("q_GSE39445", "q_GSE56931", "q_GSE113883")):
            lq = -math.log10(num(row[col]))
            color = cmap((lq - vmin) / (vmax - vmin if vmax > vmin else 1))
            axb.add_patch(plt.Rectangle((1.95 + i * 0.95, y - 0.35), 0.85, 0.7, color=color, ec="white", lw=0.4))
        fr = s4.loc[(s4["dataset"] == "L1") & (s4["clock"] == clock) & (s4["in_f3"] == "TRUE")].iloc[0]
        post = one(d11, clock=clock)
        l3 = s4.loc[(s4["dataset"] == "L3") & (s4["clock"] == clock) & (s4["in_f3"] == "TRUE")].iloc[0]
        nf = {ver: one(feat, version=ver, clock=clock).n_entered for ver in ("frozen", "D11_participant_mean")}
        agree = l3.sign_agrees if l3.sign_agrees != "" else ""
        axb.text(5.8, y, f"{fr.confirmed}\n{nf['frozen']}", ha="center", va="center", fontsize=7)
        axb.text(7.7, y, f"{post.posthoc_confirmed}\n{nf['D11_participant_mean']}", ha="center", va="center", fontsize=7)
        axb.text(10.2, y, f"{l3.predicted_sign} / {l3.observed_sign} / {agree}", ha="center", va="center", fontsize=7)
        brows.append({"clock": clock, "short": MODULE_SHORT[clock], "q_GSE39445": row.q_GSE39445, "q_GSE56931": row.q_GSE56931, "q_GSE113883": row.q_GSE113883, "frozen_confirmed": fr.confirmed, "frozen_n": nf["frozen"], "d11_confirmed": post.posthoc_confirmed, "d11_n": nf["D11_participant_mean"], "predicted_sign": l3.predicted_sign, "observed_sign": l3.observed_sign, "sign_agrees": l3.sign_agrees, "color_min_log10q": vmin, "color_max_log10q": vmax})
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin, vmax))
    cax = fig.add_axes([0.88, 0.10, 0.02, 0.22])
    fig.colorbar(sm, cax=cax).set_label("-log10 q", fontsize=7)
    write_src(pd.DataFrame(brows), "Fig5_B.tsv")
    save_fig(fig, "Fig5_modules")


def fig6():
    alias = load("results/stage8/gene_peak_differences_arntl_alias.tsv")
    wake = load("results/stage9/l2_wake_corr_genes.tsv")
    geo_p = load("results/posthoc_stage7c/l2_period1_gene_peaks.tsv")
    bmal = load("results/stage8/l2_bmal1_as_arntl.tsv")
    fig = new_fig(220)
    ax = fig.add_axes([0.10, 0.56, 0.48, 0.36])
    panel(ax, "A")
    arows = []
    for _, row in alias.iterrows():
        gene = row.gene
        x = num(row.peak_gse39445)
        y_geo = x + num(row.signed_diff_h)
        w = one(wake, gene=gene)
        y_cor = x + num(w.signed_diff_h)
        if gene == "ARNTL":
            p_geo = one(bmal, gene="BMAL1").lrt_p
            p_src = "results/stage8/l2_bmal1_as_arntl.tsv"
        else:
            p_geo = one(geo_p, gene=gene).lrt_p
            p_src = "results/posthoc_stage7c/l2_period1_gene_peaks.tsv"
        filled = num(w.lrt_p_wake_corr) < 0.05
        ax.scatter([x], [y_geo], marker="o", facecolors="none", edgecolors="0.45", s=18, zorder=2)
        ax.scatter([x], [y_cor], marker="D", facecolors=(DATA_COLOR["GSE122725"] if filled else "white"), edgecolors=DATA_COLOR["GSE122725"], s=18, zorder=3)
        ax.annotate(gene, (x, y_cor), textcoords="offset points", xytext={
            "PER1": (10, 2), "PER2": (8, -10), "PER3": (-28, 6), "NR1D1": (-36, 8), "NR1D2": (10, -4),
            "ARNTL": (6, 2), "CRY1": (6, -8), "CRY2": (6, 4), "DBP": (6, -8), "TEF": (-16, -8),
        }[gene], fontsize=7)
        arows.append({"gene": gene, "x_gse39445": row.peak_gse39445, "signed_diff_geo": row.signed_diff_h, "y_geo": y_geo, "signed_diff_corrected": w.signed_diff_h, "y_corrected": y_cor, "lrt_p_geo_axis": p_geo, "lrt_p_source": p_src, "lrt_p_wake_corr": w.lrt_p_wake_corr, "filled": filled})
    xs = [num(r["x_gse39445"]) for r in arows]
    ys = [r["y_geo"] for r in arows] + [r["y_corrected"] for r in arows]
    lo, hi = min(xs + ys) - 2, max(xs + ys) + 2
    grid = np.linspace(lo, hi, 40)
    ax.fill_between(grid, grid - 2, grid + 2, color="0.9", zorder=0)
    ax.plot(grid, grid, color="black", lw=0.6)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("GSE39445 peak (h)")
    ax.set_ylabel("GSE39445 peak + signed difference (h)")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(handles=[
        Line2D([0], [0], marker="o", markerfacecolor="none", markeredgecolor="0.45", lw=0, label="GEO axis"),
        Line2D([0], [0], marker="D", markerfacecolor=DATA_COLOR["GSE122725"], markeredgecolor=DATA_COLOR["GSE122725"], lw=0, label="relclocktime − 8 (post hoc), LRT p < 0.05"),
        Line2D([0], [0], marker="D", markerfacecolor="white", markeredgecolor=DATA_COLOR["GSE122725"], lw=0, label="relclocktime − 8 (post hoc), LRT p ≥ 0.05"),
    ], frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1))
    write_src(pd.DataFrame(arows), "Fig6_A.tsv")

    axb = fig.add_axes([0.24, 0.20, 0.30, 0.32])
    panel(axb, "B")
    tams = load("results/stage7c/tams_feature_counts.tsv")
    feat = load("results/stage8/f3_feature_counts.tsv")
    clocks = ["TAms_chrono", "TAms_mort"] + MODULES
    versions = [("frozen", "frozen", "o", "black"), ("D2_knn", "D2 kNN (post hoc)", "s", "#D55E00"), ("D11_participant_mean", "D11 participant mean (post hoc)", "D", "#CC79A7")]
    brows = []
    for i, clock in enumerate(clocks):
        for ver, label, marker, color in versions:
            if clock in ("TAms_chrono", "TAms_mort"):
                row = one(tams, version=ver, clock=clock)
                n = row.n_entered
                base = row.n_weight
            else:
                row = one(feat, version=ver, clock=clock)
                n = row.n_entered
                base = ""
            y = len(clocks) - 1 - i
            if n == "0":
                axb.annotate("0", xy=(0, y), xycoords=("axes fraction", "data"), xytext=(3, 0), textcoords="offset points", fontsize=7, color=color, va="center", ha="left")
                shown = ""
            else:
                axb.plot(num(n), y, marker=marker, color=color, ms=4)
                shown = n
            brows.append({"clock": clock, "version": ver, "n_entered": n, "n_weight": base, "plotted_x": shown, "y": y})
    axb.set_xscale("log")
    axb.set_xlim(0.4, 30000)
    axb.axvline(1839, color="0.4", lw=0.5, ls="--")
    axb.axvline(10487, color="0.4", lw=0.5, ls=":")
    axb.set_yticks(range(len(clocks)))
    axb.set_yticklabels([CLOCK_NAME[c] if c in CLOCK_NAME else MODULE_SHORT[c] for c in clocks][::-1])
    axb.set_xlabel("Features entered")
    axb.legend(handles=[Line2D([0], [0], marker=m, color=col, lw=0, label=lab) for _, lab, m, col in versions] + [
        Line2D([0], [0], color="0.4", lw=0.5, ls="--", label="1,839 nonzero weights"),
        Line2D([0], [0], color="0.4", lw=0.5, ls=":", label="10,487 nonzero weights"),
        Line2D([0], [0], color="none", lw=0, label="0 at the axis: n_entered = 0"),
    ], frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2)
    write_src(pd.DataFrame(brows), "Fig6_B.tsv")

    axc = fig.add_axes([0.68, 0.12, 0.29, 0.34])
    panel(axc, "C")
    knn = load("results/stage7b/d1_development_model_m.tsv")
    part = load("results/stage7c/d11_development_model_m.tsv")
    crows = []
    masks = ["random", "low_expression"]
    treatments = ["unmasked", "kNN", "participant mean"]
    xpos = []
    labels = []
    x = 0
    for mask in masks:
        for treat in treatments:
            xpos.append(x)
            labels.append({"unmasked": "before", "kNN": "kNN", "participant mean": "mean"}[treat])
            for k, clock in enumerate(("RA_default", "RA_Peters", "RA_GTEx")):
                if treat == "unmasked":
                    row = one(knn, scheme="unmasked", clock=clock)
                    src = "results/stage7b/d1_development_model_m.tsv"
                elif treat == "kNN":
                    row = one(knn, scheme=mask, clock=clock)
                    src = "results/stage7b/d1_development_model_m.tsv"
                else:
                    row = one(part, scheme=mask, clock=clock)
                    src = "results/stage7c/d11_development_model_m.tsv"
                lq = -math.log10(num(row.lrt_p))
                axc.plot(x + (k - 1) * 0.18, lq, marker="o", color=CLOCK_COLOR[clock], ms=4)
                crows.append({"mask": mask, "treatment": treat, "clock": clock, "lrt_p": row.lrt_p, "neg_log10_p": lq, "source": src, "x": x + (k - 1) * 0.18})
            x += 1
        x += 0.4
    axc.axhline(-math.log10(0.05), color="black", lw=0.5, ls="--")
    axc.set_xticks(xpos)
    axc.set_xticklabels(labels, fontsize=7)
    axc.text(xpos[1], -0.22, "random", transform=axc.get_xaxis_transform(), ha="center", fontsize=7)
    axc.text(xpos[4], -0.22, "low expression", transform=axc.get_xaxis_transform(), ha="center", fontsize=7)
    axc.set_ylabel("−log10 LRT p")
    axc.set_ylim(bottom=0)
    axc.legend(handles=[Line2D([0], [0], marker="o", color=CLOCK_COLOR[c], lw=0, label=CLOCK_NAME[c]) for c in ("RA_default", "RA_Peters", "RA_GTEx")], frameon=False, fontsize=7)
    write_src(pd.DataFrame(crows), "Fig6_C.tsv")
    save_fig(fig, "Fig6_public")


def wilson(k, n):
    from statistics import NormalDist
    z = NormalDist().inv_cdf(0.975)
    p = k / n
    den = 1 + z * z / n
    center = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return center - half, center + half, z


def fig7(clocks, stem, prefix):
    bias = load("results/stage9/sim_design1_bias.tsv")
    fpr1 = load("results/stage9/sim_design1_fpr.tsv")
    fpr2 = load("results/stage9/sim_design2_fpr.tsv")
    d3 = load("results/stage9/sim_design3_power.tsv")
    fig = new_fig(220)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1, 1], hspace=0.48, wspace=0.28, left=0.10, right=0.98, top=0.95, bottom=0.14)
    panels = {
        "A": fig.add_subplot(gs[0, :]),
        "B": fig.add_subplot(gs[1, 0]),
        "C": fig.add_subplot(gs[1, 1]),
        "D": fig.add_subplot(gs[2, 0]),
        "E": fig.add_subplot(gs[2, 1]),
    }
    arows, brows, crows, drows, erows = [], [], [], [], []
    for clock in clocks:
        for scenario, ls in SCENARIO_LS.items():
            block = bias.loc[(bias["clock"] == clock) & (bias["scenario"] == scenario)].copy()
            block["_s"] = block["s"].map(num)
            block = block.sort_values("_s")
            panels["A"].plot(block["_s"], block["b_over_total_sd"].map(num), color=CLOCK_COLOR[clock], ls=ls, lw=SCENARIO_LW[scenario])
            for _, row in block.iterrows():
                arows.append({"clock": clock, "scenario": scenario, "s": row.s, "b": row.b, "b_over_total_sd": row.b_over_total_sd, "unit": row.unit})
            if prefix == "Fig7" and clock == "RA_default" and scenario in ("lab", "free-living"):
                hit = one(block, s="3")
                label = "1.1 years" if scenario == "lab" else "0.6 years"
                panels["A"].annotate(label, (3, num(hit.b_over_total_sd)), textcoords="offset points", xytext=(6, 4 if scenario == "lab" else -8), fontsize=7, color=CLOCK_COLOR[clock])
                arows.append({"clock": clock, "scenario": scenario, "s": "3", "b": hit.b, "annotation": label, "record": "annotation"})
            sub = fpr1.loc[(fpr1["clock"] == clock) & (fpr1["scenario"] == scenario) & (fpr1["n"] == "100")].copy()
            sub["_s"] = sub["s"].map(num)
            sub = sub.sort_values("_s")
            panels["B"].plot(sub["_s"], sub["fpr"].map(num), color=CLOCK_COLOR[clock], ls=ls, lw=SCENARIO_LW[scenario])
            for _, row in sub.iterrows():
                brows.append({"clock": clock, "scenario": scenario, "s": row.s, "n": row.n, "fpr": row.fpr})
            sub = fpr1.loc[(fpr1["clock"] == clock) & (fpr1["scenario"] == scenario) & (fpr1["s"] == "3")].copy()
            sub["_n"] = sub["n"].map(num)
            sub = sub.sort_values("_n")
            panels["C"].plot(sub["_n"], sub["fpr"].map(num), color=CLOCK_COLOR[clock], ls=ls, lw=SCENARIO_LW[scenario])
            for _, row in sub.iterrows():
                crows.append({"clock": clock, "scenario": scenario, "s": row.s, "n": row.n, "fpr": row.fpr})
            sub = fpr2.loc[(fpr2["clock"] == clock) & (fpr2["scenario"] == scenario) & (fpr2["s"] == "4")].copy()
            sub["_n"] = sub["n_per_group"].map(num)
            sub = sub.sort_values("_n")
            panels["D"].plot(sub["_n"], sub["fpr"].map(num), color=CLOCK_COLOR[clock], ls=ls, lw=SCENARIO_LW[scenario])
            for _, row in sub.iterrows():
                drows.append({"clock": clock, "scenario": scenario, "s": row.s, "n_per_group": row.n_per_group, "fpr": row.fpr})
            if scenario == "lab":
                sub = d3.loc[(d3["clock"] == clock) & (d3["scenario"] == "lab")].copy()
                sub["_w"] = sub["W"].map(num)
                sub = sub.sort_values("_w")
                xs = sub["_w"].to_numpy()
                ys = sub["power"].map(num).to_numpy()
                lo, hi = [], []
                for _, row in sub.iterrows():
                    a, b, z = wilson(num(row.n_reject), num(row.n_rep))
                    lo.append(a)
                    hi.append(b)
                    erows.append({"clock": clock, "W": row.W, "n_reject": row.n_reject, "n_rep": row.n_rep, "n_failed": row.n_failed, "power": row.power, "wilson_low": a, "wilson_high": b, "z_0.975": z})
                panels["E"].errorbar(xs, ys, yerr=[ys - np.array(lo), np.array(hi) - ys], color=CLOCK_COLOR[clock], lw=1.0, marker="o", ms=3.5, capsize=2)
    panels["A"].set_xlabel("s (h)")
    panels["A"].set_ylabel("b(s) / total SD")
    panels["A"].set_xlim(0, 12)
    panels["B"].axhline(0.05, color="0.35", lw=0.5, ls=":")
    panels["B"].axhline(0.10, color="0.35", lw=0.5, ls="--")
    panels["B"].set_xlabel("s (h)")
    panels["B"].set_ylabel("False-positive rate")
    panels["B"].set_ylim(bottom=0)
    panels["C"].set_xlabel("n")
    panels["C"].set_ylabel("False-positive rate")
    panels["C"].set_ylim(bottom=0)
    panels["D"].set_xlabel("n per group")
    panels["D"].set_ylabel("False-positive rate")
    panels["D"].set_ylim(bottom=0)
    panels["E"].set_xlabel("W (h)")
    panels["E"].set_ylabel("Power")
    panels["E"].set_ylim(0, 1)
    for letter, ax in panels.items():
        panel(ax, letter)
    handles = [Line2D([0], [0], color=CLOCK_COLOR[c], lw=1.2, label=CLOCK_NAME[c]) for c in clocks]
    handles += [
        Line2D([0], [0], color="black", lw=1.0, ls="-", label="lab"),
        Line2D([0], [0], color="black", lw=1.0, ls="--", label="free-living"),
        Line2D([0], [0], color="black", lw=0.6, ls=":", label="free-living upper"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, title="post hoc", bbox_to_anchor=(0.5, 0.005))
    write_src(pd.DataFrame(arows), f"{prefix}_A.tsv")
    write_src(pd.DataFrame(brows), f"{prefix}_B.tsv")
    write_src(pd.DataFrame(crows), f"{prefix}_C.tsv")
    write_src(pd.DataFrame(drows), f"{prefix}_D.tsv")
    write_src(pd.DataFrame(erows), f"{prefix}_E.tsv")
    save_fig(fig, stem)


def fig8(clocks, stem, prefix):
    cos = load("results/stage9/composition_cosinor.tsv")
    l3 = load("results/stage9/composition_l3.tsv")
    fig = new_fig(150)
    gs = fig.add_gridspec(2, len(clocks), height_ratios=[1.2, 1], hspace=0.45, wspace=0.35, left=0.08, right=0.98, top=0.90, bottom=0.12)
    arows = []
    datasets = ["GSE39445", "GSE56931", "GSE113883", "L1"]
    for i, clock in enumerate(clocks):
        ax = fig.add_subplot(gs[0, i])
        if i == 0:
            panel(ax, "A")
        for dataset in datasets:
            row = one(cos, dataset=dataset, clock=clock)
            geo = geo_name(dataset)
            open_mark = num(row.lrt_p_before) >= 0.05
            ax.plot([0, 1], [num(row.two_a_before), num(row.two_a_after)], color=DATA_COLOR[geo], lw=0.8, marker=DATA_MARKER[geo], ms=4, markerfacecolor=("white" if open_mark else DATA_COLOR[geo]), markeredgecolor=DATA_COLOR[geo])
            ax.text(1.08, num(row.two_a_after), f"{num(row.two_a_ratio):.2f}", fontsize=7, va="center", color=DATA_COLOR[geo])
            arows.append({"clock": clock, "dataset": geo, "two_a_before": row.two_a_before, "two_a_after": row.two_a_after, "two_a_ratio": row.two_a_ratio, "lrt_p_before": row.lrt_p_before, "open_marker": open_mark})
        ax.set_xlim(-0.25, 1.45)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["before", "after"])
        ax.set_ylim(bottom=0)
        ax.set_title(CLOCK_NAME[clock], fontsize=8)
        if i == 0:
            ax.set_ylabel("2A")
    brows = []
    axb_parent = fig.add_subplot(gs[1, :])
    axb_parent.set_visible(False)
    fig.text(0.04, 0.42, "B", fontsize=9, fontweight="bold")
    inner = gs[1, :].subgridspec(1, len(clocks), wspace=0.45)
    for i, clock in enumerate(clocks):
        ax = fig.add_subplot(inner[0, i])
        row = one(l3, clock=clock)
        labels, xs, estimates = [], [], []
        series = [("before", row.mean_delta_before, row.ci_low_before, row.ci_high_before), ("composition", row.mean_delta_composition, row.ci_low_composition, row.ci_high_composition)]
        if row.mean_delta_time != "":
            series.append(("time", row.mean_delta_time, row.ci_low_time, row.ci_high_time))
        for j, (lab, est, lo, hi) in enumerate(series):
            ax.errorbar(j, num(est), yerr=[[num(est) - num(lo)], [num(hi) - num(est)]], color=CLOCK_COLOR[clock], marker="o", ms=4, capsize=2, lw=0.8)
            labels.append(lab)
            brows.append({"clock": clock, "adjustment": lab, "estimate": est, "low": lo, "high": hi, "time_source": row.time_source})
        ax.axhline(0, color="0.4", lw=0.5)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_title(CLOCK_NAME[clock], fontsize=8)
        if i == 0:
            ax.set_ylabel("Mean Δ")
    handles = [Line2D([0], [0], marker=DATA_MARKER[geo_name(d)], color=DATA_COLOR[geo_name(d)], lw=0.8, label=geo_name(d)) for d in datasets]
    handles.append(Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor="0.2", lw=0, label="LRT p ≥ 0.05"))
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, title="post hoc", bbox_to_anchor=(0.5, 0.0))
    write_src(pd.DataFrame(arows), f"{prefix}_A.tsv")
    write_src(pd.DataFrame(brows), f"{prefix}_B.tsv")
    save_fig(fig, stem)


def merge_review():
    stems = [
        "Fig1_design", "Fig2_trajectories", "Fig3_effects", "Fig4_composition", "Fig5_modules",
        "Fig6_public", "Fig7_bias", "Fig8_mitigation", "FigS1_effects_secondary", "FigS2_bias_secondary", "FigS3_mitigation_secondary",
    ]
    cmd = ["gs", "-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", "-sOutputFile=" + str(OUT / "all_figures_review.pdf")]
    cmd += [str(OUT / f"{stem}.pdf") for stem in stems]
    subprocess.check_call(cmd)


def main():
    fig1()
    fig2()
    fig3(PRIMARY, "3", ["A", "B", "C"])
    fig3(SECONDARY, "S1", ["A", "B"])
    fig4()
    fig5()
    fig6()
    fig7(PRIMARY, "Fig7_bias", "Fig7")
    fig7(SECONDARY, "FigS2_bias_secondary", "FigS2")
    fig8(PRIMARY, "Fig8_mitigation", "Fig8")
    fig8(SECONDARY, "FigS3_mitigation_secondary", "FigS3")
    merge_review()
    print("figures written")


if __name__ == "__main__":
    main()
