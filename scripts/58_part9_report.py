#!/usr/bin/env python3
"""Write reports/TOD_보고서_9.md from the stage-9 tables."""

from pathlib import Path

import pandas as pd

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
CLOCKS = ["TAms_chrono", "TAms_mort", "RA_default", "RA_Peters", "RA_GTEx"]


def load(path: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / path, sep="\t", dtype=str).fillna("")


def md(frame: pd.DataFrame) -> str:
    bad = frame.astype(str).apply(lambda col: col.str.contains("|", regex=False)).any()
    if bool(bad.any()):
        raise SystemExit("a table cell contains a pipe")
    cols = [c for c in frame.columns if frame[c].astype(str).str.strip().ne("").any()]
    use = frame[cols]
    header = "| " + " | ".join(cols) + " |"
    rule = "| " + " | ".join(["---"] * len(cols)) + " |"
    rows = ["| " + " | ".join(row.astype(str)) + " |" for _, row in use.iterrows()]
    return "\n".join([header, rule, *rows])


def one(frame: pd.DataFrame, **eq) -> pd.Series:
    matched = frame
    for key, value in eq.items():
        matched = matched.loc[matched[key] == str(value)]
    if len(matched) != 1:
        raise SystemExit(f"{len(matched)} rows for {eq}")
    return matched.iloc[0]


def join_clocks(frame: pd.DataFrame, text) -> str:
    return "; ".join(text(one(frame, clock=clock)) for clock in CLOCKS)


def main() -> None:
    gate = load("results/stage9/l2_wake_corr_gate.tsv")
    genes = load("results/stage9/l2_wake_corr_genes.tsv")
    circ = load("results/stage9/l2_wake_corr_gene_circular.tsv")
    clocks = load("results/stage9/l2_wake_corr_clocks.tsv")
    params = load("results/stage9/sim_parameters.tsv")
    bias = load("results/stage9/sim_design1_bias.tsv")
    fpr1 = load("results/stage9/sim_design1_fpr.tsv")
    smallest = load("results/stage9/sim_design1_smallest_s.tsv")
    mc = load("results/stage9/sim_design1_montecarlo.tsv")
    fpr2 = load("results/stage9/sim_design2_fpr.tsv")
    d3 = load("results/stage9/sim_design3_power.tsv")
    cos = load("results/stage9/composition_cosinor.tsv")
    l3 = load("results/stage9/composition_l3.tsv")
    betas = load("results/stage9/composition_l3_betas.tsv")
    age = load("results/stage9/age_correlation_l3.tsv")
    summary = load("results/stage9/mitigation_summary.tsv")
    docx = load("results/stage9/docx_check.tsv")
    sources = load("figures/figure_sources.tsv")
    sources = sources.loc[sources["figure"].isin(["Fig6", "Fig7"])]

    gate_bits = []
    for gene in ["PER1", "PER3"]:
        row = one(gate, gene=gene)
        gate_bits.append(
            f"{gene} peak_wake_corr `{row.peak_wake_corr}`, in_prereg_window `{row.in_prereg_window}`, "
            f"lrt_p_relclock `{row.lrt_p_relclock}`, lrt_p_wake_corr `{row.lrt_p_wake_corr}`, "
            f"lrt_p_difference `{row.lrt_p_difference}`, prereg_gate_pass `{row.prereg_gate_pass}`"
        )
    circ_row = circ.iloc[0]
    lines = [
        "- " + ". ".join(gate_bits) + ".",
        f"- 유전자 10개 circ_mean_h `{circ_row.circ_mean_h}`, circ_sd_h `{circ_row.circ_sd_h}`.",
    ]
    for scenario in ["lab", "free-living"]:
        bits = []
        for clock in CLOCKS:
            b = one(bias, scenario=scenario, clock=clock, s="3")
            f = one(fpr1, scenario=scenario, clock=clock, s="3", n="100")
            bits.append(f"{clock} b `{b.b}`, b/total_sd `{b.b_over_total_sd}`, fpr `{f.fpr}`")
        lines.append(f"- 설계 1, s = 3, {scenario}, n = 100: " + "; ".join(bits) + ".")
    small_bits = []
    for scenario in ["lab", "free-living", "free-living upper"]:
        for clock in CLOCKS:
            row = one(smallest, scenario=scenario, clock=clock, n="100")
            small_bits.append(f"{scenario} {clock} `{row.smallest_s}` ({row.status})")
    lines.append("- n = 100에서 위양성률이 threshold `0.1`을 넘는 가장 작은 s: " + "; ".join(small_bits) + ".")
    d2_bits = []
    for scenario in ["lab", "free-living", "free-living upper"]:
        for clock in CLOCKS:
            row = one(fpr2, scenario=scenario, clock=clock, s="4", n_per_group="100")
            d2_bits.append(f"{scenario} {clock} `{row.fpr}`")
    lines.append("- 설계 2, s = 4, 군당 n = 100: " + "; ".join(d2_bits) + ".")
    for width in ["0", "2", "4", "8"]:
        bits = []
        for clock in CLOCKS:
            row = one(d3, clock=clock, W=width)
            bits.append(f"{clock} power `{row.power}`, bias `{row.bias}`")
        lines.append(f"- 설계 3, W = {width}: " + "; ".join(bits) + ".")
    l1 = cos.loc[cos["dataset"] == "L1"]
    lines.append("- L1 two_a_ratio: " + join_clocks(l1, lambda row: f"{row.clock} `{row.two_a_ratio}`") + ".")
    lines.append(
        "- L3 mean Δ: "
        + join_clocks(
            l3,
            lambda row: (
                f"{row.clock} before `{row.mean_delta_before}`, composition `{row.mean_delta_composition}`, "
                f"time `{row.mean_delta_time}`"
            ),
        )
        + "."
    )
    lines.append(
        "- L3 나이 Pearson r: "
        + join_clocks(
            age,
            lambda row: (
                f"{row.clock} morning `{row.morning_r_before}` / `{row.morning_r_after}`, "
                f"evening `{row.evening_r_before}` / `{row.evening_r_after}`, "
                f"person mean `{row.person_mean_r}`"
            ),
        )
        + "."
    )
    if len(lines) > 14:
        raise SystemExit(f"section 1 has {len(lines)} lines")

    text = f"""# TOD 보고서 9 — L2 사후 관문, 편향 시뮬레이션, 완화

## 1. 요약 (수치만, 14줄 이내)

{chr(10).join(lines)}

## 2. 환경과 자원

- 분석: `scripts/56_part9_analysis.R`. R version 4.5.3 (2026-03-11), conda 환경 `r_env`.
- 명령: `/usr/bin/time -v Rscript scripts/56_part9_analysis.R`.
- 사용자 시간 `108.79` s, 시스템 시간 `2.18` s, CPU `102%`, 경과 `1:48.72`, 최대 RSS `2106748` kB, 종료 상태 0.
- 결정 10의 중지 기준은 최대 RSS 40 GB다. 관측 최대 RSS는 `2106748` kB다.
- 몬테카를로 시드 `20261019`. 설계 1 몬테카를로 시작에 `set.seed(20261019)`, 설계 3 시작에 다시 `set.seed(20261019)`.
- 로그 `logs/stage9_analysis.log`의 Part A 구간에 `boundary (singular) fit` 5줄, `Warning: stack imbalance` 4줄이 있다. 설계 1·2 출력 뒤에 `There were 25 warnings`가 있다.
- 그림: `/home/jkang/anaconda3/bin/python3`, `MPLBACKEND=Agg`, `scripts/57_part9_figures.py`. `savefig` dpi 300. PNG info dpi는 Fig 6과 Fig 7 모두 `(299.9994, 299.9994)`.
- 시계 표기: `TAms_chrono` tAge chronological, `TAms_mort` tAge mortality, `RA_default` RNAAgeCalc blood, `RA_Peters`, `RA_GTEx`.
- 표와 그림 제목은 `post hoc`다. `results/confirm/`은 수정하지 않았다.

## 3. L2 사후 관문 (표)

`wake_corr = (relclocktime - 8) mod 24`. period 1. LRT p는 `relclocktime` 축과 `wake_corr` 축에서 각각 계산했다.

사후 관문 (`results/stage9/l2_wake_corr_gate.tsv`):

{md(gate)}

유전자 10개 (`results/stage9/l2_wake_corr_genes.tsv`):

{md(genes)}

원형 요약 (`results/stage9/l2_wake_corr_gene_circular.tsv`):

{md(circ)}

시계 정점 (`results/stage9/l2_wake_corr_clocks.tsv`):

{md(clocks)}

## 4. 편향 시뮬레이션

### 4.1 매개변수 표

`results/stage9/sim_parameters.tsv`. free-living의 코사인 모양과 φ̄는 `assumption` 열에 있다.

{md(params)}

### 4.2 설계 1 (b(s) 표, 위양성률 표, 최소 s 표, 몬테카를로 비교)

대응 t, 양측 α = 0.05. ncp = b / SD_diff × sqrt(n), df = n − 1.

b(s) (`results/stage9/sim_design1_bias.tsv`):

{md(bias)}

위양성률 (`results/stage9/sim_design1_fpr.tsv`):

{md(fpr1)}

n = 100에서 위양성률이 0.10을 넘는 가장 작은 s (`results/stage9/sim_design1_smallest_s.tsv`):

{md(smallest)}

몬테카를로, lab, n = 50, s = 3, 2000회, 시드 20261019 (`results/stage9/sim_design1_montecarlo.tsv`):

{md(mc)}

### 4.3 설계 2

두 표본 t, 같은 분산, 양측 α = 0.05. SD는 총 SD. ncp = b / total_sd × sqrt(n/2), df = 2n − 2. `results/stage9/sim_design2_fpr.tsv`.

{md(fpr2)}

### 4.4 설계 3

lab. 사전·사후 시각은 각각 독립 Uniform(8, 8+W). δ = 0.3 × sqrt(2) × σ를 사후 표본에 더한다. n = 50. 2000회, 시드 20261019. `power` = n_reject / 2000. `bias` = mean(estimates) − δ. `results/stage9/sim_design3_power.tsv`.

{md(d3)}

## 5. 완화

### 5.1 코사이너 자료 구성 보정 표

MCPcounter 8개 점수를 참가자 안에서 중심화해 모형 M의 고정 효과로 넣었다. 보정 모형의 LRT는 cos와 sin을 빼고 구성 공변량은 남긴다. L1은 조성이 비유한인 샘플을 두 모형에서 함께 뺐다. `results/stage9/composition_cosinor.tsv`.

{md(cos)}

### 5.2 L3 구성 보정 표

Δclock을 Δcomposition 8개와 절편에 회귀했다. 절편이 구성 보정 mean Δ다. 시각 열은 H5의 mean Δ − c다. `results/stage9/composition_l3.tsv`.

{md(l3)}

3.2 회귀 계수 (`results/stage9/composition_l3_betas.tsv`):

{md(betas)}

### 5.3 나이 연관 표

아침 보정 = 점수 − β̂ᵀ(Δcomposition / 2), 저녁 보정 = 점수 + β̂ᵀ(Δcomposition / 2). `results/stage9/age_correlation_l3.tsv`.

{md(age)}

### 5.4 요약 표

코사이너 행의 효과는 2A다. L3 행의 효과는 mean Δ다. `results/stage9/mitigation_summary.tsv`.

{md(summary)}

## 6. 그림 (파일 목록, 원자료와 n), docx 확인 결과

- `figures/Fig6_simulation.pdf`, `figures/Fig6_simulation.png`. PNG 3600×4200, info dpi `(299.9994, 299.9994)`.
- `figures/Fig7_mitigation.pdf`, `figures/Fig7_mitigation.png`. PNG 3600×2400, info dpi `(299.9994, 299.9994)`.

`figures/figure_sources.tsv`의 Fig 6·Fig 7 행:

{md(sources)}

지시서 8 docx 확인. `/usr/bin/soffice` (LibreOffice 7.3) `--headless --norestore --convert-to pdf`. 변환 종료 상태 0. 경고: `failed to launch javaldx - java may not function correctly`. `results/stage9/docx_check.tsv`.

{md(docx)}

## 7. 판단이 필요한 사항

- 설계 2의 ncp 식은 지시서에 없다. 지시한 같은 분산 두 표본 t에 맞춰 ncp = b / total_sd × sqrt(n/2), df = 2n − 2로 계산했다.
- `prereg_gate_pass`는 `in_prereg_window`와 `lrt_p_below_0.01`을 함께 요구한다. 창 판정 열은 `in_prereg_window`다.

## 8. 지시서와 다르게 한 것

없음.

## 9. 산출물 목록, 커밋 해시

- Part A–C `50e5f00580c580d77fa67215b9dd097528ca8268` (2026-10-19 16:33:32 +0900). `scripts/56_part9_analysis.R`. `results/stage9/`의 L2 표 4개, 시뮬레이션 표 7개, 완화 표 5개.
- 그림 `c31ddf74158a59a282e8f3546780b50081d74d55` (2026-10-19 16:33:39 +0900). `scripts/57_part9_figures.py`, `figures/Fig6_simulation.pdf`, `figures/Fig6_simulation.png`, `figures/Fig7_mitigation.pdf`, `figures/Fig7_mitigation.png`, `figures/figure_sources.tsv`, `results/stage9/docx_check.tsv`.
- 보고서: 이 보고서를 넣는 커밋. `scripts/58_part9_report.py`, `reports/TOD_보고서_9.md`.
"""
    out = ROOT / "reports" / "TOD_보고서_9.md"
    out.write_text(text, encoding="utf-8")
    print(f"section 1 lines: {len(lines)}")
    print(f"wrote {out} ({out.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
