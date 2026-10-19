#!/usr/bin/env python3
"""Write reports/TOD_보고서_8.md from the stage-8 and confirmation tables."""

from pathlib import Path

import pandas as pd

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")


def load(path: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / path, sep="\t", dtype=str).fillna("")


def md(frame: pd.DataFrame) -> str:
    cols = [c for c in frame.columns if frame[c].astype(str).str.strip().ne("").any()]
    use = frame[cols]
    header = "| " + " | ".join(cols) + " |"
    rule = "| " + " | ".join(["---"] * len(cols)) + " |"
    rows = ["| " + " | ".join(row.astype(str)) + " |" for _, row in use.iterrows()]
    return "\n".join([header, rule, *rows])


def main() -> None:
    boot = load("results/stage8/cosinor_effects_boot.tsv")
    f1 = load("prereg/frozen/F1.tsv")
    joined = boot.merge(f1, on=["dataset", "clock"], suffixes=("_s8", "_f1"))
    rel = (joined["two_a_s8"].astype(float) - joined["two_a_f1"].astype(float)).abs() / joined["two_a_f1"].astype(float).abs()
    max_rel = format(float(rel.max()), ".16e")
    effects = load("results/stage8/effect_summary.tsv")
    cos = effects.loc[effects["kind"] == "cosinor", [
        "dataset", "clock", "role", "n_subjects", "n_samples",
        "estimate", "low", "high", "r2_within", "r2_low", "r2_high",
        "peak_h", "peak_center", "peak_halfwidth_95", "p",
    ]].rename(columns={
        "estimate": "two_a_over_sigma", "low": "ratio_p025", "high": "ratio_p975",
    })
    paired = effects.loc[effects["kind"] == "paired", [
        "dataset", "clock", "role", "n_subjects", "n_samples",
        "estimate", "low", "high", "mean_delta", "mean_low", "mean_high", "p",
    ]].rename(columns={"estimate": "d", "low": "d_p025", "high": "d_p975"})
    pred = load("results/stage8/l3_d_vs_prediction.tsv")
    d11 = load("results/stage8/l1_f3_d11_model_m.tsv")
    feat = load("results/stage8/f3_feature_counts.tsv")
    audit = load("results/stage8/d11_matrix_audit_recount.tsv")
    circ = load("results/stage8/gene_peak_circular_summary_arntl_alias.tsv")
    bmal = load("results/stage8/l2_bmal1_as_arntl.tsv")
    quotes = load("results/stage8/l2_sleep_start_search.tsv")
    sources = load("figures/figure_sources.tsv")
    table1 = load("results/stage8/tables/Table1.tsv")
    fam_p = load("results/confirm/family_P.tsv")
    fam_s = load("results/confirm/family_S1.tsv")
    rep = load("results/confirm/replication.tsv")
    h1 = load("results/confirm/h1_secondary.tsv")
    h2 = load("results/confirm/h2.tsv")
    h3 = load("results/confirm/h3.tsv")
    h4 = load("results/confirm/h4.tsv")
    h5 = load("results/confirm/h5.tsv")
    s4 = load("results/confirm/s4.tsv")
    s4 = s4.loc[s4["in_f3"] == "TRUE", [
        "dataset", "clock", "n_subjects", "n_samples", "two_a", "peak_h", "p",
        "confirmed", "circ_dist", "phi_bar", "mean_delta", "observed_sign",
        "predicted_sign", "sign_agrees",
    ]]
    s4count = load("results/confirm/s4_l1_confirmed.tsv")
    gates = []
    for dataset in ["L1", "L2", "L3"]:
        g = load(f"results/confirm/{dataset}/gate.tsv")
        g.insert(0, "dataset", dataset)
        gates.append(g)
    p_cols = [
        "dataset", "clock", "n_subjects", "n_samples", "two_a_over_sigma",
        "ratio_p025", "ratio_p975", "r2_within", "r2_p025", "r2_p975",
        "peak_h", "peak_center", "peak_halfwidth_95", "p", "q_holm", "decision",
        "d", "mean_delta", "ci_low", "ci_high",
    ]
    text = f"""# TOD 보고서 8 — 확인 결과 확정, 사후 분석, 원고 그림·표

## 1. 요약 (수치만, 12줄 이내)

- Part A `67d4c2d7b3570f2b72aba06242706057cc4e635d` (2026-10-19 11:09:06 +0900). 사전등록 확인 결과는 실행 5 `e56a1c489a74453ca33a832835fdaa312e3e8371`의 `results/confirm/`. 이탈은 D1–D12. 8.A2는 L2 기호 BMAL1을 ARNTL 별칭으로 기록. 8.A3는 원형 SD 식을 기록.
- Part B `9317da205567e399907d98a2ae1ef21e50f4aab8` (2026-10-19 12:03:27 +0900). `scripts/51_part8_posthoc.R`의 `/usr/bin/time -v`: 경과 `46:00.17`, 최대 RSS `1560928` kB, 종료 상태 0.
- 습관적 수면 시작 시각의 평균과 범위는 확인 못 함. `time_transform`은 `not performed`.
- BMAL1 별칭: n 70/11, peak_h `21.3692664519028`, lrt_p `0.00925089859203897`. 유전자 10개의 circ_mean_h `9.69765168583892`, circ_sd_h `5.49841393915787`.
- L3 검정력, 쌍 58, 양측 α=0.05: `TAms_chrono` `0.176777740279187`, `TAms_mort` `0.291920013796665`, `RA_default` `0.686538271492044`, `RA_Peters` `0.511452502493332`, `RA_GTEx` `0.108269826724681`.
- D11 모형 M `posthoc_confirmed` TRUE 6/8. 동결 S4 L1 확인 수는 `1`/`8`.
- 그림 PDF·PNG 300 dpi: `figures/Fig2_trajectories`, `figures/Fig3A_effect_forest`, `figures/Fig3B_peak_clock`, `figures/Fig4_composition`, `figures/Fig5_modules`. 표 TSV·docx: `results/stage8/tables/Table1`, `Table2`, `Table3`, `TableS1`, `TableS2`, `TableS3`.

## 2. 확정 기록

`config/decisions.tsv`의 2026-10-19 행이다.

| id | content |
|---|---|
| 8.A1 | The preregistered confirmation result is run 5, results/confirm/, commit e56a1c489a74453ca33a832835fdaa312e3e8371. The deviations that apply are D1 through D12. |
| 8.A2 | In the L2 post-hoc clock-gene table, the L2 matrix symbol BMAL1 is used as the alias of ARNTL and the circular summary is recomputed. The recomputed table is post hoc only. results/confirm/ is not changed. |
| 8.A3 | The circular SD formula in the stage-7c post-hoc table is approved: sqrt(-2*log(R))*24/(2*pi), where R is the mean resultant length of the signed differences. |

BMAL1 적합 (`results/stage8/l2_bmal1_as_arntl.tsv`):

{md(bmal)}

원형 요약 (`results/stage8/gene_peak_circular_summary_arntl_alias.tsv`):

{md(circ)}

## 3. 사후 분석

### 3.1 L2 시각 기준 보정 (인용, 변환, 사후 관문, 정점 거리)

출처는 Kervezee et al., J Biol Rhythms, PMID 30621487, DOI 10.1177/0748730418821776. 습관적 수면 시작 시각의 평균과 범위는 확인 못 함. `time_transform`은 `not performed`. 사후 관문과 보정 축의 정점 거리는 계산하지 않았다. Fig 3B에는 L2를 넣지 않았다.

원문:

- Briefly, participants maintained a regular, selfselected sleep/wake schedule with a nocturnal sleep period of 8 h for at least one week before admission to the laboratory, which was verified by daily sleep logs and actigraphy.
- Upon admission to the laboratory, participants were involved in a 24-h baseline measurement period that consisted of an 8-h sleep period during their habitual sleep times and a 16-h constant posture procedure, involving a constant semi-recumbent position, minimal activity, and hourly isocaloric snacks.
- relative clock time refers to the time of day relative to the start of each participant's habitual sleep period.

검색 표:

{md(quotes)}

### 3.2 자료 전체 효과 요약 표

개발 자료 3개와 L2의 코사이너 구간은 `scripts/51_part8_posthoc.R`가 확인 코드 `boot_effects`로 계산한 값이다. 2A/σ는 `two_a/sigma_reml`이다. 부트스트랩은 참가자 재추출 2000회, seed `20261015`, 백분위 type 7이다. 개발 15개 시계의 `two_a`와 `prereg/frozen/F1.tsv`의 `two_a` 상대 차이의 최댓값은 `{max_rel}`이다. 네 자료 모두 `n_converged`는 2000이다. L1과 L3의 점추정과 L1 구간, L3의 mean Δ 구간과 p는 `results/confirm/family_P.tsv`와 `family_S1.tsv`의 문자열이다. L3의 d 구간은 `results/stage8/l3_d_vs_prediction.tsv`이다. 표를 합치지 않았다.

코사이너. `two_a_over_sigma`의 구간은 `ratio_p025`, `ratio_p975`. `r2_within`의 구간은 `r2_p025`, `r2_p975`. 정점 구간은 `peak_center`와 `peak_halfwidth_95`이다.

{md(cos)}

L3. d의 구간은 `d_p025`, `d_p975`. mean Δ의 구간은 `mean_low`, `mean_high`.

{md(paired)}

### 3.3 L3 관측 효과와 예측

관측 d는 참가자 58명의 morning−evening에 대한 부트스트랩 2000회, seed `20261015`, quantile type 7이다. 검정력은 paired t, df = n−1, ncp = d*sqrt(n), 양측 α=0.05이다. `d_exp_median`과 `d_GSE113883`은 `prereg/frozen/F2.tsv`이다.

{md(pred)}

`RA_default`의 d는 이 표에서 `0.326704860352488`이고 `results/confirm/family_P.tsv`에서 `0.326704860352489`이다.

### 3.4 L1 모듈 결과의 D11 판

D11 점수는 `results/confirm/L1/tage_scores_d11.tsv`와 `tage_components_d11.tsv`의 period D1, 참가자 22, 샘플 147이다. 기준은 LRT p < 0.05이고 φ̄까지 원형 거리가 3시간 이내이다.

{md(d11)}

`posthoc_confirmed` TRUE는 6행, FALSE는 2행이다.

특징 수 (`results/stage8/f3_feature_counts.tsv`). 유전자 이름은 적지 않았다.

{md(feat)}

D11 행렬을 같은 규칙으로 다시 센 감사 (`results/stage8/d11_matrix_audit_recount.tsv`):

{md(audit)}

## 4. 그림 (파일 목록, 그림별 원자료와 n)

PDF와 PNG는 300 dpi이다. 제목에 post hoc, exploratory가 있다. 개발은 녹색 원, 확인은 주황 사각형, 사후는 보라 다이아몬드이다. 주 시계 이름은 tAge chronological, tAge mortality, RNAAgeCalc blood이다.

| 파일 | 내용 |
|---|---|
| figures/Fig2_trajectories.pdf, figures/Fig2_trajectories.png | GSE113883과 L1 Day 1의 참가자 궤적과 집단 코사이너. x는 벽시계 시각. |
| figures/Fig3A_effect_forest.pdf, figures/Fig3A_effect_forest.png | §3.2. 코사이너 자료는 2A/σ, L3는 d. |
| figures/Fig3B_peak_clock.pdf, figures/Fig3B_peak_clock.png | 주 시계 3개의 정점과 개발 φ̄. L2는 없다. |
| figures/Fig4_composition.pdf, figures/Fig4_composition.png | p_med와 MCPcounter 8개 집단의 정점 시각. |
| figures/Fig5_modules.pdf, figures/Fig5_modules.png | F3 8개 결과의 개발 φ̄. 바깥 원 L1 동결 확인, 가운데 사각형 L1 D11 사후 확인, 안쪽 다이아몬드 L3 부호 일치. 회색은 해당 상태가 아님. |
| figures/figure_sources.tsv | 그림별 원자료와 n. |

{md(sources)}

Fig 2의 집단 곡선은 저장된 `beta_cos`, `beta_sin`으로 계산한 파동을 점수에서 뺀 값의 평균을 세로 위치로 둔 곡선이다.

## 5. 표 (파일 목록, Table 1–3 본문 수록)

| 파일 | 내용 |
|---|---|
| results/stage8/tables/Table1.tsv, Table1.docx | 자료 8개 |
| results/stage8/tables/Table2.tsv, Table2.docx | 관문, 가족 P, 재현, H1, S1 |
| results/stage8/tables/Table3.tsv, Table3.docx | H2, H3, H4, H5, S4 |
| results/stage8/tables/TableS1.tsv, TableS1.docx | D1–D12 |
| results/stage8/tables/TableS2.tsv, TableS2.docx | 사전등록 §8과 D1, D2, D6, D9, D11 민감도 |
| results/stage8/tables/TableS3.tsv, TableS3.docx | §3.1–§3.4와 L2 H1·H4 탐색 |

### Table 1

{md(table1.drop(columns=["title"]))}

### Table 2

관문:

{md(pd.concat(gates, ignore_index=True))}

가족 P:

{md(fam_p[p_cols])}

재현:

{md(rep)}

H1:

{md(h1)}

가족 S1:

{md(fam_s[p_cols])}

### Table 3

H2:

{md(h2)}

H3:

{md(h3)}

H4:

{md(h4)}

H5:

{md(h5)}

S4, `in_f3` TRUE:

{md(s4)}

S4 L1 확인 수:

{md(s4count)}

## 6. 판단이 필요한 사항

없음.

## 7. 지시서와 다르게 한 것

- docx는 python-docx, officer, flextable이 설치되어 있지 않아 OOXML zip으로 썼다. 패키지는 설치하지 않았다.
- Fig 4의 MCPcounter 정점은 시각 0–24 가로축에 집단 이름을 세로축으로 두었다.
- L2 `RA_default` LRT p는 이번 적합에서 `0.0116625636160948`이고 `results/posthoc_stage7c/h1_exploratory.tsv`에서 `0.0116625636160945`이다.
- `scripts/52_part8_feature_counts.py`의 `/usr/bin/time -v` 기록은 없다.

## 8. 산출물 목록, 커밋 해시

- Part A: `67d4c2d7b3570f2b72aba06242706057cc4e635d`, 2026-10-19 11:09:06 +0900, `TOD8 Part A`.
- Part B: `9317da205567e399907d98a2ae1ef21e50f4aab8`, 2026-10-19 12:03:27 +0900, `TOD8 Part B`.
- 이 보고서를 넣는 커밋: 보고서 작성 후 기록.
- 사후 표: `results/stage8/l2_sleep_start_search.tsv`, `l2_bmal1_as_arntl.tsv`, `gene_peak_differences_arntl_alias.tsv`, `gene_peak_circular_summary_arntl_alias.tsv`, `cosinor_effects_boot.tsv`, `effect_summary.tsv`, `l3_d_vs_prediction.tsv`, `l1_f3_d11_model_m.tsv`, `f3_feature_counts.tsv`, `d11_matrix_audit_recount.tsv`.
- 스크립트: `scripts/51_part8_posthoc.R`, `scripts/52_part8_feature_counts.py`, `scripts/53_part8_figures.py`, `scripts/54_part8_tables.py`, `scripts/55_part8_report.py`.
- `results/confirm/`은 변경하지 않았다.
"""
    out = ROOT / "reports" / "TOD_보고서_8.md"
    out.write_text(text, encoding="utf-8")
    section = text.split("## 2.")[0]
    print("section1 lines", section.count("\n"))


if __name__ == "__main__":
    main()
