#!/usr/bin/env python3
"""Write reports/TOD_보고서_11.md from the figure files and stored tables."""

import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
MAN = ROOT / "figures" / "manuscript"


def load(path):
    return pd.read_csv(ROOT / path, sep="\t", dtype=str, keep_default_na=False)


def md(frame):
    cols = list(frame.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def page_mm(pdf):
    text = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
    for line in text.splitlines():
        if line.startswith("Page size:"):
            parts = line.split()
            return float(parts[2]) * 25.4 / 72, float(parts[4]) * 25.4 / 72
    raise SystemExit(pdf)


def main():
    stems = [
        ("Fig1_design", "Fig1"), ("Fig2_trajectories", "Fig2"), ("Fig3_effects", "Fig3"),
        ("Fig4_composition", "Fig4"), ("Fig5_modules", "Fig5"), ("Fig6_public", "Fig6"),
        ("Fig7_bias", "Fig7"), ("Fig8_mitigation", "Fig8"),
        ("FigS1_effects_secondary", "FigS1"), ("FigS2_bias_secondary", "FigS2"),
        ("FigS3_mitigation_secondary", "FigS3"),
    ]
    rows = []
    for stem, sub in stems:
        w, h = page_mm(MAN / f"{stem}.pdf")
        tiff_mb = (MAN / f"{stem}.tiff").stat().st_size / 1e6
        if sub.startswith("FigS"):
            copy = f"figures/submission/{sub}.pdf"
        else:
            copy = f"figures/submission/{sub}.eps, figures/submission/{sub}.tif"
        rows.append({
            "stem": stem, "pdf_mm": f"{w:.2f} x {h:.2f}", "tiff_MB": f"{tiff_mb:.2f}",
            "eps": f"{stem}.eps", "png": f"{stem}.png", "submission": copy,
        })
    audit = load("figures/manuscript/source_data/audit.tsv")
    audit = audit[["stem", "text_overlap", "text_marker", "legend_data", "outside", "min_fontsize", "alpha_not_1"]]
    peaks = load("figures/manuscript/source_data/Fig3_C.tsv")
    peaks = peaks[["clock", "dataset", "peak_h", "peak_center", "peak_halfwidth_95", "lrt_p", "phi_bar"]]
    genes = load("results/posthoc_stage7c/gse39445_gene_peaks.tsv")
    axes = sorted(set(genes["time_axis"]))
    knn = load("results/stage7b/d1_development_model_m.tsv")
    part = load("results/stage7c/d11_development_model_m.tsv")
    clocks = ["RA_default", "RA_Peters", "RA_GTEx"]
    crows = []
    for scheme in ("unmasked", "random", "low_expression"):
        for clock in clocks:
            a = knn.loc[(knn["scheme"] == scheme) & (knn["clock"] == clock)].iloc[0]
            b = part.loc[(part["scheme"] == scheme) & (part["clock"] == clock)].iloc[0]
            crows.append({"scheme": scheme, "clock": clock, "lrt_p_knn_file": a.lrt_p, "lrt_p_participant_mean_file": b.lrt_p})
    audit_mask = load("results/stage7b/d1_development_audit.tsv")
    bias = load("results/stage9/sim_design1_bias.tsv")
    s4 = bias.loc[(bias["clock"].isin(["RA_Peters", "RA_GTEx"])) & (bias["scenario"] == "lab") & (bias["s"] == "4"), ["clock", "scenario", "s", "b", "b_over_total_sd", "unit"]]
    text = f"""# TOD 보고서 11 — 원고 그림 최종 수정

## 1. 파일 목록 (그림별 PDF·EPS·TIFF·PNG, 크기 mm, TIFF MB, submission 사본)

PDF 쪽 크기는 `pdfinfo` Page size를 pt × 25.4 / 72로 둔 값이다. TIFF는 1200 dpi, LZW다. PNG는 300 dpi다. 글꼴을 포함한 EPS는 `figures/manuscript/`에 있고, 본문 투고용 사본은 `figures/submission/Fig1.eps`–`Fig8.eps`와 `Fig1.tif`–`Fig8.tif`다. 보충 그림 투고용 사본은 `FigS1.pdf`, `FigS2.pdf`, `FigS3.pdf`다.

{md(pd.DataFrame(rows))}

검토 묶음: `figures/manuscript/all_figures_review.pdf`.

## 2. 글꼴 (사용 글꼴, 설치 방법, pdffonts 결과, 가장 작은 글자 크기)

Arial은 없다. Liberation Sans는 패키지 `fonts-liberation` 1:1.07.4-11로 `/usr/share/fonts/truetype/liberation/`에 있다. 스크립트가 `fontManager.addfont()`로 Regular, Bold, Italic, BoldItalic TTF를 등록했다. `mathtext.fontset`은 custom이고 rm, it, bf, cal, sf, tt는 Liberation Sans다.

그림 PDF 11개와 `all_figures_review.pdf`의 `pdffonts`는 모두 emb yes다. EPS 11개를 `ps2pdf -dEPSCrop`로 바꾼 뒤의 `pdffonts`도 모두 emb yes다. 글꼴 이름은 LiberationSans, LiberationSans-Bold, LiberationSans-Italic이다.

가장 작은 `get_fontsize()`는 그림 11개 모두 8이다.

## 3. 점검 결과 (§5의 1–5 표)

겹침 네 가지와 `alpha`가 1이 아닌 artist 수는 0이다. 같은 축의 눈금 글자끼리의 겹침은 세지 않았다.

{md(audit)}

쪽 크기와 TIFF MB는 §1 표와 같다.

## 4. 확인용 값 (원고 수치 대조용, 저장값 그대로)

### 4.1 Fig. 3c: 시계 × 자료별 peak_h, peak_center, peak_halfwidth_95, LRT p, φ̄

`results/stage8/effect_summary.tsv`의 kind cosinor와 `prereg/frozen/F1.tsv`의 GSE39445 행이다. LRT p는 effect_summary의 `p`다.

{md(peaks)}

### 4.2 Fig. 6a: x에 쓴 GSE39445 시각 변수 이름과 정의

x는 `results/stage8/gene_peak_differences_arntl_alias.tsv`의 `peak_gse39445`다. 그 값은 `results/posthoc_stage7c/gse39445_gene_peaks.tsv`의 `peak_h`와 같고, 그 파일의 `time_axis`는 {", ".join(axes)}이다. 축 이름은 "Peak, GSE39445 (h since wake)"다.

### 4.3 Fig. 6c: 시계 × 가림 방식 × 처리별 LRT p

kNN 파일은 `results/stage7b/d1_development_model_m.tsv`다. 참가자 평균 파일은 `results/stage7c/d11_development_model_m.tsv`다. scheme이 unmasked일 때 두 파일의 LRT p는 세 시계 모두 같다. 그림의 Unmasked는 한 칸이다. random과 low_expression의 kNN은 7b 파일, 참가자 평균은 7c 파일이다.

가림 수는 `results/stage7b/d1_development_audit.tsv`에 있다.

{md(audit_mask)}

{md(pd.DataFrame(crows))}

### 4.4 Fig. 8a와 S3a: 채움에 쓴 파일과 열 이름

파일은 `results/stage9/composition_cosinor.tsv`다. 왼쪽 점은 `lrt_p_before`, 오른쪽 점은 `lrt_p_after`다. 각 열이 0.05 미만이면 채운 표식이다.

### 4.5 Fig. S2d: RNAAgeCalc-Peters와 RNAAgeCalc-GTEx Laboratory의 s = 4 |b|/총 SD

파일은 `results/stage9/sim_design1_bias.tsv`다. 열 `b_over_total_sd`는 부호가 있는 b / 총 SD다.

{md(s4)}

## 5. 판단이 필요한 사항

- U+0304(결합 매크론)은 Liberation Sans에 없다. φ̄는 mathtext `$\\bar{{\\phi}}$`로 그렸다.
- Fig. 3b와 Fig. S1b의 개발 범위는 관측 d보다 y로 0.18 아래에 두었다.
- Fig. 3b와 Fig. S1b의 범례는 패널 오른쪽 여백에 두었다. 오른쪽 아래 안쪽에는 관측 구간이 있다.
- Fig. 4b의 구간은 `axis(peak_center) ± peak_halfwidth_95`다. `hour mod 24`가 18 미만이면 24를 더한다. 구간 끝을 따로 옮기면 18시를 가로지르는 구간이 두 동강이 난다.
- Fig. 4b에서 같은 집단의 자료 4개는 y를 0.16 간격으로 벌렸다.
- Fig. 5a 표식 면적은 `14 * -log10(p)`다. 범례의 p = 0.05, 10⁻⁴, 10⁻⁸도 같은 식이다.
- Fig. 5b에서 Features, Confirmed, Predicted, Observed, Agree는 90도 돌렸다. 짙은 칸의 글자는 면 RGB 평균이 0.55 미만일 때 흰색이다.
- Fig. 5b의 Agree는 `sign_agrees`가 TRUE이면 Yes, FALSE이면 No다. 그 칸이 비어 있으면 "—"다. 빈 칸은 `sign_eligible`이 FALSE인 행이다.
- 특징 수가 0이면 Features는 "—"다. Confirmed는 저장된 값이다.
- Fig. 6a는 동일 축척이 아니다.
- Fig. 8b와 Fig. S3b의 범주 글자는 30도 돌렸다.
- Fig. 1의 GSE295758는 `n_subjects` 58을 "58 pairs"로 적었다. 관문 파일에는 PER1만 있다.
- Fig. 2의 모형 곡선은 저장된 `curve_t`, `curve_y`다.

## 6. 지시서와 다르게 한 것 (빠짐없이)

- Fig. 3b와 Fig. S1b 범례는 패널 안 오른쪽 아래가 아니라 패널 오른쪽 여백에 있다.
- Fig. 5b의 Agree는 `sign_agrees`가 비어 있으면 Yes나 No가 아니라 "—"다.
- Fig. 5b의 Features, Confirmed, Predicted, Observed, Agree는 가로가 아니라 90도 돌린 글자다.
- Fig. 8b와 Fig. S3b의 범주 글자는 30도 돌렸다.
- φ̄는 U+03C6에 U+0304를 붙이지 않고 mathtext 가로선으로 그렸다.

## 7. 커밋 해시

이 보고서를 넣는 커밋. 스크립트 `scripts/62_part11_figures.py`, `scripts/63_part11_report.py`. 그림 `figures/manuscript/`, `figures/submission/`.
"""
    out = ROOT / "reports" / "TOD_보고서_11.md"
    out.write_text(text, encoding="utf-8")
    print(out, out.stat().st_size)


if __name__ == "__main__":
    main()
