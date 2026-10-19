#!/usr/bin/env python3
"""Write reports/TOD_보고서_10.md from the manuscript figure files."""

import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
MAN = ROOT / "figures" / "manuscript"


def load(path):
    return pd.read_csv(ROOT / path, sep="\t", dtype=str, keep_default_na=False)


def md(frame):
    cols = list(frame.columns)
    header = "| " + " | ".join(cols) + " |"
    rule = "| " + " | ".join(["---"] * len(cols)) + " |"
    rows = ["| " + " | ".join(row.astype(str)) + " |" for _, row in frame.iterrows()]
    return "\n".join([header, rule, *rows])


def page_mm(pdf):
    text = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
    width = height = None
    for line in text.splitlines():
        if line.startswith("Page size:"):
            parts = line.split()
            width = float(parts[2])
            height = float(parts[4])
    return width, height, width * 25.4 / 72, height * 25.4 / 72


def font_lines(pdf):
    text = subprocess.check_output(["pdffonts", str(pdf)], text=True)
    rows = []
    for line in text.splitlines()[2:]:
        if not line.strip():
            continue
        emb = "yes" if " yes " in f" {line} " else "no"
        rows.append((line.split()[0], emb))
    return rows


def main():
    stems = [
        "Fig1_design", "Fig2_trajectories", "Fig3_effects", "Fig4_composition", "Fig5_modules",
        "Fig6_public", "Fig7_bias", "Fig8_mitigation", "FigS1_effects_secondary",
        "FigS2_bias_secondary", "FigS3_mitigation_secondary",
    ]
    file_rows = []
    font_rows = []
    for stem in stems:
        pdf = MAN / f"{stem}.pdf"
        wpt, hpt, wmm, hmm = page_mm(pdf)
        file_rows.append({
            "stem": stem,
            "pdf_pt": f"{wpt:.2f} x {hpt:.2f}",
            "pdf_mm": f"{wmm:.2f} x {hmm:.2f}",
            "png": f"{stem}.png",
            "tiff": f"{stem}.tiff",
        })
        for name, emb in font_lines(pdf):
            font_rows.append({"file": pdf.name, "font": name.split("+")[-1], "embedded": emb})
    fonts = pd.DataFrame(font_rows)
    if (fonts["embedded"] != "yes").any():
        raise SystemExit("a font is not embedded")
    review_fonts = pd.DataFrame([{"file": "all_figures_review.pdf", "font": name.split("+")[-1], "embedded": emb} for name, emb in font_lines(MAN / "all_figures_review.pdf")])
    if (review_fonts["embedded"] != "yes").any():
        raise SystemExit("review pdf font is not embedded")
    peaks = load("results/stage10/fig4b_composition_peaks.tsv")
    show = peaks[[
        "dataset", "population", "n_subjects", "n_samples", "n_dropped", "two_a", "peak_h",
        "peak_center", "peak_halfwidth_95", "lrt_p", "n_boot", "n_boot_finite",
    ]]
    summary = load("figures/manuscript/source_data/Fig2_summary.tsv")
    sources = sorted(p.name for p in (MAN / "source_data").glob("*.tsv"))
    text = f"""# TOD 보고서 10 — 원고 그림

## 1. 파일 목록 (그림별 PDF·TIFF·PNG, 크기 mm, 원자료 표)

PDF 쪽 크기는 `pdfinfo`의 Page size다. mm는 pt × 25.4 / 72이다. 너비는 174.00 mm이다. TIFF는 600 dpi, LZW다. PNG는 300 dpi다.

{md(pd.DataFrame(file_rows))}

검토용 묶음: `figures/manuscript/all_figures_review.pdf`. 순서: {", ".join(stems)}.

원자료 표 (`figures/manuscript/source_data/`): {", ".join(sources)}.

## 2. 글꼴 확인 (pdffonts 결과: 글꼴 포함 여부)

Arial은 설치되어 있지 않다. 그림 글꼴은 DejaVu Sans다. 그림 PDF 11개와 `all_figures_review.pdf`의 `pdffonts` 행은 모두 `emb` yes다. 글꼴 이름은 DejaVuSans, DejaVuSans-Bold이고 형식은 CID TrueType이다.

## 3. Fig. 4B 새 계산 표 (집단 × 자료: n, 2A, 정점 [95% 구간], LRT p)

`results/stage10/fig4b_composition_peaks.tsv`. 시드 20261019. 참가자 부트스트랩 2000회. 정점 구간은 `peak_center`와 `peak_halfwidth_95`다. GSE48113은 조성 점수가 유한하지 않은 샘플 1개를 뺐다. 나머지 자료의 `n_dropped`는 0이다. `n_boot_finite`는 모두 2000이다.

명령 `/usr/bin/time -v Rscript scripts/59_part10_fig4b.R`. 사용자 시간 758.03 s, 경과 12:37.60, 최대 RSS 412268 kB, 종료 상태 0.

{md(show)}

## 4. Fig. 2 구간 평균 표 요약 (패널별 구간 수, 참가자 수)

구간은 벽시계 시각 [0, 2), [2, 4), …, [22, 24)이다. 참가자 안 평균을 낸 뒤 참가자 사이 평균과 SEM을 냈다. `n_bins`는 그 패널에서 참가자가 하나 이상인 구간 수다.

{md(summary)}

## 5. 판단이 필요한 사항

- Fig. 4B 부트스트랩의 `lmer` 제어는 `calc.derivs=FALSE`, `check.conv.singular='ignore'`다. 점추정은 `fit_cosinor`의 기본 제어다. 수렴은 `conv$opt`다.
- Fig. 4B에서 같은 집단의 자료 4개는 y를 0.16 간격으로 벌렸다. x는 정점 시각이다. 12–36 밖으로 나간 구간은 화살표로 잘랐다.
- Fig. 3C의 예측 정점 반지름은 지시서에 없다. 1.15에 두었다.
- Fig. 3B의 개발 예측 범위와 빈 마름모는 관측 d보다 y로 0.18 아래에 두었다.
- Fig. 5A 표식 면적은 `14 * -log10(p)`다.
- Fig. 6B에서 `n_entered`가 0인 점은 로그 축에 두지 않고 축 왼쪽 끝에 숫자 0으로 적었다.
- Fig. 7A의 `1.1 years`와 `0.6 years`는 지시서의 표기다. 저장된 b는 lab `1.09380254924844`, free-living `0.582482849833315`이다.
- Fig. 8A와 Fig. S3A의 비율 글자는 저장된 `two_a_ratio`를 소수 둘째 자리로 적었다. 원자료 표에는 저장값을 두었다.
- Fig. 1의 GSE295758 관문 파일 `results/confirm/L3/gate.tsv`에는 PER1만 있다. `dataset_pass`는 TRUE다.
- Fig. S3B에는 시각 보정 열이 비어 있는 보조 시계의 time 점을 그리지 않았다.

## 6. 지시서와 다르게 한 것

없음.

## 7. 커밋 해시

이 보고서를 넣는 커밋. 스크립트 `scripts/59_part10_fig4b.R`, `scripts/60_part10_figures.py`, `scripts/61_part10_report.py`. 표 `results/stage10/fig4b_composition_peaks.tsv`. 그림 `figures/manuscript/`.
"""
    out = ROOT / "reports" / "TOD_보고서_10.md"
    out.write_text(text, encoding="utf-8")
    print(f"wrote {out} ({out.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
