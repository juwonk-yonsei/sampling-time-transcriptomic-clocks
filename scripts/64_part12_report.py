#!/usr/bin/env python3
"""Write reports/TOD_보고서_12.md from stored audit tables and file measurements."""

import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path("/HDD8T1/WORK/GeroSci/TOD")
MAN = ROOT / "figures" / "manuscript"
AUDIT = pd.read_csv(MAN / "source_data" / "audit.tsv", sep="\t")

STEMS = [
    "Fig1_design", "Fig2_trajectories", "Fig3_effects", "Fig4_composition", "Fig5_modules",
    "Fig6_public", "Fig7_bias", "Fig8_mitigation", "FigS1_effects_secondary", "FigS2_bias_secondary",
    "FigS3_mitigation_secondary",
]
SUB = {
    "Fig1_design": "figures/submission/Fig1.eps, figures/submission/Fig1.tif",
    "Fig2_trajectories": "figures/submission/Fig2.eps, figures/submission/Fig2.tif",
    "Fig3_effects": "figures/submission/Fig3.eps, figures/submission/Fig3.tif",
    "Fig4_composition": "figures/submission/Fig4.eps, figures/submission/Fig4.tif",
    "Fig5_modules": "figures/submission/Fig5.eps, figures/submission/Fig5.tif",
    "Fig6_public": "figures/submission/Fig6.eps, figures/submission/Fig6.tif",
    "Fig7_bias": "figures/submission/Fig7.eps, figures/submission/Fig7.tif",
    "Fig8_mitigation": "figures/submission/Fig8.eps, figures/submission/Fig8.tif",
    "FigS1_effects_secondary": "figures/submission/FigS1.pdf",
    "FigS2_bias_secondary": "figures/submission/FigS2.pdf",
    "FigS3_mitigation_secondary": "figures/submission/FigS3.pdf",
}


def page_mm(path):
    text = subprocess.check_output(["pdfinfo", str(path)], text=True)
    for line in text.splitlines():
        if line.startswith("Page size:"):
            parts = line.split()
            width, height = float(parts[2]), float(parts[4])
            return width * 25.4 / 72, height * 25.4 / 72
    raise SystemExit(f"no page size {path}")


def font_names(path):
    text = subprocess.check_output(["pdffonts", str(path)], text=True)
    rows = text.splitlines()
    emb_at = rows[0].index("emb")
    names, emb = set(), set()
    for line in rows[2:]:
        if not line.strip():
            continue
        names.add(line.split()[0].split("+", 1)[-1])
        emb.add(line[emb_at:emb_at + 3].strip())
    return sorted(names), sorted(emb)


rows = []
for stem in STEMS:
    width, height = page_mm(MAN / f"{stem}.pdf")
    tiff_mb = (MAN / f"{stem}.tiff").stat().st_size / 1e6
    rows.append((stem, width, height, tiff_mb))

pdf_fonts = {}
for stem in STEMS:
    pdf_fonts[stem] = font_names(MAN / f"{stem}.pdf")
review_fonts = font_names(MAN / "all_figures_review.pdf")

eps_fonts = {}
tmp = Path("/tmp/tod12_eps")
tmp.mkdir(exist_ok=True)
for stem in STEMS:
    out = tmp / f"{stem}.pdf"
    subprocess.check_call(["ps2pdf", "-dEPSCrop", str(MAN / f"{stem}.eps"), str(out)])
    eps_fonts[stem] = font_names(out)

lines = []
lines.append("# TOD 보고서 12 — 그림 마지막 손질")
lines.append("")
lines.append("## 1. 바꾼 그림과 파일 (크기 mm, TIFF MB)")
lines.append("")
lines.append("PDF 쪽 크기는 `pdfinfo` Page size를 pt × 25.4 / 72로 둔 값이다. TIFF MB는 파일 바이트 수를 1e6으로 나눈 값이다. TIFF는 1200 dpi, LZW다. PNG는 300 dpi다. 11개 그림을 다시 썼다.")
lines.append("")
lines.append("| stem | pdf_mm | tiff_MB | eps | png | submission |")
lines.append("| --- | --- | --- | --- | --- | --- |")
for stem, width, height, tiff_mb in rows:
    lines.append(f"| {stem} | {width:.2f} x {height:.2f} | {tiff_mb:.2f} | {stem}.eps | {stem}.png | {SUB[stem]} |")
lines.append("")
lines.append("검토 묶음: `figures/manuscript/all_figures_review.pdf`.")
lines.append("")
lines.append("## 2. 점검 결과 (§3의 1–6 표)")
lines.append("")
lines.append("PDF 잘림은 `pdftotext -bbox` 단어 상자가 쪽 경계에서 1 mm 안쪽 선을 넘는 단어 수다. 1 mm는 72/25.4 pt다. 글자 밖은 `get_window_extent()`가 `fig.bbox` 밖인 글자 수다. 글자 겹침은 눈금 글자를 포함한 모든 글자 쌍이다. 지시선–글자는 annotate 화살표가 다른 글자 상자를 지나는 수다. 글자가 비어 있는 annotate는 빼었다. 축 밖 자료는 `transData`인 선과 점 가운데 축 범위 밖인 점 수다. 의도 자름은 Fig. 3c, 4a, 4b, S1c에서 화살표로 자른 구간 수다. 최소 글자는 `get_fontsize()`다. alpha는 None도 1도 아닌 artist 수다.")
lines.append("")
lines.append("| stem | pdf_clip | n_text | outside | text_overlap | leader_text | data_outside | intentional_clip | min_fontsize | alpha_not_1 |")
lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
for stem in STEMS:
    row = AUDIT.loc[AUDIT["stem"] == stem].iloc[0]
    lines.append(
        f"| {stem} | {int(row.pdf_clip)} | {int(row.n_text)} | {int(row.outside)} | {int(row.text_overlap)} | {int(row.leader_text)} | {int(row.data_outside)} | {int(row.intentional_clip)} | {row.min_fontsize} | {int(row.alpha_not_1)} |"
    )
lines.append("")
lines.append("자료 선이 글자 상자를 지나는 수(line_text), 글자 상자와 표식 중심이 겹치는 수(text_marker), 범례 상자에 들어간 자료 점 수(legend_data)는 위 표의 항목이 아니다. 11개 그림 모두 0이다.")
lines.append("")
pdf_name_set = sorted({name for names, _ in pdf_fonts.values() for name in names})
pdf_emb_set = sorted({emb for _, embs in pdf_fonts.values() for emb in embs})
eps_name_set = sorted({name for names, _ in eps_fonts.values() for name in names})
eps_emb_set = sorted({emb for _, embs in eps_fonts.values() for emb in embs})
lines.append(f"그림 PDF 11개의 `pdffonts` 글꼴은 {', '.join(pdf_name_set)}이고 emb는 {', '.join(pdf_emb_set)}이다. `all_figures_review.pdf`의 글꼴은 {', '.join(review_fonts[0])}이고 emb는 {', '.join(review_fonts[1])}이다. EPS 11개를 `ps2pdf -dEPSCrop`로 바꾼 뒤의 글꼴은 {', '.join(eps_name_set)}이고 emb는 {', '.join(eps_emb_set)}이다.")
lines.append("")
lines.append("## 3. 판단이 필요한 사항")
lines.append("")
lines.append("Fig. 6c에서 세 시계의 x 이동은 RA_default −0.12, RA_Peters 0, RA_GTEx +0.12이다. 지시서는 세 값만 적고 어느 시계에 어느 값을 둘지는 적지 않았다. 범례에 적은 순서대로 두었다.")
lines.append("")
lines.append("Fig. 8a와 S3a에서 같은 x의 자료 이동은 GSE39445 −0.06, GSE56931 −0.02, GSE113883 +0.02, GSE48113 +0.06이다. 지시서의 \"자료 순서\"는 이 네 자료의 기존 나열 순서로 읽었다.")
lines.append("")
lines.append("Fig. 3과 Fig. S1은 패널 a, b, c가 한 gridspec이다. b의 범례가 쪽 안에 들어가도록 right를 0.56으로 두었고, a와 c도 같은 너비다.")
lines.append("")
lines.append("Fig. 5a에서 옮긴 글자 좌표(자료 좌표, 왼쪽 정렬)는 chrono-red 이름 (6.2, 5.4), chrono-orange 이름 (6.2, 1.0), chrono-orange의 n = 1 (6.2, 8.2), chrono-brown 이름 (8.6, −1.6), chrono-yellow의 n = 1 (16.8, −0.6)이다. 나머지 이름과 n 좌표는 지시서 11과 같다. 좌표는 `figures/manuscript/source_data/Fig5_A_labels.tsv`에 있다.")
lines.append("")
lines.append("Fig. 6a에서 PER2 지시선과 CRY2 지시선은 다른 유전자 표식을 지난다. 이 수는 2이다. 지시서 12는 이 좌표를 바꾸지 않아서 지시서 11의 좌표를 유지했다.")
lines.append("")
lines.append("Fig. 6b 범례는 1열로 그림 왼쪽 아래 (0.02, 0.012)에 두고, Fig. 6c 범례는 그림 오른쪽 아래 (0.62, 0.012)에 두었다. 두 범례를 한 행에 두면 글자가 만난다.")
lines.append("")
lines.append("Fig. 7a와 S2a의 y 범위는 그 패널에 그린 값의 최솟값, 최댓값, 0을 포함한 구간에 위아래 5%를 더한 것이다.")
lines.append("")
lines.append("## 4. 지시서와 다르게 한 것 (빠짐없이)")
lines.append("")
lines.append("- Fig. 2 높이는 172 mm이다. 지시서 12는 234 mm 이하만 적는다. 지시서 11의 Fig. 2 높이는 156 mm였다.")
lines.append("- Fig. 2 패널 글자 a–f는 축 안 왼쪽 위(axes 좌표 0.03, 0.97, 세로 맞춤 top)에 두었다. 축 위 제목과 같은 높이에 두면 p 줄과 만난다.")
lines.append("- Fig. 6c 눈금 \"Participant mean\"은 \"Participant\"와 \"mean\" 두 줄이다.")
lines.append("- Fig. 7a와 S2a의 y = 0 선은 #BFBFBF, 0.5 pt이다.")
lines.append("- Fig. 8a와 S3a의 비율 지시선은 옮긴 점 (x, y)에서 (1.28, 글자를 놓은 y)로 가는 annotate이다. 비율 글자의 x는 1.32이다.")
lines.append("- 지시선–글자 수에서 글자가 비어 있는 annotate를 빼었다. 자른 구간은 intentional_clip으로 세었다. Fig. 4가 5이고 나머지 그림은 0이다.")
lines.append("- 자료 선이 글자 상자를 지나는 수, 글자–표식 수, 범례 안의 자료 점 수를 추가로 세었다. 11개 그림 모두 0이다. 지시서 §3의 1–6에는 이 세 수가 없다.")
lines.append("- Fig. 3과 Fig. S1의 패널 a와 c 축 너비도 b와 같다. §3에 적었다.")
lines.append("- Fig. 6c 시계별 x 이동의 시계 배정, Fig. 8a와 S3a의 자료별 x 이동 배정, Fig. 5a에서 옮긴 다섯 좌표, Fig. 6a PER2와 CRY2 지시선을 그대로 둔 것, Fig. 6 범례를 좌우로 나눈 것, Fig. 7a와 S2a의 y 범위에 0을 넣은 것은 §3에 적었다.")
lines.append("")
lines.append("## 5. 커밋 해시")
lines.append("")
lines.append("이 보고서를 넣는 커밋")
lines.append("")

out = ROOT / "reports" / "TOD_보고서_12.md"
out.write_text("\n".join(lines), encoding="utf-8")
print(out)
print("pdf fonts", pdf_name_set, pdf_emb_set)
print("eps fonts", eps_name_set, eps_emb_set)
print("review", review_fonts)
