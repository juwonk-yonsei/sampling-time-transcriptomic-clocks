# TOD 보고서 13 — 그림 미세 수정

## 1. 바꾼 그림 (크기 mm, TIFF MB)

PDF 쪽 크기는 `pdfinfo` Page size를 pt × 25.4 / 72로 둔 값이다. TIFF MB는 파일 바이트 수를 1e6으로 나눈 값이다. 바꾼 그림은 Fig. 2, Fig. 3, Fig. 6, Fig. S1이다.

| stem | pdf_mm | tiff_MB | submission |
| --- | --- | --- | --- |
| Fig2_trajectories | 174.00 x 172.00 | 6.29 | figures/submission/Fig2.eps, figures/submission/Fig2.tif |
| Fig3_effects | 174.00 x 228.00 | 5.72 | figures/submission/Fig3.eps, figures/submission/Fig3.tif |
| Fig6_public | 174.00 x 220.00 | 5.87 | figures/submission/Fig6.eps, figures/submission/Fig6.tif |
| FigS1_effects_secondary | 174.00 x 188.00 | 3.94 | figures/submission/FigS1.pdf |

검토 묶음 `figures/manuscript/all_figures_review.pdf`와 위 submission 사본을 다시 썼다.

## 2. 점검 결과 (지시서 12 §3과 이번 §2의 표)

지시서 12 §3의 열은 pdf_clip, n_text, outside, text_overlap, leader_text, data_outside, intentional_clip, min_fontsize, alpha_not_1이다. 이번 §2의 열은 leader_cross, leader_foreign, axes_text, ylabel_gap_mm이다. ylabel_gap_mm는 y축 이름과 y 눈금 글자의 왼쪽 끝이 왼쪽 이웃 축 테두리에서 떨어진 최소 거리(mm)다. 왼쪽 이웃 축이 없으면 없음이다.

| stem | pdf_clip | n_text | outside | text_overlap | leader_text | data_outside | intentional_clip | min_fontsize | alpha_not_1 | leader_cross | leader_foreign | axes_text | ylabel_gap_mm |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Fig1_design | 0 | 15 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 |
| Fig2_trajectories | 0 | 100 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 9.92 |
| Fig3_effects | 0 | 67 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 |
| Fig4_composition | 0 | 43 | 0 | 0 | 0 | 0 | 5 | 8.0 | 0 | 0 | 0 | 0 | 없음 |
| Fig5_modules | 0 | 143 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 |
| Fig6_public | 0 | 86 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 10.08 |
| Fig7_bias | 0 | 98 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 9.42 |
| Fig8_mitigation | 0 | 96 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 5.60 |
| FigS1_effects_secondary | 0 | 50 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 |
| FigS2_bias_secondary | 0 | 95 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 9.42 |
| FigS3_mitigation_secondary | 0 | 65 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 21.39 |

Fig. S1 패널 기호 a의 단어 상자 위쪽 끝은 `pdftotext -bbox`의 yMin 24.263128 pt이고, pt × 25.4 / 72로 8.56 mm다.

## 3. 지시서와 다르게 한 것 (빠짐없이)

- Fig. 2의 wspace는 0.85다. left, right, top, bottom, hspace는 지시서 12의 값이다. e의 y 눈금은 바꾸지 않았다. y축 글자와 왼쪽 이웃 축 테두리의 최소 거리는 9.92 mm다.
- Fig. 6a 오른쪽 열의 글자 y는 대상 점의 y 순서를 유지하고, 이웃 글자 중심의 간격이 2.3 h보다 작으면 아래 글자를 내렸다. PER1의 대상 y는 24.627072384054568이고 PER3의 대상 y는 22.3757375270708이다. PER3 글자 y는 22.327072384054567이다.
- PER2 글자는 대상 점의 위에 두지 않고 왼쪽에 두었다. 글자 좌표는 (20.3330572724929, 22.7975702798119)이고 ha는 right다. 대상 x는 20.7330572724929다.
- TEF, ARNTL, DBP의 글자 좌표는 지시서 12와 같다.
- Fig. 3에서 a와 c의 축은 left 0.16, right 0.75다. b의 축은 left 0.16, right 0.56이다.
- Fig. S1에서 a와 c의 축은 left 0.18, right 0.75다. b의 축은 left 0.18, right 0.56이다. gridspec top은 0.93이다.

## 4. 커밋 해시

이 보고서를 넣는 커밋
