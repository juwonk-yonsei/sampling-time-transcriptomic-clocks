# TOD 보고서 14 — 그림 마지막 두 가지

## 1. 바꾼 그림 (크기 mm, TIFF MB)

PDF 쪽 크기는 `pdfinfo` Page size를 pt × 25.4 / 72로 둔 값이다. TIFF MB는 파일 바이트 수를 1e6으로 나눈 값이다. 바꾼 그림은 Fig. 2, Fig. 6이다.

| stem | pdf_mm | tiff_MB | submission |
| --- | --- | --- | --- |
| Fig2_trajectories | 174.00 x 142.00 | 6.26 | figures/submission/Fig2.eps, figures/submission/Fig2.tif |
| Fig6_public | 174.00 x 220.00 | 5.89 | figures/submission/Fig6.eps, figures/submission/Fig6.tif |

Fig. 2 TIFF 바이트는 6258328이다. Fig. 6 TIFF 바이트는 5889520이다. 검토 묶음 `figures/manuscript/all_figures_review.pdf`와 위 submission 사본을 다시 썼다.

## 2. 점검 결과 (지시서 12 §3, 지시서 13 §2, 글자-선 겹침)

지시서 12 §3의 열은 pdf_clip, n_text, outside, text_overlap, leader_text, data_outside, intentional_clip, min_fontsize, alpha_not_1이다. 지시서 13 §2의 열은 leader_cross, leader_foreign, axes_text, ylabel_gap_mm이다. ylabel_gap_mm는 y축 이름과 y 눈금 글자의 왼쪽 끝이 왼쪽 이웃 축 테두리에서 떨어진 최소 거리(mm)다. 왼쪽 이웃 축이 없으면 없음이다. text_geom은 글자 상자가 연결선, 다른 지시선, 대각선, 표식과 겹친 수다. 그 글자 자신의 지시선은 세지 않는다.

| stem | pdf_clip | n_text | outside | text_overlap | leader_text | data_outside | intentional_clip | min_fontsize | alpha_not_1 | leader_cross | leader_foreign | axes_text | ylabel_gap_mm | text_geom |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Fig1_design | 0 | 15 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 | 0 |
| Fig2_trajectories | 0 | 100 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 2.58 | 0 |
| Fig3_effects | 0 | 67 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 | 0 |
| Fig4_composition | 0 | 43 | 0 | 0 | 0 | 0 | 5 | 8.0 | 0 | 0 | 0 | 0 | 없음 | 0 |
| Fig5_modules | 0 | 143 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 | 0 |
| Fig6_public | 0 | 86 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 10.08 | 0 |
| Fig7_bias | 0 | 98 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 9.42 | 0 |
| Fig8_mitigation | 0 | 96 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 5.60 | 0 |
| FigS1_effects_secondary | 0 | 50 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 | 0 |
| FigS2_bias_secondary | 0 | 95 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 9.42 | 0 |
| FigS3_mitigation_secondary | 0 | 65 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 21.39 | 0 |

Fig. 6a의 글자 14개 중 text_geom에 든 수는 0이다. Fig. 2와 Fig. 6의 PDF와 EPS에서 `pdffonts` 글꼴은 LiberationSans와 LiberationSans-Bold이고 emb는 yes다.

## 3. 지시서와 다르게 한 것 (빠짐없이)

- Fig. 2의 wspace는 0.52다. left 0.16, right 0.97, top 0.892, bottom 0.088, hspace 0.38이다. 그림 높이는 142 mm다. 행 이름 위치는 축 위쪽 끝 + 0.058이다. y축 글자와 왼쪽 이웃 축 테두리의 최소 거리는 2.58 mm다.
- Fig. 2a의 y 눈금은 −0.075, −0.050, −0.025, 0.000, 0.025, 0.050, 0.075, 0.100이다. b, c, d, e, f의 y 눈금은 지정하지 않았다.
- Fig. 6a의 PER2는 오른쪽 이름 열의 PER1과 PER3 사이에 두었다. ha는 left다. 글자 좌표는 (28.0, 23.501404955562684)다. 지시선은 `angle,angleA=180,angleB=90,rad=0`이다.
- PER2 글자 y는 PER3 대상 y 22.3757375270708과 PER1 대상 y 24.627072384054568의 평균이다.
- PER1 글자 좌표는 (28.0, 25.801404955562685)다. CRY2 글자 좌표는 (28.0, 27.67074663411381)다. PER3 글자 좌표는 (28.0, 21.201404955562683)다. NR1D2 글자 좌표는 (28.0, 18.901404955562683)다. NR1D1 글자 좌표는 (28.0, 16.601404955562682)다.
- CRY1 글자 좌표는 (28.0, 8.9670687601885)다. TEF, ARNTL, DBP의 글자 좌표는 지시서 13과 같다.
- Fig. 6의 b, c와 Fig. 2, Fig. 6 이외의 그림 파일은 다시 쓰지 않았다.

## 4. 커밋 해시

이 보고서를 넣는 커밋
