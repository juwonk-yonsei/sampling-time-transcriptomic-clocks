# TOD 보고서 10 — 원고 그림

## 1. 파일 목록 (그림별 PDF·TIFF·PNG, 크기 mm, 원자료 표)

PDF 쪽 크기는 `pdfinfo`의 Page size다. mm는 pt × 25.4 / 72이다. 너비는 174.00 mm이다. TIFF는 600 dpi, LZW다. PNG는 300 dpi다.

| stem | pdf_pt | pdf_mm | png | tiff |
| --- | --- | --- | --- | --- |
| Fig1_design | 493.23 x 476.22 | 174.00 x 168.00 | Fig1_design.png | Fig1_design.tiff |
| Fig2_trajectories | 493.23 x 425.20 | 174.00 x 150.00 | Fig2_trajectories.png | Fig2_trajectories.tiff |
| Fig3_effects | 493.23 x 623.62 | 174.00 x 220.00 | Fig3_effects.png | Fig3_effects.tiff |
| Fig4_composition | 493.23 x 566.93 | 174.00 x 200.00 | Fig4_composition.png | Fig4_composition.tiff |
| Fig5_modules | 493.23 x 595.28 | 174.00 x 210.00 | Fig5_modules.png | Fig5_modules.tiff |
| Fig6_public | 493.23 x 623.62 | 174.00 x 220.00 | Fig6_public.png | Fig6_public.tiff |
| Fig7_bias | 493.23 x 623.62 | 174.00 x 220.00 | Fig7_bias.png | Fig7_bias.tiff |
| Fig8_mitigation | 493.23 x 425.20 | 174.00 x 150.00 | Fig8_mitigation.png | Fig8_mitigation.tiff |
| FigS1_effects_secondary | 493.23 x 481.89 | 174.00 x 170.00 | FigS1_effects_secondary.png | FigS1_effects_secondary.tiff |
| FigS2_bias_secondary | 493.23 x 623.62 | 174.00 x 220.00 | FigS2_bias_secondary.png | FigS2_bias_secondary.tiff |
| FigS3_mitigation_secondary | 493.23 x 425.20 | 174.00 x 150.00 | FigS3_mitigation_secondary.png | FigS3_mitigation_secondary.tiff |

검토용 묶음: `figures/manuscript/all_figures_review.pdf`. 순서: Fig1_design, Fig2_trajectories, Fig3_effects, Fig4_composition, Fig5_modules, Fig6_public, Fig7_bias, Fig8_mitigation, FigS1_effects_secondary, FigS2_bias_secondary, FigS3_mitigation_secondary.

원자료 표 (`figures/manuscript/source_data/`): Fig1_A.tsv, Fig2_A.tsv, Fig2_B.tsv, Fig2_C.tsv, Fig2_D.tsv, Fig2_E.tsv, Fig2_F.tsv, Fig2_summary.tsv, Fig3_A.tsv, Fig3_B.tsv, Fig3_C.tsv, Fig4_A.tsv, Fig4_B.tsv, Fig5_A.tsv, Fig5_B.tsv, Fig6_A.tsv, Fig6_B.tsv, Fig6_C.tsv, Fig7_A.tsv, Fig7_B.tsv, Fig7_C.tsv, Fig7_D.tsv, Fig7_E.tsv, Fig8_A.tsv, Fig8_B.tsv, FigS1_A.tsv, FigS1_B.tsv, FigS2_A.tsv, FigS2_B.tsv, FigS2_C.tsv, FigS2_D.tsv, FigS2_E.tsv, FigS3_A.tsv, FigS3_B.tsv.

## 2. 글꼴 확인 (pdffonts 결과: 글꼴 포함 여부)

Arial은 설치되어 있지 않다. 그림 글꼴은 DejaVu Sans다. 그림 PDF 11개와 `all_figures_review.pdf`의 `pdffonts` 행은 모두 `emb` yes다. 글꼴 이름은 DejaVuSans, DejaVuSans-Bold이고 형식은 CID TrueType이다.

## 3. Fig. 4B 새 계산 표 (집단 × 자료: n, 2A, 정점 [95% 구간], LRT p)

`results/stage10/fig4b_composition_peaks.tsv`. 시드 20261019. 참가자 부트스트랩 2000회. 정점 구간은 `peak_center`와 `peak_halfwidth_95`다. GSE48113은 조성 점수가 유한하지 않은 샘플 1개를 뺐다. 나머지 자료의 `n_dropped`는 0이다. `n_boot_finite`는 모두 2000이다.

명령 `/usr/bin/time -v Rscript scripts/59_part10_fig4b.R`. 사용자 시간 758.03 s, 경과 12:37.60, 최대 RSS 412268 kB, 종료 상태 0.

| dataset | population | n_subjects | n_samples | n_dropped | two_a | peak_h | peak_center | peak_halfwidth_95 | lrt_p | n_boot | n_boot_finite |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | T cells | 24 | 221 | 0 | 0.37047951805148 | 2.01178453285483 | 2.03099048028316 | 1.15405702544208 | 4.6785097735e-20 | 2000 | 2000 |
| GSE39445 | CD8 T cells | 24 | 221 | 0 | 0.456416154866155 | 1.46898840392427 | 1.46172072434736 | 1.02141417438358 | 8.86356306008214e-18 | 2000 | 2000 |
| GSE39445 | Cytotoxic lymphocytes | 24 | 221 | 0 | 0.231865236884741 | 13.4282218619436 | 13.4408762843507 | 1.61864135317054 | 1.92009443964517e-06 | 2000 | 2000 |
| GSE39445 | NK cells | 24 | 221 | 0 | 0.429616717109436 | 13.6363947099391 | 13.6328421808056 | 0.935929380938144 | 3.1680572247023e-19 | 2000 | 2000 |
| GSE39445 | B lineage | 24 | 221 | 0 | 0.0796877116596472 | 0.543554859815641 | 0.563950669268381 | 4.55583118373969 | 0.0944215101801147 | 2000 | 2000 |
| GSE39445 | Monocytic lineage | 24 | 221 | 0 | 0.0943620945050906 | 1.48106304939847 | 1.55194370008111 | 1.95426989424538 | 5.45963753263961e-05 | 2000 | 2000 |
| GSE39445 | Myeloid dendritic cells | 24 | 221 | 0 | 0.125350227474236 | 23.7180924279838 | 23.732718767764 | 3.32358109791331 | 0.0376237493997132 | 2000 | 2000 |
| GSE39445 | Neutrophils | 24 | 221 | 0 | 0.491398502744687 | 14.7133857365347 | 14.7045276304474 | 1.15972749133112 | 1.72783871706906e-17 | 2000 | 2000 |
| GSE56931 | T cells | 14 | 130 | 0 | 0.252352434185952 | 0.268450187540798 | 0.284454997231921 | 1.56420662298931 | 1.02272255210566e-09 | 2000 | 2000 |
| GSE56931 | CD8 T cells | 14 | 130 | 0 | 0.365374367890902 | 0.821668717174426 | 0.809063274792159 | 1.34923953469051 | 7.44808106876266e-10 | 2000 | 2000 |
| GSE56931 | Cytotoxic lymphocytes | 14 | 130 | 0 | 0.274071086486198 | 10.888296991903 | 10.8657372371338 | 1.35133122550403 | 3.92173321072081e-08 | 2000 | 2000 |
| GSE56931 | NK cells | 14 | 130 | 0 | 0.508119137506881 | 10.9420393769587 | 10.9029287017295 | 1.04914483639715 | 1.21575312140795e-11 | 2000 | 2000 |
| GSE56931 | B lineage | 14 | 130 | 0 | 0.196618375878883 | 22.5809101432576 | 22.6014645438817 | 1.6119672900677 | 2.98864519498149e-07 | 2000 | 2000 |
| GSE56931 | Monocytic lineage | 14 | 130 | 0 | 0.115108886651412 | 13.2327615839505 | 13.3212146736105 | 2.26310164220032 | 0.00388916842751795 | 2000 | 2000 |
| GSE56931 | Myeloid dendritic cells | 14 | 130 | 0 | 0.0676595364044592 | 17.8752932282875 | 17.8277300027295 | 3.17166102122748 | 0.0597202722447278 | 2000 | 2000 |
| GSE56931 | Neutrophils | 14 | 130 | 0 | 0.184671515887009 | 14.2483972010018 | 14.1965158005652 | 1.47460169400801 | 3.3700439462561e-07 | 2000 | 2000 |
| GSE113883 | T cells | 11 | 153 | 0 | 0.501333019807356 | 1.74942881938888 | 1.71733596691049 | 1.30330434643716 | 4.99931417782842e-11 | 2000 | 2000 |
| GSE113883 | CD8 T cells | 11 | 153 | 0 | 0.419891128542506 | 0.619346944281585 | 0.599665593645633 | 1.49102364177043 | 1.15565194309793e-08 | 2000 | 2000 |
| GSE113883 | Cytotoxic lymphocytes | 11 | 153 | 0 | 0.243090121688901 | 9.50166976042667 | 9.56253755659584 | 2.76903957924238 | 0.000732999198460937 | 2000 | 2000 |
| GSE113883 | NK cells | 11 | 153 | 0 | 0.4851262151375 | 11.3747771559794 | 11.3713622026429 | 1.49038082714599 | 5.23804105419159e-12 | 2000 | 2000 |
| GSE113883 | B lineage | 11 | 153 | 0 | 0.332784618573013 | 23.9510317757198 | 23.9229433568439 | 1.68613079099825 | 1.42327204249322e-08 | 2000 | 2000 |
| GSE113883 | Monocytic lineage | 11 | 153 | 0 | 0.123338435985523 | 3.56329592771249 | 3.40114334504273 | 3.97471451640269 | 0.0315096138649216 | 2000 | 2000 |
| GSE113883 | Myeloid dendritic cells | 11 | 153 | 0 | 0.125222679215967 | 0.754093779900894 | 0.702356724057426 | 3.11049824833534 | 0.0141134915362991 | 2000 | 2000 |
| GSE113883 | Neutrophils | 11 | 153 | 0 | 0.283981120369161 | 13.8007121554087 | 13.7938673616655 | 2.7377556800773 | 0.00187867425470749 | 2000 | 2000 |
| GSE48113 | T cells | 22 | 146 | 1 | 0.436217526043816 | 1.71982939550351 | 1.67472317385882 | 1.21127151857267 | 9.46728158859098e-09 | 2000 | 2000 |
| GSE48113 | CD8 T cells | 22 | 146 | 1 | 0.540261495642269 | 1.04961295081186 | 1.042880018829 | 1.1287136681649 | 1.29738358924692e-13 | 2000 | 2000 |
| GSE48113 | Cytotoxic lymphocytes | 22 | 146 | 1 | 0.126499636459599 | 9.42427095572418 | 9.14589501575537 | 4.53655030030652 | 0.0645790604933873 | 2000 | 2000 |
| GSE48113 | NK cells | 22 | 146 | 1 | 0.409266554512416 | 12.2169080385766 | 12.2213182738741 | 1.40742476781953 | 1.96647316419963e-07 | 2000 | 2000 |
| GSE48113 | B lineage | 22 | 146 | 1 | 0.0753524731147373 | 2.52326772409224 | 2.66993955104818 | 10.6570974575096 | 0.541521651002368 | 2000 | 2000 |
| GSE48113 | Monocytic lineage | 22 | 146 | 1 | 0.118784207153735 | 4.97374251824763 | 5.241774874122 | 3.75554586587364 | 0.0363735503261206 | 2000 | 2000 |
| GSE48113 | Myeloid dendritic cells | 22 | 146 | 1 | 0.15476726112969 | 1.29508684878054 | 1.26895788359483 | 3.2578121422068 | 0.0509949056397604 | 2000 | 2000 |
| GSE48113 | Neutrophils | 22 | 146 | 1 | 0.504889525181411 | 14.6352979255312 | 14.6736307660879 | 1.03078612067456 | 2.52153749159091e-10 | 2000 | 2000 |

## 4. Fig. 2 구간 평균 표 요약 (패널별 구간 수, 참가자 수)

구간은 벽시계 시각 [0, 2), [2, 4), …, [22, 24)이다. 참가자 안 평균을 낸 뒤 참가자 사이 평균과 SEM을 냈다. `n_bins`는 그 패널에서 참가자가 하나 이상인 구간 수다.

| panel | dataset | clock | n_subjects | n_samples | n_bins | p | peak_h | beta_source |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | GSE113883 | TAms_chrono | 11 | 153 | 12 | 2.96936394371969e-10 | 12.188050742382 | results/stage5/h1_point_estimates.tsv |
| B | GSE113883 | TAms_mort | 11 | 153 | 12 | 4.44466652817908e-06 | 10.0600990696032 | results/stage5/h1_point_estimates.tsv |
| C | GSE113883 | RA_default | 11 | 153 | 12 | 3.24675832907807e-05 | 11.8126319723118 | results/stage4/rnaseq/h1_point_estimates.tsv |
| D | GSE48113 | TAms_chrono | 22 | 147 | 12 | 2.70474257565244e-08 | 10.879869555482 | results/confirm/family_P.tsv |
| E | GSE48113 | TAms_mort | 22 | 147 | 12 | 2.67642593880552e-05 | 11.4858304321163 | results/confirm/family_P.tsv |
| F | GSE48113 | RA_default | 22 | 147 | 12 | 0.0125285842066838 | 11.5918483633937 | results/confirm/family_P.tsv |

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
