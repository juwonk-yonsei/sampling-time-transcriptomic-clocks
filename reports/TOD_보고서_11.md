# TOD 보고서 11 — 원고 그림 최종 수정

## 1. 파일 목록 (그림별 PDF·EPS·TIFF·PNG, 크기 mm, TIFF MB, submission 사본)

PDF 쪽 크기는 `pdfinfo` Page size를 pt × 25.4 / 72로 둔 값이다. TIFF는 1200 dpi, LZW다. PNG는 300 dpi다. 글꼴을 포함한 EPS는 `figures/manuscript/`에 있고, 본문 투고용 사본은 `figures/submission/Fig1.eps`–`Fig8.eps`와 `Fig1.tif`–`Fig8.tif`다. 보충 그림 투고용 사본은 `FigS1.pdf`, `FigS2.pdf`, `FigS3.pdf`다.

| stem | pdf_mm | tiff_MB | eps | png | submission |
| --- | --- | --- | --- | --- | --- |
| Fig1_design | 174.00 x 120.00 | 5.05 | Fig1_design.eps | Fig1_design.png | figures/submission/Fig1.eps, figures/submission/Fig1.tif |
| Fig2_trajectories | 174.00 x 156.00 | 6.88 | Fig2_trajectories.eps | Fig2_trajectories.png | figures/submission/Fig2.eps, figures/submission/Fig2.tif |
| Fig3_effects | 174.00 x 228.00 | 5.72 | Fig3_effects.eps | Fig3_effects.png | figures/submission/Fig3.eps, figures/submission/Fig3.tif |
| Fig4_composition | 174.00 x 214.00 | 4.56 | Fig4_composition.eps | Fig4_composition.png | figures/submission/Fig4.eps, figures/submission/Fig4.tif |
| Fig5_modules | 174.00 x 214.00 | 6.28 | Fig5_modules.eps | Fig5_modules.png | figures/submission/Fig5.eps, figures/submission/Fig5.tif |
| Fig6_public | 174.00 x 220.00 | 6.47 | Fig6_public.eps | Fig6_public.png | figures/submission/Fig6.eps, figures/submission/Fig6.tif |
| Fig7_bias | 174.00 x 220.00 | 6.19 | Fig7_bias.eps | Fig7_bias.png | figures/submission/Fig7.eps, figures/submission/Fig7.tif |
| Fig8_mitigation | 174.00 x 168.00 | 4.96 | Fig8_mitigation.eps | Fig8_mitigation.png | figures/submission/Fig8.eps, figures/submission/Fig8.tif |
| FigS1_effects_secondary | 174.00 x 188.00 | 3.94 | FigS1_effects_secondary.eps | FigS1_effects_secondary.png | figures/submission/FigS1.pdf |
| FigS2_bias_secondary | 174.00 x 220.00 | 5.37 | FigS2_bias_secondary.eps | FigS2_bias_secondary.png | figures/submission/FigS2.pdf |
| FigS3_mitigation_secondary | 174.00 x 168.00 | 4.12 | FigS3_mitigation_secondary.eps | FigS3_mitigation_secondary.png | figures/submission/FigS3.pdf |

검토 묶음: `figures/manuscript/all_figures_review.pdf`.

## 2. 글꼴 (사용 글꼴, 설치 방법, pdffonts 결과, 가장 작은 글자 크기)

Arial은 없다. Liberation Sans는 패키지 `fonts-liberation` 1:1.07.4-11로 `/usr/share/fonts/truetype/liberation/`에 있다. 스크립트가 `fontManager.addfont()`로 Regular, Bold, Italic, BoldItalic TTF를 등록했다. `mathtext.fontset`은 custom이고 rm, it, bf, cal, sf, tt는 Liberation Sans다.

그림 PDF 11개와 `all_figures_review.pdf`의 `pdffonts`는 모두 emb yes다. EPS 11개를 `ps2pdf -dEPSCrop`로 바꾼 뒤의 `pdffonts`도 모두 emb yes다. 글꼴 이름은 LiberationSans, LiberationSans-Bold, LiberationSans-Italic이다.

가장 작은 `get_fontsize()`는 그림 11개 모두 8이다.

## 3. 점검 결과 (§5의 1–5 표)

겹침 네 가지와 `alpha`가 1이 아닌 artist 수는 0이다. 같은 축의 눈금 글자끼리의 겹침은 세지 않았다.

| stem | text_overlap | text_marker | legend_data | outside | min_fontsize | alpha_not_1 |
| --- | --- | --- | --- | --- | --- | --- |
| Fig1_design | 0 | 0 | 0 | 0 | 8.0 | 0 |
| Fig2_trajectories | 0 | 0 | 0 | 0 | 8.0 | 0 |
| Fig3_effects | 0 | 0 | 0 | 0 | 8.0 | 0 |
| Fig4_composition | 0 | 0 | 0 | 0 | 8.0 | 0 |
| Fig5_modules | 0 | 0 | 0 | 0 | 8.0 | 0 |
| Fig6_public | 0 | 0 | 0 | 0 | 8.0 | 0 |
| Fig7_bias | 0 | 0 | 0 | 0 | 8.0 | 0 |
| Fig8_mitigation | 0 | 0 | 0 | 0 | 8.0 | 0 |
| FigS1_effects_secondary | 0 | 0 | 0 | 0 | 8.0 | 0 |
| FigS2_bias_secondary | 0 | 0 | 0 | 0 | 8.0 | 0 |
| FigS3_mitigation_secondary | 0 | 0 | 0 | 0 | 8.0 | 0 |

쪽 크기와 TIFF MB는 §1 표와 같다.

## 4. 확인용 값 (원고 수치 대조용, 저장값 그대로)

### 4.1 Fig. 3c: 시계 × 자료별 peak_h, peak_center, peak_halfwidth_95, LRT p, φ̄

`results/stage8/effect_summary.tsv`의 kind cosinor와 `prereg/frozen/F1.tsv`의 GSE39445 행이다. LRT p는 effect_summary의 `p`다.

| clock | dataset | peak_h | peak_center | peak_halfwidth_95 | lrt_p | phi_bar |
| --- | --- | --- | --- | --- | --- | --- |
| TAms_chrono | GSE39445 | 14.2335397784943 | 14.2186064321564 | 0.820399743143368 | 5.49745798912395e-20 | 12.616768998085 |
| TAms_chrono | GSE56931 | 11.4579688514277 | 11.2197113952556 | 7.28242031911783 | 0.181165188487945 | 12.616768998085 |
| TAms_chrono | GSE113883 | 12.188050742382 | 12.2120885405226 | 1.30484190182366 | 2.96936394371969e-10 | 12.616768998085 |
| TAms_chrono | GSE48113 | 10.879869555482 | 10.855130663467 | 1.31149483591321 | 2.70474257565244e-08 | 12.616768998085 |
| TAms_mort | GSE39445 | 13.4288594848278 | 13.4014342140422 | 1.16438541795076 | 3.45851149483147e-15 | 11.8599407935094 |
| TAms_mort | GSE56931 | 12.0688456403715 | 12.1544296381706 | 3.06206316312207 | 0.00673707816031792 | 11.8599407935094 |
| TAms_mort | GSE113883 | 10.0600990696032 | 10.0668531335917 | 2.22504509606368 | 4.44466652817908e-06 | 11.8599407935094 |
| TAms_mort | GSE48113 | 11.4858304321163 | 11.4467106578545 | 2.03534016526948 | 2.67642593880552e-05 | 11.8599407935094 |
| RA_default | GSE39445 | 12.0064473932843 | 12.4584652559485 | 8.83488247674881 | 0.571779550797631 | 11.2423863477089 |
| RA_default | GSE56931 | 9.88701727265618 | 9.77493370481417 | 2.411697769131 | 0.0187399469939367 | 11.2423863477089 |
| RA_default | GSE113883 | 11.812631972311 | 11.8212982069481 | 1.44667958366095 | 3.24675832907807e-05 | 11.2423863477089 |
| RA_default | GSE48113 | 11.5918483633937 | 11.5611538419167 | 2.79211694956184 | 0.0125285842066838 | 11.2423863477089 |

### 4.2 Fig. 6a: x에 쓴 GSE39445 시각 변수 이름과 정의

x는 `results/stage8/gene_peak_differences_arntl_alias.tsv`의 `peak_gse39445`다. 그 값은 `results/posthoc_stage7c/gse39445_gene_peaks.tsv`의 `peak_h`와 같고, 그 파일의 `time_axis`는 hours_awake이다. 축 이름은 "Peak, GSE39445 (h since wake)"다.

### 4.3 Fig. 6c: 시계 × 가림 방식 × 처리별 LRT p

kNN 파일은 `results/stage7b/d1_development_model_m.tsv`다. 참가자 평균 파일은 `results/stage7c/d11_development_model_m.tsv`다. scheme이 unmasked일 때 두 파일의 LRT p는 세 시계 모두 같다. 그림의 Unmasked는 한 칸이다. random과 low_expression의 kNN은 7b 파일, 참가자 평균은 7c 파일이다.

가림 수는 `results/stage7b/d1_development_audit.tsv`에 있다.

| scheme | n_masked | n_finite_before_mask | n_genes_dropped | n_genes_kept | n_values_imputed |
| --- | --- | --- | --- | --- | --- |
| random | 863712 | 4318561 | 0 | 19541 | 863712 |
| low_expression | 863712 | 4318561 | 0 | 19541 | 863712 |

| scheme | clock | lrt_p_knn_file | lrt_p_participant_mean_file |
| --- | --- | --- | --- |
| unmasked | RA_default | 0.571779550797631 | 0.571779550797631 |
| unmasked | RA_Peters | 0.0822398090507771 | 0.0822398090507771 |
| unmasked | RA_GTEx | 0.643663244864561 | 0.643663244864561 |
| random | RA_default | 0.131795197933585 | 0.493775609974818 |
| random | RA_Peters | 0.396462971000013 | 0.233868759239569 |
| random | RA_GTEx | 0.638021860434852 | 0.651673618269485 |
| low_expression | RA_default | 0.0312868740107463 | 0.0125311384531322 |
| low_expression | RA_Peters | 0.175697860237701 | 0.0962455882817528 |
| low_expression | RA_GTEx | 0.223295276970678 | 0.290132197293332 |

### 4.4 Fig. 8a와 S3a: 채움에 쓴 파일과 열 이름

파일은 `results/stage9/composition_cosinor.tsv`다. 왼쪽 점은 `lrt_p_before`, 오른쪽 점은 `lrt_p_after`다. 각 열이 0.05 미만이면 채운 표식이다.

### 4.5 Fig. S2d: RNAAgeCalc-Peters와 RNAAgeCalc-GTEx Laboratory의 s = 4 |b|/총 SD

파일은 `results/stage9/sim_design1_bias.tsv`다. 열 `b_over_total_sd`는 부호가 있는 b / 총 SD다.

| clock | scenario | s | b | b_over_total_sd | unit |
| --- | --- | --- | --- | --- | --- |
| RA_Peters | lab | 4 | -3.12267692924951 | -0.187553096932267 | year |
| RA_GTEx | lab | 4 | 5.68776312390269 | 0.18804909857676 | year |

## 5. 판단이 필요한 사항

- U+0304(결합 매크론)은 Liberation Sans에 없다. φ̄는 mathtext `$\bar{\phi}$`로 그렸다.
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
