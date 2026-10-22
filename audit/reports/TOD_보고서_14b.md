# TOD 보고서 14b — 그림 용어 통일과 보충 자료 묶음

## 1. 바꾼 그림과 "reregist" 검색 결과

검색한 파일은 `scripts/62_part11_figures.py`와 `scripts/60_part10_figures.py`다. 바꾼 곳은 다섯 줄이다.

| 파일 | 줄 | 바꾸기 전 | 바꾸기 후 |
| --- | --- | --- | --- |
| scripts/62_part11_figures.py | 368 | Preregistration and locked confirmation | Prespecified plan and locked confirmation |
| scripts/62_part11_figures.py | 656 | Preregistered scoring | Prespecified scoring |
| scripts/62_part11_figures.py | 670 | Preregistered / scoring | Prespecified / scoring |
| scripts/62_part11_figures.py | 789 | Preregistered scoring | Prespecified scoring |
| scripts/60_part10_figures.py | 216 | Preregistration and locked confirmation | Prespecified plan and locked confirmation |

바꾼 뒤 두 파일에서 `reregist` 검색 수는 0이다. `pdftotext`에서 `prespecified` 수는 Fig. 1이 1, Fig. 5가 2, Fig. 6이 1이고, `reregist` 수는 0이다.

다시 쓴 그림 파일은 Fig. 1, Fig. 5, Fig. 6이다. PDF 쪽 크기는 `pdfinfo` Page size를 pt × 25.4 / 72로 둔 값이다. TIFF MB는 바이트 수를 1e6으로 나눈 값이다.

| stem | pdf_mm | tiff_bytes | tiff_MB |
| --- | --- | --- | --- |
| Fig1_design | 174.00 x 120.00 | 5053190 | 5.05 |
| Fig5_modules | 174.00 x 214.00 | 6300304 | 6.30 |
| Fig6_public | 174.00 x 220.00 | 5887900 | 5.89 |

`figures/submission/`의 Fig1, Fig5, Fig6 사본과 `figures/manuscript/all_figures_review.pdf`를 다시 썼다.

## 2. 그림 점검 (지시서 12 §3, 13 §2 표)

ylabel_gap_mm는 y축 이름과 y 눈금 글자의 왼쪽 끝이 왼쪽 이웃 축 테두리에서 떨어진 최소 거리(mm)다. 왼쪽 이웃 축이 없으면 없음이다.

| stem | pdf_clip | n_text | outside | text_overlap | leader_text | data_outside | intentional_clip | min_fontsize | alpha_not_1 | leader_cross | leader_foreign | axes_text | ylabel_gap_mm |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Fig1_design | 0 | 15 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 |
| Fig2_trajectories | 0 | 100 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 2.58 |
| Fig3_effects | 0 | 67 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 |
| Fig4_composition | 0 | 43 | 0 | 0 | 0 | 0 | 5 | 8.0 | 0 | 0 | 0 | 0 | 없음 |
| Fig5_modules | 0 | 143 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 |
| Fig6_public | 0 | 86 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 10.08 |
| Fig7_bias | 0 | 98 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 9.42 |
| Fig8_mitigation | 0 | 96 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 5.60 |
| FigS1_effects_secondary | 0 | 50 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 없음 |
| FigS2_bias_secondary | 0 | 95 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 9.42 |
| FigS3_mitigation_secondary | 0 | 65 | 0 | 0 | 0 | 0 | 0 | 8.0 | 0 | 0 | 0 | 0 | 21.39 |

Fig. 1, Fig. 5, Fig. 6의 PDF 글꼴은 LiberationSans와 LiberationSans-Bold이고 `pdffonts` emb는 yes다.

## 3. ESM 파일 (경로, 크기, 쪽·시트 구성, 표별 출처 파일)

| 파일 | 바이트 | MB |
| --- | --- | --- |
| figures/submission/esm/ESM_1.pdf | 197050 | 0.20 |
| figures/submission/esm/ESM_2.xlsx | 221449 | 0.22 |

ESM_1.pdf는 4쪽이다. 쪽 크기는 pt × 25.4 / 72이다.

| 쪽 | mm | 내용 |
| --- | --- | --- |
| 1 | 174.00 x 160.00 | Cover |
| 2 | 174.00 x 199.99 | Fig. S1 |
| 3 | 174.00 x 231.99 | Fig. S2 |
| 4 | 174.00 x 176.11 | Fig. S3 |

그림 PDF는 `figures/manuscript/FigS1_effects_secondary.pdf`, `FigS2_bias_secondary.pdf`, `FigS3_mitigation_secondary.pdf`이다. 설명문은 그 쪽 아래에 있다.

ESM_2.xlsx의 행·열 수는 `openpyxl`의 max_row, max_column이다. 데이터 행은 4행부터다.

| 시트 | 행 | 열 | 데이터 행 | 1행 |
| --- | --- | --- | --- | --- |
| Cover | 4 | 1 | 4 | Title: Time of blood sampling biases transcriptomic aging clocks: development, prespecified confirmation, and practical guidance |
| S1 | 15 | 6 | 12 | Table S1. Deviations D1–D12: content, timing relative to results, and reason. |
| S2 | 40 | 36 | 37 | Table S2. Sensitivity analyses: prespecified (flagged-sample exclusion, random slopes, batch-corrected counts, comprehensive annotation) and deviation-defined (complete-case, kNN and participant-mean versions; imputed composition). |
| S3a | 58 | 57 | 55 | Table S3a. GSE122725 time-axis correction and exploratory H1 and H4 |
| S3b | 35 | 18 | 32 | Table S3b. Module outcomes under participant-mean scoring |
| S3c | 158 | 16 | 155 | Table S3c. Masking experiment |
| S3d | 18 | 17 | 15 | Table S3d. Simulation parameters |
| S3e | 1123 | 34 | 1120 | Table S3e. Simulation results, designs 1–3 |
| S3f | 28 | 8 | 25 | Table S3f. Mitigation |
| S4 | 31 | 25 | 28 | Table S4. Frozen development parameters: peak times, amplitudes, variance components, expected morning−evening effects and module outcomes. |
| S5 | 110 | 13 | 107 | Table S5. Software versions, model-file checksums and the frozen manifest. |
| S6 | 33 | 14 | 30 | Table S6. Prospective power of the confirmatory tests. |

Cover의 2–4행은 Journal: GeroScience, Authors: [Author names], Corresponding author: [Affiliation], [e-mail]이다.

| 시트 | 출처 | 행 |
| --- | --- | --- |
| S1 | results/stage8/tables/TableS1.tsv | 12 |
| S2 | results/stage8/tables/TableS2.tsv | 37 |
| S3a | results/stage8/tables/TableS3.tsv block l2_sleep_start | 3 |
| S3a | results/stage8/tables/TableS3.tsv block arntl_alias | 1 |
| S3a | results/stage8/tables/TableS3.tsv block l2_h1 | 13 |
| S3a | results/stage8/tables/TableS3.tsv block l2_h4 | 3 |
| S3a | results/stage8/tables/TableS3.tsv block effects | 30 |
| S3a | results/stage8/tables/TableS3.tsv block l3_prediction | 5 |
| S3b | results/stage8/tables/TableS3.tsv block d11_modules | 8 |
| S3b | results/stage8/tables/TableS3.tsv block f3_features | 24 |
| S3c | results/stage7b/d1_development_audit.tsv | 2 |
| S3c | results/stage7b/d1_development_model_m.tsv | 9 |
| S3c | results/stage7b/d1_development_within_pearson.tsv | 144 |
| S3d | results/stage9/sim_parameters.tsv | 15 |
| S3e | results/stage9/sim_design1_bias.tsv | 120 |
| S3e | results/stage9/sim_design1_fpr.tsv | 480 |
| S3e | results/stage9/sim_design1_smallest_s.tsv | 15 |
| S3e | results/stage9/sim_design1_montecarlo.tsv | 5 |
| S3e | results/stage9/sim_design2_fpr.tsv | 480 |
| S3e | results/stage9/sim_design3_power.tsv | 20 |
| S3f | results/stage9/mitigation_summary.tsv | 25 |
| S4 | prereg/frozen/F1.tsv | 15 |
| S4 | prereg/frozen/F2.tsv | 5 |
| S4 | prereg/frozen/F3.tsv | 8 |
| S5 | prereg/frozen/environment.tsv | 17 |
| S5 | results/stage5/code_vs_zenodo_md5.tsv | 64 |
| S5 | prereg/manifest.tsv | 26 |
| S6 | results/stage5/locked_design_power.tsv | 30 |

S5의 64행 중 md5가 있는 행은 4이고 출처는 `results/stage4/module_clocks/zenodo_v2_files.tsv`다. sha256이 있는 행은 60이고 출처는 `results/stage4/module_clocks/module_archive_members.tsv`다. md5와 sha256이 둘 다 있는 행은 0이다.

`(h)`를 붙인 열 이름은 peak (h), circ_mean (h), circ_sd (h), peak_center (h), peak_halfwidth_95 (h), phi_bar (h), max_circ_dist (h), peak_GSE39445 (h), peak_GSE56931 (h), peak_GSE113883 (h), s (h), smallest_s (h), W (h), phi (h)다.

## 4. ESM 점검 (§4의 1–5)

1. `pdffonts figures/submission/esm/ESM_1.pdf`의 글꼴은 LiberationSans와 LiberationSans-Bold이고 emb는 yes다.
2. 쪽 수는 4다. 1쪽은 Cover, 2쪽은 Fig. S1, 3쪽은 Fig. S2, 4쪽은 Fig. S3다.
3. 시트 이름, 행·열 수, 1행 제목은 §3의 표와 같다.
4. 검색 수.

| 파일 | reregist | TAms_ | RA_ | CMP_ | L1 | L2 | L3 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ESM_1.pdf | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ESM_2.xlsx | 0 | 0 | 0 | 0 | 1 | 0 | 0 |

ESM_2.xlsx의 L1 1건은 시트 S3a 7행의 저장 문장 `GSE122725 symbol BMAL1; GSE39445 symbol ARNTL`이다.

5. 모듈 구성 유전자와 계수. pkl 파일은 열지 않았다. 유전자 목록 파일은 넣지 않았다. F1, F2, F3의 열은 유전자 이름 열이나 계수 열이 아니다. ESM_2의 열 이름 중 `gene`이 들어간 이름은 S3a의 n_genes, S3c의 n_genes_dropped, n_genes_kept, S3e의 generator다. `allmodulegenes` 문자열 수는 0이다. S5의 체크섬 열은 filename, md5, sha256, code_bytes, zenodo_bytes, md5_equal, class, used_for_scoring다.

## 5. 지시서와 다르게 한 것 (빠짐없이)

- `scripts/60_part10_figures.py` 216행도 바꿨다. 이 파일은 이번 그림 PDF를 쓰지 않는다.
- Table S3는 S3a–S3f로 나눴다. S3a에는 block effects 30행과 l3_prediction 5행도 있다. S3a의 1행은 지시서의 S3a 문장이다.
- S3b, S3c, S3d, S3e, S3f의 1행 문장은 §3에 적은 문장이다.
- 한 시트에 표를 여러 개 넣을 때 stored_file 열을 앞에 더했다.
- 값이 모두 빈 열은 뺐다.
- 셀 값은 문자열이다.
- `preregistration`은 `prespecified plan`, `preregistered`는 `prespecified`로 바꿨다.
- 셀 전체가 `frozen`이면 `prespecified scoring`, 셀 전체가 `D11` 또는 `D11_participant_mean`이면 `participant-mean imputation`으로 바꿨다.
- `L1`, `L2`, `L3`는 앞뒤가 영숫자가 아닐 때만 `GSE48113`, `GSE122725`, `GSE295758`로 바꿨다.
- BMAL1은 바꾸지 않았다.
- 체크섬 다이제스트는 다시 계산하지 않았다.
- `results/stage9/composition_cosinor.tsv`, `composition_l3.tsv`, `composition_l3_betas.tsv`, `age_correlation_l3.tsv`는 넣지 않았다.
- ESM_1의 쪽 너비는 174.00 mm다. 보충 그림 PDF는 배율을 바꾸지 않고 설명문 위에 두었다.
- Fig. 2, Fig. 3, Fig. 4, Fig. 7, Fig. 8, Fig. S1, Fig. S2, Fig. S3 파일은 다시 쓰지 않았다.

## 6. 커밋 해시

이 보고서를 넣는 커밋
