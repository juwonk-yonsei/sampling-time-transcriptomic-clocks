# TOD 보고서 5 — 판단 사항 처리, tAge 재채점, 복합 시계 모듈 분해, 사전등록 준비 계산

## 1. 요약 (수치만, 16줄 이내)
- 리비전 코드: `array_to_diff`, `clocks.composite`, `clocks.modules`, `scripts/10_score_human.py`. 사람 자료 복합 시계 2개 파일의 md5는 Zenodo 사본과 같다. 사람 자료 모듈은 multispecies mortality 14개. 채점은 0이 아닌 선택 가중치와 scaled-diff의 내적이다.
- 다종 복합 시계 2개: 특징 10487, 설치류 23개 모듈과 겹침 1292, 밖 9195. 논문 값과 같다. `TAms_chrono`, `TAms_mort`로 넘겼다. 설치류 복합 시계 특징 수는 16517, 16760이라 넘기지 않았다.
- 보고서 4 점수와 짝이 있는 재채점의 개인 내 Pearson 최솟값 0.3525667717160795.
- 성분 합 + 절편과 `predict`의 최대 절대차 3.9968028886505635e-15.
- H1 TAms: GSE39445 TAms_chrono R2 0.362193176953973 p 5.49745798912395e-20 peak 14.2335397784943; GSE39445 TAms_mort R2 0.28699925739554 p 3.45851149483147e-15 peak 13.4288594848278; GSE56931 TAms_chrono R2 0.0286364463384637 p 0.181165188487945 peak 11.4579688514277; GSE56931 TAms_mort R2 0.0821416660164724 p 0.00673707816031792 peak 12.0688456403715; GSE113883 TAms_chrono R2 0.265399481367344 p 2.96936394371969e-10 peak 12.188050742382; GSE113883 TAms_mort R2 0.158916853730887 p 4.44466652817908e-06 peak 10.0600990696032.
- BH q < 0.05: TM GSE39445 7/14; GSE56931 2/14; GSE113883 9/14. CMP GSE39445 TAms_chrono 17/24; GSE56931 TAms_chrono 11/24; GSE113883 TAms_chrono 19/24; GSE39445 TAms_mort 8/24; GSE56931 TAms_mort 13/24; GSE113883 TAms_mort 16/24.
- H2 p_med: GSE39445 TAms_chrono 0.604554553943655 [0.360490943625388, 0.881407966272375]; GSE39445 TAms_mort 0.648619163577675 [0.328462743626483, 1.00587179922653]; GSE56931 TAms_chrono 1.23230529037825 [-2.39490072930965, 4.49887335706176]; GSE56931 TAms_mort 0.937587420102395 [0.0854764204114879, 2.35128746858658]; GSE113883 TAms_chrono 0.827434926978766 [0.449125356663948, 1.15576326319054]; GSE113883 TAms_mort 0.668143587639813 [0.433397727889269, 1.11458876791915].
- φ̄, 최대 원형 거리: RA_GTEx 11.518486825228 max 5.35019887986941; RA_Peters 7.56604107260516 max 4.27017978445731; RA_default 11.2423863477089 max 1.35536907505411; TAms_chrono 12.616768998085 max 1.61677078040933; TAms_mort 11.8599407935094 max 1.79984172390618.
- 재현 리듬 모듈 8/62: TM_turquoise_mortality, CMP_TAms_chrono_turquoise, CMP_TAms_chrono_yellow, CMP_TAms_chrono_orange, CMP_TAms_chrono_brown, CMP_TAms_chrono_red, CMP_TAms_mort_turquoise, CMP_TAms_mort_darkred.
- GSE295758 d_exp 중앙값, 검정력: RA_GTEx 0.0799546478017167 power0.05 0.0919710023400413 power0.01 0.0243768273749048; RA_Peters 0.286305297017999 power0.05 0.572873321105643 power0.01 0.323267563613161; RA_default 0.412448894924069 power0.05 0.870300205086949 power0.01 0.682114842365414; TAms_chrono 0.208443813199012 power0.05 0.345077590733327 power0.01 0.150545763753431; TAms_mort 0.291250817026961 power0.05 0.587339593367421 power0.01 0.336502748601378.
- 잠금 설계 검정력 범위: GSE48113 RA_GTEx p<0.05 0.0935–0.579 p<0.01 0.027–0.3245; GSE48113 RA_Peters p<0.05 0.361–0.998 p<0.01 0.1545–0.968; GSE48113 RA_default p<0.05 0.112–0.9795 p<0.01 0.0375–0.931; GSE48113 TAms_chrono p<0.05 0.431–1 p<0.01 0.201–1; GSE48113 TAms_mort p<0.05 0.8865–1 p<0.01 0.7175–1; GSE122725 RA_GTEx p<0.05 0.0755–0.3395 p<0.01 0.0225–0.1505; GSE122725 RA_Peters p<0.05 0.1635–0.83 p<0.01 0.059–0.6185; GSE122725 RA_default p<0.05 0.0975–0.842 p<0.01 0.03–0.644; GSE122725 TAms_chrono p<0.05 0.237–1 p<0.01 0.092–0.9985; GSE122725 TAms_mort p<0.05 0.6175–0.9975 p<0.01 0.3815–0.988.
- 마지막 분석 커밋 037236f77ac8d4d7efed92761f375e68caee32cf.

## 2. 환경과 자원 (바뀐 것, /usr/bin/time 기록)
- 시계를 읽는 Python은 `/HDD8T1/WORK/GeroSci/tod4_clock_env`다. R은 conda `r_env`다.
- Part B 재채점: elapsed 2:56.93, maximum resident set size 644760 kB.
- Part D cosinor: elapsed 26:17.35, maximum resident set size 302776 kB.
- Part E 시뮬레이션: elapsed 29:00.35, maximum resident set size 237092 kB.
- 세 기록 모두 40 GB 아래다. `ulimit`은 두지 않았다.

## 3. 결정 기록 (config/decisions.tsv에 이번에 더한 행)
| date | id | content | reason | timing |
| --- | --- | --- | --- | --- |
| 2026-10-16 | 5.1.1 | Zenodo record 22166800 is approved as version 2. | Instruction 5 section 1.1 item 1. Report 4 section 10. The metadata version field is empty. Publication order identifies 22166800, published 2026-08-29, as version 2. | prior |
| 2026-10-16 | 5.1.2 | The downloaded clock-file scope is approved. The two multispecies composite files already matched md5 and are not downloaded again. | Instruction 5 section 1.1 item 2. Report 4 section 10. Files: EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl and EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl. | prior |
| 2026-10-16 | 5.1.3 | The module-archive filename decision is approved. | Instruction 5 section 1.1 item 3. Report 4 section 10. Local names differed from the Zenodo keys, so the official names were downloaded and the md5 matched. | prior |
| 2026-10-16 | 5.1.4 | The scoring path changes in Part B. Preprocessing uses the revision-code function array_to_diff. | Instruction 5 section 1.1 item 4. Report 4 section 10 item 4 is replaced by this Part B. Instruction 3 section 7.1 omitted per-sample scaling. | prior |
| 2026-10-16 | 5.1.5 | The rodent composite feature counts disagreed with the paper value 10487. That was a breach of rule 1.4. Those clocks are sensitivity only. | Instruction 5 section 1.1 item 5. Report 4 recorded Chronoage 16517 and Mortality 16760 features. Part B names the rescored rodent composites TArd_chrono and TArd_mort and does not pass them as the primary clocks. | 사후 |
| 2026-10-16 | 5.1.6 | The GSE113883 QC scale is approved. Sample medians use log2(TPM+1) on genes with mean TPM at least 1. | Instruction 5 section 1.1 item 6. Report 4 section 10. | prior |
| 2026-10-16 | 5.1.7 | Every GSE113883 subject woke at 06:00. T_wake is the same model as T_clock and adds no information. The record is kept. | Instruction 5 section 1.1 item 7. Report 4 section 10. clock_time_h minus hours_since_wake is 6 for 87 samples and -18 for 66 samples. The difference modulo 24 is 6 for all 153 samples. | prior |
| 2026-10-16 | 5.1.8 | Report 4 section 11 items 1 to 3 are approved. Item 4 is replaced by Part B. | Instruction 5 section 1.1. Item 1: clock files outside the placed set were not downloaded. Item 2: requirements.txt packages that the clock loader does not import were not installed. Item 3: the stage-4 scoring script was not wrapped in /usr/bin/time -v. | prior |

## 4. tAge 재채점
### 4.1 리비전 코드 확인 (함수, 규칙, 시계 파일, 모듈 목록, 방식)
| item | value |
| --- | --- |
| human_script | clock-specificity/scripts/10_score_human.py |
| human_composite_files | EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl;EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl |
| human_modules | multispecies mortality modules returned by clocks.modules |
| human_module_count | 14 |
| scoring_method | dot product of effective weights of selected nonzero genes with scaled-diff columns present in every sample |
| intercept_in_weight_score | false |
| module_score_definition | weights of that module clock, not a split of the composite coefficient vector |
| predict_also_used_by_human_script | false |
| gene_filter | drop genes whose mean is below the 0.25 quantile of gene means |
| standardization_gene_set | genes remaining after the quantile filter and the ortholog map, before alignment to the clock gene list |
| standardization | z-score each sample across those genes, ddof=1 |
| reference | subtract the per-gene median of the reference samples after alignment |
| log_transform | none inside array_to_diff; the input is already on a log scale |
| ortholog | human Entrez to mouse Entrez through the ortholog map; unmapped rows dropped; duplicate mouse ids keep the first row; values are not summed |
| wrapper | scripts/31_partB_score.py symbols_to_human_entrez |
| gse113883_input | log2(TPM+1) because counts are not stored; counts_to_diff was not used |
| code_files_md5 | equal to the stage-4 Zenodo copies |

`scripts/10_score_human.py`의 목록은 multispecies mortality 복합 시계, multispecies chronological 복합 시계, `clocks.modules("multispecies", "mortality")` 14개다. 그 스크립트의 점수는 `one_clock`의 가중치 내적이고, 절편은 넣지 않는다. `predict`는 그 스크립트가 호출하지 않는다.

코드가 가리키는 시계 파일 64개의 md5는 `results/stage5/code_vs_zenodo_md5.tsv`에서 Zenodo 22166800 사본과 모두 같다. 논문 값 확인은 그 파일로 했다.

### 4.2 복합 시계 확인 표 (논문 값과 나란히)
| clock_set | outcome | n_features | n_rodent_modules | n_rodent_module_genes | n_modules_with_shared_gene | n_in_rodent_modules | n_outside_rodent_modules | paper_n_features | paper_n_in_module | paper_n_outside | all_three_match | passed_to_next_stage |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| multispecies | chronological | 10487 | 23 | 2007 | 0 | 1292 | 9195 | 10487 | 1292 | 9195 | True | True |
| multispecies | mortality | 10487 | 23 | 2007 | 0 | 1292 | 9195 | 10487 | 1292 | 9195 | True | True |
| rodent | chronological | 16517 | 23 | 2007 | 0 | 2007 | 14510 | 10487 | 1292 | 9195 | False | False |
| rodent | mortality | 16760 | 23 | 2007 | 0 | 2007 | 14753 | 10487 | 1292 | 9195 | False | False |

### 4.3 입력 준비 (래퍼 경로, 자료별 유전자 수, 기준 샘플 수)
래퍼는 `scripts/31_partB_score.py`의 `symbols_to_human_entrez`다. 기호를 사람 Entrez로 바꾼 뒤 `array_to_diff`를 수정 없이 호출한다. 같은 사람 Entrez가 여러 행이면 첫 행을 남기고, 값은 합하지 않는다.

GSE39445와 GSE56931의 입력은 1단계 log2 행렬이다. 참조는 그 자료의 대조 표본 전체다. GSE113883은 계수 행렬이 없어 `counts_to_diff`를 쓰지 않고, `log2(TPM+1)`에 `array_to_diff`를 적용했다. 참조는 153개 표본 전체이고, 전부 대조다.

| dataset | input_scale | reference | n_samples | n_reference | n_input_rows | n_symbol_mapped | n_unique_human_entrez | n_after_quantile | n_with_mouse_ortholog | n_unique_mouse_entrez | n_clock_gene_list | n_clock_genes_finite_in_all_samples |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | stage1_log2 | all_control_samples | 438 | 221 | 19541 | 16185 | 16185 | 12139 | 10775 | 10775 | 18696 | 10698 |
| GSE56931 | stage1_log2 | all_control_samples | 249 | 130 | 20755 | 18432 | 18432 | 13824 | 12054 | 12054 | 18696 | 11955 |
| GSE113883 | log2_TPM_plus_1 | all_153_samples_all_control | 153 | 153 | 35581 | 32238 | 32238 | 32238 | 16063 | 16063 | 18696 | 15681 |

표준화는 분위 필터와 ortholog 대응 뒤에 남은 유전자에 대해 표본마다 z-score(ddof=1)를 낸 다음, 시계 유전자 목록에 맞춘다. 그 다음 참조 표본의 유전자별 중앙값을 뺀다.

### 4.4 재채점과 보고서 4 점수의 개인 내 상관
상관은 대조 표본에서 참가자 평균을 뺀 뒤의 Pearson이다. `TAms_*`와 사람 자료 모듈 14개에는 보고서 4 점수가 없다.

| dataset | clock | stage4_clock | n_control | within_pearson | paired |
| --- | --- | --- | --- | --- | --- |
| GSE113883 | TAms_chrono |  | 153 |  | False |
| GSE113883 | TAms_mort |  | 153 |  | False |
| GSE113883 | TArd_chrono | TA_chrono | 153 | 0.9273743967986662 | True |
| GSE113883 | TArd_mort | TA_mort | 153 | 0.9648995200183595 | True |
| GSE113883 | TM_blue_mortality |  | 153 |  | False |
| GSE113883 | TM_brown4_mortality |  | 153 |  | False |
| GSE113883 | TM_darkgreen_mortality |  | 153 |  | False |
| GSE113883 | TM_darkmagenta_mortality |  | 153 |  | False |
| GSE113883 | TM_darkred_mortality |  | 153 |  | False |
| GSE113883 | TM_darkslateblue_mortality |  | 153 |  | False |
| GSE113883 | TM_green_mortality |  | 153 |  | False |
| GSE113883 | TM_ivory_mortality |  | 153 |  | False |
| GSE113883 | TM_orange_mortality |  | 153 |  | False |
| GSE113883 | TM_pink_mortality |  | 153 |  | False |
| GSE113883 | TM_plum1_mortality |  | 153 |  | False |
| GSE113883 | TM_sienna3_mortality |  | 153 |  | False |
| GSE113883 | TM_turquoise_mortality |  | 153 |  | False |
| GSE113883 | TM_white_mortality |  | 153 |  | False |
| GSE113883 | TMrd_black_chronological | TM_black_chronological | 153 | 0.9800600518454348 | True |
| GSE113883 | TMrd_black_mortality | TM_black_mortality | 153 | 0.9880241565526882 | True |
| GSE113883 | TMrd_blue_chronological | TM_blue_chronological | 153 | 0.9743576290533579 | True |
| GSE113883 | TMrd_blue_mortality | TM_blue_mortality | 153 | 0.99692298168825 | True |
| GSE113883 | TMrd_brown_chronological | TM_brown_chronological | 153 | 0.9942616017174875 | True |
| GSE113883 | TMrd_brown_mortality | TM_brown_mortality | 153 | 0.9959842938618159 | True |
| GSE113883 | TMrd_cyan_chronological | TM_cyan_chronological | 153 | 0.9934402760793947 | True |
| GSE113883 | TMrd_cyan_mortality | TM_cyan_mortality | 153 | 0.8933738227908817 | True |
| GSE113883 | TMrd_darkgreen_chronological | TM_darkgreen_chronological | 153 | 0.985450267290755 | True |
| GSE113883 | TMrd_darkgreen_mortality | TM_darkgreen_mortality | 153 | 0.9781804600952898 | True |
| GSE113883 | TMrd_darkgrey_chronological | TM_darkgrey_chronological | 153 | 0.9969030727147566 | True |
| GSE113883 | TMrd_darkgrey_mortality | TM_darkgrey_mortality | 153 | 0.9953310537614303 | True |
| GSE113883 | TMrd_darkorange_chronological | TM_darkorange_chronological | 153 | 0.9974253496591118 | True |
| GSE113883 | TMrd_darkorange_mortality | TM_darkorange_mortality | 153 | 0.9797011079915096 | True |
| GSE113883 | TMrd_darkred_chronological | TM_darkred_chronological | 153 | 0.9130604113958338 | True |
| GSE113883 | TMrd_darkred_mortality | TM_darkred_mortality | 153 | 0.9510122634068934 | True |
| GSE113883 | TMrd_green_chronological | TM_green_chronological | 153 | 0.9967601613881646 | True |
| GSE113883 | TMrd_green_mortality | TM_green_mortality | 153 | 0.9962446960459188 | True |
| GSE113883 | TMrd_greenyellow_chronological | TM_greenyellow_chronological | 153 | 0.9933709156128822 | True |
| GSE113883 | TMrd_greenyellow_mortality | TM_greenyellow_mortality | 153 | 0.9986120139141015 | True |
| GSE113883 | TMrd_grey60_chronological | TM_grey60_chronological | 153 | 0.9919142063350096 | True |
| GSE113883 | TMrd_grey60_mortality | TM_grey60_mortality | 153 | 0.9578102854157081 | True |
| GSE113883 | TMrd_lightgreen_chronological | TM_lightgreen_chronological | 153 | 0.8011096102752467 | True |
| GSE113883 | TMrd_lightgreen_mortality | TM_lightgreen_mortality | 153 | 0.7692852863271314 | True |
| GSE113883 | TMrd_lightyellow_chronological | TM_lightyellow_chronological | 153 | 0.9611976215928508 | True |
| GSE113883 | TMrd_lightyellow_mortality | TM_lightyellow_mortality | 153 | 0.9394621347438846 | True |
| GSE113883 | TMrd_magenta_chronological | TM_magenta_chronological | 153 | 0.9914365079239943 | True |
| GSE113883 | TMrd_magenta_mortality | TM_magenta_mortality | 153 | 0.9816705661586161 | True |
| GSE113883 | TMrd_midnightblue_chronological | TM_midnightblue_chronological | 153 | 0.9899846895662666 | True |
| GSE113883 | TMrd_midnightblue_mortality | TM_midnightblue_mortality | 153 | 0.9940291573792055 | True |
| GSE113883 | TMrd_orange_chronological | TM_orange_chronological | 153 | 0.9993880745861151 | True |
| GSE113883 | TMrd_orange_mortality | TM_orange_mortality | 153 | 0.9834667654535715 | True |
| GSE113883 | TMrd_pink_chronological | TM_pink_chronological | 153 | 0.9275087330493074 | True |
| GSE113883 | TMrd_pink_mortality | TM_pink_mortality | 153 | 0.9734411738395223 | True |
| GSE113883 | TMrd_red_chronological | TM_red_chronological | 153 | 0.9135207746329882 | True |
| GSE113883 | TMrd_red_mortality | TM_red_mortality | 153 | 0.9092568696572723 | True |
| GSE113883 | TMrd_salmon_chronological | TM_salmon_chronological | 153 | 0.7345392469299148 | True |
| GSE113883 | TMrd_salmon_mortality | TM_salmon_mortality | 153 | 0.7721477391522547 | True |
| GSE113883 | TMrd_skyblue_chronological | TM_skyblue_chronological | 153 | 0.9956094963997038 | True |
| GSE113883 | TMrd_skyblue_mortality | TM_skyblue_mortality | 153 | 0.9777240903103189 | True |
| GSE113883 | TMrd_tan_chronological | TM_tan_chronological | 153 | 0.9944919883449045 | True |
| GSE113883 | TMrd_tan_mortality | TM_tan_mortality | 153 | 0.9969855122124783 | True |
| GSE113883 | TMrd_turquoise_chronological | TM_turquoise_chronological | 153 | 0.9848767250534688 | True |
| GSE113883 | TMrd_turquoise_mortality | TM_turquoise_mortality | 153 | 0.709240083079747 | True |
| GSE113883 | TMrd_yellow_chronological | TM_yellow_chronological | 153 | 0.9915398637010883 | True |
| GSE113883 | TMrd_yellow_mortality | TM_yellow_mortality | 153 | 0.9864458640309185 | True |
| GSE39445 | TAms_chrono |  | 221 |  | False |
| GSE39445 | TAms_mort |  | 221 |  | False |
| GSE39445 | TArd_chrono | TA_chrono | 221 | 0.7168729049592001 | True |
| GSE39445 | TArd_mort | TA_mort | 221 | 0.604143486046446 | True |
| GSE39445 | TM_blue_mortality |  | 221 |  | False |
| GSE39445 | TM_brown4_mortality |  | 221 |  | False |
| GSE39445 | TM_darkgreen_mortality |  | 221 |  | False |
| GSE39445 | TM_darkmagenta_mortality |  | 221 |  | False |
| GSE39445 | TM_darkred_mortality |  | 221 |  | False |
| GSE39445 | TM_darkslateblue_mortality |  | 221 |  | False |
| GSE39445 | TM_green_mortality |  | 221 |  | False |
| GSE39445 | TM_ivory_mortality |  | 221 |  | False |
| GSE39445 | TM_orange_mortality |  | 221 |  | False |
| GSE39445 | TM_pink_mortality |  | 221 |  | False |
| GSE39445 | TM_plum1_mortality |  | 221 |  | False |
| GSE39445 | TM_sienna3_mortality |  | 221 |  | False |
| GSE39445 | TM_turquoise_mortality |  | 221 |  | False |
| GSE39445 | TM_white_mortality |  | 221 |  | False |
| GSE39445 | TMrd_black_chronological | TM_black_chronological | 221 | 0.6156716575491316 | True |
| GSE39445 | TMrd_black_mortality | TM_black_mortality | 221 | 0.5623415806536157 | True |
| GSE39445 | TMrd_blue_chronological | TM_blue_chronological | 221 | 0.9398712601679585 | True |
| GSE39445 | TMrd_blue_mortality | TM_blue_mortality | 221 | 0.9209845447127848 | True |
| GSE39445 | TMrd_brown_chronological | TM_brown_chronological | 221 | 0.9812093579975115 | True |
| GSE39445 | TMrd_brown_mortality | TM_brown_mortality | 221 | 0.9691094301491392 | True |
| GSE39445 | TMrd_cyan_chronological | TM_cyan_chronological | 221 | 0.9961761163688575 | True |
| GSE39445 | TMrd_cyan_mortality | TM_cyan_mortality | 221 | 0.971482364556249 | True |
| GSE39445 | TMrd_darkgreen_chronological | TM_darkgreen_chronological | 221 | 0.3525667717160795 | True |
| GSE39445 | TMrd_darkgreen_mortality | TM_darkgreen_mortality | 221 | 0.44579757341960596 | True |
| GSE39445 | TMrd_darkgrey_chronological | TM_darkgrey_chronological | 221 | 0.9995573550172681 | True |
| GSE39445 | TMrd_darkgrey_mortality | TM_darkgrey_mortality | 221 | 0.9975641967666337 | True |
| GSE39445 | TMrd_darkorange_chronological | TM_darkorange_chronological | 221 | 0.9119034054728778 | True |
| GSE39445 | TMrd_darkorange_mortality | TM_darkorange_mortality | 221 | 0.785374244844864 | True |
| GSE39445 | TMrd_darkred_chronological | TM_darkred_chronological | 221 | 0.9749944994817658 | True |
| GSE39445 | TMrd_darkred_mortality | TM_darkred_mortality | 221 | 0.9581828092551625 | True |
| GSE39445 | TMrd_green_chronological | TM_green_chronological | 221 | 0.4064609656658057 | True |
| GSE39445 | TMrd_green_mortality | TM_green_mortality | 221 | 0.4200221249991296 | True |
| GSE39445 | TMrd_greenyellow_chronological | TM_greenyellow_chronological | 221 | 0.4883709702082612 | True |
| GSE39445 | TMrd_greenyellow_mortality | TM_greenyellow_mortality | 221 | 0.49800425310188895 | True |
| GSE39445 | TMrd_grey60_chronological | TM_grey60_chronological | 221 | 0.9942081653164923 | True |
| GSE39445 | TMrd_grey60_mortality | TM_grey60_mortality | 221 | 0.9958297207707947 | True |
| GSE39445 | TMrd_lightgreen_chronological | TM_lightgreen_chronological | 221 | 0.8150769921274621 | True |
| GSE39445 | TMrd_lightgreen_mortality | TM_lightgreen_mortality | 221 | 0.8455143141744835 | True |
| GSE39445 | TMrd_lightyellow_chronological | TM_lightyellow_chronological | 221 | 0.9948932818724635 | True |
| GSE39445 | TMrd_lightyellow_mortality | TM_lightyellow_mortality | 221 | 0.9954120816615634 | True |
| GSE39445 | TMrd_magenta_chronological | TM_magenta_chronological | 221 | 0.8085282110073224 | True |
| GSE39445 | TMrd_magenta_mortality | TM_magenta_mortality | 221 | 0.6348671949227344 | True |
| GSE39445 | TMrd_midnightblue_chronological | TM_midnightblue_chronological | 221 | 0.9971563437854497 | True |
| GSE39445 | TMrd_midnightblue_mortality | TM_midnightblue_mortality | 221 | 0.9950175249212819 | True |
| GSE39445 | TMrd_orange_chronological | TM_orange_chronological | 221 | 0.6914705020226076 | True |
| GSE39445 | TMrd_orange_mortality | TM_orange_mortality | 221 | 0.8086918691824693 | True |
| GSE39445 | TMrd_pink_chronological | TM_pink_chronological | 221 | 0.9733867593024046 | True |
| GSE39445 | TMrd_pink_mortality | TM_pink_mortality | 221 | 0.9563999903879549 | True |
| GSE39445 | TMrd_red_chronological | TM_red_chronological | 221 | 0.9969619528313719 | True |
| GSE39445 | TMrd_red_mortality | TM_red_mortality | 221 | 0.9931028633816468 | True |
| GSE39445 | TMrd_salmon_chronological | TM_salmon_chronological | 221 | 0.8785429207204127 | True |
| GSE39445 | TMrd_salmon_mortality | TM_salmon_mortality | 221 | 0.8419400722125159 | True |
| GSE39445 | TMrd_skyblue_chronological | TM_skyblue_chronological | 221 | 0.9760309709575761 | True |
| GSE39445 | TMrd_skyblue_mortality | TM_skyblue_mortality | 221 | 0.98870002507293 | True |
| GSE39445 | TMrd_tan_chronological | TM_tan_chronological | 221 | 0.9939683871101468 | True |
| GSE39445 | TMrd_tan_mortality | TM_tan_mortality | 221 | 0.9979011746024823 | True |
| GSE39445 | TMrd_turquoise_chronological | TM_turquoise_chronological | 221 | 0.848698196093259 | True |
| GSE39445 | TMrd_turquoise_mortality | TM_turquoise_mortality | 221 | 0.8003235214254996 | True |
| GSE39445 | TMrd_yellow_chronological | TM_yellow_chronological | 221 | 0.8819165252524093 | True |
| GSE39445 | TMrd_yellow_mortality | TM_yellow_mortality | 221 | 0.8400942630024313 | True |
| GSE56931 | TAms_chrono |  | 130 |  | False |
| GSE56931 | TAms_mort |  | 130 |  | False |
| GSE56931 | TArd_chrono | TA_chrono | 130 | 0.7450980858928713 | True |
| GSE56931 | TArd_mort | TA_mort | 130 | 0.8050100676413562 | True |
| GSE56931 | TM_blue_mortality |  | 130 |  | False |
| GSE56931 | TM_brown4_mortality |  | 130 |  | False |
| GSE56931 | TM_darkgreen_mortality |  | 130 |  | False |
| GSE56931 | TM_darkmagenta_mortality |  | 130 |  | False |
| GSE56931 | TM_darkred_mortality |  | 130 |  | False |
| GSE56931 | TM_darkslateblue_mortality |  | 130 |  | False |
| GSE56931 | TM_green_mortality |  | 130 |  | False |
| GSE56931 | TM_ivory_mortality |  | 130 |  | False |
| GSE56931 | TM_orange_mortality |  | 130 |  | False |
| GSE56931 | TM_pink_mortality |  | 130 |  | False |
| GSE56931 | TM_plum1_mortality |  | 130 |  | False |
| GSE56931 | TM_sienna3_mortality |  | 130 |  | False |
| GSE56931 | TM_turquoise_mortality |  | 130 |  | False |
| GSE56931 | TM_white_mortality |  | 130 |  | False |
| GSE56931 | TMrd_black_chronological | TM_black_chronological | 130 | 0.7856103504166516 | True |
| GSE56931 | TMrd_black_mortality | TM_black_mortality | 130 | 0.8316839917085791 | True |
| GSE56931 | TMrd_blue_chronological | TM_blue_chronological | 130 | 0.9898740739221221 | True |
| GSE56931 | TMrd_blue_mortality | TM_blue_mortality | 130 | 0.9867674467234642 | True |
| GSE56931 | TMrd_brown_chronological | TM_brown_chronological | 130 | 0.9905076163970152 | True |
| GSE56931 | TMrd_brown_mortality | TM_brown_mortality | 130 | 0.9670686128819599 | True |
| GSE56931 | TMrd_cyan_chronological | TM_cyan_chronological | 130 | 0.97181714709758 | True |
| GSE56931 | TMrd_cyan_mortality | TM_cyan_mortality | 130 | 0.9483787348247747 | True |
| GSE56931 | TMrd_darkgreen_chronological | TM_darkgreen_chronological | 130 | 0.7118663952043937 | True |
| GSE56931 | TMrd_darkgreen_mortality | TM_darkgreen_mortality | 130 | 0.8289160087315823 | True |
| GSE56931 | TMrd_darkgrey_chronological | TM_darkgrey_chronological | 130 | 0.995176199108177 | True |
| GSE56931 | TMrd_darkgrey_mortality | TM_darkgrey_mortality | 130 | 0.9978659578592979 | True |
| GSE56931 | TMrd_darkorange_chronological | TM_darkorange_chronological | 130 | 0.991738309069352 | True |
| GSE56931 | TMrd_darkorange_mortality | TM_darkorange_mortality | 130 | 0.9492728020110949 | True |
| GSE56931 | TMrd_darkred_chronological | TM_darkred_chronological | 130 | 0.9788087925204904 | True |
| GSE56931 | TMrd_darkred_mortality | TM_darkred_mortality | 130 | 0.9894182188732517 | True |
| GSE56931 | TMrd_green_chronological | TM_green_chronological | 130 | 0.9483111504918079 | True |
| GSE56931 | TMrd_green_mortality | TM_green_mortality | 130 | 0.9567898709909836 | True |
| GSE56931 | TMrd_greenyellow_chronological | TM_greenyellow_chronological | 130 | 0.7596470817790817 | True |
| GSE56931 | TMrd_greenyellow_mortality | TM_greenyellow_mortality | 130 | 0.8423416225668665 | True |
| GSE56931 | TMrd_grey60_chronological | TM_grey60_chronological | 130 | 0.9984000381713707 | True |
| GSE56931 | TMrd_grey60_mortality | TM_grey60_mortality | 130 | 0.9957683239953156 | True |
| GSE56931 | TMrd_lightgreen_chronological | TM_lightgreen_chronological | 130 | 0.9573828590626116 | True |
| GSE56931 | TMrd_lightgreen_mortality | TM_lightgreen_mortality | 130 | 0.9532986119488821 | True |
| GSE56931 | TMrd_lightyellow_chronological | TM_lightyellow_chronological | 130 | 0.9926511764321824 | True |
| GSE56931 | TMrd_lightyellow_mortality | TM_lightyellow_mortality | 130 | 0.9920288346095986 | True |
| GSE56931 | TMrd_magenta_chronological | TM_magenta_chronological | 130 | 0.6504538567169743 | True |
| GSE56931 | TMrd_magenta_mortality | TM_magenta_mortality | 130 | 0.701119359944372 | True |
| GSE56931 | TMrd_midnightblue_chronological | TM_midnightblue_chronological | 130 | 0.9963127587783001 | True |
| GSE56931 | TMrd_midnightblue_mortality | TM_midnightblue_mortality | 130 | 0.9912427750370514 | True |
| GSE56931 | TMrd_orange_chronological | TM_orange_chronological | 130 | 0.8006072425594125 | True |
| GSE56931 | TMrd_orange_mortality | TM_orange_mortality | 130 | 0.9034585530713036 | True |
| GSE56931 | TMrd_pink_chronological | TM_pink_chronological | 130 | 0.9449039003333325 | True |
| GSE56931 | TMrd_pink_mortality | TM_pink_mortality | 130 | 0.9792063675778028 | True |
| GSE56931 | TMrd_red_chronological | TM_red_chronological | 130 | 0.9754819508269903 | True |
| GSE56931 | TMrd_red_mortality | TM_red_mortality | 130 | 0.9708538757903807 | True |
| GSE56931 | TMrd_salmon_chronological | TM_salmon_chronological | 130 | 0.9943032849091368 | True |
| GSE56931 | TMrd_salmon_mortality | TM_salmon_mortality | 130 | 0.9929822539605522 | True |
| GSE56931 | TMrd_skyblue_chronological | TM_skyblue_chronological | 130 | 0.9250708861289182 | True |
| GSE56931 | TMrd_skyblue_mortality | TM_skyblue_mortality | 130 | 0.9654663197500384 | True |
| GSE56931 | TMrd_tan_chronological | TM_tan_chronological | 130 | 0.9905783309282685 | True |
| GSE56931 | TMrd_tan_mortality | TM_tan_mortality | 130 | 0.9871237652408624 | True |
| GSE56931 | TMrd_turquoise_chronological | TM_turquoise_chronological | 130 | 0.9126819956345401 | True |
| GSE56931 | TMrd_turquoise_mortality | TM_turquoise_mortality | 130 | 0.9131844816596119 | True |
| GSE56931 | TMrd_yellow_chronological | TM_yellow_chronological | 130 | 0.9875160307799327 | True |
| GSE56931 | TMrd_yellow_mortality | TM_yellow_mortality | 130 | 0.9891394799797621 | True |

가중치 내적과 `predict`의 개인 내 Pearson 최솟값은 0.9999999999999994이다.

## 5. 모듈 분해 (방법, 성분 합 확인)
리비전 코드에서 복합 시계를 모듈 성분으로 나누는 함수는 찾지 못했다. 성분은 파이프라인의 대체, 표준화, 특징 선택을 거친 항에 계수를 곱한 값의 합이다. 같은 outcome의 설치류 23개 모듈 유전자 집합으로 나누었고, 그 집합들 사이의 겹침은 0이다. 어느 모듈에도 없는 항은 `outside`다. 이름은 `CMP_{시계}_{모듈}`, `CMP_{시계}_outside`다. 유전자와 계수는 적지 않았다.

| dataset | clock | n_samples | max_abs_component_sum_plus_intercept_minus_predict | n_components |
| --- | --- | --- | --- | --- |
| GSE39445 | TAms_chrono | 438 | 2.636779683484747e-16 | 24 |
| GSE39445 | TAms_mort | 438 | 3.9968028886505635e-15 | 24 |
| GSE56931 | TAms_chrono | 249 | 2.498001805406602e-16 | 24 |
| GSE56931 | TAms_mort | 249 | 2.942091015256665e-15 | 24 |
| GSE113883 | TAms_chrono | 153 | 2.220446049250313e-16 | 24 |
| GSE113883 | TAms_mort | 153 | 3.774758283725532e-15 | 24 |

## 6. H1·H2
주 모형은 `score ~ cos_t + sin_t + (1 | subject_id)`, t = `clock_time_h`다. 추정은 REML, LRT와 R2_within은 ML이다. σ는 주 모형 REML 잔차 SD다. 부트스트랩은 기존 추출 목록 2000회다. 점추정 336개 모두 REML, ML, 귀무 모형이 수렴했다.

### 6.1 복합 시계 H1 표 (보고서 2 §5.2와 같은 열 + 잔차 SD)
| dataset | clock | role | n_subjects | n_samples | beta_cos | beta_sin | two_a | two_a_ci | r2_within | r2_ci | two_a_over_sd_between | two_a_over_sd_between_ci | two_a_over_sd_resid | two_a_over_sd_resid_ci | peak_h | peak_center | peak_halfwidth_95 | lrt_p | icc_null | sigma_reml | n_boot | n_converged |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | TAms_chrono | primary | 24 | 221 | -0.0143927860487163 | -0.00952749598885632 | 0.0345210353298852 | [0.0267674222392142, 0.0430883577560432] | 0.362193176953973 | [0.23947737962447, 0.501480783091199] | 1.20795424395622 | [0.849253884882144, 1.85139008290745] | 1.61793845569666 | [1.31990237340981, 1.90868192983506] | 14.2335397784943 | 14.2186064321564 | 0.820399743143368 | 5.49745798912395e-20 | 0.615942340893745 | 0.0171269287358502 | 2000 | 2000 |
| GSE39445 | TAms_mort | primary | 24 | 221 | -0.06923427137566 | -0.027178474750457 | 0.148755555495313 | [0.115382654405193, 0.189326598312899] | 0.28699925739554 | [0.1997711867053, 0.40680716136472] | 0.898158153237692 | [0.630953225101292, 1.476290396064] | 1.46659660556601 | [1.23481846668175, 1.73283417786612] | 13.4288594848278 | 13.4014342140422 | 1.16438541795076 | 3.45851149483147e-15 | 0.70812475944782 | 0.0860831456543238 | 2000 | 2000 |
| GSE39445 | TArd_chrono | sensitivity | 24 | 221 | -0.00682260136129383 | -0.0050400285417378 | 0.0169646428829681 | [0.0101424608835498, 0.0265451135455378] | 0.0772645457809408 | [0.0263382625451298, 0.184336938073032] | 0.547230524933817 | [0.318500455628711, 0.995040245812401] | 0.741563361235455 | [0.44059346084177, 1.14645019034183] | 14.4302814414026 | 14.380290740902 | 2.58000634033399 | 0.000381339040015034 | 0.621204123212522 | 0.0220874342286457 | 2000 | 2000 |
| GSE39445 | TArd_mort | sensitivity | 24 | 221 | -0.0400057912279214 | -0.0294787623901695 | 0.0993873385060312 | [0.0493460307243264, 0.165822593798907] | 0.0689139005346088 | [0.0165955063641151, 0.189138051986657] | 0.49738665876213 | [0.247036265052592, 0.922280138007249] | 0.705234253102531 | [0.350204296627395, 1.16074209625344] | 14.4256736372303 | 14.3092444125853 | 2.44539257988713 | 0.000830944712054089 | 0.644543763166638 | 0.136680292294924 | 2000 | 2000 |
| GSE56931 | TAms_chrono | primary | 14 | 130 | -0.00466793532259587 | 0.000666878236547854 | 0.00943066206760042 | [0.0021951992943962, 0.0236400908585612] | 0.0286364463384637 | [0.00160660509292813, 0.168067851985678] | 0.279780359653385 | [0.0643450057738143, 0.923811029379264] | 0.474645399901586 | [0.11166630561518, 1.16159862863218] | 11.4579688514277 | 11.1841054814353 | 7.33573818796443 | 0.181165188487945 | 0.720030290747521 | 0.0197534498990921 | 2000 | 2000 |
| GSE56931 | TAms_mort | primary | 14 | 130 | -0.0314683059679105 | -0.000567238194411096 | 0.0629468359700225 | [0.0286807480528984, 0.10877472086395] | 0.0821416660164724 | [0.0175483148051841, 0.248795554972632] | 0.351471851364103 | [0.152472909328563, 0.843843528981526] | 0.821911051303492 | [0.378799914214617, 1.39571761865303] | 12.0688456403715 | 12.1939802509022 | 3.06033805667255 | 0.00673707816031792 | 0.832803971543668 | 0.0740143176007128 | 2000 | 2000 |
| GSE56931 | TArd_chrono | sensitivity | 14 | 130 | 0.0045007109554629 | -0.000382135634490581 | 0.00903380910751856 | [0.00309170167390487, 0.0160645695251975] | 0.0538147539888595 | [0.00658059390390159, 0.180974069075688] | 0.427470224851876 | [0.130340724563592, 1.24530885886232] | 0.651576340234743 | [0.217557796323542, 1.20534215756649] | 23.6764604105288 | 23.6502748463113 | 3.90872812549738 | 0.0408948152699152 | 0.672675117824225 | 0.0136041736057695 | 2000 | 2000 |
| GSE56931 | TArd_mort | sensitivity | 14 | 130 | 0.0187874873638303 | -0.000500196659950122 | 0.0375882895670823 | [0.00865410311307256, 0.0935952081248861] | 0.0180930649584713 | [0.000981908722421493, 0.108530349641043] | 0.408733131188537 | [0.112863882584578, 1.18423349719421] | 0.37105773052954 | [0.086995597739064, 0.944851317809185] | 23.8983281147265 | 23.804279114472 | 7.92413124959414 | 0.369352805274575 | 0.391987878622619 | 0.101250000347587 | 2000 | 2000 |
| GSE113883 | TAms_chrono | primary | 11 | 153 | -0.0209097494754014 | -0.00103025227051267 | 0.0418702301338288 | [0.0312530124473112, 0.0542300889232854] | 0.265399481367344 | [0.165601193416533, 0.404691023457529] | 0.989182246734361 | [0.642898439783664, 2.18666354352763] | 1.42787363803614 | [1.14547826431237, 1.76569323523371] | 12.188050742382 | 12.2120885405226 | 1.30484190182366 | 2.96936394371969e-10 | 0.646170962385484 | 0.0253117943835499 | 2000 | 2000 |
| GSE113883 | TAms_mort | primary | 11 | 153 | -0.0448086053054778 | 0.0249386275840794 | 0.102562103238955 | [0.063425659634912, 0.157393544808045] | 0.158916853730887 | [0.0696695591092995, 0.312364575276871] | 0.550424138429307 | [0.305408325003736, 1.28501828612413] | 1.05508637345719 | [0.709685078753047, 1.48381598827964] | 10.0600990696032 | 10.0668531335917 | 2.22504509606368 | 4.44466652817908e-06 | 0.765881495710005 | 0.0897840656082817 | 2000 | 2000 |
| GSE113883 | TArd_chrono | sensitivity | 11 | 153 | -0.00377949274533191 | -0.00598199653903614 | 0.014151868852566 | [0.00665885141367422, 0.02329878517137] | 0.0608651139287235 | [0.0125536250920475, 0.180982878204077] | 0.511325197197603 | [0.308550897157257, 1.05216556477667] | 0.716979586324932 | [0.330251153681065, 1.20498873494956] | 15.8476561201441 | 15.7501376541605 | 2.26085084006661 | 0.0116422321333177 | 0.631802587615795 | 0.0192642028494265 | 2000 | 2000 |
| GSE113883 | TArd_mort | sensitivity | 11 | 153 | 0.0132930132957737 | 0.00362344278433063 | 0.0275560185870843 | [0.00669288768637212, 0.0758098095587496] | 0.0092404532676319 | [0.000571305222578605, 0.0691850338061225] | 0.210310163316293 | [0.0494811518706968, 1.13994495495581] | 0.273123665012873 | [0.0664178666002545, 0.743216510684308] | 1.01649318423095 | 0.82267772778087 | 8.52496047006783 | 0.515022143606158 | 0.593709509318723 | 0.101139550636572 | 2000 | 2000 |

### 6.2 모듈 시계 H1 표 (자료 × 목표 × 모듈)
BH q는 자료 × outcome 안에서 `TM_*` 14개에 붙였다. `TMrd_*`는 민감도 재채점이고 이 표에 넣지 않았다.

| dataset | clock | module | outcome | n_subjects | n_samples | two_a | r2_within | peak_h | lrt_p | sigma_reml | two_a_over_sigma | q_bh |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | TM_blue_mortality | blue | mortality | 24 | 221 | 0.028423015342962 | 0.0099227945476269 | 13.5885034887055 | 0.373270670646011 | 0.104438790992741 | 0.272149984433827 | 0.515875649687256 |
| GSE39445 | TM_brown4_mortality | brown4 | mortality | 24 | 221 | 0.096055795987165 | 0.0378072194300597 | 16.8843725314732 | 0.0239273631834682 | 0.188617838706357 | 0.509261460347374 | 0.0478547263669363 |
| GSE39445 | TM_darkgreen_mortality | darkgreen | mortality | 24 | 221 | 0.102756085856872 | 0.0668626567123483 | 21.6270807263401 | 0.00114365878272909 | 0.137004173261766 | 0.750021575332171 | 0.00320224459164145 |
| GSE39445 | TM_darkmagenta_mortality | darkmagenta | mortality | 24 | 221 | 0.013835797480418 | 0.00105697435235164 | 4.74106845127198 | 0.904290769415569 | 0.167065745529944 | 0.0828164830350462 | 0.973851597832151 |
| GSE39445 | TM_darkred_mortality | darkred | mortality | 24 | 221 | 0.29509567836589 | 0.178270296305496 | 13.7348579810592 | 3.90155877588307e-09 | 0.234111634241769 | 1.26049130074904 | 2.73109114311815e-08 |
| GSE39445 | TM_darkslateblue_mortality | darkslateblue | mortality | 24 | 221 | 0.0726208045227383 | 0.0191910312723096 | 8.87263237795086 | 0.15234097503287 | 0.188418574747401 | 0.385422746245139 | 0.266596706307523 |
| GSE39445 | TM_green_mortality | green | mortality | 24 | 221 | 0.00643537131801889 | 8.468053293198e-05 | 22.6214954477059 | 0.984700668978151 | 0.182744499264945 | 0.0352151301073573 | 0.984700668978151 |
| GSE39445 | TM_ivory_mortality | ivory | mortality | 24 | 221 | 0.15416647734129 | 0.0796437356489054 | 12.4966279007714 | 0.000299421858173475 | 0.189064090215468 | 0.815419137317899 | 0.00104797650360716 |
| GSE39445 | TM_orange_mortality | orange | mortality | 24 | 221 | 0.149801158685789 | 0.0444142706863354 | 11.822008947681 | 0.0139049002532302 | 0.252340846963059 | 0.59364609609842 | 0.0324447672575371 |
| GSE39445 | TM_pink_mortality | pink | mortality | 24 | 221 | 0.0592331501447375 | 0.0137330246116952 | 21.9073270105632 | 0.239909592126786 | 0.173899463980596 | 0.340617209443192 | 0.373192698863889 |
| GSE39445 | TM_plum1_mortality | plum1 | mortality | 24 | 221 | 0.053249192488396 | 0.00927404394092191 | 17.4525490054337 | 0.405330867611416 | 0.213607366424732 | 0.249285375217428 | 0.515875649687256 |
| GSE39445 | TM_sienna3_mortality | sienna3 | mortality | 24 | 221 | 0.0431372587217732 | 0.00539529073773515 | 17.2441489753633 | 0.563852843716794 | 0.218153909833279 | 0.197737729086497 | 0.657828317669594 |
| GSE39445 | TM_turquoise_mortality | turquoise | mortality | 24 | 221 | 0.202066580091206 | 0.161203377324141 | 12.6215119266462 | 3.6492902252445e-08 | 0.166914692967102 | 1.21059791980705 | 1.7030021051141e-07 |
| GSE39445 | TM_white_mortality | white | mortality | 24 | 221 | 0.460384717293802 | 0.191346047227201 | 12.4295874818769 | 9.87748381754577e-10 | 0.341043749396897 | 1.34992861797921 | 1.38284773445641e-08 |
| GSE56931 | TM_blue_mortality | blue | mortality | 14 | 130 | 0.0808339449425606 | 0.0497422574676456 | 13.6791142412563 | 0.051072178545924 | 0.119159347573391 | 0.678368475396146 | 0.238336833214312 |
| GSE56931 | TM_brown4_mortality | brown4 | mortality | 14 | 130 | 0.0287502309565396 | 0.00709118066786285 | 20.4557398041171 | 0.663715274652139 | 0.135842927841505 | 0.211643192718019 | 0.774334487094162 |
| GSE56931 | TM_darkgreen_mortality | darkgreen | mortality | 14 | 130 | 0.0523828216038218 | 0.0375017581937666 | 21.9619381872075 | 0.115676315329639 | 0.103403969957197 | 0.506584240677654 | 0.269911402435825 |
| GSE56931 | TM_darkmagenta_mortality | darkmagenta | mortality | 14 | 130 | 0.035081075302834 | 0.0199994893739058 | 9.27084430301793 | 0.313742149841025 | 0.0968830718140463 | 0.36209705829897 | 0.48804334419715 |
| GSE56931 | TM_darkred_mortality | darkred | mortality | 14 | 130 | 0.0402668379033227 | 0.00820772603884301 | 11.7734322180973 | 0.611871009796173 | 0.15622295049316 | 0.257752383860435 | 0.774334487094162 |
| GSE56931 | TM_darkslateblue_mortality | darkslateblue | mortality | 14 | 130 | 0.0731912455962926 | 0.0237167195111732 | 2.87285023224577 | 0.248179182541772 | 0.160130617277015 | 0.457072150478737 | 0.435304924885386 |
| GSE56931 | TM_green_mortality | green | mortality | 14 | 130 | 0.0794279551644897 | 0.024876843927335 | 5.05587535385136 | 0.248745671363078 | 0.188547708744067 | 0.421261842392921 | 0.435304924885386 |
| GSE56931 | TM_ivory_mortality | ivory | mortality | 14 | 130 | 0.0322554999554806 | 0.00997665573305462 | 7.87063690852299 | 0.564224124867114 | 0.12907929180513 | 0.24988903722974 | 0.774334487094162 |
| GSE56931 | TM_orange_mortality | orange | mortality | 14 | 130 | 0.1010908931547 | 0.0408645878620076 | 14.108225165769 | 0.0887718636698445 | 0.165314926836461 | 0.611504932368897 | 0.250767680024482 |
| GSE56931 | TM_pink_mortality | pink | mortality | 14 | 130 | 0.0853250534395423 | 0.040705865606894 | 14.677684396896 | 0.0895598857230293 | 0.140632016872788 | 0.606725661317401 | 0.250767680024482 |
| GSE56931 | TM_plum1_mortality | plum1 | mortality | 14 | 130 | 0.113753052017128 | 0.18349718838146 | 8.37750606840522 | 7.80676070385971e-06 | 0.0955192037703387 | 1.19089196231817 | 0.000109294649854036 |
| GSE56931 | TM_sienna3_mortality | sienna3 | mortality | 14 | 130 | 0.0139511697174672 | 0.00109021573396884 | 4.81468470019662 | 0.945816074805526 | 0.16432438884991 | 0.0849001771137568 | 0.945816074805526 |
| GSE56931 | TM_turquoise_mortality | turquoise | mortality | 14 | 130 | 0.117090918717615 | 0.113992458020858 | 11.1568084467969 | 0.000938968376308667 | 0.120114281244289 | 0.974829283451104 | 0.00657277863416067 |
| GSE56931 | TM_white_mortality | white | mortality | 14 | 130 | 0.0315788803234449 | 0.00407500508827618 | 7.32458484251542 | 0.805171045341596 | 0.204325101407252 | 0.154552133369572 | 0.867107279598642 |
| GSE113883 | TM_blue_mortality | blue | mortality | 11 | 153 | 0.173527240755599 | 0.214396832503267 | 12.8022407140148 | 2.99240551604033e-08 | 0.118224761050428 | 1.46777408737229 | 4.18936772245646e-07 |
| GSE113883 | TM_brown4_mortality | brown4 | mortality | 11 | 153 | 0.178674571814849 | 0.157724635491724 | 15.7126208094821 | 4.94780975518217e-06 | 0.142656440746966 | 1.25248163265036 | 2.30897788575168e-05 |
| GSE113883 | TM_darkgreen_mortality | darkgreen | mortality | 11 | 153 | 0.0745648132983122 | 0.0672849101782151 | 2.20521547998296 | 0.00716681997896647 | 0.0963839693256292 | 0.77362256213373 | 0.0125419349631913 |
| GSE113883 | TM_darkmagenta_mortality | darkmagenta | mortality | 11 | 153 | 0.0615364120567157 | 0.0480979286241355 | 23.9756265260507 | 0.0311852781764266 | 0.100560340029116 | 0.61193520267432 | 0.0485104327188859 |
| GSE113883 | TM_darkred_mortality | darkred | mortality | 11 | 153 | 0.0496435030252984 | 0.00562455198119649 | 11.5236749070229 | 0.670475135944095 | 0.244613209250725 | 0.202946942960936 | 0.719630378574929 |
| GSE113883 | TM_darkslateblue_mortality | darkslateblue | mortality | 11 | 153 | 0.128910944461665 | 0.106020153430254 | 14.2275237641524 | 0.000350399669600654 | 0.129857630593891 | 0.992709815142199 | 0.000817599229068192 |
| GSE113883 | TM_green_mortality | green | mortality | 11 | 153 | 0.0364319508430803 | 0.0046307945028452 | 14.6287593458878 | 0.719630378574929 | 0.184746658013381 | 0.19719951221224 | 0.719630378574929 |
| GSE113883 | TM_ivory_mortality | ivory | mortality | 11 | 153 | 0.0876105034842862 | 0.0701260496254157 | 7.5008481617592 | 0.00569706454852678 | 0.120288638544536 | 0.72833564785796 | 0.0113941290970536 |
| GSE113883 | TM_orange_mortality | orange | mortality | 11 | 153 | 0.279901830228446 | 0.107020523084554 | 14.234500535859 | 0.000275684661996804 | 0.278026238496061 | 1.00674609613298 | 0.000771917053591052 |
| GSE113883 | TM_pink_mortality | pink | mortality | 11 | 153 | 0.098663384479717 | 0.0294133848960331 | 19.6610431539467 | 0.117640340371181 | 0.213509521786821 | 0.462102971586567 | 0.164696476519653 |
| GSE113883 | TM_plum1_mortality | plum1 | mortality | 11 | 153 | 0.0862195342400483 | 0.0232137265431132 | 5.68063447147643 | 0.189923440275026 | 0.201917726166924 | 0.427003294246544 | 0.241720742168214 |
| GSE113883 | TM_sienna3_mortality | sienna3 | mortality | 11 | 153 | 0.207536886361561 | 0.110658443147596 | 14.0211721741388 | 0.000221635449347993 | 0.204013577753964 | 1.01726997117734 | 0.000771917053591052 |
| GSE113883 | TM_turquoise_mortality | turquoise | mortality | 11 | 153 | 0.159616127666438 | 0.184278022966686 | 10.9484328837053 | 5.27669282979285e-07 | 0.126043076003201 | 1.26636172908367 | 3.693684980855e-06 |
| GSE113883 | TM_white_mortality | white | mortality | 11 | 153 | 0.0558821201943435 | 0.0116043427581355 | 8.13285702602638 | 0.443152816507367 | 0.198156566593241 | 0.282009933635222 | 0.517011619258595 |

### 6.3 복합 시계 성분 H1 표 (자료 × 시계 × 성분)
BH q는 자료 × 부모 시계 안에서 성분 24개에 붙였다.

| dataset | parent_clock | module | clock | n_subjects | n_samples | two_a | r2_within | peak_h | lrt_p | sigma_reml | two_a_over_sigma | q_bh |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | TAms_chrono | black | CMP_TAms_chrono_black | 24 | 221 | 0.000385169382820019 | 0.0351770721952061 | 22.184882672488 | 0.0305697894001868 | 0.00071712016842318 | 0.537105773592923 | 0.0431573497414402 |
| GSE39445 | TAms_chrono | blue | CMP_TAms_chrono_blue | 24 | 221 | 0.0046716963880829 | 0.17078334093659 | 13.3428866744713 | 1.12197109042399e-08 | 0.00378441131923289 | 1.23445788367155 | 4.48788436169598e-08 |
| GSE39445 | TAms_chrono | brown | CMP_TAms_chrono_brown | 24 | 221 | 0.00279405597143679 | 0.177477524922973 | 0.91282893412168 | 5.08260476981387e-09 | 0.00218952594143507 | 1.276100875793 | 2.43965028951066e-08 |
| GSE39445 | TAms_chrono | cyan | CMP_TAms_chrono_cyan | 24 | 221 | 7.04599793322103e-05 | 0.00725314478971895 | 16.3009015659973 | 0.487407594178916 | 0.000317718241859036 | 0.221768756241172 | 0.584889113014699 |
| GSE39445 | TAms_chrono | darkgreen | CMP_TAms_chrono_darkgreen | 24 | 221 | 0.000379525465303403 | 0.18514240834314 | 0.953567339921794 | 2.00928004513153e-09 | 0.000289908012206364 | 1.30912375416947 | 1.60742403610522e-08 |
| GSE39445 | TAms_chrono | darkgrey | CMP_TAms_chrono_darkgrey | 24 | 221 | 0.000172930073635673 | 0.055096759084984 | 6.66692463256511 | 0.00370347117557397 | 0.000270800415016274 | 0.638588658090798 | 0.00592555388091835 |
| GSE39445 | TAms_chrono | darkorange | CMP_TAms_chrono_darkorange | 24 | 221 | 2.8143396789956e-21 | -0.11969300040921 | 17.8385981962938 | 1 | 5.15796942353567e-21 | 0.545629383949786 | 1 |
| GSE39445 | TAms_chrono | darkred | CMP_TAms_chrono_darkred | 24 | 221 | 0.0013272510034408 | 0.0744979974100127 | 0.0273931320867367 | 0.000567829369035122 | 0.00168346467497796 | 0.78840442758811 | 0.00113565873807024 |
| GSE39445 | TAms_chrono | green | CMP_TAms_chrono_green | 24 | 221 | 0.000230688577514093 | 0.00122196516546524 | 2.31220113470282 | 0.878528838821772 | 0.00238573797570703 | 0.0966948507602674 | 0.916725744857501 |
| GSE39445 | TAms_chrono | greenyellow | CMP_TAms_chrono_greenyellow | 24 | 221 | 0.00297838964147632 | 0.0418463080639903 | 12.4067813552267 | 0.0132218132724629 | 0.00504964050525418 | 0.589822114738127 | 0.0198327199086943 |
| GSE39445 | TAms_chrono | grey60 | CMP_TAms_chrono_grey60 | 24 | 221 | 0.000133784408547701 | 0.00314133783078596 | 2.22409414111676 | 0.747069432748383 | 0.000918053637023927 | 0.145726135328425 | 0.814984835725509 |
| GSE39445 | TAms_chrono | lightgreen | CMP_TAms_chrono_lightgreen | 24 | 221 | 0.000255736953538969 | 0.0144076652804407 | 23.9770924259699 | 0.234146616117882 | 0.000747660750801775 | 0.342049456608124 | 0.299419434845185 |
| GSE39445 | TAms_chrono | lightyellow | CMP_TAms_chrono_lightyellow | 24 | 221 | 0.000626257741048089 | 0.124763990385334 | 12.5223151879221 | 2.16885919316036e-06 | 0.000598391552134617 | 1.04656848649361 | 6.00047155091236e-06 |
| GSE39445 | TAms_chrono | magenta | CMP_TAms_chrono_magenta | 24 | 221 | 0.000123814810420957 | 0.0145164744269971 | 13.6170581919721 | 0.237040385919105 | 0.000376225164770737 | 0.329097630926434 | 0.299419434845185 |
| GSE39445 | TAms_chrono | midnightblue | CMP_TAms_chrono_midnightblue | 24 | 221 | 0.000191191521795514 | 0.0949682899792058 | 1.49895476665748 | 6.04871110150721e-05 | 0.00021823201457688 | 0.876092914993277 | 0.000145169066436173 |
| GSE39445 | TAms_chrono | orange | CMP_TAms_chrono_orange | 24 | 221 | 0.00273977200755254 | 0.199927892518597 | 2.44005820482784 | 3.47997379391709e-10 | 0.002067441273765 | 1.32519943483723 | 4.17596855270051e-09 |
| GSE39445 | TAms_chrono | outside | CMP_TAms_chrono_outside | 24 | 221 | 0.0248257844623488 | 0.2714245765219 | 14.5464678599487 | 2.73026438742455e-14 | 0.015307414993185 | 1.62181429545103 | 6.55263452981891e-13 |
| GSE39445 | TAms_chrono | pink | CMP_TAms_chrono_pink | 24 | 221 | 0.00308338842326159 | 0.08254146025789 | 15.1251668500058 | 0.000179836301483848 | 0.00388469500114373 | 0.793727286789253 | 0.000392370112328395 |
| GSE39445 | TAms_chrono | red | CMP_TAms_chrono_red | 24 | 221 | 0.000319266696021939 | 0.123681119898952 | 2.08110231751181 | 2.25017683159214e-06 | 0.000317035128253838 | 1.00703886594647 | 6.00047155091236e-06 |
| GSE39445 | TAms_chrono | salmon | CMP_TAms_chrono_salmon | 24 | 221 | 0.000866469684001931 | 0.151370720346455 | 0.220172586946842 | 1.09725422252285e-07 | 0.000736511623295524 | 1.17645079398056 | 3.7620144772212e-07 |
| GSE39445 | TAms_chrono | skyblue | CMP_TAms_chrono_skyblue | 24 | 221 | 0.00141293675210103 | 0.0683337458832771 | 12.080368522175 | 0.00146545399052842 | 0.00192201005283178 | 0.735134943763319 | 0.00251220684090586 |
| GSE39445 | TAms_chrono | tan | CMP_TAms_chrono_tan | 24 | 221 | 0.000260039406690764 | 0.0026938555288839 | 19.6145710918719 | 0.742431224685133 | 0.00175737067725301 | 0.147970721292129 | 0.814984835725509 |
| GSE39445 | TAms_chrono | turquoise | CMP_TAms_chrono_turquoise | 24 | 221 | 0.00582473108169271 | 0.178263183429586 | 13.1432263174321 | 4.04519641858719e-09 | 0.00456077363979793 | 1.27713663113322 | 2.42711785115231e-08 |
| GSE39445 | TAms_chrono | yellow | CMP_TAms_chrono_yellow | 24 | 221 | 0.000985373040293526 | 0.0697413041380407 | 13.045845902978 | 0.000842131738830101 | 0.00131280060907163 | 0.750588500252409 | 0.00155470474860942 |
| GSE39445 | TAms_mort | black | CMP_TAms_mort_black | 24 | 221 | 0.00346504710347034 | 0.0555859536381143 | 21.9721917086195 | 0.00379945291118488 | 0.00508219453593923 | 0.681801351555303 | 0.0113983587335546 |
| GSE39445 | TAms_mort | blue | CMP_TAms_mort_blue | 24 | 221 | 0.00755216423815708 | 0.0962757490358684 | 13.5990903771125 | 4.42502439411487e-05 | 0.00850732541210431 | 0.887724857381356 | 0.000212401170917514 |
| GSE39445 | TAms_mort | brown | CMP_TAms_mort_brown | 24 | 221 | 0.00579936089754586 | 0.0343048751806052 | 18.2199094109322 | 0.033054519880319 | 0.0117962511928895 | 0.491627450341307 | 0.0650647454648595 |
| GSE39445 | TAms_mort | cyan | CMP_TAms_mort_cyan | 24 | 221 | 0.00176086319032717 | 0.0354051501693395 | 15.191930954719 | 0.0274976721731882 | 0.00348354081701319 | 0.505480854918459 | 0.0634393514012682 |
| GSE39445 | TAms_mort | darkgreen | CMP_TAms_mort_darkgreen | 24 | 221 | 0.00148510941397869 | 0.0349415133862532 | 22.0466118861608 | 0.0289575929697768 | 0.00274798725956552 | 0.540435334555918 | 0.0634393514012682 |
| GSE39445 | TAms_mort | darkgrey | CMP_TAms_mort_darkgrey | 24 | 221 | 0.000284328589232324 | 0.00107507590516864 | 3.6918066626364 | 0.903602200825879 | 0.00339904116761961 | 0.0836496456532897 | 0.903602200825879 |
| GSE39445 | TAms_mort | darkorange | CMP_TAms_mort_darkorange | 24 | 221 | 0.000730940511697513 | 0.0221388345831234 | 5.94070336203217 | 0.112084999270766 | 0.00186808431284697 | 0.391278116662495 | 0.141581051710441 |
| GSE39445 | TAms_mort | darkred | CMP_TAms_mort_darkred | 24 | 221 | 0.00293675296326274 | 0.0686029149245168 | 17.0651005937284 | 0.000874831345854333 | 0.00416526780795402 | 0.705057417353741 | 0.00299942175721486 |
| GSE39445 | TAms_mort | green | CMP_TAms_mort_green | 24 | 221 | 0.00117521009651132 | 0.00361701830260175 | 21.9975484531286 | 0.678699836798665 | 0.006636472312622 | 0.177083553000918 | 0.74039982196218 |
| GSE39445 | TAms_mort | greenyellow | CMP_TAms_mort_greenyellow | 24 | 221 | 0.00791897747233188 | 0.0262992487528657 | 15.2244616199829 | 0.0734486697644695 | 0.0184220480987746 | 0.42986411879245 | 0.11047008286511 |
| GSE39445 | TAms_mort | grey60 | CMP_TAms_mort_grey60 | 24 | 221 | 0.00520829528542398 | 0.154615629594791 | 3.21032356030966 | 6.58333236000225e-08 | 0.00464594332733108 | 1.12104150190225 | 5.2666658880018e-07 |
| GSE39445 | TAms_mort | lightgreen | CMP_TAms_mort_lightgreen | 24 | 221 | 0.00198579317992431 | 0.0187608156813303 | 20.6147851893821 | 0.147648803103419 | 0.00514112169165349 | 0.386256793560871 | 0.177178563724102 |
| GSE39445 | TAms_mort | lightyellow | CMP_TAms_mort_lightyellow | 24 | 221 | 0.000486552980879079 | 0.0329064071660083 | 7.09835660357778 | 0.0386877941924554 | 0.000996488048267813 | 0.488267753662325 | 0.0663219329013521 |
| GSE39445 | TAms_mort | magenta | CMP_TAms_mort_magenta | 24 | 221 | 0.000702607544500171 | 0.00289683075405223 | 12.0631943087535 | 0.739888101016421 | 0.00453156982438536 | 0.155047273180982 | 0.772057148886701 |
| GSE39445 | TAms_mort | midnightblue | CMP_TAms_mort_midnightblue | 24 | 221 | 0.000776093588366589 | 0.0227465361096185 | 7.42351295848493 | 0.105596805133434 | 0.00190307729133655 | 0.407809809879834 | 0.140795740177912 |
| GSE39445 | TAms_mort | orange | CMP_TAms_mort_orange | 24 | 221 | 0.00122904134242102 | 0.0246101607548138 | 21.1496828705772 | 0.0768612488835501 | 0.0027158776839414 | 0.452539283962663 | 0.11047008286511 |
| GSE39445 | TAms_mort | outside | CMP_TAms_mort_outside | 24 | 221 | 0.129911986217037 | 0.237096332526042 | 12.884283653621 | 2.89045569401259e-12 | 0.0845788028708217 | 1.53598752651362 | 6.93709366563023e-11 |
| GSE39445 | TAms_mort | pink | CMP_TAms_mort_pink | 24 | 221 | 0.0130425611205924 | 0.0873380565104375 | 23.847585220868 | 0.00012905502608409 | 0.0150244120573155 | 0.868091281764462 | 0.00051622010433636 |
| GSE39445 | TAms_mort | red | CMP_TAms_mort_red | 24 | 221 | 0.00312046619974404 | 0.148963096843102 | 1.93054574874185 | 1.37683666459772e-07 | 0.00277775800597622 | 1.12337582792688 | 8.26101998758634e-07 |
| GSE39445 | TAms_mort | salmon | CMP_TAms_mort_salmon | 24 | 221 | 0.00159277869908285 | 0.0255557295882924 | 23.4050115532476 | 0.078249642029453 | 0.00348011164858107 | 0.457680344747636 | 0.11047008286511 |
| GSE39445 | TAms_mort | skyblue | CMP_TAms_mort_skyblue | 24 | 221 | 0.00170528157945959 | 0.0358922480027858 | 10.0013177262601 | 0.0290763693922479 | 0.00315541442486254 | 0.540430304819272 | 0.0634393514012682 |
| GSE39445 | TAms_mort | tan | CMP_TAms_mort_tan | 24 | 221 | 0.00135660873890544 | 0.013289025818377 | 16.2931458868751 | 0.270058847822125 | 0.00452409246388841 | 0.299863176920891 | 0.308638683225286 |
| GSE39445 | TAms_mort | turquoise | CMP_TAms_mort_turquoise | 24 | 221 | 0.021004897906594 | 0.194311362220051 | 13.2360241903645 | 5.64830096146801e-10 | 0.0156225623512403 | 1.34452322444572 | 6.77796115376162e-09 |
| GSE39445 | TAms_mort | yellow | CMP_TAms_mort_yellow | 24 | 221 | 0.00360689859335469 | 0.034718154551784 | 11.8688708007029 | 0.0352434037934656 | 0.00689146787355445 | 0.523386114472931 | 0.0650647454648595 |
| GSE56931 | TAms_chrono | black | CMP_TAms_chrono_black | 14 | 130 | 0.000392226799696898 | 0.0372431827055033 | 15.7909220347951 | 0.11419983941984 | 0.000704703079408971 | 0.556584483816724 | 0.171299759129761 |
| GSE56931 | TAms_chrono | blue | CMP_TAms_chrono_blue | 14 | 130 | 0.0008882891088178 | 0.0433549367161582 | 7.47099497282694 | 0.0672849430631925 | 0.00162324255677055 | 0.547231284143422 | 0.115345616679759 |
| GSE56931 | TAms_chrono | brown | CMP_TAms_chrono_brown | 14 | 130 | 0.00111228343107183 | 0.0704861736594903 | 0.84038477850989 | 0.0149948103149181 | 0.00139264153389192 | 0.798686096890568 | 0.0359875447558034 |
| GSE56931 | TAms_chrono | cyan | CMP_TAms_chrono_cyan | 14 | 130 | 9.80369923598424e-06 | 0.0838566150669738 | 14.641458287959 | 0.00689362837685512 | 1.11034326994849e-05 | 0.882943095286111 | 0.0206808851305654 |
| GSE56931 | TAms_chrono | darkgreen | CMP_TAms_chrono_darkgreen | 14 | 130 | 8.91959461601454e-05 | 0.0490980390761122 | 22.1952154007943 | 0.0507730048344837 | 0.000148349461262685 | 0.601255612261404 | 0.101546009668967 |
| GSE56931 | TAms_chrono | darkgrey | CMP_TAms_chrono_darkgrey | 14 | 130 | 6.69999759573e-05 | 0.0272965184045829 | 23.3051532450979 | 0.185498331210528 | 0.000142615588482263 | 0.469794197607175 | 0.234313681529088 |
| GSE56931 | TAms_chrono | darkorange | CMP_TAms_chrono_darkorange | 14 | 130 | 3.41965687789239e-06 | 0.0225016387636331 | 18.2149911380555 | 0.247006490554749 | 8.50322854640472e-06 | 0.402159821911204 | 0.296407788665699 |
| GSE56931 | TAms_chrono | darkred | CMP_TAms_chrono_darkred | 14 | 130 | 0.001415502198911 | 0.141248858555156 | 16.6156267629132 | 0.000140441433517734 | 0.00126613623368026 | 1.11796990028204 | 0.000842648601106405 |
| GSE56931 | TAms_chrono | green | CMP_TAms_chrono_green | 14 | 130 | 0.00104658360551093 | 0.042815093590585 | 7.90937139661737 | 0.0741951373568655 | 0.00195083541383421 | 0.536479704074035 | 0.118712219770985 |
| GSE56931 | TAms_chrono | greenyellow | CMP_TAms_chrono_greenyellow | 14 | 130 | 0.00190434199803476 | 0.0344067738797542 | 11.8973116427241 | 0.127392605730811 | 0.0035612716869422 | 0.534736511403284 | 0.179848384561145 |
| GSE56931 | TAms_chrono | grey60 | CMP_TAms_chrono_grey60 | 14 | 130 | 0.000446196970418535 | 0.0772890261098173 | 3.32251958105996 | 0.0104212032910459 | 0.000537625109514048 | 0.82994071988536 | 0.0277898754427891 |
| GSE56931 | TAms_chrono | lightgreen | CMP_TAms_chrono_lightgreen | 14 | 130 | 0.000449236733119421 | 0.032346665130795 | 1.73520703694657 | 0.146976760525907 | 0.000827776274269605 | 0.54270307942301 | 0.195969014034543 |
| GSE56931 | TAms_chrono | lightyellow | CMP_TAms_chrono_lightyellow | 14 | 130 | 0.000188934039327846 | 0.0118697785415157 | 12.9002100073182 | 0.473593151727197 | 0.00057040097226672 | 0.331230219641877 | 0.541249316259654 |
| GSE56931 | TAms_chrono | magenta | CMP_TAms_chrono_magenta | 14 | 130 | 0.00049326571830635 | 0.0672551851282774 | 6.22299583533311 | 0.0214821143869443 | 0.000728083309996789 | 0.677485270618998 | 0.0468700677533329 |
| GSE56931 | TAms_chrono | midnightblue | CMP_TAms_chrono_midnightblue | 14 | 130 | 0.000173302407683553 | 0.0489476840708664 | 2.79385213790353 | 0.0563936703288362 | 0.000261762183822056 | 0.662060520557707 | 0.104111391376313 |
| GSE56931 | TAms_chrono | orange | CMP_TAms_chrono_orange | 14 | 130 | 0.00112988806312844 | 0.101344814081915 | 0.0303979642477727 | 0.00154166408634827 | 0.00116608442002731 | 0.96895905967252 | 0.00616665634539307 |
| GSE56931 | TAms_chrono | outside | CMP_TAms_chrono_outside | 14 | 130 | 0.00396618419794712 | 0.00519885425286193 | 14.0313775069439 | 0.744157321537192 | 0.0187168331408668 | 0.211904661867571 | 0.776511987690983 |
| GSE56931 | TAms_chrono | pink | CMP_TAms_chrono_pink | 14 | 130 | 0.000332531984845647 | 0.00356403674202022 | 8.07688889414835 | 0.804544452826182 | 0.00216565274312831 | 0.1535481558162 | 0.804544452826182 |
| GSE56931 | TAms_chrono | red | CMP_TAms_chrono_red | 14 | 130 | 7.45817119379737e-05 | 0.12606738198406 | 23.2474694413108 | 0.000430880231917367 | 7.20342771826523e-05 | 1.03536420236247 | 0.00206822511320336 |
| GSE56931 | TAms_chrono | salmon | CMP_TAms_chrono_salmon | 14 | 130 | 0.000835903959848489 | 0.274438462566386 | 22.9454850095866 | 7.8991980142257e-09 | 0.000502498512517678 | 1.66349539157906 | 6.31935841138056e-08 |
| GSE56931 | TAms_chrono | skyblue | CMP_TAms_chrono_skyblue | 14 | 130 | 0.00152390463176805 | 0.318397229009213 | 5.58846464034619 | 1.7585000358238e-10 | 0.000840706461597774 | 1.81264769735664 | 4.22040008597713e-09 |
| GSE56931 | TAms_chrono | tan | CMP_TAms_chrono_tan | 14 | 130 | 0.000321797600853565 | 0.0131949244135086 | 3.76942463243976 | 0.520552974834114 | 0.00105444537777899 | 0.305181859236157 | 0.567875972546307 |
| GSE56931 | TAms_chrono | turquoise | CMP_TAms_chrono_turquoise | 14 | 130 | 0.0053191322688683 | 0.101993915383935 | 12.1422987334514 | 0.00182870083132063 | 0.00552528154633841 | 0.962689814855366 | 0.00626983142167074 |
| GSE56931 | TAms_chrono | yellow | CMP_TAms_chrono_yellow | 14 | 130 | 0.00121942710904472 | 0.280452248149977 | 11.2503142784199 | 3.91448877682289e-09 | 0.000709510352395572 | 1.71868825441021 | 4.69738653218747e-08 |
| GSE56931 | TAms_mort | black | CMP_TAms_mort_black | 14 | 130 | 0.00326826909655884 | 0.0665924521539014 | 12.6800110610643 | 0.0176351203278607 | 0.00420167581873584 | 0.777848943506108 | 0.0384766261698779 |
| GSE56931 | TAms_mort | blue | CMP_TAms_mort_blue | 14 | 130 | 0.00385975547050352 | 0.0675332781223587 | 7.78311507522425 | 0.0156068270976116 | 0.0056496055476001 | 0.683190257794742 | 0.0374563850342679 |
| GSE56931 | TAms_mort | brown | CMP_TAms_mort_brown | 14 | 130 | 0.00252369114184631 | 0.00708191498340383 | 10.2358127038323 | 0.628675997649183 | 0.0107691754373877 | 0.23434395293485 | 0.656009736677408 |
| GSE56931 | TAms_mort | cyan | CMP_TAms_mort_cyan | 14 | 130 | 0.00133327325079187 | 0.052684800478908 | 12.7449978091296 | 0.0382647061012226 | 0.00191248133424931 | 0.697143144309541 | 0.0612235297619561 |
| GSE56931 | TAms_mort | darkgreen | CMP_TAms_mort_darkgreen | 14 | 130 | 0.00181508107037263 | 0.0234113750175663 | 13.3106109192703 | 0.254796212185218 | 0.00399126035359492 | 0.454763886484577 | 0.305755454622262 |
| GSE56931 | TAms_mort | darkgrey | CMP_TAms_mort_darkgrey | 14 | 130 | 0.00170821021298905 | 0.0794053311380034 | 1.82771852633565 | 0.00718988983536723 | 0.00193933267815723 | 0.88082371437798 | 0.0275912785277717 |
| GSE56931 | TAms_mort | darkorange | CMP_TAms_mort_darkorange | 14 | 130 | 0.000794398336686693 | 0.0252001954706249 | 22.4159354140251 | 0.196990785459946 | 0.00179018122208359 | 0.443753027284072 | 0.259750335803272 |
| GSE56931 | TAms_mort | darkred | CMP_TAms_mort_darkred | 14 | 130 | 0.00297772779631072 | 0.0781747725487902 | 13.327427235357 | 0.00804745623726676 | 0.00343952068724361 | 0.865739173296565 | 0.0275912785277717 |
| GSE56931 | TAms_mort | green | CMP_TAms_mort_green | 14 | 130 | 0.00685462858870916 | 0.0713690635674427 | 1.35169912239713 | 0.0117555948291992 | 0.00827499936304637 | 0.828353971762203 | 0.0352667844875976 |
| GSE56931 | TAms_mort | greenyellow | CMP_TAms_mort_greenyellow | 14 | 130 | 0.0103127141459329 | 0.0729520739815737 | 10.5690522983707 | 0.0139180011563526 | 0.0140056223201804 | 0.736326734376772 | 0.0371146697502737 |
| GSE56931 | TAms_mort | grey60 | CMP_TAms_mort_grey60 | 14 | 130 | 0.00109556555008661 | 0.0153003256462636 | 22.0560982138305 | 0.419546979135674 | 0.00341937177252509 | 0.32039965905128 | 0.479482261869342 |
| GSE56931 | TAms_mort | lightgreen | CMP_TAms_mort_lightgreen | 14 | 130 | 0.00240229752076174 | 0.0268651016071796 | 4.55451090711936 | 0.205635682510924 | 0.00523480178650526 | 0.458908974730351 | 0.259750335803272 |
| GSE56931 | TAms_mort | lightyellow | CMP_TAms_mort_lightyellow | 14 | 130 | 0.00074136847425276 | 0.0391759865145055 | 11.3642799856295 | 0.0913435280768493 | 0.00131498358642629 | 0.563785344475336 | 0.137015292115274 |
| GSE56931 | TAms_mort | magenta | CMP_TAms_mort_magenta | 14 | 130 | 0.0029584980043267 | 0.0615931830928793 | 11.8366943498171 | 0.0251836588484475 | 0.00411627298492145 | 0.718732216051786 | 0.0464929086432877 |
| GSE56931 | TAms_mort | midnightblue | CMP_TAms_mort_midnightblue | 14 | 130 | 0.000333414625153434 | 0.00598018231557751 | 4.42947728221442 | 0.726961128372138 | 0.00161452004103084 | 0.206510056660898 | 0.726961128372138 |
| GSE56931 | TAms_mort | orange | CMP_TAms_mort_orange | 14 | 130 | 0.0019468178741336 | 0.0621672377405967 | 19.9022042211049 | 0.0240133404981232 | 0.00301202787934283 | 0.646347893220154 | 0.0464929086432877 |
| GSE56931 | TAms_mort | outside | CMP_TAms_mort_outside | 14 | 130 | 0.0406417324647928 | 0.0363460322409275 | 12.9853507178368 | 0.116101863386523 | 0.0714734212350075 | 0.568627215019709 | 0.163908513016268 |
| GSE56931 | TAms_mort | pink | CMP_TAms_mort_pink | 14 | 130 | 0.0121356168681622 | 0.239616612648394 | 22.3083584264042 | 1.28990489819454e-07 | 0.00821688014024966 | 1.47691297195842 | 1.54788587783345e-06 |
| GSE56931 | TAms_mort | red | CMP_TAms_mort_red | 14 | 130 | 0.00108538038341601 | 0.0904374128684702 | 18.0588953196784 | 0.0041092724355883 | 0.0013264395599432 | 0.818265992807454 | 0.0197245076908238 |
| GSE56931 | TAms_mort | salmon | CMP_TAms_mort_salmon | 14 | 130 | 0.00323562935156911 | 0.138573731004823 | 23.5368039424505 | 0.000160119510629775 | 0.00289691825532671 | 1.11692117843491 | 0.000960717063778652 |
| GSE56931 | TAms_mort | skyblue | CMP_TAms_mort_skyblue | 14 | 130 | 0.00342476373105239 | 0.259948194608667 | 4.17254798822611 | 1.82401287212257e-08 | 0.0020451385328198 | 1.67458765071058 | 4.37763089309416e-07 |
| GSE56931 | TAms_mort | tan | CMP_TAms_mort_tan | 14 | 130 | 0.000828850797095842 | 0.0114943879113093 | 8.20660259390749 | 0.500787242448664 | 0.00301928633496445 | 0.274518778658865 | 0.546313355398543 |
| GSE56931 | TAms_mort | turquoise | CMP_TAms_mort_turquoise | 14 | 130 | 0.0236853751787755 | 0.173230353075117 | 11.9113647277402 | 1.56295095012086e-05 | 0.0183433697616038 | 1.29122268626747 | 0.000125036076009668 |
| GSE56931 | TAms_mort | yellow | CMP_TAms_mort_yellow | 14 | 130 | 0.00247859258589921 | 0.0560636681904658 | 6.89753361717946 | 0.0364899150549081 | 0.00402858652762719 | 0.615251172812485 | 0.0612235297619561 |
| GSE113883 | TAms_chrono | black | CMP_TAms_chrono_black | 11 | 153 | 0.000657957431221466 | 0.0623752994193139 | 13.5436274515614 | 0.00983624998753281 | 0.000891644941156301 | 0.737914163869101 | 0.0124247368263572 |
| GSE113883 | TAms_chrono | blue | CMP_TAms_chrono_blue | 11 | 153 | 0.00260496949725263 | 0.0646940580324065 | 5.14107004723602 | 0.00899167684965024 | 0.00353072371321684 | 0.737800436636047 | 0.0119889024662003 |
| GSE113883 | TAms_chrono | brown | CMP_TAms_chrono_brown | 11 | 153 | 0.00230936710972663 | 0.153108476410065 | 2.22823394895771 | 7.80140348046347e-06 | 0.00188691269607259 | 1.2238865711876 | 2.67476690758748e-05 |
| GSE113883 | TAms_chrono | cyan | CMP_TAms_chrono_cyan | 11 | 153 | 0.000152050708791115 | 0.0277164916638674 | 12.2379423959593 | 0.132120550797719 | 0.000325047046244856 | 0.467780619906251 | 0.144131509961148 |
| GSE113883 | TAms_chrono | darkgreen | CMP_TAms_chrono_darkgreen | 11 | 153 | 0.000202384455416656 | 0.0329587118067981 | 17.7987677620622 | 0.0912796579958413 | 0.000395287980903796 | 0.511992433855235 | 0.10953558959501 |
| GSE113883 | TAms_chrono | darkgrey | CMP_TAms_chrono_darkgrey | 11 | 153 | 3.59185560500513e-05 | 0.00106732748373606 | 3.35360966440041 | 0.922870284067173 | 0.000369240627780886 | 0.097276825320983 | 0.922870284067173 |
| GSE113883 | TAms_chrono | darkorange | CMP_TAms_chrono_darkorange | 11 | 153 | 8.29726850739398e-06 | 0.106710641958543 | 2.32475896664902 | 0.000369067064852841 | 8.36952583599756e-06 | 0.991366616219429 | 0.000738134129705682 |
| GSE113883 | TAms_chrono | darkred | CMP_TAms_chrono_darkred | 11 | 153 | 0.00221571996190656 | 0.108201670720904 | 14.9745384915279 | 0.000321036443880568 | 0.00220417357916201 | 1.00523841808726 | 0.000700443150284876 |
| GSE113883 | TAms_chrono | green | CMP_TAms_chrono_green | 11 | 153 | 0.00673498782664778 | 0.11806711547286 | 14.7923040912005 | 0.000141421173117746 | 0.0063752011303465 | 1.05643534830433 | 0.00037712312831399 |
| GSE113883 | TAms_chrono | greenyellow | CMP_TAms_chrono_greenyellow | 11 | 153 | 0.00587824566835351 | 0.0955287130475538 | 3.09831492514547 | 0.000818580234656384 | 0.00624587880625929 | 0.941139886105802 | 0.00130972837545021 |
| GSE113883 | TAms_chrono | grey60 | CMP_TAms_chrono_grey60 | 11 | 153 | 0.000349553223656221 | 0.0125626553116532 | 9.05432753476736 | 0.424556156235223 | 0.00121359930121276 | 0.288030178747558 | 0.443015119549798 |
| GSE113883 | TAms_chrono | lightgreen | CMP_TAms_chrono_lightgreen | 11 | 153 | 0.00129757549155914 | 0.171562610464398 | 2.78585455953181 | 1.60178210460723e-06 | 0.00098445179349066 | 1.31806910215299 | 7.68855410211471e-06 |
| GSE113883 | TAms_chrono | lightyellow | CMP_TAms_chrono_lightyellow | 11 | 153 | 0.000306642329836029 | 0.0287326568918248 | 4.56659393564997 | 0.128395540265793 | 0.000627240414901801 | 0.488875274218477 | 0.144131509961148 |
| GSE113883 | TAms_chrono | magenta | CMP_TAms_chrono_magenta | 11 | 153 | 0.00068157798236862 | 0.124594097586489 | 13.0667941717723 | 6.91812597711582e-05 | 0.000637826048425182 | 1.06859540160121 | 0.000207543779313475 |
| GSE113883 | TAms_chrono | midnightblue | CMP_TAms_chrono_midnightblue | 11 | 153 | 0.000435849468662386 | 0.219931048846424 | 3.91450566828495 | 1.90950094670547e-08 | 0.000283652568910438 | 1.53656097787714 | 1.14570056802328e-07 |
| GSE113883 | TAms_chrono | orange | CMP_TAms_chrono_orange | 11 | 153 | 0.00346256188964057 | 0.165758178115896 | 0.765259718398529 | 2.46533427328184e-06 | 0.00277626557822185 | 1.24720124645218 | 9.86133709312738e-06 |
| GSE113883 | TAms_chrono | outside | CMP_TAms_chrono_outside | 11 | 153 | 0.0328298556060485 | 0.236551735229012 | 11.9465235982556 | 4.57504735372346e-09 | 0.0215717073863694 | 1.52189416526171 | 5.49005682446815e-08 |
| GSE113883 | TAms_chrono | pink | CMP_TAms_chrono_pink | 11 | 153 | 0.00680815200931047 | 0.230040314309368 | 14.314100268862 | 9.15856845832467e-09 | 0.00432115768038648 | 1.57553889787738 | 7.32685476665974e-08 |
| GSE113883 | TAms_chrono | red | CMP_TAms_chrono_red | 11 | 153 | 0.00042149831634286 | 0.0977528909986353 | 1.00061811727099 | 0.000662615378892636 | 0.000455047632044851 | 0.926272958390685 | 0.00113591207810166 |
| GSE113883 | TAms_chrono | salmon | CMP_TAms_chrono_salmon | 11 | 153 | 0.000765574703980403 | 0.112122432247006 | 16.3306801444531 | 0.000231142744731363 | 0.000754400986837482 | 1.0148113766258 | 0.000554742587355272 |
| GSE113883 | TAms_chrono | skyblue | CMP_TAms_chrono_skyblue | 11 | 153 | 0.00188372849808141 | 0.0984016743638867 | 11.7951001049552 | 0.000594767508286905 | 0.00208778908865599 | 0.902259959263444 | 0.00109803232299121 |
| GSE113883 | TAms_chrono | tan | CMP_TAms_chrono_tan | 11 | 153 | 0.00144932804800611 | 0.0688926264195883 | 5.10456064392476 | 0.00586060395985471 | 0.00187962293830651 | 0.77107382468524 | 0.00879090593978207 |
| GSE113883 | TAms_chrono | turquoise | CMP_TAms_chrono_turquoise | 11 | 153 | 0.00453731354114232 | 0.0663095830202147 | 14.5730749724508 | 0.00762018516796083 | 0.00588020002013836 | 0.771625714363974 | 0.0107579084724153 |
| GSE113883 | TAms_chrono | yellow | CMP_TAms_chrono_yellow | 11 | 153 | 0.00148510108931472 | 0.261653517520032 | 13.523129022035 | 3.75331743604448e-10 | 0.000873732057528858 | 1.69972141518416 | 9.00796184650676e-09 |
| GSE113883 | TAms_mort | black | CMP_TAms_mort_black | 11 | 153 | 0.00454345704042036 | 0.07710590128205 | 12.3891954147455 | 0.00322941848813597 | 0.00566712482680632 | 0.801721715909477 | 0.00553614597966166 |
| GSE113883 | TAms_mort | blue | CMP_TAms_mort_blue | 11 | 153 | 0.0144536916852471 | 0.128641093388183 | 4.12966405311176 | 5.91191299150702e-05 | 0.0131029952097838 | 1.10308303207307 | 0.000283771823592337 |
| GSE113883 | TAms_mort | brown | CMP_TAms_mort_brown | 11 | 153 | 0.0047891060829532 | 0.0181079780728034 | 6.90433520830122 | 0.279448438076593 | 0.0132368500982061 | 0.361801036305625 | 0.319369643516106 |
| GSE113883 | TAms_mort | cyan | CMP_TAms_mort_cyan | 11 | 153 | 0.00321385465249041 | 0.0862829098225801 | 12.4889961052185 | 0.00166177113484966 | 0.0037740629271863 | 0.851563610489785 | 0.00362568247603562 |
| GSE113883 | TAms_mort | darkgreen | CMP_TAms_mort_darkgreen | 11 | 153 | 0.000792005507731334 | 0.00604253043836689 | 1.84759302418195 | 0.656620461570209 | 0.00358248929680315 | 0.221076866423044 | 0.685169177290653 |
| GSE113883 | TAms_mort | darkgrey | CMP_TAms_mort_darkgrey | 11 | 153 | 0.000713965944546712 | 0.00356007715542794 | 15.2942140319912 | 0.769065322219441 | 0.00405174142361658 | 0.176212119654325 | 0.769065322219441 |
| GSE113883 | TAms_mort | darkorange | CMP_TAms_mort_darkorange | 11 | 153 | 0.000654816804172815 | 0.0265258266938334 | 4.68589564712116 | 0.152109515630903 | 0.00140270370867264 | 0.466824747182326 | 0.182531418757083 |
| GSE113883 | TAms_mort | darkred | CMP_TAms_mort_darkred | 11 | 153 | 0.0126268193669204 | 0.237154586759835 | 14.5608059570561 | 3.95006739375895e-09 | 0.00780601440572888 | 1.61757571926251 | 9.48016174502148e-08 |
| GSE113883 | TAms_mort | green | CMP_TAms_mort_green | 11 | 153 | 0.00314274789970999 | 0.0154428149654439 | 14.0149420982083 | 0.320564721582854 | 0.00861906194581616 | 0.364627603266679 | 0.349706968999477 |
| GSE113883 | TAms_mort | greenyellow | CMP_TAms_mort_greenyellow | 11 | 153 | 0.0161528540751019 | 0.0820706118511145 | 3.94310460902604 | 0.00233299965027511 | 0.0187609587020545 | 0.860982337396917 | 0.00466599930055021 |
| GSE113883 | TAms_mort | grey60 | CMP_TAms_mort_grey60 | 11 | 153 | 0.008567822881796 | 0.154047973090427 | 1.22310276370287 | 6.58355129617387e-06 | 0.00708959166494144 | 1.20850724367731 | 3.95013077770432e-05 |
| GSE113883 | TAms_mort | lightgreen | CMP_TAms_mort_lightgreen | 11 | 153 | 0.0119131362113981 | 0.124283681740643 | 14.4944913242334 | 7.31856736361922e-05 | 0.0108867189545389 | 1.09428159771051 | 0.000292742694544769 |
| GSE113883 | TAms_mort | lightyellow | CMP_TAms_mort_lightyellow | 11 | 153 | 0.00292913154685448 | 0.116564095474992 | 14.297277528273 | 0.000147524201118047 | 0.00279248957333188 | 1.04893195477882 | 0.00044257260335414 |
| GSE113883 | TAms_mort | magenta | CMP_TAms_mort_magenta | 11 | 153 | 0.00404347687142 | 0.0864865548167676 | 13.9322339553813 | 0.001583297094354 | 0.00457128330547321 | 0.884538673544633 | 0.00362568247603562 |
| GSE113883 | TAms_mort | midnightblue | CMP_TAms_mort_midnightblue | 11 | 153 | 0.00197255004492776 | 0.120113028607869 | 13.9036780354862 | 0.000133894459928286 | 0.0018768173057816 | 1.05100802238516 | 0.00044257260335414 |
| GSE113883 | TAms_mort | orange | CMP_TAms_mort_orange | 11 | 153 | 0.00376540072470335 | 0.171200374984839 | 21.1655553261607 | 1.66222942699139e-06 | 0.00317793313930011 | 1.18485838425557 | 1.32978354159311e-05 |
| GSE113883 | TAms_mort | outside | CMP_TAms_mort_outside | 11 | 153 | 0.0757424094937537 | 0.111065104733955 | 8.13039494852243 | 0.00022854370077958 | 0.0815348294772949 | 0.928957722476697 | 0.000609449868745547 |
| GSE113883 | TAms_mort | pink | CMP_TAms_mort_pink | 11 | 153 | 0.0061695144917665 | 0.0330900739715059 | 20.9264739942792 | 0.0914731928017524 | 0.0127741875846652 | 0.482967269024037 | 0.115545085644319 |
| GSE113883 | TAms_mort | red | CMP_TAms_mort_red | 11 | 153 | 0.00252331573424274 | 0.0428348963836099 | 20.4884484002574 | 0.0451196830698779 | 0.00456745484950153 | 0.552455539766993 | 0.0601595774265038 |
| GSE113883 | TAms_mort | salmon | CMP_TAms_mort_salmon | 11 | 153 | 0.00242970021760854 | 0.0699471466957966 | 16.7932299360392 | 0.00601852409727282 | 0.0031302633235597 | 0.776196749749956 | 0.00962963855563651 |
| GSE113883 | TAms_mort | skyblue | CMP_TAms_mort_skyblue | 11 | 153 | 0.00225236069271712 | 0.0767874125658462 | 9.97548166942692 | 0.00307669059807747 | 0.00295357715170901 | 0.762587390484739 | 0.00553614597966166 |
| GSE113883 | TAms_mort | tan | CMP_TAms_mort_tan | 11 | 153 | 0.00283673936694391 | 0.0424939914907068 | 8.49901194071154 | 0.0449791855441583 | 0.00513581176799237 | 0.552344886279352 | 0.0601595774265038 |
| GSE113883 | TAms_mort | turquoise | CMP_TAms_mort_turquoise | 11 | 153 | 0.0381879157893178 | 0.215655110275658 | 13.8564175468774 | 3.24144361130252e-08 | 0.0254113333151108 | 1.5027907160861 | 3.88973233356303e-07 |
| GSE113883 | TAms_mort | yellow | CMP_TAms_mort_yellow | 11 | 153 | 0.00534046591063136 | 0.0645610235911281 | 12.7739149905998 | 0.00818671073149388 | 0.00722955989528177 | 0.738698619001235 | 0.0122800660972408 |

### 6.4 TAms_* H2
투영은 기존 8개 혈액 구성 점수를 넣은 REML 모형이다. 부트스트랩 2000회에서 두 REML이 함께 수렴한 수는 6개 행 모두 2000이다.

| dataset | clock | n_boot | n_converged_both_reml | p_med | p_med_p025 | p_med_p975 | two_a_orth | two_a_orth_p025 | two_a_orth_p975 | lrt_p_unadj |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | TAms_chrono | 2000 | 2000 | 0.604554553943655 | 0.360490943625388 | 0.881407966272375 | 0.000122429695499625 | 6.25500536820803e-05 | 0.00444575632084338 | 5.49745798912395e-20 |
| GSE39445 | TAms_mort | 2000 | 2000 | 0.648619163577675 | 0.328462743626483 | 1.00587179922653 | 0.00038278778315292 | 0.000261464410423898 | 0.0236878820150735 | 3.45851149483147e-15 |
| GSE56931 | TAms_chrono | 2000 | 2000 | 1.23230529037825 | -2.39490072930965 | 4.49887335706176 | 0.00578403104994379 | 0.000275022165738306 | 0.016229170326604 | 0.181165188487945 |
| GSE56931 | TAms_mort | 2000 | 2000 | 0.937587420102395 | 0.0854764204114879 | 2.35128746858658 | 0.0176994744478302 | 0.0010003640951316 | 0.0577097795356693 | 0.00673707816031792 |
| GSE113883 | TAms_chrono | 2000 | 2000 | 0.827434926978766 | 0.449125356663948 | 1.15576326319054 | 0.0107288701925842 | 0.00118910766937035 | 0.0151767679565033 | 2.96936394371969e-10 |
| GSE113883 | TAms_mort | 2000 | 2000 | 0.668143587639813 | 0.433397727889269 | 1.11458876791915 | 0.0305438078227955 | 0.00209175804837061 | 0.0732161132561262 | 4.44466652817908e-06 |

### 6.5 전처리 오류판과의 비교표
`TArd_*`는 보고서 4의 `TA_*`와, `TMrd_*`는 보고서 4의 `TM_*`와 같은 시계 파일이다. `r2_within`, `lrt_p`, `peak_h`는 이번 재채점이고 `stage4_r2`, `stage4_lrt_p`, `stage4_peak`는 보고서 4다.

| dataset | stage4_clock | clock | r2_within | lrt_p | peak_h | stage4_r2 | stage4_lrt_p | stage4_peak |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE113883 | TA_chrono | TArd_chrono | 0.0608651139287235 | 0.0116422321333177 | 15.8476561201441 | 0.0966791813140306 | 0.000751939790616891 | 15.667602590042 |
| GSE113883 | TA_mort | TArd_mort | 0.0092404532676319 | 0.515022143606158 | 1.01649318423095 | 0.00416979906927029 | 0.745144998027147 | 17.6460030950774 |
| GSE113883 | TM_black_chronological | TMrd_black_chronological | 0.0960376545899501 | 0.00075776232951682 | 10.3694799349063 | 0.0920649868165196 | 0.00104272227494776 | 9.43426408390714 |
| GSE113883 | TM_black_mortality | TMrd_black_mortality | 0.18744582962441 | 3.93365845618274e-07 | 11.335743620785 | 0.180051879317388 | 7.51473389146692e-07 | 10.9097923409538 |
| GSE113883 | TM_blue_chronological | TMrd_blue_chronological | 0.0705054273065931 | 0.00584760939885993 | 12.7086652463249 | 0.0949335321810294 | 0.000889833155799814 | 13.343483035576 |
| GSE113883 | TM_blue_mortality | TMrd_blue_mortality | 0.0106656153710409 | 0.474033941479901 | 14.5304209830112 | 0.0204800727249376 | 0.234320942124992 | 14.5296601202447 |
| GSE113883 | TM_brown_chronological | TMrd_brown_chronological | 0.0124314604842792 | 0.396745687238909 | 22.8276951816844 | 0.00824164069931765 | 0.54229698366573 | 20.9352127938297 |
| GSE113883 | TM_brown_mortality | TMrd_brown_mortality | 0.0430380467254647 | 0.044643480033738 | 13.5821394745719 | 0.0601982394593815 | 0.0123276056507036 | 13.752169203839 |
| GSE113883 | TM_cyan_chronological | TMrd_cyan_chronological | 0.107228306262268 | 0.00025756746016359 | 15.4482379729396 | 0.123362172471242 | 6.94656665720149e-05 | 15.4118534698436 |
| GSE113883 | TM_cyan_mortality | TMrd_cyan_mortality | 0.0923144421906369 | 0.000956515682018712 | 16.2643140842568 | 0.12634953899778 | 7.2150772242239e-05 | 15.9715619563109 |
| GSE113883 | TM_darkgreen_chronological | TMrd_darkgreen_chronological | 0.0187005438816327 | 0.26072454009064 | 13.3673188016985 | 0.0443816487496821 | 0.0394479954252844 | 13.8172747297159 |
| GSE113883 | TM_darkgreen_mortality | TMrd_darkgreen_mortality | 0.018327345322235 | 0.270465620682247 | 5.69382504888998 | 0.0113963005443716 | 0.448160024993602 | 9.01603504720726 |
| GSE113883 | TM_darkgrey_chronological | TMrd_darkgrey_chronological | 0.19185268416746 | 3.19313381359815e-07 | 13.80658498837 | 0.193471250411986 | 2.6468084551431e-07 | 13.7807786019332 |
| GSE113883 | TM_darkgrey_mortality | TMrd_darkgrey_mortality | 0.152911101944683 | 8.2961132347846e-06 | 14.0311008024758 | 0.165501059230859 | 2.81529071490669e-06 | 14.0262131726109 |
| GSE113883 | TM_darkorange_chronological | TMrd_darkorange_chronological | 0.195095130666254 | 1.38463688522557e-07 | 12.808509553616 | 0.202292621573607 | 7.38058900210152e-08 | 13.0407547772519 |
| GSE113883 | TM_darkorange_mortality | TMrd_darkorange_mortality | 0.0888787415927354 | 0.00112623648161785 | 11.6867217890285 | 0.0974056889863783 | 0.000606000741228307 | 12.5957894777605 |
| GSE113883 | TM_darkred_chronological | TMrd_darkred_chronological | 0.0747419615657492 | 0.00389104967457445 | 12.8492271985134 | 0.0489837158397306 | 0.0281393678234745 | 11.6969906951804 |
| GSE113883 | TM_darkred_mortality | TMrd_darkred_mortality | 0.0764691879237553 | 0.00340774232575982 | 12.5759528837204 | 0.0593356600743192 | 0.0129193968690727 | 11.7418291295997 |
| GSE113883 | TM_green_chronological | TMrd_green_chronological | 0.0238842191202923 | 0.197439286916834 | 12.7449981885838 | 0.0283018515741664 | 0.140448938560477 | 12.7952611234869 |
| GSE113883 | TM_green_mortality | TMrd_green_mortality | 0.00550180865425376 | 0.698550816118086 | 9.87268133497187 | 0.00727744333767699 | 0.614711725512297 | 10.1493538043249 |
| GSE113883 | TM_greenyellow_chronological | TMrd_greenyellow_chronological | 0.0895400114743112 | 0.00126931088042712 | 11.2624919286433 | 0.0848094442133197 | 0.00184236749904285 | 10.8979759962828 |
| GSE113883 | TM_greenyellow_mortality | TMrd_greenyellow_mortality | 0.204147543594578 | 8.88928727100064e-08 | 11.4839256146908 | 0.206825201478375 | 7.00536210296014e-08 | 11.5406136291133 |
| GSE113883 | TM_grey60_chronological | TMrd_grey60_chronological | 0.108668091426533 | 0.000266615604969793 | 13.7163187881964 | 0.114840502179783 | 0.000154386719262304 | 13.6397119948119 |
| GSE113883 | TM_grey60_mortality | TMrd_grey60_mortality | 0.00958115524821091 | 0.514803201820023 | 23.5298143909164 | 0.0186371248605184 | 0.274871285203076 | 1.253999661326 |
| GSE113883 | TM_lightgreen_chronological | TMrd_lightgreen_chronological | 0.0222275410703552 | 0.202984485152965 | 21.6783265740308 | 0.0217209816336884 | 0.211500951209 | 0.610438851335317 |
| GSE113883 | TM_lightgreen_mortality | TMrd_lightgreen_mortality | 0.0211830489291146 | 0.218288463366695 | 21.8918729680255 | 0.0254412903465814 | 0.161792310366885 | 1.00127115461664 |
| GSE113883 | TM_lightyellow_chronological | TMrd_lightyellow_chronological | 0.0183441445452177 | 0.258505002136834 | 11.8469537262561 | 0.0360141107383472 | 0.0613543948049248 | 13.4227516081739 |
| GSE113883 | TM_lightyellow_mortality | TMrd_lightyellow_mortality | 0.00875529158587762 | 0.53288795183278 | 9.67816622136416 | 0.0198339626171844 | 0.206866312809622 | 13.2076443735556 |
| GSE113883 | TM_magenta_chronological | TMrd_magenta_chronological | 0.150994282935349 | 7.78269559182711e-06 | 3.17837769638435 | 0.126198905064618 | 6.18663055924296e-05 | 3.24652801917783 |
| GSE113883 | TM_magenta_mortality | TMrd_magenta_mortality | 0.0989578582224564 | 0.000509218507569465 | 3.642235751741 | 0.0664444723016787 | 0.00675632441920848 | 4.0424626850677 |
| GSE113883 | TM_midnightblue_chronological | TMrd_midnightblue_chronological | 0.0208721362672132 | 0.232273274054151 | 23.2935360808353 | 0.0146870597153587 | 0.359631823335839 | 21.4272941407535 |
| GSE113883 | TM_midnightblue_mortality | TMrd_midnightblue_mortality | 0.0522675494297801 | 0.0222457238648075 | 20.1828963870539 | 0.0562925968431228 | 0.0163194315335967 | 19.2183714558209 |
| GSE113883 | TM_orange_chronological | TMrd_orange_chronological | 0.239553228606538 | 3.3154463060905e-09 | 14.1338847766488 | 0.237566290178325 | 4.02441765520324e-09 | 14.1583138939012 |
| GSE113883 | TM_orange_mortality | TMrd_orange_mortality | 0.148293212143321 | 9.49892551511532e-06 | 14.0958629166956 | 0.135688293041572 | 2.8459735767537e-05 | 13.9279380615639 |
| GSE113883 | TM_pink_chronological | TMrd_pink_chronological | 0.295072247159039 | 1.88604117468432e-11 | 12.1916860039504 | 0.212319385973019 | 4.80387011478345e-08 | 11.4550545518287 |
| GSE113883 | TM_pink_mortality | TMrd_pink_mortality | 0.038809426888497 | 0.0616201739984647 | 8.86911046448394 | 0.0399901537273766 | 0.0553988457722534 | 7.39782198995527 |
| GSE113883 | TM_red_chronological | TMrd_red_chronological | 0.0687806224199643 | 0.00653517123436425 | 1.04220719640953 | 0.0763074511089803 | 0.00377567349522923 | 1.70607891543789 |
| GSE113883 | TM_red_mortality | TMrd_red_mortality | 0.00705726510049709 | 0.610185047132444 | 0.909265338487589 | 0.0107498435749711 | 0.469208443664722 | 2.85945718876094 |
| GSE113883 | TM_salmon_chronological | TMrd_salmon_chronological | 0.130446757417008 | 5.48603853733598e-05 | 17.1046499569921 | 0.0377226599091116 | 0.0655159201986105 | 17.9390325766023 |
| GSE113883 | TM_salmon_mortality | TMrd_salmon_mortality | 0.134312951689649 | 4.00930940838677e-05 | 15.1902682755822 | 0.0421977270777493 | 0.0471365433288297 | 14.7114513603138 |
| GSE113883 | TM_skyblue_chronological | TMrd_skyblue_chronological | 0.0246334114865252 | 0.161910488712392 | 10.4424217521067 | 0.0282175443143632 | 0.122490396691927 | 11.5618933160481 |
| GSE113883 | TM_skyblue_mortality | TMrd_skyblue_mortality | 0.0126506482918325 | 0.403598902646806 | 5.39017084704905 | 0.0204953298620161 | 0.230454562251633 | 5.24138228441545 |
| GSE113883 | TM_tan_chronological | TMrd_tan_chronological | 0.0240137821060807 | 0.177422872310181 | 12.148541586731 | 0.0372348524161215 | 0.0669911700771828 | 12.752035911046 |
| GSE113883 | TM_tan_mortality | TMrd_tan_mortality | 0.0454012595663962 | 0.0364999832211194 | 10.6950213526642 | 0.050643987897681 | 0.0246558077612224 | 11.0836100014759 |
| GSE113883 | TM_turquoise_chronological | TMrd_turquoise_chronological | 0.0361397552841218 | 0.0743883535245499 | 0.678019413592925 | 0.0642818018652043 | 0.00909627741769173 | 1.13969806986858 |
| GSE113883 | TM_turquoise_mortality | TMrd_turquoise_mortality | 0.0364471774041928 | 0.0717570281955081 | 13.2353766671555 | 0.01958278310867 | 0.244652824553788 | 8.39017465316885 |
| GSE113883 | TM_yellow_chronological | TMrd_yellow_chronological | 0.0799254501719494 | 0.00275299651503155 | 15.8227551886139 | 0.110856561888506 | 0.000240544783693973 | 15.5766940777842 |
| GSE113883 | TM_yellow_mortality | TMrd_yellow_mortality | 0.0337861647728142 | 0.0897691248823284 | 18.7689498889017 | 0.048984193053558 | 0.0291533549700461 | 17.4615677991321 |
| GSE39445 | TA_chrono | TArd_chrono | 0.0772645457809408 | 0.000381339040015034 | 14.4302814414026 | 0.0818707783732363 | 0.000235567727394529 | 15.3360608213662 |
| GSE39445 | TA_mort | TArd_mort | 0.0689139005346088 | 0.000830944712054089 | 14.4256736372303 | 0.0719801227245767 | 0.000615890141275604 | 15.5538826108025 |
| GSE39445 | TM_black_chronological | TMrd_black_chronological | 0.00965626708127676 | 0.383212948552562 | 16.7465478529784 | 0.00419662297051537 | 0.649316597852252 | 15.5430089372762 |
| GSE39445 | TM_black_mortality | TMrd_black_mortality | 0.021405920372352 | 0.118456944664595 | 14.6913085707251 | 0.0334487137985466 | 0.0347355763656406 | 15.2074353819786 |
| GSE39445 | TM_blue_chronological | TMrd_blue_chronological | 0.273212935390491 | 3.23194334314227e-14 | 12.6246949272105 | 0.237829014095836 | 3.15097496339781e-12 | 12.8858696782194 |
| GSE39445 | TM_blue_mortality | TMrd_blue_mortality | 0.194734730342118 | 7.85827959476316e-10 | 13.0254761199196 | 0.173384151988194 | 9.86053798875224e-09 | 13.3350666192169 |
| GSE39445 | TM_brown_chronological | TMrd_brown_chronological | 0.130158699758699 | 9.86993595041124e-07 | 2.51855625767681 | 0.120012185415971 | 3.1110076635881e-06 | 2.41251767514426 |
| GSE39445 | TM_brown_mortality | TMrd_brown_mortality | 0.0159206975938524 | 0.226673497507782 | 11.9237484132225 | 0.0180956475419376 | 0.179471889493929 | 12.7854394600048 |
| GSE39445 | TM_cyan_chronological | TMrd_cyan_chronological | 0.08779457870228 | 0.000112459224737045 | 16.9480659090376 | 0.0784202942552863 | 0.000310206618827339 | 17.0881605720823 |
| GSE39445 | TM_cyan_mortality | TMrd_cyan_mortality | 0.0539058662302986 | 0.00395408060153505 | 14.6031991885237 | 0.0387155662466691 | 0.0194919313586336 | 14.963377925844 |
| GSE39445 | TM_darkgreen_chronological | TMrd_darkgreen_chronological | 0.246852248319311 | 1.01368160741315e-12 | 13.5793636871899 | 0.0438650714444653 | 0.0141317219871215 | 12.2297235510656 |
| GSE39445 | TM_darkgreen_mortality | TMrd_darkgreen_mortality | 0.155736246692013 | 6.99005187392206e-08 | 13.5566918472623 | 0.0253755609023505 | 0.085967523581588 | 12.2070204131449 |
| GSE39445 | TM_darkgrey_chronological | TMrd_darkgrey_chronological | 0.0224804469401431 | 0.104245358231263 | 13.3620485646281 | 0.0228754471660846 | 0.0995629987206397 | 13.4425874433925 |
| GSE39445 | TM_darkgrey_mortality | TMrd_darkgrey_mortality | 0.0174130925431409 | 0.184807247099146 | 10.3613547975916 | 0.0157714699716435 | 0.215049197736021 | 10.5754611115216 |
| GSE39445 | TM_darkorange_chronological | TMrd_darkorange_chronological | 0.290858953865026 | 3.56018412547148e-15 | 13.2581749430489 | 0.252125912766848 | 6.43156705950182e-13 | 13.1994196787556 |
| GSE39445 | TM_darkorange_mortality | TMrd_darkorange_mortality | 0.244474001877664 | 1.49364628139268e-12 | 12.5753919828086 | 0.222732128205742 | 2.5440199674316e-11 | 11.990507591373 |
| GSE39445 | TM_darkred_chronological | TMrd_darkred_chronological | 0.15377771007466 | 7.73579797817402e-08 | 14.5335008844385 | 0.170268168064267 | 1.11038075040807e-08 | 14.6035753971536 |
| GSE39445 | TM_darkred_mortality | TMrd_darkred_mortality | 0.119837772739331 | 3.98454737587873e-06 | 13.5735536410406 | 0.126281955122939 | 1.93892619480763e-06 | 13.6277744098779 |
| GSE39445 | TM_green_chronological | TMrd_green_chronological | 0.124601943077417 | 1.73763791764011e-06 | 1.51678056266431 | 0.0142119820256309 | 0.192140166425451 | 23.1118641649033 |
| GSE39445 | TM_green_mortality | TMrd_green_mortality | 0.137097045603345 | 5.60568666194975e-07 | 1.36384408508586 | 0.0169112079772625 | 0.156286236015497 | 21.0399406078005 |
| GSE39445 | TM_greenyellow_chronological | TMrd_greenyellow_chronological | 0.260134461342361 | 1.65875077718807e-13 | 14.1154000919707 | 0.0916708400447678 | 8.55966574357814e-05 | 14.0126566688726 |
| GSE39445 | TM_greenyellow_mortality | TMrd_greenyellow_mortality | 0.308729757324834 | 2.25652006599657e-16 | 14.3958370908396 | 0.121801675571493 | 3.18247726036681e-06 | 14.1088658354686 |
| GSE39445 | TM_grey60_chronological | TMrd_grey60_chronological | 0.0495578888919218 | 0.00771377440696271 | 1.85725060675388 | 0.0452800114168385 | 0.0118589267499855 | 1.61122714677413 |
| GSE39445 | TM_grey60_mortality | TMrd_grey60_mortality | 0.0950034675336462 | 6.34923891838673e-05 | 3.43382687592683 | 0.0884199606988083 | 0.000127459660491946 | 3.47713463219173 |
| GSE39445 | TM_lightgreen_chronological | TMrd_lightgreen_chronological | 0.00458956507932473 | 0.655412947059859 | 17.6104018194813 | 0.0113722292224852 | 0.337849225785314 | 19.5440607248531 |
| GSE39445 | TM_lightgreen_mortality | TMrd_lightgreen_mortality | 0.000617517341806595 | 0.942535203050797 | 10.9040970066057 | 0.00200577581082695 | 0.834398256894445 | 19.6939855837568 |
| GSE39445 | TM_lightyellow_chronological | TMrd_lightyellow_chronological | 0.0196186467178187 | 0.141984800951545 | 16.0083362409656 | 0.016958944730711 | 0.188011532026613 | 16.3703421898219 |
| GSE39445 | TM_lightyellow_mortality | TMrd_lightyellow_mortality | 0.00473031699076632 | 0.618656931830644 | 15.4028342604518 | 0.00352594888358004 | 0.703834950188363 | 16.0614428390137 |
| GSE39445 | TM_magenta_chronological | TMrd_magenta_chronological | 0.0182057018172028 | 0.165107255559739 | 5.01859088316123 | 0.031841709089661 | 0.0418042393016934 | 4.83609712470271 |
| GSE39445 | TM_magenta_mortality | TMrd_magenta_mortality | 0.0682755373164626 | 0.0010159654280022 | 3.8393171644644 | 0.0169441690686331 | 0.191846089160145 | 4.65555674736115 |
| GSE39445 | TM_midnightblue_chronological | TMrd_midnightblue_chronological | 0.035434862404592 | 0.0283381556827208 | 2.68889503736959 | 0.0360812244223455 | 0.0263980738213452 | 2.48971135004274 |
| GSE39445 | TM_midnightblue_mortality | TMrd_midnightblue_mortality | 0.0147086914393125 | 0.232566265257363 | 2.51203106403811 | 0.0143688103662002 | 0.240032836220979 | 2.00118991393135 |
| GSE39445 | TM_orange_chronological | TMrd_orange_chronological | 0.151154578410192 | 1.1587071471288e-07 | 13.7026792528042 | 0.0623826475609638 | 0.00194298914813481 | 15.5211451674636 |
| GSE39445 | TM_orange_mortality | TMrd_orange_mortality | 0.20621999840851 | 1.76028532228109e-10 | 14.2350078503596 | 0.139130759332057 | 5.11987945961834e-07 | 14.9948982783575 |
| GSE39445 | TM_pink_chronological | TMrd_pink_chronological | 0.196839332958786 | 3.55422012948682e-10 | 14.1201847015606 | 0.202568746239115 | 1.81628517628059e-10 | 13.8622823951712 |
| GSE39445 | TM_pink_mortality | TMrd_pink_mortality | 0.0686632536438139 | 0.000777434866548605 | 15.2921982325918 | 0.0863849327369744 | 0.00011647346563911 | 14.6593137052705 |
| GSE39445 | TM_red_chronological | TMrd_red_chronological | 0.0231123111712674 | 0.112388447528015 | 4.107008351832 | 0.0166988451905946 | 0.210605494933635 | 4.09149417114038 |
| GSE39445 | TM_red_mortality | TMrd_red_mortality | 0.0301233524737148 | 0.049589071780739 | 3.39850642987535 | 0.0214812002813103 | 0.118609431495671 | 3.15560654057404 |
| GSE39445 | TM_salmon_chronological | TMrd_salmon_chronological | 0.0131173168707371 | 0.278747908971104 | 15.0328830600364 | 0.0207582463074085 | 0.12989270482601 | 14.2862780888814 |
| GSE39445 | TM_salmon_mortality | TMrd_salmon_mortality | 0.0458513358488147 | 0.0104320800213753 | 14.9376126728354 | 0.0535088229874735 | 0.00468818860828374 | 14.3737687655268 |
| GSE39445 | TM_skyblue_chronological | TMrd_skyblue_chronological | 0.0420096720151547 | 0.0157783945728216 | 2.3855767682922 | 0.0367524785410869 | 0.0269056225941273 | 2.12416317227273 |
| GSE39445 | TM_skyblue_mortality | TMrd_skyblue_mortality | 0.0110221031447092 | 0.3418090390457 | 2.29399160853936 | 0.00897552098664223 | 0.418185721734746 | 2.03710442275896 |
| GSE39445 | TM_tan_chronological | TMrd_tan_chronological | 0.0156356520787937 | 0.208937213638317 | 16.209114706835 | 0.0197098625397636 | 0.139182162275866 | 16.4114414075102 |
| GSE39445 | TM_tan_mortality | TMrd_tan_mortality | 0.000715478030452355 | 0.94218317521691 | 2.63750073093247 | 0.00106264679483215 | 0.912885926819834 | 1.85297786563406 |
| GSE39445 | TM_turquoise_chronological | TMrd_turquoise_chronological | 0.039993461160248 | 0.0201307917526974 | 12.633562450264 | 0.0482271919172995 | 0.00919790868805224 | 11.6135068062167 |
| GSE39445 | TM_turquoise_mortality | TMrd_turquoise_mortality | 0.130347538836632 | 1.34568663538648e-06 | 13.4078416799291 | 0.144979647914123 | 2.59158198342643e-07 | 13.079788251819 |
| GSE39445 | TM_yellow_chronological | TMrd_yellow_chronological | 0.121453391777851 | 2.94019724233245e-06 | 4.20967119862452 | 0.0788763500807268 | 0.00030816196430124 | 3.98006739489746 |
| GSE39445 | TM_yellow_mortality | TMrd_yellow_mortality | 0.0347001190696372 | 0.0314132006327156 | 4.95910402269568 | 0.00586286783638201 | 0.580147184298968 | 6.16149376718439 |
| GSE56931 | TA_chrono | TArd_chrono | 0.0538147539888595 | 0.0408948152699152 | 23.6764604105288 | 0.0589265110850511 | 0.029920138819241 | 0.827909618359159 |
| GSE56931 | TA_mort | TArd_mort | 0.0180930649584713 | 0.369352805274575 | 23.8983281147265 | 0.032665947996622 | 0.157968351876906 | 1.28654792915375 |
| GSE56931 | TM_black_chronological | TMrd_black_chronological | 0.00808850351420909 | 0.623350364950938 | 18.2859272060427 | 0.00689307540415207 | 0.675530523241086 | 18.8430210798695 |
| GSE56931 | TM_black_mortality | TMrd_black_mortality | 0.0183850363172609 | 0.338209526108064 | 13.3407508725362 | 0.0172923543267151 | 0.364508273835368 | 15.0010467690942 |
| GSE56931 | TM_blue_chronological | TMrd_blue_chronological | 0.0878035246902666 | 0.00420264739126932 | 8.44292735140082 | 0.0835709307842948 | 0.00551597823577023 | 8.49374611650761 |
| GSE56931 | TM_blue_mortality | TMrd_blue_mortality | 0.0577995847834829 | 0.0292200551611706 | 6.28425639397042 | 0.0556424862973739 | 0.0333910809856194 | 6.38160062769946 |
| GSE56931 | TM_brown_chronological | TMrd_brown_chronological | 0.109738119793836 | 0.00124431319146678 | 23.3344947162811 | 0.103992356891917 | 0.00180315277649665 | 23.3615266605128 |
| GSE56931 | TM_brown_mortality | TMrd_brown_mortality | 0.077212694616785 | 0.00971900812660163 | 15.3478728541987 | 0.0779639907500801 | 0.00936479530210405 | 15.3192656248541 |
| GSE56931 | TM_cyan_chronological | TMrd_cyan_chronological | 0.0941383073965043 | 0.00341274003243299 | 18.9978533208491 | 0.0931102048235051 | 0.00367203908014258 | 18.4887569220951 |
| GSE56931 | TM_cyan_mortality | TMrd_cyan_mortality | 0.0557565936856742 | 0.0396270306607154 | 18.2345664587987 | 0.0530021929012933 | 0.046891518458832 | 17.4680181184238 |
| GSE56931 | TM_darkgreen_chronological | TMrd_darkgreen_chronological | 0.0283435950940413 | 0.183465554279949 | 9.40156128921537 | 0.0467864370573649 | 0.0613364693244821 | 5.47809913544504 |
| GSE56931 | TM_darkgreen_mortality | TMrd_darkgreen_mortality | 0.0282687188027784 | 0.196795909791021 | 4.26950285952008 | 0.0453470709027217 | 0.0713540512810864 | 4.53105837810403 |
| GSE56931 | TM_darkgrey_chronological | TMrd_darkgrey_chronological | 0.000687888287785987 | 0.960573996289304 | 1.49527285708347 | 0.000406752578562819 | 0.975456154614006 | 0.151978490634241 |
| GSE56931 | TM_darkgrey_mortality | TMrd_darkgrey_mortality | 0.0124105956399536 | 0.473303852911435 | 3.73546250266167 | 0.0107799719222351 | 0.522396032647298 | 3.67332251022731 |
| GSE56931 | TM_darkorange_chronological | TMrd_darkorange_chronological | 0.218786825450032 | 7.31578430610281e-07 | 12.9431586792082 | 0.223188558409583 | 5.24850078796061e-07 | 13.0257575359056 |
| GSE56931 | TM_darkorange_mortality | TMrd_darkorange_mortality | 0.0786312063175245 | 0.00965498315408789 | 13.7948056715957 | 0.0895951442676672 | 0.00477623855766983 | 13.9271338124809 |
| GSE56931 | TM_darkred_chronological | TMrd_darkred_chronological | 0.0484550837159842 | 0.0668219896129173 | 13.0107335212523 | 0.0523315873682693 | 0.0517382344463975 | 12.9319319521483 |
| GSE56931 | TM_darkred_mortality | TMrd_darkred_mortality | 0.104714543259832 | 0.00166188175754768 | 11.7888486583853 | 0.108126184963575 | 0.00134898352061123 | 11.7921623328549 |
| GSE56931 | TM_green_chronological | TMrd_green_chronological | 0.0232660231492419 | 0.243603405861654 | 7.61136591764577 | 0.0225949280167521 | 0.253978923752743 | 8.61391421556577 |
| GSE56931 | TM_green_mortality | TMrd_green_mortality | 0.0139901064272934 | 0.429643706131985 | 6.17764788108991 | 0.0183216635047599 | 0.32919155148734 | 7.43524324525559 |
| GSE56931 | TM_greenyellow_chronological | TMrd_greenyellow_chronological | 0.132268423624651 | 0.000293238583900257 | 8.41039799994602 | 0.0900155127962415 | 0.00454601902938541 | 8.25986430634499 |
| GSE56931 | TM_greenyellow_mortality | TMrd_greenyellow_mortality | 0.216568091909472 | 8.30842365278521e-07 | 9.10693438096604 | 0.103721383650612 | 0.00215346916179761 | 8.77025120655723 |
| GSE56931 | TM_grey60_chronological | TMrd_grey60_chronological | 0.0310837823301093 | 0.187078945444375 | 14.735358823595 | 0.0320803366959598 | 0.176594826357309 | 14.7965669723146 |
| GSE56931 | TM_grey60_mortality | TMrd_grey60_mortality | 0.0608266883249575 | 0.0296509923515952 | 16.4898330720896 | 0.0576479338343247 | 0.0358886599528894 | 16.3703934608109 |
| GSE56931 | TM_lightgreen_chronological | TMrd_lightgreen_chronological | 0.0285764426987187 | 0.189622886681096 | 12.5676484687013 | 0.0282868105573336 | 0.192796504568752 | 12.2750491530141 |
| GSE56931 | TM_lightgreen_mortality | TMrd_lightgreen_mortality | 0.0427669073533289 | 0.0803044557139391 | 11.1388157049831 | 0.046798175391677 | 0.0628541571264112 | 10.7480267161231 |
| GSE56931 | TM_lightyellow_chronological | TMrd_lightyellow_chronological | 0.109121342528465 | 0.00112929360414981 | 5.2015047540198 | 0.0919596079778938 | 0.00343578887956388 | 5.21997300629944 |
| GSE56931 | TM_lightyellow_mortality | TMrd_lightyellow_mortality | 0.139500890953261 | 0.000163690329178288 | 5.16663462198874 | 0.119040291733443 | 0.000633691830898987 | 5.18623631730539 |
| GSE56931 | TM_magenta_chronological | TMrd_magenta_chronological | 0.117420333748274 | 0.000727118254561163 | 5.894872553388 | 0.0180737310695097 | 0.339025317633109 | 9.36134792638769 |
| GSE56931 | TM_magenta_mortality | TMrd_magenta_mortality | 0.123398322636573 | 0.000476251480956685 | 7.78282950239244 | 0.0462426563030093 | 0.0599008637003445 | 10.1520816986573 |
| GSE56931 | TM_midnightblue_chronological | TMrd_midnightblue_chronological | 0.0789888276112155 | 0.00892494745813125 | 2.33734142799907 | 0.0712450857755046 | 0.0143823762637541 | 2.31207348973244 |
| GSE56931 | TM_midnightblue_mortality | TMrd_midnightblue_mortality | 0.0256938136813356 | 0.224216853342828 | 2.35235061939007 | 0.0228218027858318 | 0.263429132501458 | 2.10206579880247 |
| GSE56931 | TM_orange_chronological | TMrd_orange_chronological | 0.243896724649109 | 7.63175475924698e-08 | 15.0332126419389 | 0.206661212076165 | 1.59780340143946e-06 | 15.2637127251035 |
| GSE56931 | TM_orange_mortality | TMrd_orange_mortality | 0.0206741031439754 | 0.286370897310261 | 14.5577382283612 | 0.0356302204809017 | 0.123510868170419 | 15.0857874262409 |
| GSE56931 | TM_pink_chronological | TMrd_pink_chronological | 0.114313821575229 | 0.00080086481809504 | 10.6609779874152 | 0.126342340445824 | 0.00036335103065202 | 11.1199417302676 |
| GSE56931 | TM_pink_mortality | TMrd_pink_mortality | 0.152421763986453 | 6.2069672108746e-05 | 11.6144088757819 | 0.17016687754802 | 1.82336978045962e-05 | 11.6840808585176 |
| GSE56931 | TM_red_chronological | TMrd_red_chronological | 0.0641390993113912 | 0.0201410751245809 | 22.0655872610843 | 0.0868918422466177 | 0.00480536104685813 | 22.3663160298919 |
| GSE56931 | TM_red_mortality | TMrd_red_mortality | 0.0487431445499853 | 0.0541508383756257 | 0.241835704794313 | 0.0695703562698264 | 0.0149398846306077 | 0.162919136687949 |
| GSE56931 | TM_salmon_chronological | TMrd_salmon_chronological | 0.0736177023211254 | 0.0104884563624888 | 19.3278639496425 | 0.0657790403304638 | 0.017234347310599 | 19.4527771849152 |
| GSE56931 | TM_salmon_mortality | TMrd_salmon_mortality | 0.025813056214241 | 0.206165723790726 | 18.9096264515715 | 0.0213253750617666 | 0.270314337390179 | 19.0996621216644 |
| GSE56931 | TM_skyblue_chronological | TMrd_skyblue_chronological | 0.265059116633377 | 1.74040884746183e-08 | 4.12110817901442 | 0.166661129139752 | 2.50696964183316e-05 | 4.23405108385446 |
| GSE56931 | TM_skyblue_mortality | TMrd_skyblue_mortality | 0.13033004138163 | 0.000284024723587035 | 4.4786857829388 | 0.0866878748408164 | 0.0049047617602511 | 4.64568014131419 |
| GSE56931 | TM_tan_chronological | TMrd_tan_chronological | 0.00869784693183284 | 0.588651263000799 | 4.6873810150619 | 0.00644165551802955 | 0.673349937764662 | 4.62056644992187 |
| GSE56931 | TM_tan_mortality | TMrd_tan_mortality | 0.0695438901161031 | 0.0159271223539135 | 2.61329547743164 | 0.0592320049884755 | 0.0298849040758074 | 2.51201116672037 |
| GSE56931 | TM_turquoise_chronological | TMrd_turquoise_chronological | 0.0786075729524733 | 0.0090410631402461 | 3.38579486814741 | 0.091941341024543 | 0.00388621528189854 | 3.51438834939511 |
| GSE56931 | TM_turquoise_mortality | TMrd_turquoise_mortality | 0.0897575591067797 | 0.00419034037330313 | 8.67491640119687 | 0.0864204056444807 | 0.00515054296721403 | 8.61329057080909 |
| GSE56931 | TM_yellow_chronological | TMrd_yellow_chronological | 0.0361118252698642 | 0.116533953246999 | 7.27712786073087 | 0.0330910772045983 | 0.139519326354172 | 7.42135124811049 |
| GSE56931 | TM_yellow_mortality | TMrd_yellow_mortality | 0.0777363447864636 | 0.00916088388292554 | 6.844138751368 | 0.0671274627141432 | 0.0177643086151436 | 6.89325179690463 |

## 7. 사전등록 준비 계산
### 7.1 시계 요약 표
RA 세 시계는 저장된 lin 점수를 같은 주 모형으로 다시 적합해 σ와 σ_u를 냈다. 저장된 LRT p, 2A와의 절대차 최댓값은 각각 1.1529861093651e-09, 3.35997896172557e-12이다. σ_u는 귀무 모형 REML의 참가자 절편 SD이고, 그 SD의 길이가 1일 때만 기록했다.

| dataset | clock | lrt_p | peak_h | two_a | sigma_reml | two_a_over_sigma | sigma_u |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | 0.643663244864414 | 13.7137336093391 | 4.5622147045564 | 25.1603887545092 | 0.181325286706421 | 16.2536695399875 |
| GSE39445 | RA_Peters | 0.0822398090507864 | 6.84503499428439 | 3.86882060477993 | 9.12159876450608 | 0.424138432818846 | 14.613080804435 |
| GSE39445 | RA_default | 0.571779550797501 | 12.0064473932836 | 1.98498899424841 | 9.34019861854301 | 0.212521068910422 | 8.81869510215541 |
| GSE56931 | RA_GTEx | 0.0977660921594727 | 13.2671292133393 | 8.41295616809185 | 14.1878999046968 | 0.592966980638678 | 11.0606886507247 |
| GSE56931 | RA_Peters | 0.000105116530090887 | 4.54192615314062 | 6.37997357324785 | 5.58777765048951 | 1.14177298602584 | 13.9973116863226 |
| GSE56931 | RA_default | 0.0187399469939356 | 9.88701727265479 | 4.61694071236019 | 6.68803391836383 | 0.690328543293286 | 6.90600678789726 |
| GSE113883 | RA_GTEx | 0.104384139695734 | 6.16828794535857 | 4.34369415553422 | 8.81987355742966 | 0.492489390834316 | 10.4317015658613 |
| GSE113883 | RA_Peters | 0.0773443137299791 | 11.8362208570625 | 3.1551856815567 | 6.04743769494992 | 0.521739262264997 | 13.6399297750127 |
| GSE113883 | RA_default | 3.24675832907807e-05 | 11.812631972311 | 6.770165944749 | 6.28931435261543 | 1.07645532806507 | 6.86875908783177 |
| GSE39445 | TAms_chrono | 5.49745798912395e-20 | 14.2335397784943 | 0.0345210353298852 | 0.0171269287358502 | 2.01559986979017 | 0.0276429862743754 |
| GSE39445 | TAms_mort | 3.45851149483147e-15 | 13.4288594848278 | 0.148755555495313 | 0.0860831456543238 | 1.72804507043291 | 0.161562476919705 |
| GSE56931 | TAms_chrono | 0.181165188487945 | 11.4579688514277 | 0.00943066206760042 | 0.0197534498990921 | 0.477418482127208 | 0.0331124542924689 |
| GSE56931 | TAms_mort | 0.00673707816031792 | 12.0688456403715 | 0.0629468359700225 | 0.0740143176007128 | 0.850468368966173 | 0.177501208368076 |
| GSE113883 | TAms_chrono | 2.96936394371969e-10 | 12.188050742382 | 0.0418702301338288 | 0.0253117943835499 | 1.65417866072112 | 0.0416315613478055 |
| GSE113883 | TAms_mort | 4.44466652817908e-06 | 10.0600990696032 | 0.102562103238955 | 0.0897840656082817 | 1.14231965933045 | 0.184575482929547 |

| clock | phi_bar | max_circ_dist | dist_GSE39445 | dist_GSE56931 | dist_GSE113883 |
| --- | --- | --- | --- | --- | --- |
| RA_GTEx | 11.518486825228 | 5.35019887986941 | 2.19524678411108 | 1.74864238811129 | 5.35019887986941 |
| RA_Peters | 7.56604107260516 | 4.27017978445731 | 0.72100607832077 | 3.02411491946453 | 4.27017978445731 |
| RA_default | 11.2423863477089 | 1.35536907505411 | 0.764061045574719 | 1.35536907505411 | 0.570245624602082 |
| TAms_chrono | 12.616768998085 | 1.61677078040933 | 1.61677078040933 | 1.15880014665727 | 0.428718255702973 |
| TAms_mort | 11.8599407935094 | 1.79984172390618 | 1.56891869131842 | 0.208904846862117 | 1.79984172390618 |

### 7.2 재현 리듬 모듈 표
규칙: 세 자료 모두 BH q < 0.05이고, 세 정점의 원형 거리가 비가중 원형 평균 φ̄에서 3시간 이내(`circ_dist_0_12` <= 3)다. 대상은 `TM_*` 14개와 `CMP_*` 48개, 합 62개다.

| clock | family | outcome | module | parent_clock | phi_bar | max_circ_dist | q_GSE39445 | q_GSE56931 | q_GSE113883 | dist_GSE39445 | dist_GSE56931 | dist_GSE113883 | two_a_over_sigma_GSE39445 | two_a_over_sigma_GSE56931 | two_a_over_sigma_GSE113883 | mean_two_a_over_sigma | reproducible |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TM_blue_mortality | TM | mortality | blue |  | 13.3570949474744 | 0.554854233459569 | 0.515875649687256 | 0.238336833214312 | 4.18936772245646e-07 | 0.231408541231131 | 0.322019293781931 | 0.554854233459569 | 0.272149984433827 | 0.678368475396146 | 1.46777408737229 | 0.806097515734088 | FALSE |
| TM_brown4_mortality | TM | mortality | brown4 |  | 17.6285337838569 | 2.82720602026017 | 0.0478547263669363 | 0.774334487094162 | 2.30897788575168e-05 | 0.744161252383726 | 2.82720602026017 | 1.91591297437483 | 0.509261460347374 | 0.211643192718019 | 1.25248163265036 | 0.657795428571918 | FALSE |
| TM_darkgreen_mortality | TM | mortality | darkgreen |  | 23.1841646637323 | 3.02105081625064 | 0.00320224459164145 | 0.269911402435825 | 0.0125419349631913 | 1.55708393739222 | 1.22222647652482 | 3.02105081625064 | 0.750021575332171 | 0.506584240677654 | 0.77362256213373 | 0.676742792714518 | FALSE |
| TM_darkmagenta_mortality | TM | mortality | darkmagenta |  | 4.69281957981297 | 4.71719305376227 | 0.973851597832151 | 0.48804334419715 | 0.0485104327188859 | 0.0482488714590135 | 4.57802472320496 | 4.71719305376227 | 0.0828164830350462 | 0.36209705829897 | 0.61193520267432 | 0.352282914669445 | FALSE |
| TM_darkred_mortality | TM | mortality | darkred |  | 12.3363608436555 | 1.39849713740366 | 2.73109114311815e-08 | 0.774334487094162 | 0.719630378574929 | 1.39849713740366 | 0.562928625558238 | 0.812685936632638 | 1.26049130074904 | 0.257752383860435 | 0.202946942960936 | 0.573730209190137 | FALSE |
| TM_darkslateblue_mortality | TM | mortality | darkslateblue |  | 8.82611087508097 | 5.9532606428352 | 0.266596706307523 | 0.435304924885386 | 0.000817599229068192 | 0.0465215028698935 | 5.9532606428352 | 5.40141288907143 | 0.385422746245139 | 0.457072150478737 | 0.992709815142199 | 0.611734903955358 | FALSE |
| TM_green_mortality | TM | mortality | green |  | 23.8254217382045 | 9.19666239231668 | 0.984700668978151 | 0.435304924885386 | 0.719630378574929 | 1.20392629049858 | 5.23045361564688 | 9.19666239231668 | 0.0352151301073573 | 0.421261842392921 | 0.19719951221224 | 0.217892161570839 | FALSE |
| TM_ivory_mortality | TM | mortality | ivory |  | 9.1823459451001 | 3.3142819556713 | 0.00104797650360716 | 0.774334487094162 | 0.0113941290970536 | 3.3142819556713 | 1.31170903657711 | 1.6814977833409 | 0.815419137317899 | 0.24988903722974 | 0.72833564785796 | 0.5978812741352 | FALSE |
| TM_orange_mortality | TM | mortality | orange |  | 13.3995006852737 | 1.57749173759269 | 0.0324447672575371 | 0.250767680024482 | 0.000771917053591052 | 1.57749173759269 | 0.708724480495304 | 0.834999850585305 | 0.59364609609842 | 0.611504932368897 | 1.00674609613298 | 0.737299041533432 | FALSE |
| TM_pink_mortality | TM | mortality | pink |  | 18.9228014231657 | 4.2451170262697 | 0.373192698863889 | 0.250767680024482 | 0.164696476519653 | 2.9845255873975 | 4.2451170262697 | 0.738241730780999 | 0.340617209443192 | 0.606725661317401 | 0.462102971586567 | 0.469815280782387 | FALSE |
| TM_plum1_mortality | TM | mortality | plum1 |  | 8.53993333470402 | 8.91261567072968 | 0.515875649687256 | 0.000109294649854036 | 0.241720742168214 | 8.91261567072968 | 0.162427266298799 | 2.85929886322759 | 0.249285375217428 | 1.19089196231817 | 0.427003294246544 | 0.622393543927381 | FALSE |
| TM_sienna3_mortality | TM | mortality | sienna3 |  | 14.3494271689576 | 9.534742468761 | 0.657828317669594 | 0.945816074805526 | 0.000771917053591052 | 2.89472180640568 | 9.534742468761 | 0.328254994818819 | 0.197737729086497 | 0.0849001771137568 | 1.01726997117734 | 0.433302625792531 | FALSE |
| TM_turquoise_mortality | TM | mortality | turquoise |  | 11.5724011994242 | 1.04911072722202 | 1.7030021051141e-07 | 0.00657277863416067 | 3.693684980855e-06 | 1.04911072722202 | 0.415592752627278 | 0.623968315718878 | 1.21059791980705 | 0.974829283451104 | 1.26636172908367 | 1.15059631078061 | TRUE |
| TM_white_mortality | TM | mortality | white |  | 9.2016089095212 | 3.2279785723557 | 1.38284773445641e-08 | 0.867107279598642 | 0.517011619258595 | 3.2279785723557 | 1.87702406700578 | 1.06875188349482 | 1.34992861797921 | 0.154552133369572 | 0.282009933635222 | 0.595496894994668 | FALSE |
| CMP_TAms_chrono_outside | CMP | chronological | outside | TAms_chrono | 13.5181394764404 | 1.57161587818483 | 6.55263452981891e-13 | 0.776511987690983 | 5.49005682446815e-08 | 1.02832838350827 | 0.513238030503469 | 1.57161587818483 | 1.62181429545103 | 0.211904661867571 | 1.52189416526171 | 1.11853770752677 | FALSE |
| CMP_TAms_chrono_turquoise | CMP | chronological | turquoise | TAms_chrono | 13.2837328415609 | 1.2893421308899 | 2.42711785115231e-08 | 0.00626983142167074 | 0.0107579084724153 | 0.140506524128796 | 1.1414341081095 | 1.2893421308899 | 1.27713663113322 | 0.962689814855366 | 0.771625714363974 | 1.00381738678419 | TRUE |
| CMP_TAms_chrono_yellow | CMP | chronological | yellow | TAms_chrono | 12.6128269162576 | 1.36251263783774 | 0.00155470474860942 | 4.69738653218747e-08 | 9.00796184650676e-09 | 0.433018986720361 | 1.36251263783774 | 0.910302105777362 | 0.750588500252409 | 1.71868825441021 | 1.69972141518416 | 1.38966605661559 | TRUE |
| CMP_TAms_chrono_orange | CMP | chronological | orange | TAms_chrono | 1.07332762742524 | 1.3667305774026 | 4.17596855270051e-09 | 0.00616665634539307 | 9.86133709312738e-06 | 1.3667305774026 | 1.04292966317747 | 0.308067909026716 | 1.32519943483723 | 0.96895905967252 | 1.24720124645218 | 1.18045324698731 | TRUE |
| CMP_TAms_chrono_brown | CMP | chronological | brown | TAms_chrono | 1.32505139804729 | 0.903182550910418 | 2.43965028951066e-08 | 0.0359875447558034 | 2.67476690758748e-05 | 0.412222463925612 | 0.484666619537402 | 0.903182550910418 | 1.276100875793 | 0.798686096890568 | 1.2238865711876 | 1.09955784795706 | TRUE |
| CMP_TAms_chrono_pink | CMP | chronological | pink | TAms_chrono | 12.8244775263056 | 4.74758863215725 | 0.000392370112328395 | 0.804544452826182 | 7.32685476665974e-08 | 2.3006893237002 | 4.74758863215725 | 1.4896227425564 | 0.793727286789253 | 0.1535481558162 | 1.57553889787738 | 0.840938113494278 | FALSE |
| CMP_TAms_chrono_greenyellow | CMP | chronological | greenyellow | TAms_chrono | 10.2454008309947 | 7.14708590584921 | 0.0198327199086943 | 0.179848384561145 | 0.00130972837545021 | 2.16138052423202 | 1.65191081172942 | 7.14708590584921 | 0.589822114738127 | 0.534736511403284 | 0.941139886105802 | 0.688566170749071 | FALSE |
| CMP_TAms_chrono_green | CMP | chronological | green | TAms_chrono | 7.81765959791766 | 6.97464449328284 | 0.916725744857501 | 0.118712219770985 | 0.00037712312831399 | 5.50545846321484 | 0.0917117986997091 | 6.97464449328284 | 0.0966948507602674 | 0.536479704074035 | 1.05643534830433 | 0.563203301046211 | FALSE |
| CMP_TAms_chrono_grey60 | CMP | chronological | grey60 | TAms_chrono | 4.6147811120749 | 4.43954642269246 | 0.814984835725509 | 0.0277898754427891 | 0.443015119549798 | 2.39068697095814 | 1.29226153101494 | 4.43954642269246 | 0.145726135328425 | 0.82994071988536 | 0.288030178747558 | 0.421232344653781 | FALSE |
| CMP_TAms_chrono_salmon | CMP | chronological | salmon | TAms_chrono | 21.5933027478127 | 5.26262260335956 | 3.7620144772212e-07 | 6.31935841138056e-08 | 0.000554742587355272 | 2.62686983913418 | 1.35218226177394 | 5.26262260335956 | 1.17645079398056 | 1.66349539157906 | 1.0148113766258 | 1.28491918739514 | FALSE |
| CMP_TAms_chrono_black | CMP | chronological | black | TAms_chrono | 16.7422051719165 | 5.44267750057146 | 0.0431573497414402 | 0.171299759129761 | 0.0124247368263572 | 5.44267750057146 | 0.951283137121438 | 3.19857772035514 | 0.537105773592923 | 0.556584483816724 | 0.737914163869101 | 0.610534807092916 | FALSE |
| CMP_TAms_chrono_red | CMP | chronological | red | TAms_chrono | 0.781691026210673 | 1.53422158489987 | 6.00047155091236e-06 | 0.00206822511320336 | 0.00113591207810166 | 1.29941129130114 | 1.53422158489987 | 0.218927091060317 | 1.00703886594647 | 1.03536420236247 | 0.926272958390685 | 0.989558675566542 | TRUE |
| CMP_TAms_chrono_darkred | CMP | chronological | darkred | TAms_chrono | 17.8445651208081 | 6.18282801127864 | 0.00113565873807024 | 0.000842648601106405 | 0.000700443150284876 | 6.18282801127864 | 1.22893835789489 | 2.87002662928019 | 0.78840442758811 | 1.11796990028204 | 1.00523841808726 | 0.970537581985803 | FALSE |
| CMP_TAms_chrono_blue | CMP | chronological | blue | TAms_chrono | 8.33521204030182 | 5.00767463416949 | 4.48788436169598e-08 | 0.115345616679759 | 0.0119889024662003 | 5.00767463416949 | 0.864217067474876 | 3.1941419930658 | 1.23445788367155 | 0.547231284143422 | 0.737800436636047 | 0.83982986815034 | FALSE |
| CMP_TAms_chrono_magenta | CMP | chronological | magenta | TAms_chrono | 11.388576471974 | 5.16558063664093 | 0.299419434845185 | 0.0468700677533329 | 0.000207543779313475 | 2.22848171999807 | 5.16558063664093 | 1.67821769979826 | 0.329097630926434 | 0.677485270618998 | 1.06859540160121 | 0.691726101048881 | FALSE |
| CMP_TAms_chrono_darkgreen | CMP | chronological | darkgreen | TAms_chrono | 21.7508640851462 | 3.95209632308396 | 1.60742403610522e-08 | 0.101546009668967 | 0.10953558959501 | 3.20270325477563 | 0.444351315648134 | 3.95209632308396 | 1.30912375416947 | 0.601255612261404 | 0.511992433855235 | 0.807457266762036 | FALSE |
| CMP_TAms_chrono_tan | CMP | chronological | tan | TAms_chrono | 2.45710572487218 | 6.84253463300028 | 0.814984835725509 | 0.567875972546307 | 0.00879090593978207 | 6.84253463300028 | 1.31231890756758 | 2.64745491905258 | 0.147970721292129 | 0.305181859236157 | 0.77107382468524 | 0.408075468404509 | FALSE |
| CMP_TAms_chrono_cyan | CMP | chronological | cyan | TAms_chrono | 14.4059753373161 | 2.16803294135685 | 0.584889113014699 | 0.0206808851305654 | 0.144131509961148 | 1.89492622868115 | 0.235482950642854 | 2.16803294135685 | 0.221768756241172 | 0.882943095286111 | 0.467780619906251 | 0.524164157144511 | FALSE |
| CMP_TAms_chrono_midnightblue | CMP | chronological | midnightblue | TAms_chrono | 2.73676291569705 | 1.23780814903957 | 0.000145169066436173 | 0.104111391376313 | 1.14570056802328e-07 | 1.23780814903957 | 0.0570892222064812 | 1.1777427525879 | 0.876092914993277 | 0.662060520557707 | 1.53656097787714 | 1.02490480447604 | FALSE |
| CMP_TAms_chrono_darkorange | CMP | chronological | darkorange | TAms_chrono | 20.0225692488611 | 6.30218971778788 | 1 | 0.296407788665699 | 0.000738134129705682 | 2.18397105256734 | 1.80757811080564 | 6.30218971778788 | 0.545629383949786 | 0.402159821911204 | 0.991366616219429 | 0.646385274026806 | FALSE |
| CMP_TAms_chrono_darkgrey | CMP | chronological | darkgrey | TAms_chrono | 3.15770771154075 | 3.85255446644285 | 0.00592555388091835 | 0.234313681529088 | 0.922870284067173 | 3.50921692102436 | 3.85255446644285 | 0.195901952859661 | 0.638588658090798 | 0.469794197607175 | 0.097276825320983 | 0.401886560339652 | FALSE |
| CMP_TAms_chrono_lightgreen | CMP | chronological | lightgreen | TAms_chrono | 1.50484834311761 | 1.5277559171477 | 0.299419434845185 | 0.195969014034543 | 7.68855410211471e-06 | 1.5277559171477 | 0.230358693828964 | 1.2810062164142 | 0.342049456608124 | 0.54270307942301 | 1.31806910215299 | 0.734273879394708 | FALSE |
| CMP_TAms_chrono_skyblue | CMP | chronological | skyblue | TAms_chrono | 10.0997592433309 | 4.51129460298475 | 0.00251220684090586 | 4.22040008597713e-09 | 0.00109803232299121 | 1.98060927884406 | 4.51129460298475 | 1.69534086162426 | 0.735134943763319 | 1.81264769735664 | 0.902259959263444 | 1.1500142001278 | FALSE |
| CMP_TAms_chrono_lightyellow | CMP | chronological | lightyellow | TAms_chrono | 10.710120185328 | 6.14352624967807 | 6.00047155091236e-06 | 0.541249316259654 | 0.144131509961148 | 1.81219500259406 | 2.19008982199016 | 6.14352624967807 | 1.04656848649361 | 0.331230219641877 | 0.488875274218477 | 0.622224660117988 | FALSE |
| CMP_TAms_mort_outside | CMP | mortality | outside | TAms_mort | 11.4411562397544 | 3.31076129123194 | 6.93709366563023e-11 | 0.163908513016268 | 0.000609449868745547 | 1.44312741386663 | 1.54419447808243 | 3.31076129123194 | 1.53598752651362 | 0.568627215019709 | 0.928957722476697 | 1.01119082133668 | FALSE |
| CMP_TAms_mort_turquoise | CMP | mortality | turquoise | TAms_mort | 13.003811149926 | 1.09244642218583 | 6.77796115376162e-09 | 0.000125036076009668 | 3.88973233356303e-07 | 0.232213040438472 | 1.09244642218583 | 0.852606396951373 | 1.34452322444572 | 1.29122268626747 | 1.5027907160861 | 1.3795122089331 | TRUE |
| CMP_TAms_mort_yellow | CMP | mortality | yellow | TAms_mort | 10.665756829595 | 3.76822321241556 | 0.0650647454648595 | 0.0612235297619561 | 0.0122800660972408 | 1.20311397110788 | 3.76822321241556 | 2.10815816100478 | 0.523386114472931 | 0.615251172812485 | 0.738698619001235 | 0.625778635428884 | FALSE |
| CMP_TAms_mort_orange | CMP | mortality | orange | TAms_mort | 20.7408362913636 | 0.838632070258711 | 0.11047008286511 | 0.0464929086432877 | 1.32978354159311e-05 | 0.408846579213591 | 0.838632070258711 | 0.424719034797089 | 0.452539283962663 | 0.646347893220154 | 1.18485838425557 | 0.761248520479462 | FALSE |
| CMP_TAms_mort_brown | CMP | mortality | brown | TAms_mort | 10.5758824538491 | 7.64402695708308 | 0.0650647454648595 | 0.656009736677408 | 0.319369643516106 | 7.64402695708308 | 0.34006975001682 | 3.6715472455479 | 0.491627450341307 | 0.23434395293485 | 0.361801036305625 | 0.362590813193927 | FALSE |
| CMP_TAms_mort_pink | CMP | mortality | pink | TAms_mort | 22.3594798489862 | 1.48810537188183 | 0.00051622010433636 | 1.54788587783345e-06 | 0.115545085644319 | 1.48810537188183 | 0.0511214225819714 | 1.43300585470697 | 0.868091281764462 | 1.47691297195842 | 0.482967269024037 | 0.942657174248973 | FALSE |
| CMP_TAms_mort_greenyellow | CMP | mortality | greenyellow | TAms_mort | 10.4142264899409 | 6.47112188091483 | 0.11047008286511 | 0.0371146697502737 | 0.00466599930055021 | 4.81023513004203 | 0.154825808429827 | 6.47112188091483 | 0.42986411879245 | 0.736326734376772 | 0.860982337396917 | 0.67572439685538 | FALSE |
| CMP_TAms_mort_green | CMP | mortality | green | TAms_mort | 21.668207798873 | 7.65326570066466 | 0.74039982196218 | 0.0352667844875976 | 0.349706968999477 | 0.329340654255638 | 3.68349132352417 | 7.65326570066466 | 0.177083553000918 | 0.828353971762203 | 0.364627603266679 | 0.456688376009933 | FALSE |
| CMP_TAms_mort_grey60 | CMP | mortality | grey60 | TAms_mort | 0.863352350058284 | 2.80725413622778 | 5.2666658880018e-07 | 0.479482261869342 | 3.95013077770432e-05 | 2.34697121025138 | 2.80725413622778 | 0.359750413644586 | 1.12104150190225 | 0.32039965905128 | 1.20850724367731 | 0.883316134876947 | FALSE |
| CMP_TAms_mort_salmon | CMP | mortality | salmon | TAms_mort | 21.5800168742457 | 4.78678693820649 | 0.11047008286511 | 0.000960717063778652 | 0.00962963855563651 | 1.82499467900191 | 1.95678706820481 | 4.78678693820649 | 0.457680344747636 | 1.11692117843491 | 0.776196749749956 | 0.783599424310834 | FALSE |
| CMP_TAms_mort_black | CMP | mortality | black | TAms_mort | 14.3406911014335 | 7.63150060718604 | 0.0113983587335546 | 0.0384766261698779 | 0.00553614597966166 | 7.63150060718604 | 1.66068004036916 | 1.95149568668796 | 0.681801351555303 | 0.777848943506108 | 0.801721715909477 | 0.753790670323629 | FALSE |
| CMP_TAms_mort_red | CMP | mortality | red | TAms_mort | 21.2524747856653 | 4.67807096307654 | 8.26101998758634e-07 | 0.0197245076908238 | 0.0601595774265038 | 4.67807096307654 | 3.19357946598691 | 0.764026385407906 | 1.12337582792688 | 0.818265992807454 | 0.552455539766993 | 0.831365786833776 | FALSE |
| CMP_TAms_mort_darkred | CMP | mortality | darkred | TAms_mort | 14.9666596497573 | 2.0984409439711 | 0.00299942175721486 | 0.0275912785277717 | 9.48016174502148e-08 | 2.0984409439711 | 1.6392324144003 | 0.405853692701202 | 0.705057417353741 | 0.865739173296565 | 1.61757571926251 | 1.06279076997094 | TRUE |
| CMP_TAms_mort_blue | CMP | mortality | blue | TAms_mort | 8.20851096959007 | 5.39057940752243 | 0.000212401170917514 | 0.0374563850342679 | 0.000283771823592337 | 5.39057940752243 | 0.425395894365824 | 4.07884691647831 | 0.887724857381356 | 0.683190257794742 | 1.10308303207307 | 0.891332715749723 | FALSE |
| CMP_TAms_mort_magenta | CMP | mortality | magenta | TAms_mort | 12.6041620057806 | 1.32807194960072 | 0.772057148886701 | 0.0464929086432877 | 0.00362568247603562 | 0.54096769702708 | 0.767467655963481 | 1.32807194960072 | 0.155047273180982 | 0.718732216051786 | 0.884538673544633 | 0.586106054259134 | FALSE |
| CMP_TAms_mort_darkgreen | CMP | mortality | darkgreen | TAms_mort | 22.4093145702652 | 9.09870365099487 | 0.0634393514012682 | 0.305755454622262 | 0.685169177290653 | 0.362702684104367 | 9.09870365099487 | 3.43827845391678 | 0.540435334555918 | 0.454763886484577 | 0.221076866423044 | 0.405425362487846 | FALSE |
| CMP_TAms_mort_tan | CMP | mortality | tan | TAms_mort | 10.3541424816313 | 5.93900340524383 | 0.308638683225286 | 0.546313355398543 | 0.0601595774265038 | 5.93900340524383 | 2.14753988772378 | 1.85513054091973 | 0.299863176920891 | 0.274518778658865 | 0.552344886279352 | 0.375575613953036 | FALSE |
| CMP_TAms_mort_cyan | CMP | mortality | cyan | TAms_mort | 13.4606247901352 | 1.73130616458376 | 0.0634393514012682 | 0.0612235297619561 | 0.00362568247603562 | 1.73130616458376 | 0.71562698100564 | 0.971628684916739 | 0.505480854918459 | 0.697143144309541 | 0.851563610489785 | 0.684729203239262 | FALSE |
| CMP_TAms_mort_midnightblue | CMP | mortality | midnightblue | TAms_mort | 8.10657634227367 | 5.79710169321253 | 0.140795740177912 | 0.726961128372138 | 0.00044257260335414 | 0.683063383788743 | 3.67709906005925 | 5.79710169321253 | 0.407809809879834 | 0.206510056660898 | 1.05100802238516 | 0.555109296308631 | FALSE |
| CMP_TAms_mort_darkorange | CMP | mortality | darkorange | TAms_mort | 3.36678987490333 | 4.95085446087823 | 0.141581051710441 | 0.259750335803272 | 0.182531418757083 | 2.57391348712884 | 4.95085446087823 | 1.31910577221783 | 0.391278116662495 | 0.443753027284072 | 0.466824747182326 | 0.433951963709631 | FALSE |
| CMP_TAms_mort_darkgrey | CMP | mortality | darkgrey | TAms_mort | 2.2032974731237 | 10.9090834411325 | 0.903602200825879 | 0.0275912785277717 | 0.769065322219441 | 1.4885091895127 | 0.375578946788052 | 10.9090834411325 | 0.0836496456532897 | 0.88082371437798 | 0.176212119654325 | 0.380228493228532 | FALSE |
| CMP_TAms_mort_lightgreen | CMP | mortality | lightgreen | TAms_mort | 19.6408365070096 | 8.91367440010976 | 0.177178563724102 | 0.259750335803272 | 0.000292742694544769 | 0.973948682372502 | 8.91367440010976 | 5.1463451827762 | 0.386256793560871 | 0.458908974730351 | 1.09428159771051 | 0.646482455333911 | FALSE |
| CMP_TAms_mort_skyblue | CMP | mortality | skyblue | TAms_mort | 8.25525605601567 | 4.08270806778956 | 0.0634393514012682 | 4.37763089309416e-07 | 0.00553614597966166 | 1.74606167024443 | 4.08270806778956 | 1.72022561341125 | 0.540430304819272 | 1.67458765071058 | 0.762587390484739 | 0.992535115338197 | FALSE |
| CMP_TAms_mort_lightyellow | CMP | mortality | lightyellow | TAms_mort | 11.0040562052861 | 3.90569960170828 | 0.0663219329013521 | 0.137015292115274 | 0.00044257260335414 | 3.90569960170828 | 0.360223780343439 | 3.29322132298694 | 0.488267753662325 | 0.563785344475336 | 1.04893195477882 | 0.70032835097216 | FALSE |

### 7.3 GSE295758 기대 차이와 검정력 표
아침 8시와 저녁 20시의 기대 차이는 `2A · cos(2π(8 − φ)/24)`다. `d_exp = 그 값 / (√2 · σ)`다. 검정력은 쌍 58개, 양측 대응 t, 자유도 57, 비중심 모수 `d_exp_median · √58`이다. RA 세 시계의 관측 아침−저녁은 `results/stage4/rnaseq/g2_clocks.tsv`다. `TAms_*`는 같은 창(07:00–09:00과 19:00–21:00, 경계 포함, 아침−저녁, 참가자 11명)을 이번 점수에 적용했다. 관측 d는 `mean_diff / sd_diff`다.

| clock | d_GSE39445 | d_GSE56931 | d_GSE113883 | d_exp_median | g2_mean_diff | g2_sd_diff | g2_observed_d | power_alpha_0.05 | power_alpha_0.01 | n_pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | 0.00960009963837758 | 0.0799546478017167 | 0.308963171982023 | 0.0799546478017167 | 3.38675095929399 | 7.6092984988986 | 0.44508057605891 | 0.0919710023400413 | 0.0243768273749048 | 58 |
| RA_Peters | 0.286305297017999 | 0.498487573681539 | 0.197988187071947 | 0.286305297017999 | 0.956543945170376 | 7.27976563446149 | 0.131397629154738 | 0.572873321105643 | 0.323267563613161 | 58 |
| RA_default | 0.0749177677231186 | 0.429771444888433 | 0.412448894924069 | 0.412448894924069 | 2.33773972761062 | 7.16235159340991 | 0.326392763203859 | 0.870300205086949 | 0.682114842365414 | 58 |
| TAms_chrono | -0.087085971883512 | 0.208443813199012 | 0.534281736707176 | 0.208443813199012 | 0.019393521149103 | 0.0332147523276872 | 0.583882756606827 | 0.345077590733327 | 0.150545763753431 | 58 |
| TAms_mort | 0.182025471981377 | 0.291250817026961 | 0.693084284817919 | 0.291250817026961 | 0.104008595904116 | 0.11903765298797 | 0.873745351100187 | 0.587339593367421 | 0.336502748601378 | 58 |
| TM_turquoise_mortality | 0.302265260190575 | 0.467000081226947 | 0.641670998412943 | 0.467000081226947 |  |  |  | 0.937738727319687 | 0.809657830878481 | 58 |
| CMP_TAms_chrono_turquoise | 0.20086732450169 | 0.318169181792475 | -0.0815532491415976 | 0.20086732450169 |  |  |  | 0.32447915139951 | 0.137838068985723 | 58 |
| CMP_TAms_chrono_yellow | 0.131204299668764 | 0.801225229773768 | 0.1496592336641 | 0.1496592336641 |  |  |  | 0.201673482334222 | 0.0710530468829523 | 58 |
| CMP_TAms_chrono_orange | 0.1077169039087 | -0.337845873402953 | -0.280140424913859 | -0.280140424913859 |  |  |  | 0.554704876198954 | 0.307078529849453 | 58 |
| CMP_TAms_chrono_brown | -0.253371044612015 | -0.168830899708354 | 0.0516792960726569 | -0.168830899708354 |  |  |  | 0.24390970952334 | 0.0923075897685969 | 58 |
| CMP_TAms_chrono_red | 0.0151182141701988 | -0.483080229057626 | -0.169417336437209 | -0.169417336437209 |  |  |  | 0.245277194009682 | 0.0930256680718282 | 58 |
| CMP_TAms_mort_turquoise | 0.1888870282562 | 0.474739748425362 | 0.0399347730444082 | 0.1888870282562 |  |  |  | 0.293026258353957 | 0.119279433884641 | 58 |
| CMP_TAms_mort_darkred | -0.358485469831977 | 0.107234249480259 | -0.167328334678284 | -0.167328334678284 |  |  |  | 0.24042577405142 | 0.0904866240983017 | 58 |

### 7.4 GSE48113·GSE122725 검정력 표
GSE48113은 in-phase D1 control에서 GSM1168837, GSM1168735를 뺀 147표본, 22명, t = `clock_time_h`다. GSE122725는 period1 control 70표본, 11명, t = `wake_ref_h`다. 표현 행렬은 열지 않았다. 모형은 `score = u_i + A·cos(2π(t − φ)/24) + ε`, A = 2A/2, u_i ~ N(0, σ_u²), ε ~ N(0, σ²)다. 시계마다 개발 자료 3개의 (A, φ, σ_u, σ)를 따로 썼다. 2000회, 시드 20261015. 비율의 분모는 2000이고, 수렴 수는 30행 모두 2000이다.

| design | clock | param_dataset | n_subjects | n_samples | A | phi | sigma_u | sigma | n_rep | n_converged | frac_p_lt_0.05 | frac_p_lt_0.01 | seed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE48113 | RA_GTEx | GSE39445 | 22 | 147 | 2.2811073522782 | 13.7137336093391 | 16.2536695399875 | 25.1603887545092 | 2000 | 2000 | 0.0935 | 0.027 | 20261015 |
| GSE48113 | RA_GTEx | GSE56931 | 22 | 147 | 4.20647808404592 | 13.2671292133393 | 11.0606886507247 | 14.1878999046968 | 2000 | 2000 | 0.579 | 0.3245 | 20261015 |
| GSE48113 | RA_GTEx | GSE113883 | 22 | 147 | 2.17184707776711 | 6.16828794535857 | 10.4317015658613 | 8.81987355742966 | 2000 | 2000 | 0.468 | 0.2485 | 20261015 |
| GSE48113 | RA_Peters | GSE39445 | 22 | 147 | 1.93441030238997 | 6.84503499428439 | 14.613080804435 | 9.12159876450608 | 2000 | 2000 | 0.361 | 0.1545 | 20261015 |
| GSE48113 | RA_Peters | GSE56931 | 22 | 147 | 3.18998678662392 | 4.54192615314062 | 13.9973116863226 | 5.58777765048951 | 2000 | 2000 | 0.998 | 0.968 | 20261015 |
| GSE48113 | RA_Peters | GSE113883 | 22 | 147 | 1.57759284077835 | 11.8362208570625 | 13.6399297750127 | 6.04743769494992 | 2000 | 2000 | 0.477 | 0.244 | 20261015 |
| GSE48113 | RA_default | GSE39445 | 22 | 147 | 0.992494497124205 | 12.0064473932836 | 8.81869510215541 | 9.34019861854301 | 2000 | 2000 | 0.112 | 0.0375 | 20261015 |
| GSE48113 | RA_default | GSE56931 | 22 | 147 | 2.30847035618009 | 9.88701727265479 | 6.90600678789726 | 6.68803391836383 | 2000 | 2000 | 0.773 | 0.5335 | 20261015 |
| GSE48113 | RA_default | GSE113883 | 22 | 147 | 3.3850829723745 | 11.812631972311 | 6.86875908783177 | 6.28931435261543 | 2000 | 2000 | 0.9795 | 0.931 | 20261015 |
| GSE48113 | TAms_chrono | GSE39445 | 22 | 147 | 0.0172605176649426 | 14.2335397784943 | 0.0276429862743754 | 0.0171269287358502 | 2000 | 2000 | 1 | 1 | 20261015 |
| GSE48113 | TAms_chrono | GSE56931 | 22 | 147 | 0.00471533103380021 | 11.4579688514277 | 0.0331124542924689 | 0.0197534498990921 | 2000 | 2000 | 0.431 | 0.201 | 20261015 |
| GSE48113 | TAms_chrono | GSE113883 | 22 | 147 | 0.0209351150669144 | 12.188050742382 | 0.0416315613478055 | 0.0253117943835499 | 2000 | 2000 | 1 | 1 | 20261015 |
| GSE48113 | TAms_mort | GSE39445 | 22 | 147 | 0.0743777777476565 | 13.4288594848278 | 0.161562476919705 | 0.0860831456543238 | 2000 | 2000 | 1 | 1 | 20261015 |
| GSE48113 | TAms_mort | GSE56931 | 22 | 147 | 0.0314734179850113 | 12.0688456403715 | 0.177501208368076 | 0.0740143176007128 | 2000 | 2000 | 0.8865 | 0.7175 | 20261015 |
| GSE48113 | TAms_mort | GSE113883 | 22 | 147 | 0.0512810516194775 | 10.0600990696032 | 0.184575482929547 | 0.0897840656082817 | 2000 | 2000 | 0.9945 | 0.9775 | 20261015 |
| GSE122725 | RA_GTEx | GSE39445 | 11 | 70 | 2.2811073522782 | 13.7137336093391 | 16.2536695399875 | 25.1603887545092 | 2000 | 2000 | 0.0755 | 0.0225 | 20261015 |
| GSE122725 | RA_GTEx | GSE56931 | 11 | 70 | 4.20647808404592 | 13.2671292133393 | 11.0606886507247 | 14.1878999046968 | 2000 | 2000 | 0.3395 | 0.1505 | 20261015 |
| GSE122725 | RA_GTEx | GSE113883 | 11 | 70 | 2.17184707776711 | 6.16828794535857 | 10.4317015658613 | 8.81987355742966 | 2000 | 2000 | 0.1995 | 0.067 | 20261015 |
| GSE122725 | RA_Peters | GSE39445 | 11 | 70 | 1.93441030238997 | 6.84503499428439 | 14.613080804435 | 9.12159876450608 | 2000 | 2000 | 0.1635 | 0.059 | 20261015 |
| GSE122725 | RA_Peters | GSE56931 | 11 | 70 | 3.18998678662392 | 4.54192615314062 | 13.9973116863226 | 5.58777765048951 | 2000 | 2000 | 0.83 | 0.6185 | 20261015 |
| GSE122725 | RA_Peters | GSE113883 | 11 | 70 | 1.57759284077835 | 11.8362208570625 | 13.6399297750127 | 6.04743769494992 | 2000 | 2000 | 0.2995 | 0.1155 | 20261015 |
| GSE122725 | RA_default | GSE39445 | 11 | 70 | 0.992494497124205 | 12.0064473932836 | 8.81869510215541 | 9.34019861854301 | 2000 | 2000 | 0.0975 | 0.03 | 20261015 |
| GSE122725 | RA_default | GSE56931 | 11 | 70 | 2.30847035618009 | 9.88701727265479 | 6.90600678789726 | 6.68803391836383 | 2000 | 2000 | 0.456 | 0.22 | 20261015 |
| GSE122725 | RA_default | GSE113883 | 11 | 70 | 3.3850829723745 | 11.812631972311 | 6.86875908783177 | 6.28931435261543 | 2000 | 2000 | 0.842 | 0.644 | 20261015 |
| GSE122725 | TAms_chrono | GSE39445 | 11 | 70 | 0.0172605176649426 | 14.2335397784943 | 0.0276429862743754 | 0.0171269287358502 | 2000 | 2000 | 1 | 0.9985 | 20261015 |
| GSE122725 | TAms_chrono | GSE56931 | 11 | 70 | 0.00471533103380021 | 11.4579688514277 | 0.0331124542924689 | 0.0197534498990921 | 2000 | 2000 | 0.237 | 0.092 | 20261015 |
| GSE122725 | TAms_chrono | GSE113883 | 11 | 70 | 0.0209351150669144 | 12.188050742382 | 0.0416315613478055 | 0.0253117943835499 | 2000 | 2000 | 0.997 | 0.976 | 20261015 |
| GSE122725 | TAms_mort | GSE39445 | 11 | 70 | 0.0743777777476565 | 13.4288594848278 | 0.161562476919705 | 0.0860831456543238 | 2000 | 2000 | 0.9975 | 0.988 | 20261015 |
| GSE122725 | TAms_mort | GSE56931 | 11 | 70 | 0.0314734179850113 | 12.0688456403715 | 0.177501208368076 | 0.0740143176007128 | 2000 | 2000 | 0.6175 | 0.3815 | 20261015 |
| GSE122725 | TAms_mort | GSE113883 | 11 | 70 | 0.0512810516194775 | 10.0600990696032 | 0.184575482929547 | 0.0897840656082817 | 2000 | 2000 | 0.8635 | 0.6735 | 20261015 |

### 7.5 GSE295758 메타데이터 보충
값 행렬은 열지 않았다. 표본이 4개인 참가자는 subject8, subject16이다. 제목과 characteristics는 원문이다.

```
gsm	subject_id	sample_title	raw_characteristics	age	sex	char_opioids	char_timepoint
GSM8956774	subject8	subject8_batch1_timeD_Homosapiens_RNAseq	tissue: Whole peripheral blood || age: 44 || Sex: male || timepoint: D || batch: 1 || opioids: 0 || timepoint: D	44	male	0	D
GSM8956782	subject16	subject16_batch1_timeD_Homosapiens_RNAseq	tissue: Whole peripheral blood || age: 73 || Sex: male || timepoint: D || batch: 1 || opioids: 0 || timepoint: D	73	male	0	D
GSM8956806	subject8	subject8_batch1_timeN_Homosapiens_RNAseq	tissue: Whole peripheral blood || age: 44 || Sex: male || timepoint: N || batch: 1 || opioids: 0 || timepoint: N	44	male	0	N
GSM8956814	subject16	subject16_batch1_timeN_Homosapiens_RNAseq	tissue: Whole peripheral blood || age: 73 || Sex: male || timepoint: N || batch: 1 || opioids: 0 || timepoint: N	73	male	0	N
GSM8956835	subject8	subject8_batch2_timeD_Homosapiens_RNAseq	tissue: Whole peripheral blood || age: 44 || Sex: male || timepoint: D || batch: 2 || opioids: 0 || timepoint: D	44	male	0	D
GSM8956843	subject16	subject16_batch2_timeD_Homosapiens_RNAseq	tissue: Whole peripheral blood || age: 73 || Sex: male || timepoint: D || batch: 2 || opioids: 0 || timepoint: D	73	male	0	D
GSM8956863	subject8	subject8_batch2_timeN_Homosapiens_RNAseq	tissue: Whole peripheral blood || age: 44 || Sex: male || timepoint: N || batch: 2 || opioids: 0 || timepoint: N	44	male	0	N
GSM8956871	subject16	subject16_batch2_timeN_Homosapiens_RNAseq	tissue: Whole peripheral blood || age: 73 || Sex: male || timepoint: N || batch: 2 || opioids: 0 || timepoint: N	73	male	0	N
```

| subject_id | age | sex | opioid_group | n_samples |
| --- | --- | --- | --- | --- |
| subject1 | 45 | female | 0 | 2 |
| subject2 | 30 | female | 2 | 2 |
| subject3 | 48 | female | 0 | 2 |
| subject4 | 29 | female | 0 | 2 |
| subject5 | 52 | female | 1 | 2 |
| subject6 | 52 | male | 1 | 2 |
| subject7 | 59 | male | 1 | 2 |
| subject8 | 44 | male | 0 | 4 |
| subject9 | 51 | male | 1 | 2 |
| subject10 | 64 | male | 2 | 2 |
| subject11 | 58 | female | 0 | 2 |
| subject12 | 52 | female | 0 | 2 |
| subject13 | 63 | female | 1 | 2 |
| subject14 | 47 | male | 1 | 2 |
| subject15 | 48 | male | 1 | 2 |
| subject16 | 73 | male | 0 | 4 |
| subject17 | 58 | male | 0 | 2 |
| subject18 | 20 | male | 1 | 2 |
| subject19 | 40 | female | 1 | 2 |
| subject20 | 75 | female | 0 | 2 |
| subject21 | 60 | female | 1 | 2 |
| subject22 | 38 | male | 0 | 2 |
| subject23 | 75 | male | 0 | 2 |
| subject24 | 72 | female | 0 | 2 |
| subject25 | 63 | female | 0 | 2 |
| subject26 | 33 | female | 0 | 2 |
| subject27 | 72 | female | 0 | 2 |
| subject28 | 31 | male | 0 | 2 |
| subject29 | 69 | male | 0 | 2 |
| subject30 | 29 | male | 0 | 2 |
| subject31 | 64 | male | 0 | 2 |
| subject32 | 49 | male | 0 | 2 |
| subject33 | 56 | female | 0 | 2 |
| subject34 | 61 | female | 1 | 2 |
| subject35 | 26 | female | 1 | 2 |
| subject36 | 38 | male | 1 | 2 |
| subject37 | 37 | female | 0 | 2 |
| subject38 | 65 | female | 0 | 2 |
| subject39 | 25 | female | 0 | 2 |
| subject40 | 57 | male | 0 | 2 |
| subject41 | 52 | male | 0 | 2 |
| subject42 | 74 | male | 0 | 2 |
| subject43 | 61 | male | 0 | 2 |
| subject44 | 27 | male | 0 | 2 |
| subject45 | 65 | female | 1 | 2 |
| subject46 | 57 | female | 0 | 2 |
| subject47 | 42 | female | 0 | 2 |
| subject48 | 47 | female | 0 | 2 |
| subject49 | 40 | female | 2 | 2 |
| subject50 | 50 | female | 1 | 2 |
| subject51 | 49 | female | 1 | 2 |
| subject52 | 30 | female | 2 | 2 |
| subject53 | 58 | male | 2 | 2 |
| subject54 | 41 | male | 0 | 2 |
| subject55 | 57 | male | 0 | 2 |
| subject56 | 57 | male | 0 | 2 |
| subject57 | 39 | male | 0 | 2 |
| subject58 | 64 | male | 0 | 2 |

## 8. 버그·재실행 기록 (없으면 "없음")
- Part B 첫 실행은 점수 표에서 `gsm` 열이 빠져 `KeyError`로 끝났다. elapsed 0:05.12, maximum resident set size 649448 kB. 열을 넣은 뒤 같은 스크립트를 다시 실행했고, §2의 2:56.93 기록이 그 실행이다.
- Part D 출력에 singular fit 메시지가 있다. 점추정 336개와 부트스트랩 12개 묶음은 모두 수렴했다.
- `TAms_*`의 GSE113883 아침−저녁은 본 Part E 실행이 `g2_clocks.tsv`의 RA 3개만 붙인 뒤에, 같은 창 정의로 따로 계산해 표에 넣었다. 그 계산은 `scripts/33_partE_prereg.R`에 들어 있다. 시뮬레이션은 다시 실행하지 않았다.

## 9. 판단이 필요한 사항
- 사람 자료 스크립트의 모듈은 multispecies mortality 14개다. 보고서 4의 설치류 모듈과 색 이름 7개가 같다. 비교표를 만들기 위해 설치류 모듈 46개를 다시 채점하고 이름을 `TMrd_`로 두었다. 재현 리듬 규칙의 `TM_*`에는 넣지 않았다.
- 분해에 쓸 함수가 코드에 없다. 논문 확인에 쓴 설치류 23개 모듈 유전자 집합으로 항을 나누었다.
- σ_u의 추정량은 귀무 모형 REML이다. 주 모형 잔차 SD도 REML이다.
- 재현 리듬의 3시간 조건은 원형 거리 3 이하로 두었다.

## 10. 지시서와 다르게 한 것
- 설치류 모듈 재채점의 시계 이름은 `TM_{모듈}_{목표}`가 아니라 `TMrd_{모듈}_{목표}`다. 사람 자료 모듈 14개와 색 이름이 겹치기 때문이다.

## 11. 산출물 목록, Part별 커밋 해시
- Part A `0ff7151906c1d6b096e1798f2e2e4fd308d9e880`. `scripts/30_partA_decisions.py`, `results/stage5/partA_decisions.tsv`, `config/decisions.tsv`.
- Part B `adc4488798dda354b854c8985070d7f848114ff9`. `scripts/31_partB_score.py`, `clock_scores.tsv`, `clock_catalog.tsv`, `paper_feature_check.tsv`, `preprocess_spec.tsv`, `preprocess_counts.tsv`, `code_vs_zenodo_md5.tsv`, `weight_vs_predict.tsv`, `within_person_vs_stage4.tsv`.
- Part C `b3459f97a9ab547d6ebdf2f7dc3b4a2f0bcec356`. `component_scores.tsv`, `decomposition_check.tsv`.
- Part D `94cbe4af89821ca15f1b69588876f0172bc32bb7`. `scripts/32_partD_analyze.R`, `h1_point_estimates.tsv`, `h1_bootstrap_ci.tsv`, `h1_vs_stage4.tsv`, `projection_point.tsv`, `projection_ci.tsv`, `results/stage5/boot/`.
- Part E `037236f77ac8d4d7efed92761f375e68caee32cf`. `scripts/33_partE_prereg.R`, `ra_sigma_refit.tsv`, `clock_summary_long.tsv`, `clock_summary_phase.tsv`, `reproducible_modules.tsv`, `gse295758_expected_d.tsv`, `g2_tams.tsv`, `locked_design_power.tsv`, `gse295758_four_sample_subjects.tsv`, `gse295758_subject_covariates.tsv`.
