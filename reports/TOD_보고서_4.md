# TOD 보고서 4 — 판단 사항 처리, 관문 재정정, GSE113883 RNA-seq 개발, 모듈 시계, H5 개발

## 1. 요약 (수치만, 16줄 이내)

- 관문 재판정(refit False, 기준 (b) removed): GSE39445 pass, GSE56931 pass, GSE113883 pass. GSE196728은 결정 4.3으로 이후 분석에서 제외.
- DLMO 값이 있는 참가자 7, 18–26 안 0, `use_t_dlmo` False. `hour`−`hours_since_wake` 최소=중앙값=최대=−1, n=153, 고유 차이 1.
- Zenodo 22166800, `publication_date` 2026-08-29, `version_label` 빈칸, 개념 18763484의 버전 2개 중 나중. `license_field` 빈칸. 배치한 시계 파일 md5 일치 7/7 (로컬 하드링크 5, 다운로드 2). 채점 코드 `/HDD8T1/WORK/GeroSci/clock-specificity`. Part D 수행.
- 모듈 23 (논문 23), 고유 유전자 2007 (논문 2007), 둘 이상 모듈에 들어간 유전자 0. 복합 특징 Chronoage 16517, Mortality 16760 (논문 10487). 교집합 16441, 합집합 16836.
- GSE113883 QC: `log2(TPM+1)`, mean TPM ≥ 1 유전자 12390, 샘플 153, 샘플 중앙값의 중앙값 3.69376571221778, MAD 0.18397853773122, 3 MAD 표시 6.
- GSE113883 lin. RA_GTEx 평균 62.7599299733162, 2A 4.34369415553687 [2.14402772768952, 10.7945031669817], R2_within 0.0315491727978717 [0.00905836854543264, 0.14463204052845], 정점 6.1682879453552 (center 6.7371459498579, half 4.82078466749953), LRT p 0.104384139695722, 무작위 백분위 20.3.
- GSE113883 lin. RA_Peters 평균 39.6065298952733, 2A 3.15518568155639 [0.724064511975616, 6.78626343037065], R2_within 0.03527416607544 [0.00202148111734601, 0.138905271852363], 정점 11.8362208570621 (center 11.6741391650562, half 4.62490864529145), LRT p 0.077344312576993, 무작위 백분위 22.9.
- GSE113883 lin. RA_default 평균 37.8218395219532, 2A 6.77016594474564 [4.54731764000349, 9.17576797367988], R2_within 0.134425171769798 [0.0636149012292298, 0.247718423922361], 정점 11.8126319723118 (center 11.8212982069323, half 1.44667958367679), LRT p 3.24675832907807e-05, 무작위 백분위 72.9.
- official–lin 개인 내 Pearson: RA_GTEx 0.997677855578136, RA_Peters 0.998982885935051, RA_default 0.999014038243109. 세 시계 `use_lin_only` TRUE. official 부트스트랩 없음.
- `p_med` [구간]: RA_GTEx -0.149683650058583 [-0.549633569275945, 1.49242245851642], RA_Peters 0.812219272485387 [-0.47918439806038, 3.41052546250761], RA_default 0.870626030579166 [0.464188530693954, 1.54833140363586]. 두 REML 수렴 2000/2000.
- 아침−저녁 평균 차이 [95% 구간], 참가자 11: RA_GTEx 3.38675095929399 [-1.72524533014459, 8.49874724873257], RA_Peters 0.956543945170376 [-3.93406915362397, 5.84715704396472], RA_default 2.33773972761062 [-2.47399354510375, 7.14947300032499].
- ΔAIC (T_wake AIC − T_clock AIC): RA_GTEx 0, RA_Peters 2.98141458188184e-08, RA_default 0, PER1 5.6843418860808e-14, PER3 -2.8421709430404e-14, BMAL1 -2.8421709430404e-14. T_dlmo 미수행.
- 모듈 BH q < 0.05: GSE113883 chronological 14/23, GSE113883 mortality 10/23, GSE39445 chronological 14/23, GSE39445 mortality 10/23, GSE56931 chronological 13/23, GSE56931 mortality 10/23. 복합 시계 부트스트랩 수렴 2000/2000, 비수렴 0.
- H5 리듬 유전자 q<0.05: GSE39445 10041/19541, GSE56931 4087/20755, GSE113883 8755/12390. GSE39445 RA_GTEx 남음 107/229, R2 0.00480935568643748→0.00839588294303228, 참가자평균 r 0.952675076502066; GSE39445 RA_Peters 남음 58/280, R2 0.0237318052396572→0.0472106995689958, 참가자평균 r 0.320866373280244; GSE39445 RA_default 남음 92/214, R2 0.00680181474094321→0.00817794646007386, 참가자평균 r 0.700138202872276; GSE56931 RA_GTEx 남음 92/255, R2 0.0366339113390549→0.0292048630098719, 참가자평균 r 0.702451762976866; GSE56931 RA_Peters 남음 37/278, R2 0.145491804603583→0.0640713388376701, 참가자평균 r 0.698616403013853; GSE56931 RA_default 남음 78/235, R2 0.0658778147622046→0.0211941126866814, 참가자평균 r 0.464065753956846; GSE113883 RA_GTEx 남음 127/271, R2 0.0315491727978717→0.003459203222695, 참가자평균 r 0.358379331746022; GSE113883 RA_Peters 남음 36/276, R2 0.03527416607544→0.0383896778926445, 참가자평균 r 0.446910683302218; GSE113883 RA_default 남음 97/242, R2 0.134425171769798→0.026148980350796, 참가자평균 r 0.57106656133164.
- 마지막 분석 커밋 `963e75aee3e3b97870c0a537ae3e9d048f76b45c`.

## 2. 환경과 자원 (바뀐 것, 새 가상환경, /usr/bin/time 기록)

시계를 읽는 가상환경을 새로 만들었다. 경로 `/HDD8T1/WORK/GeroSci/tod4_clock_env`. 기준 Python은 `/home/jkang/anaconda3/bin/python3` 3.12.7. 고정 파일은 `/HDD8T1/WORK/GeroSci/sc-module-aging/requirements.txt`. 설치한 고정 버전은 joblib 1.6.0, numpy 1.26.4, pandas 2.2.3, scikit-learn 1.3.2, scipy 1.13.1. 이어서 scmodage 0.2.0, clockspec 0.1.0을 그 환경에 설치했다. `clock-specificity`에는 requirements 또는 lock 파일이 없다. requirements.txt의 anndata, h5py, matplotlib, openpyxl, pyarrow, pytest, PyYAML, scanpy, seaborn, statsmodels, tqdm은 시계를 읽는 함수가 불러오지 않아 설치하지 않았다. 직렬화 파일을 읽을 때 버전 경고는 없었다 (`version_warnings.tsv`의 message는 none).

R은 기존 `r_env`다. RNAAgeCalc 1.22.0, MCPcounter 1.2.0, lme4 2.0-6, limma 3.66.0. GPU를 쓰지 않았다. `ulimit`을 두지 않았다. OMP/MKL 스레드 상한을 두지 않았다. `config/paths.yaml`의 `module_clock_dir`는 `/HDD8T1/WORK/GeroSci/zenodo_22166800_v2`, `module_clock_code_dir`는 `/HDD8T1/WORK/GeroSci/clock-specificity`.

| 실행 | elapsed | max RSS (kB) | exit |
| --- | --- | --- | --- |
| 26_partC_rnaseq.R, data.table 필터 수정 전 1회 | 0:25.85 | time -v 파일이 다음 실행으로 덮임 | 1 |
| 26_partC_rnaseq.R, data.table 필터 수정 전 2회 | 0:25.97 | time -v 파일이 다음 실행으로 덮임 | 1 |
| 26_partC_rnaseq.R | 11:39.96 | 2195848 | 0 |
| 27_partD_score.py | 0:06.29 | time -v로 재지 않음 | 0 |
| 28_partD_analyze.R | 18:10.03 | 1533696 | 0 |
| 29_partE_h5.R | 0:26.99 | 1710216 | 0 |

기록된 최대 RSS는 40 GB 아래였다.

## 3. 결정 기록 (config/decisions.tsv에 이번에 더한 행)

| date | id | content | reason | timing |
| --- | --- | --- | --- | --- |
| 2026-10-16 | 4.1.1 | Three-platform coverage is approved. The four-platform sentence in instruction 3 section 1.3 was an error. The three coverages in decision 3.3 stand. | Instruction 4 section 1.1 item 1. Report 3 section 11. | prior |
| 2026-10-16 | 4.1.2 | GSE113883 condition=control and period=time_course for every sample is approved. | Instruction 4 section 1.1 item 2. Report 3 section 11. | prior |
| 2026-10-16 | 4.1.3 | The GSE196728 sea-level control label is not applicable. The dataset is excluded by decision 4.3. | Instruction 4 section 1.1 item 3. | prior |
| 2026-10-16 | 4.1.4 | The species-call thresholds are not applicable. The record is kept. GSE196728 is excluded. | Instruction 4 section 1.1 item 4. | prior |
| 2026-10-16 | 4.1.5 | The GSE196728 count-column case match is not applicable. The dataset is excluded. | Instruction 4 section 1.1 item 5. | prior |
| 2026-10-16 | 4.1.6 | Dropping 583 ENSEMBL keys that mapped to more than one SYMBOL is approved. | Instruction 4 section 1.1 item 6. Report 3 section 11. | prior |
| 2026-10-16 | 4.1.7 | Using circadianphase as an hour-unit axis is approved. The stored values are 3-hour bin labels. State that on every result that uses them. | Instruction 4 section 1.1 item 7. Report 3 section 11. | prior |
| 2026-10-16 | 4.1.8 | The unresolved zero point of salivary_dlmo_circadian_time is handled by instruction 4 section 1.4. Stage-3 T_dlmo rows are not used. | Instruction 4 section 1.1 item 8. | prior |
| 2026-10-16 | 4.1.9 | Keeping expression matrices out of git is approved. | Instruction 4 section 1.1 item 9. Report 3 section 11. | prior |
| 2026-10-16 | 4.1.10 | The module-clock path was not found in stage 3. Part B repeats the search with Zenodo checksums. | Instruction 4 section 1.1 item 10. | prior |
| 2026-10-16 | 4.1.11 | GSE113883 sample-median MAD of 0 is not used as the QC flag. Recompute sample medians on genes with mean TPM at least 1. | Instruction 4 section 1.1 item 11 and section 3.1. | prior |
| 2026-10-16 | 4.1.12 | Leaving sections 3.2-3.5 unrun after the stage-3 positive-control failure is approved. After decision 4.2, Part C runs those sections for GSE113883. | Instruction 4 section 1.1 item 12. Report 3 section 12 is also approved. | prior |
| 2026-10-16 | 4.2 | Development positive-control criteria are (a') PER1 and PER3 peaks inside 04-12 inclusive, and PER1 and PER3 LRT p both < 0.01. Criterion (b) is removed. Locked-data gates are set at preregistration and do not use criterion (b). Rejudge uses stored peaks and does not refit. GSE39445 pass, GSE56931 pass, GSE113883 pass. | Instruction 4 section 1.2. Post-hoc correction. GSE113883 clock results had not been calculated. | 사후 |
| 2026-10-16 | 4.3 | GSE196728 is excluded from every later analysis. Prior outputs are not deleted. | Instruction 4 section 1.3. Sea-level positive control was weak (PER1 LRT p 0.187843896897575, PER3 LRT p 0.0573355762517646, 8 subjects). Processing text maps reads to mm10. ARNTL FPKM is missing for every sample because gene length was unavailable. | 사후 |

관문 재판정 수치 (`results/stage4/gate_rejudge.tsv`, 저장 적합, refit False):

| dataset | per1_peak_h | per1_lrt_p | per3_peak_h | per3_lrt_p | criterion_a_prime | both_lrt_p_lt_0.01 | criterion_b | refit | pass |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | 6.90324597433964 | 9.73597605831301e-18 | 6.08583286807854 | 1.177374716677e-06 | True | True | removed | False | True |
| GSE56931 | 5.79155389251083 | 1.58787990764694e-24 | 4.8945603273957 | 3.48666322673004e-14 | True | True | removed | False | True |
| GSE113883 | 6.88888783794359 | 2.75027621424795e-09 | 4.3428731541087 | 7.582950678639909e-11 | True | True | removed | False | True |

## 4. GSE113883 시각 변수

### 4.1 DLMO 원값 표와 판정

규칙: 숫자 시간이 [18, 26] 안. 24 이상은 다음 날 새벽 표기. 한 명이라도 밖이면 T_dlmo를 하지 않는다.

| subject_id | salivary_dlmo_circadian_time | n_samples | n_unique | hour_decimal | in_18_26 | use_as_dlmo_clock_h |
| --- | --- | --- | --- | --- | --- | --- |
| 1741_001 | 12.5 | 14 | 1 | 12.5 | False | False |
| 1741_005 | 13.0 | 13 | 1 | 13.0 | False | False |
| 1741_006 | 11.0 | 13 | 1 | 11.0 | False | False |
| 1741_009 | 15.0 | 14 | 1 | 15.0 | False | False |
| 1741_010 | 12.0 | 14 | 1 | 12.0 | False | False |
| 1741_011 | 12.0 | 14 | 1 | 12.0 | False | False |
| 1741_014 | 12.5 | 13 | 1 | 12.5 | False | False |

| n_subjects_with_value | all_in_18_26 | use_t_dlmo | rule | t_dlmo_definition |
| --- | --- | --- | --- | --- |
| 7 | False | False | numeric hour in [18, 26]; values at or above 24 are next-day morning labels |  |

7명 모두 18–26 밖이다. T_dlmo는 수행하지 않았다. 단계 3의 T_dlmo 행은 쓰지 않았다.

### 4.2 hours_since_wake와 hour 필드

`hours_since_wake` = 2 + (참가자 첫 채혈부터 지난 시간). 지난 시간은 `local_time`을 `hour` 순으로 두고, 시계가 뒤로 가면 24를 더해 누적했다. 프로토콜 문장: "starting 2 h after wake on on Day 2".

| n_samples | n_subjects | diff_min | diff_median | diff_max | n_unique_diff | definition |
| --- | --- | --- | --- | --- | --- | --- |
| 153 | 11 | -1.0 | -1.0 | -1.0 | 1 | hours_since_wake = 2 + cumulative local_time steps within subject, ordered by hour, with a negative step wrapped by +24 |

`clock_time_h` − `hours_since_wake`는 6이 87샘플, −18이 66샘플이다. 그 차이를 24로 나눈 나머지는 153샘플 모두 6이다.

## 5. 모듈 시계 확보

### 5.1 Zenodo 기록 (버전, 공개일, 라이선스, 파일 구분 표)

API `conceptrecid:18763484`. 메타데이터 `version` 필드는 두 기록 모두 빈칸이다. 공개일 순서로 버전 1은 18763485 (2026-02-25), 버전 2는 22166800 (2026-08-29)이다. 2026-08-29보다 새로운 공개일은 없다.

| order_by_publication_date | record_id | doi | publication_date | created | version_label | is_2026_08_29 |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 18763485 | 10.5281/zenodo.18763485 | 2026-02-25 | 2026-03-01T20:37:34.520921+00:00 |  | False |
| 2 | 22166800 | 10.5281/zenodo.22166800 | 2026-08-29 | 2026-08-29T22:41:13.647515+00:00 |  | True |

| record_id | conceptrecid | publication_date | version_label | access_right | license_field | description_repository | local_license_title | n_versions | newer_than_2026_08_29 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 22166800 | 18763484 | 2026-08-29 |  | open |  | https://github.com/Gladyshev-Lab/tAge | MGB OPEN ACCESS LICENSE 1.0 | 2 | False |

파일 구분. `module clocks`로 끝나고 `.zip`이면 module_clock. `.pkl`이면 composite_clock (tAge README의 복합 시계 60개와 개수가 같다). `Expression_data`와 `Data_annotation`은 other. `used_for_scoring`는 이번 채점에 쓴 EN multitissue scaleddiff 5개와 모듈 압축 2개다.

| name | size | md5 | class | used_for_scoring |
| --- | --- | --- | --- | --- |
| BR_NormalizedAge_Mouse_Multitissue_scaleddiff.pkl | 2393070870 | 8d1d0ac05cd33945666b7b5149341703 | composite_clock | False |
| BR_NormalizedAge_Multispecies_Multitissue_yugenediff.pkl | 880570274 | a09fbc41ccc94202ed280a194e3447cf | composite_clock | False |
| BR_NormalizedAge_Mouse_Multitissue_yugenediff.pkl | 2419988892 | 61c07c0d58f0e05de4bd639e8ced88bb | composite_clock | False |
| BR_NormalizedAge_Multispecies_Multitissue_scaleddiff.pkl | 880570274 | e7b8d25fb30d82aec0a9e0c071ad9a83 | composite_clock | False |
| BR_Chronoage_Mouse_Multitissue_scaleddiff.pkl | 2351467264 | aa9f753dc279eac9688f08ef0de1e569 | composite_clock | False |
| BR_Chronoage_Mouse_Multitissue_yugenediff.pkl | 2374290449 | c77c24d194a983cece0a7e1ce3cbcebd | composite_clock | False |
| BR_NormalizedAge_Rodents_Multitissue_yugenediff.pkl | 2265314557 | 67bdc7efb474f2b69e7d9bcc17338530 | composite_clock | False |
| BR_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl | 2248384016 | 6f9063381b97f56d1932aeecb7f91ed1 | composite_clock | False |
| EN_NormalizedAge_Rodents_Liver_scaleddiff.pkl | 913774 | 791f866582741402cabb6eed4c0bde30 | composite_clock | False |
| BR_Chronoage_Rodents_Brain_WT_scaleddiff.pkl | 2182090561 | 3f6921476fb19693ac84146635373817 | composite_clock | False |
| EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl | 584873 | 87307d7f5572d9932382912d3fd9571c | composite_clock | True |
| EN_NormalizedAge_Multispecies_Multitissue_yugenediff.pkl | 584873 | 634b7ef57c1ce535cce1172023f444b3 | composite_clock | False |
| EN_NormalizedAge_Mouse_Multitissue_scaleddiff.pkl | 964509 | fc3ddb5a235db2d280f1fd37284f9f73 | composite_clock | False |
| EN_Mortality_Rodents_Multitissue_yugenediff.pkl | 938484 | 27e6fcd83926eb3fc3680fa0ab23c5ff | composite_clock | False |
| Expression_data_relative_rodents_Scaled.csv | 721701059 | c283869cf69e35afa200c6a6fc281519 | other | False |
| BR_Chronoage_Multispecies_Multitissue_yugenediff.pkl | 880570274 | c33d7b7ec750249f55967079b27bc37d | composite_clock | False |
| BR_Chronoage_Rodents_Skeletal_muscle_WT_scaleddiff.pkl | 1910481982 | aaf40e3c69042ff8f0c0cad8c432d25c | composite_clock | False |
| BR_Chronoage_Rodents_Kidney_WT_yugenediff.pkl | 1968526857 | 6f87c9be1ba1cb7d3efb45f9f8efbf34 | composite_clock | False |
| BR_Chronoage_Multispecies_Multitissue_scaleddiff.pkl | 880570274 | efa60dab8b6160bfd39381be73acee9f | composite_clock | False |
| BR_Mortality_Mouse_Multitissue_yugenediff.pkl | 2419988892 | daea37146186fda8b5a7f258ac3b8b9b | composite_clock | False |
| BR_Mortality_Multispecies_Multitissue_yugenediff.pkl | 880570274 | 8a123360449cd022d695c783a4b13d6c | composite_clock | False |
| BR_NormalizedAge_Rodents_Liver_scaleddiff.pkl | 2148922007 | fb1da22da9a43058d49c655c53a3e02d | composite_clock | False |
| BR_Mortality_Rodents_Multitissue_yugenediff.pkl | 2265314557 | 4b316c2ba5397a3ba3db57c159cf68aa | composite_clock | False |
| BR_Chronoage_Rodents_Multitissue_scaleddiff.pkl | 2183676137 | eddcba9e9a2cbff21f01ce23d27cd73e | composite_clock | False |
| BR_NormalizedAge_Rodents_Liver_yugenediff.pkl | 2170478953 | 53f38dbc98ffa9df5ce532a5e42ea0a2 | composite_clock | False |
| BR_Mortality_Rodents_Liver_scaleddiff.pkl | 2148921991 | 6e378ec38438dc4919125c4dcb5b18c3 | composite_clock | False |
| BR_Mortality_Rodents_Multitissue_scaleddiff.pkl | 2248384016 | da7dc56fada1ddb89bb4c6a0b947faca | composite_clock | False |
| EN_Chronoage_Rodents_Skeletal_muscle_yugenediff.pkl | 864063 | 30f9e2cc01a982bbf2e7deaee4aab62f | composite_clock | False |
| EN_Chronoage_Rodents_Brain_yugenediff.pkl | 925925 | 4de57fb362eb885013bb0edd0150a3d7 | composite_clock | False |
| EN_Chronoage_Rodents_Brain_scaleddiff.pkl | 921144 | 035013a679acd5ddb4ecb5f092626b6d | composite_clock | False |
| EN_NormalizedAge_Multispecies_Multitissue_scaleddiff.pkl | 584873 | a009e1cc7049c7cfc7efd41e1c348476 | composite_clock | False |
| EN_Mortality_Mouse_Multitissue_yugenediff.pkl | 970099 | a96f0068691f853a89e9c28935d49561 | composite_clock | False |
| BR_Mortality_Multispecies_Multitissue_scaleddiff.pkl | 880570274 | 20fedb08b78828b8a355c5a53c9da369 | composite_clock | False |
| EN_Mortality_Rodents_Multitissue_scaleddiff.pkl | 934823 | 26cc9e089e45bf978205738e56f190fb | composite_clock | True |
| EN_Mortality_Rodents_Liver_yugenediff.pkl | 918512 | d4d0f33b64f23b734c3e95d06e78eba3 | composite_clock | False |
| EN_Mortality_Rodents_Liver_scaleddiff.pkl | 913774 | cd76d203c42f88a7580707c36a7e4a0b | composite_clock | False |
| EN_Chronoage_Rodents_Kidney_scaleddiff.pkl | 871690 | c58a1430f36c738cc3d72825f28d44ec | composite_clock | False |
| EN_Mortality_Multispecies_Multitissue_yugenediff.pkl | 584873 | bd371c5a7c6706c363dcac9662ceab00 | composite_clock | False |
| EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl | 921296 | 19fee9e782b4eda36ab13a459fcdb3db | composite_clock | True |
| EN_Chronoage_Rodents_Liver_yugenediff.pkl | 912387 | 1967d3ea992f697ae7c1761ac9bf1aea | composite_clock | False |
| EN_Chronoage_Multispecies_Multitissue_yugenediff.pkl | 584873 | ad893ce7c008094b9f36450257d1934a | composite_clock | False |
| Expression_data_relative_rodents.rds | 891751277 | 16af0a906aa243a7066894242fb6f3c9 | other | False |
| Expression_data_relative_rodents_YuGene.csv | 1181727132 | 26a44ce691c5781bcbd4f3c2423b198b | other | False |
| Expression_data_absolute_rodents.rds | 615602456 | 9d4d8627f548a73e8cb851ea2d2f1d13 | other | False |
| Expression_data_absolute_rodents_YuGene.csv | 1112885301 | c822d9bcf3cb05f516ce8408fe5a08f3 | other | False |
| Expression_data_absolute_rodents_Scaled.csv | 1119842142 | 3b182622809b1361a6266f119102d017 | other | False |
| BR_Chronoage_Rodents_Brain_WT_yugenediff.pkl | 2204078350 | 68e03da91cf58a53c202ad1b223baaab | composite_clock | False |
| BR_Chronoage_Rodents_Kidney_WT_scaleddiff.pkl | 1955247155 | 065b5f041f73e26d7f9719dd24e569dc | composite_clock | False |
| BR_Chronoage_Rodents_Multitissue_yugenediff.pkl | 2199035132 | 509e45c1d11e7278066ae1e7ffa4b3fc | composite_clock | False |
| BR_Chronoage_Rodents_Liver_yugenediff.pkl | 2143156620 | a9f74d3b86c0b5e29b0d56caab184fec | composite_clock | False |
| EN_Chronoage_Rodents_Multitissue_yugenediff.pkl | 924659 | 2e19c9ef45b45ffdebd0b2882722041a | composite_clock | False |
| BR_Chronoage_Rodents_Skeletal muscle_WT_yugenediff.pkl | 1920385128 | c6ebdc57c4ef19fbce4a21da6678147d | composite_clock | False |
| EN_Chronoage_Rodents_Skeletal_muscle_scaleddiff.pkl | 861717 | 5640befd2990e1601e86b0ccfcbaf771 | composite_clock | False |
| EN_Chronoage_Rodents_Kidney_yugenediff.pkl | 874720 | 5a465266c64f29187fff656231d5a8d6 | composite_clock | False |
| EN_Mortality_Mouse_Multitissue_scaleddiff.pkl | 964509 | 43ee221e8cbe8577107f1a74230891e0 | composite_clock | False |
| EN_NormalizedAge_Rodents_Multitissue_yugenediff.pkl | 938484 | b37545ab983594d4eae123ce7947c880 | composite_clock | False |
| EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl | 934823 | bbf2b0f8feecaafb29456e3403616144 | composite_clock | True |
| EN_Chronoage_Mouse_Multitissue_scaleddiff.pkl | 955991 | c31d3b5ea8af7aa0891aee12d95109b5 | composite_clock | False |
| EN_NormalizedAge_Mouse_Multitissue_yugenediff.pkl | 970099 | c3ef18b0b8a8ae5609f9b41cd87234af | composite_clock | False |
| EN_NormalizedAge_Rodents_Liver_yugenediff.pkl | 918512 | 8dabb0f4f2354d9e0eaff2b01e3113db | composite_clock | False |
| EN_Chronoage_Rodents_Liver_scaleddiff.pkl | 909931 | a73012b22b5f4bc9881edd1ac2b0cb72 | composite_clock | False |
| EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl | 584873 | 7dc8a470f5db20a75d7052829efbd37a | composite_clock | True |
| EN_Chronoage_Mouse_Multitissue_yugenediff.pkl | 960792 | 7dff7ed53329448476836fa052ed6097 | composite_clock | False |
| BR_Mortality_Rodents_Liver_yugenediff.pkl | 2170478969 | de581da49e81703429d07109293fdb39 | composite_clock | False |
| BR_Mortality_Mouse_Multitissue_scaleddiff.pkl | 2393070886 | 8ff83db036636b11a8da86b5adaecabc | composite_clock | False |
| BR_Chronoage_Rodents_Liver_scaleddiff.pkl | 2132171652 | 8c9bf298ab0783145f8721bcb4fc32b2 | composite_clock | False |
| Rodent module clocks.zip | 394795 | be9234547aacfbb2e83f0f1c5e44b3bd | module_clock | True |
| Data_annotation_relative_rodents.xlsx | 495148 | d9024c6f91e89e170b9c0a869aded48b | other | False |
| Multispecies module clocks.zip | 252009 | f2a01d5d0fc07745f78cf704e2dd3e98 | module_clock | True |
| Data_annotation_absolute_rodents.xlsx | 392619 | c4f51562c38ab278dbb3b4b8fb7f64cc | other | False |

| item | n_composite_pkl_on_record | n_composite_placed | n_composite_not_downloaded | bytes_composite_not_downloaded | n_training_files | bytes_training_not_downloaded | training_files_downloaded |
| --- | --- | --- | --- | --- | --- | --- | --- |
| clock_file_set | 60 | 5 | 55 | 58122080804 | 8 | 5644397134 | 0 |

### 5.2 로컬 탐색 (찾아본 위치, 일치한 파일, md5 표)

| location | what | result |
| --- | --- | --- |
| /HDD8T1/WORK/GeroSci/clock-specificity | repository tree | no Zenodo clock pkl or module zip |
| /HDD8T1/WORK/GeroSci/sc-module-aging | data/external | EN multitissue scaleddiff pkl files and renamed module zips; see local_md5.tsv |
| /home/jkang/.cache | name search for 22166800, scaleddiff, module clock | no match |
| tAge download_clocks | default destination in the package | dest_dir is an argument; vignette uses tempdir(); no fixed cache path |

이름·크기·md5가 모두 같은 파일만 `accepted`다. 모듈 압축은 로컬 파일명 `rodent_modules.zip`, `multispecies_modules.zip`이라 이름 일치가 False다. 그 두 파일의 크기와 md5는 Zenodo와 같다.

| zenodo_name | local_path | exists | size | md5 | name_match | size_match | md5_match | accepted |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl | /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_small/EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl | True | 921296 | 19fee9e782b4eda36ab13a459fcdb3db | True | True | True | True |
| EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl | /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_zenodo/EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl | True | 921296 | 19fee9e782b4eda36ab13a459fcdb3db | True | True | True | True |
| EN_Mortality_Rodents_Multitissue_scaleddiff.pkl | /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_small/EN_Mortality_Rodents_Multitissue_scaleddiff.pkl | True | 934823 | 26cc9e089e45bf978205738e56f190fb | True | True | True | True |
| EN_Mortality_Rodents_Multitissue_scaleddiff.pkl | /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_zenodo/EN_Mortality_Rodents_Multitissue_scaleddiff.pkl | True | 934823 | 26cc9e089e45bf978205738e56f190fb | True | True | True | True |
| EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl | /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_small/EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl | True | 934823 | bbf2b0f8feecaafb29456e3403616144 | True | True | True | True |
| EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl | /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_zenodo/EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl | True | 934823 | bbf2b0f8feecaafb29456e3403616144 | True | True | True | True |
| EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl | /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_small/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl | True | 584873 | 7dc8a470f5db20a75d7052829efbd37a | True | True | True | True |
| EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl | /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_small/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl | True | 584873 | 87307d7f5572d9932382912d3fd9571c | True | True | True | True |
| EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl | /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_zenodo/EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl | True | 584873 | 87307d7f5572d9932382912d3fd9571c | True | True | True | True |
| Rodent module clocks.zip | /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_modules/rodent_modules.zip | True | 394795 | be9234547aacfbb2e83f0f1c5e44b3bd | False | True | True | False |
| Multispecies module clocks.zip | /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tage_modules/multispecies_modules.zip | True | 252009 | f2a01d5d0fc07745f78cf704e2dd3e98 | False | True | True | False |

### 5.3 내려받기 (했을 때: 위치, 파일, md5 표)

위치 `/HDD8T1/WORK/GeroSci/zenodo_22166800_v2/` (TOD 바깥). pkl 5개는 이름·크기·md5가 맞은 로컬 파일의 하드링크다. zip 2개는 Zenodo에서 받아 md5를 다시 계산했다. 학습 파일은 받지 않았다. 압축은 md5 확인 뒤 같은 폴더에 풀었다. 풀린 파일의 이름·바이트·SHA-256은 `module_archive_members.tsv`에 있다. 그 표는 유전자나 계수를 담지 않는다.

| name | size | md5 | zenodo_md5 | md5_match | source |
| --- | --- | --- | --- | --- | --- |
| EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl | 921296 | 19fee9e782b4eda36ab13a459fcdb3db | 19fee9e782b4eda36ab13a459fcdb3db | True | local_hardlink |
| EN_Mortality_Rodents_Multitissue_scaleddiff.pkl | 934823 | 26cc9e089e45bf978205738e56f190fb | 26cc9e089e45bf978205738e56f190fb | True | local_hardlink |
| EN_NormalizedAge_Rodents_Multitissue_scaleddiff.pkl | 934823 | bbf2b0f8feecaafb29456e3403616144 | bbf2b0f8feecaafb29456e3403616144 | True | local_hardlink |
| EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl | 584873 | 7dc8a470f5db20a75d7052829efbd37a | 7dc8a470f5db20a75d7052829efbd37a | True | local_hardlink |
| EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl | 584873 | 87307d7f5572d9932382912d3fd9571c | 87307d7f5572d9932382912d3fd9571c | True | local_hardlink |
| Rodent module clocks.zip | 394795 | be9234547aacfbb2e83f0f1c5e44b3bd | be9234547aacfbb2e83f0f1c5e44b3bd | True | downloaded |
| Multispecies module clocks.zip | 252009 | f2a01d5d0fc07745f78cf704e2dd3e98 | f2a01d5d0fc07745f78cf704e2dd3e98 | True | downloaded |

### 5.4 채점 코드 저장소 (판정 근거)

`module_clock_code_dir`는 비어 있어서 §2.4를 수행했다. 두 저장소 모두 사람 식별자를 설치류 오솔로그로 바꾸고 예측하는 코드가 있다. 사람 자료를 채점하는 파일의 커밋터 날짜가 늦은 쪽은 `clock-specificity`의 `scripts/10_score_human.py` (2026-10-14T19:57:59+09:00)다. `sc-module-aging`의 `scripts/16_score_bulk.py`는 2026-10-12T14:05:01+09:00이다.

| repo | file | function | commit | committer_date |
| --- | --- | --- | --- | --- |
| /HDD8T1/WORK/GeroSci/clock-specificity | scripts/10_score_human.py | main | b85cbd75b332ee2fc3f833c28c900ae6278edee3 | 2026-10-14T19:57:59+09:00 |
| /HDD8T1/WORK/GeroSci/clock-specificity | src/clockspec/preprocess.py | array_to_diff;counts_to_diff | bc1265385d3fabf6639e8a70f7da95cd44a21f0f | 2026-10-14T14:39:06+09:00 |
| /HDD8T1/WORK/GeroSci/clock-specificity | src/clockspec/clocks.py | composite;modules | de9f09961e7e44962a6b671816854b8871a52ae8 | 2026-10-14T14:02:38+09:00 |
| /HDD8T1/WORK/GeroSci/sc-module-aging | scripts/16_score_bulk.py | module top level | 1af38ed8d782e4af95e8a68622897902f79feb58 | 2026-10-12T14:05:01+09:00 |
| /HDD8T1/WORK/GeroSci/sc-module-aging | src/scmodage/preprocess.py | map_to_mouse_entrez;preprocess_stratum | b35a3f035987e7cdb1a9e357c418a08dea071d3e | 2026-10-12T12:13:57+09:00 |
| /HDD8T1/WORK/GeroSci/sc-module-aging | src/scmodage/clocks.py | predict;predict_all | 985c18d3f3ac0f54f342e97fd86825eaa2f4d24e | 2026-10-12T14:48:49+09:00 |

| chosen_repo | reason | ortholog_source | ortholog_table | ortholog_rule | code_relative_expression | code_filters | part_d_input |
| --- | --- | --- | --- | --- | --- | --- | --- |
| /HDD8T1/WORK/GeroSci/clock-specificity | both repositories contain human-to-rodent ortholog mapping and predict; scripts/10_score_human.py committer date 2026-10-14 is later than scripts/16_score_bulk.py 2026-10-12 | scmodage.preprocess._ortholog_map | tAge tag 1.5.0 metadata/Table_of_orthologs.csv | Entrez human to Entrez mouse; drop_duplicates keep first | per-sample z-score then subtract reference-sample median; microarray log base treated as irrelevant because of that scaling; RNA-seq uses log10(RLE+1) | counts: at least 10 counts in at least 20 percent of samples; microarray: drop genes below the 0.25 quantile of mean expression | instruction 3 section 7.1: log2 microarray or log2(TPM+1); same ortholog method; relative expression = log minus control-sample gene mean; then predict |

| repo | commit | describe |
| --- | --- | --- |
| /HDD8T1/WORK/GeroSci/clock-specificity | 9d109499f0bf47a7869cd7b4bcffe311b8be9b75 | 9d10949 |
| /HDD8T1/WORK/GeroSci/sc-module-aging | a0e9d090cbabcdde441a599bb6443ed3f5ba0350 | a0e9d09 |
| /HDD8T1/WORK/GeroSci/sc-module-aging/data/external/tAge | d048ba6016a7e5a114f45195e0b11071e2b8c424 | 1.5.0 |

### 5.5 읽기 환경, 구조 확인 표, 포함률

파이프라인 형식은 sklearn `Pipeline`이다. 단계는 imputation `SimpleImputer`, scaler `StandardScaler`, featureselection `SelectKBest`, estimator `ElasticNet`. 특징 ID는 십진 문자열이고, 코드 설명은 mouse Entrez다. `allmodulegenes` 파일은 23개 모듈에 넣지 않았다. 출력 단위는 chronological이 fraction of maximum lifespan, mortality가 log10 hazard ratio다. 버전 경고 원문은 없다.

| outcome | n_modules | paper_n_modules | n_unique_genes | paper_n_unique_genes | n_genes_in_more_than_one_module | paper_overlap | output_unit |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Chronoage | 23 | 23 | 2007 | 2007 | 0 | 0 | fraction of maximum lifespan |
| Mortality | 23 | 23 | 2007 | 2007 | 0 | 0 | log10 hazard ratio |

| outcome | n_features | paper_n_features | file |
| --- | --- | --- | --- |
| Chronoage | 16517 | 10487 | EN_Chronoage_Rodents_Multitissue_scaleddiff.pkl |
| Mortality | 16760 | 10487 | EN_Mortality_Rodents_Multitissue_scaleddiff.pkl |

| n_chronoage_features | n_mortality_features | n_intersection | n_union | paper_n_features |
| --- | --- | --- | --- | --- |
| 16517 | 16760 | 16441 | 16836 | 10487 |

플랫폼 포함률은 오솔로그 대응 뒤 파이프라인 특징이 유전자 목록에 있는 수다. 유전자 이름은 적지 않는다.

| outcome | module | gene_list | n_list_symbols | n_list_mouse_entrez | n_features | n_features_present | fraction_present |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Chronoage | composite | GPL15331 | 19541 | 14120 | 16517 | 13280 | 0.8040201005025126 |
| Chronoage | composite | GPL10379 | 21003 | 15190 | 16517 | 14010 | 0.848216988557244 |
| Chronoage | composite | GPL24324 | 20179 | 16115 | 16517 | 14756 | 0.8933825755282436 |
| Chronoage | composite | GSE113883 | 35581 | 16063 | 16517 | 14715 | 0.8909002845553067 |
| Mortality | composite | GPL15331 | 19541 | 14120 | 16760 | 13368 | 0.7976133651551313 |
| Mortality | composite | GPL10379 | 21003 | 15190 | 16760 | 14100 | 0.8412887828162291 |
| Mortality | composite | GPL24324 | 20179 | 16115 | 16760 | 14864 | 0.886873508353222 |
| Mortality | composite | GSE113883 | 35581 | 16063 | 16760 | 14822 | 0.8843675417661098 |

| outcome | gene_list | n_modules | fraction_min | fraction_median | fraction_max | n_features_present_min | n_features_present_max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Chronoage | GPL10379 | 23 | 0.6779661016949152 | 0.9032258064516129 | 1.0 | 9 | 278 |
| Chronoage | GPL15331 | 23 | 0.6440677966101694 | 0.8809523809523809 | 1.0 | 10 | 276 |
| Chronoage | GPL24324 | 23 | 0.6779661016949152 | 0.9259259259259259 | 1.0 | 10 | 319 |
| Chronoage | GSE113883 | 23 | 0.6779661016949152 | 0.9259259259259259 | 1.0 | 10 | 318 |
| Mortality | GPL10379 | 23 | 0.6779661016949152 | 0.9032258064516129 | 1.0 | 9 | 278 |
| Mortality | GPL15331 | 23 | 0.6440677966101694 | 0.8809523809523809 | 1.0 | 10 | 276 |
| Mortality | GPL24324 | 23 | 0.6779661016949152 | 0.9259259259259259 | 1.0 | 10 | 319 |
| Mortality | GSE113883 | 23 | 0.6779661016949152 | 0.9259259259259259 | 1.0 | 10 | 318 |

모듈별 포함 수:

| outcome | module | gene_list | n_features | n_features_present | fraction_present |
| --- | --- | --- | --- | --- | --- |
| Chronoage | black | GPL10379 | 42 | 38 | 0.9047619047619048 |
| Chronoage | black | GPL15331 | 42 | 37 | 0.8809523809523809 |
| Chronoage | black | GPL24324 | 42 | 38 | 0.9047619047619048 |
| Chronoage | black | GSE113883 | 42 | 38 | 0.9047619047619048 |
| Chronoage | blue | GPL10379 | 357 | 278 | 0.7787114845938375 |
| Chronoage | blue | GPL15331 | 357 | 268 | 0.7507002801120448 |
| Chronoage | blue | GPL24324 | 357 | 319 | 0.8935574229691877 |
| Chronoage | blue | GSE113883 | 357 | 318 | 0.8907563025210085 |
| Chronoage | brown | GPL10379 | 255 | 239 | 0.9372549019607843 |
| Chronoage | brown | GPL15331 | 255 | 224 | 0.8784313725490196 |
| Chronoage | brown | GPL24324 | 255 | 245 | 0.9607843137254902 |
| Chronoage | brown | GSE113883 | 255 | 245 | 0.9607843137254902 |
| Chronoage | cyan | GPL10379 | 20 | 20 | 1.0 |
| Chronoage | cyan | GPL15331 | 20 | 20 | 1.0 |
| Chronoage | cyan | GPL24324 | 20 | 20 | 1.0 |
| Chronoage | cyan | GSE113883 | 20 | 20 | 1.0 |
| Chronoage | darkgreen | GPL10379 | 36 | 35 | 0.9722222222222222 |
| Chronoage | darkgreen | GPL15331 | 36 | 33 | 0.9166666666666666 |
| Chronoage | darkgreen | GPL24324 | 36 | 35 | 0.9722222222222222 |
| Chronoage | darkgreen | GSE113883 | 36 | 35 | 0.9722222222222222 |
| Chronoage | darkgrey | GPL10379 | 31 | 28 | 0.9032258064516129 |
| Chronoage | darkgrey | GPL15331 | 31 | 31 | 1.0 |
| Chronoage | darkgrey | GPL24324 | 31 | 31 | 1.0 |
| Chronoage | darkgrey | GSE113883 | 31 | 31 | 1.0 |
| Chronoage | darkorange | GPL10379 | 17 | 15 | 0.8823529411764706 |
| Chronoage | darkorange | GPL15331 | 17 | 15 | 0.8823529411764706 |
| Chronoage | darkorange | GPL24324 | 17 | 16 | 0.9411764705882353 |
| Chronoage | darkorange | GSE113883 | 17 | 16 | 0.9411764705882353 |
| Chronoage | darkred | GPL10379 | 26 | 24 | 0.9230769230769231 |
| Chronoage | darkred | GPL15331 | 26 | 22 | 0.8461538461538461 |
| Chronoage | darkred | GPL24324 | 26 | 24 | 0.9230769230769231 |
| Chronoage | darkred | GSE113883 | 26 | 24 | 0.9230769230769231 |
| Chronoage | green | GPL10379 | 97 | 94 | 0.9690721649484536 |
| Chronoage | green | GPL15331 | 97 | 93 | 0.9587628865979382 |
| Chronoage | green | GPL24324 | 97 | 96 | 0.9896907216494846 |
| Chronoage | green | GSE113883 | 97 | 96 | 0.9896907216494846 |
| Chronoage | greenyellow | GPL10379 | 59 | 57 | 0.9661016949152542 |
| Chronoage | greenyellow | GPL15331 | 59 | 58 | 0.9830508474576272 |
| Chronoage | greenyellow | GPL24324 | 59 | 59 | 1.0 |
| Chronoage | greenyellow | GSE113883 | 59 | 59 | 1.0 |
| Chronoage | grey60 | GPL10379 | 52 | 45 | 0.8653846153846154 |
| Chronoage | grey60 | GPL15331 | 52 | 46 | 0.8846153846153846 |
| Chronoage | grey60 | GPL24324 | 52 | 48 | 0.9230769230769231 |
| Chronoage | grey60 | GSE113883 | 52 | 48 | 0.9230769230769231 |
| Chronoage | lightgreen | GPL10379 | 59 | 40 | 0.6779661016949152 |
| Chronoage | lightgreen | GPL15331 | 59 | 38 | 0.6440677966101694 |
| Chronoage | lightgreen | GPL24324 | 59 | 40 | 0.6779661016949152 |
| Chronoage | lightgreen | GSE113883 | 59 | 40 | 0.6779661016949152 |
| Chronoage | lightyellow | GPL10379 | 16 | 12 | 0.75 |
| Chronoage | lightyellow | GPL15331 | 16 | 11 | 0.6875 |
| Chronoage | lightyellow | GPL24324 | 16 | 13 | 0.8125 |
| Chronoage | lightyellow | GSE113883 | 16 | 13 | 0.8125 |
| Chronoage | magenta | GPL10379 | 54 | 49 | 0.9074074074074074 |
| Chronoage | magenta | GPL15331 | 54 | 49 | 0.9074074074074074 |
| Chronoage | magenta | GPL24324 | 54 | 50 | 0.9259259259259259 |
| Chronoage | magenta | GSE113883 | 54 | 50 | 0.9259259259259259 |
| Chronoage | midnightblue | GPL10379 | 52 | 42 | 0.8076923076923077 |
| Chronoage | midnightblue | GPL15331 | 52 | 39 | 0.75 |
| Chronoage | midnightblue | GPL24324 | 52 | 42 | 0.8076923076923077 |
| Chronoage | midnightblue | GSE113883 | 52 | 42 | 0.8076923076923077 |
| Chronoage | orange | GPL10379 | 30 | 28 | 0.9333333333333333 |
| Chronoage | orange | GPL15331 | 30 | 27 | 0.9 |
| Chronoage | orange | GPL24324 | 30 | 28 | 0.9333333333333333 |
| Chronoage | orange | GSE113883 | 30 | 28 | 0.9333333333333333 |
| Chronoage | pink | GPL10379 | 93 | 88 | 0.946236559139785 |
| Chronoage | pink | GPL15331 | 93 | 86 | 0.9247311827956989 |
| Chronoage | pink | GPL24324 | 93 | 91 | 0.978494623655914 |
| Chronoage | pink | GSE113883 | 93 | 91 | 0.978494623655914 |
| Chronoage | red | GPL10379 | 117 | 91 | 0.7777777777777778 |
| Chronoage | red | GPL15331 | 117 | 96 | 0.8205128205128205 |
| Chronoage | red | GPL24324 | 117 | 100 | 0.8547008547008547 |
| Chronoage | red | GSE113883 | 117 | 101 | 0.8632478632478633 |
| Chronoage | salmon | GPL10379 | 12 | 9 | 0.75 |
| Chronoage | salmon | GPL15331 | 12 | 10 | 0.8333333333333334 |
| Chronoage | salmon | GPL24324 | 12 | 10 | 0.8333333333333334 |
| Chronoage | salmon | GSE113883 | 12 | 10 | 0.8333333333333334 |
| Chronoage | skyblue | GPL10379 | 20 | 18 | 0.9 |
| Chronoage | skyblue | GPL15331 | 20 | 18 | 0.9 |
| Chronoage | skyblue | GPL24324 | 20 | 18 | 0.9 |
| Chronoage | skyblue | GSE113883 | 20 | 18 | 0.9 |
| Chronoage | tan | GPL10379 | 61 | 51 | 0.8360655737704918 |
| Chronoage | tan | GPL15331 | 61 | 52 | 0.8524590163934426 |
| Chronoage | tan | GPL24324 | 61 | 59 | 0.9672131147540983 |
| Chronoage | tan | GSE113883 | 61 | 58 | 0.9508196721311475 |
| Chronoage | turquoise | GPL10379 | 332 | 275 | 0.8283132530120482 |
| Chronoage | turquoise | GPL15331 | 332 | 276 | 0.8313253012048193 |
| Chronoage | turquoise | GPL24324 | 332 | 293 | 0.8825301204819277 |
| Chronoage | turquoise | GSE113883 | 332 | 292 | 0.8795180722891566 |
| Chronoage | yellow | GPL10379 | 169 | 153 | 0.9053254437869822 |
| Chronoage | yellow | GPL15331 | 169 | 142 | 0.8402366863905325 |
| Chronoage | yellow | GPL24324 | 169 | 158 | 0.9349112426035503 |
| Chronoage | yellow | GSE113883 | 169 | 157 | 0.9289940828402367 |
| Mortality | black | GPL10379 | 42 | 38 | 0.9047619047619048 |
| Mortality | black | GPL15331 | 42 | 37 | 0.8809523809523809 |
| Mortality | black | GPL24324 | 42 | 38 | 0.9047619047619048 |
| Mortality | black | GSE113883 | 42 | 38 | 0.9047619047619048 |
| Mortality | blue | GPL10379 | 357 | 278 | 0.7787114845938375 |
| Mortality | blue | GPL15331 | 357 | 268 | 0.7507002801120448 |
| Mortality | blue | GPL24324 | 357 | 319 | 0.8935574229691877 |
| Mortality | blue | GSE113883 | 357 | 318 | 0.8907563025210085 |
| Mortality | brown | GPL10379 | 255 | 239 | 0.9372549019607843 |
| Mortality | brown | GPL15331 | 255 | 224 | 0.8784313725490196 |
| Mortality | brown | GPL24324 | 255 | 245 | 0.9607843137254902 |
| Mortality | brown | GSE113883 | 255 | 245 | 0.9607843137254902 |
| Mortality | cyan | GPL10379 | 20 | 20 | 1.0 |
| Mortality | cyan | GPL15331 | 20 | 20 | 1.0 |
| Mortality | cyan | GPL24324 | 20 | 20 | 1.0 |
| Mortality | cyan | GSE113883 | 20 | 20 | 1.0 |
| Mortality | darkgreen | GPL10379 | 36 | 35 | 0.9722222222222222 |
| Mortality | darkgreen | GPL15331 | 36 | 33 | 0.9166666666666666 |
| Mortality | darkgreen | GPL24324 | 36 | 35 | 0.9722222222222222 |
| Mortality | darkgreen | GSE113883 | 36 | 35 | 0.9722222222222222 |
| Mortality | darkgrey | GPL10379 | 31 | 28 | 0.9032258064516129 |
| Mortality | darkgrey | GPL15331 | 31 | 31 | 1.0 |
| Mortality | darkgrey | GPL24324 | 31 | 31 | 1.0 |
| Mortality | darkgrey | GSE113883 | 31 | 31 | 1.0 |
| Mortality | darkorange | GPL10379 | 17 | 15 | 0.8823529411764706 |
| Mortality | darkorange | GPL15331 | 17 | 15 | 0.8823529411764706 |
| Mortality | darkorange | GPL24324 | 17 | 16 | 0.9411764705882353 |
| Mortality | darkorange | GSE113883 | 17 | 16 | 0.9411764705882353 |
| Mortality | darkred | GPL10379 | 26 | 24 | 0.9230769230769231 |
| Mortality | darkred | GPL15331 | 26 | 22 | 0.8461538461538461 |
| Mortality | darkred | GPL24324 | 26 | 24 | 0.9230769230769231 |
| Mortality | darkred | GSE113883 | 26 | 24 | 0.9230769230769231 |
| Mortality | green | GPL10379 | 97 | 94 | 0.9690721649484536 |
| Mortality | green | GPL15331 | 97 | 93 | 0.9587628865979382 |
| Mortality | green | GPL24324 | 97 | 96 | 0.9896907216494846 |
| Mortality | green | GSE113883 | 97 | 96 | 0.9896907216494846 |
| Mortality | greenyellow | GPL10379 | 59 | 57 | 0.9661016949152542 |
| Mortality | greenyellow | GPL15331 | 59 | 58 | 0.9830508474576272 |
| Mortality | greenyellow | GPL24324 | 59 | 59 | 1.0 |
| Mortality | greenyellow | GSE113883 | 59 | 59 | 1.0 |
| Mortality | grey60 | GPL10379 | 52 | 45 | 0.8653846153846154 |
| Mortality | grey60 | GPL15331 | 52 | 46 | 0.8846153846153846 |
| Mortality | grey60 | GPL24324 | 52 | 48 | 0.9230769230769231 |
| Mortality | grey60 | GSE113883 | 52 | 48 | 0.9230769230769231 |
| Mortality | lightgreen | GPL10379 | 59 | 40 | 0.6779661016949152 |
| Mortality | lightgreen | GPL15331 | 59 | 38 | 0.6440677966101694 |
| Mortality | lightgreen | GPL24324 | 59 | 40 | 0.6779661016949152 |
| Mortality | lightgreen | GSE113883 | 59 | 40 | 0.6779661016949152 |
| Mortality | lightyellow | GPL10379 | 16 | 12 | 0.75 |
| Mortality | lightyellow | GPL15331 | 16 | 11 | 0.6875 |
| Mortality | lightyellow | GPL24324 | 16 | 13 | 0.8125 |
| Mortality | lightyellow | GSE113883 | 16 | 13 | 0.8125 |
| Mortality | magenta | GPL10379 | 54 | 49 | 0.9074074074074074 |
| Mortality | magenta | GPL15331 | 54 | 49 | 0.9074074074074074 |
| Mortality | magenta | GPL24324 | 54 | 50 | 0.9259259259259259 |
| Mortality | magenta | GSE113883 | 54 | 50 | 0.9259259259259259 |
| Mortality | midnightblue | GPL10379 | 52 | 42 | 0.8076923076923077 |
| Mortality | midnightblue | GPL15331 | 52 | 39 | 0.75 |
| Mortality | midnightblue | GPL24324 | 52 | 42 | 0.8076923076923077 |
| Mortality | midnightblue | GSE113883 | 52 | 42 | 0.8076923076923077 |
| Mortality | orange | GPL10379 | 30 | 28 | 0.9333333333333333 |
| Mortality | orange | GPL15331 | 30 | 27 | 0.9 |
| Mortality | orange | GPL24324 | 30 | 28 | 0.9333333333333333 |
| Mortality | orange | GSE113883 | 30 | 28 | 0.9333333333333333 |
| Mortality | pink | GPL10379 | 93 | 88 | 0.946236559139785 |
| Mortality | pink | GPL15331 | 93 | 86 | 0.9247311827956989 |
| Mortality | pink | GPL24324 | 93 | 91 | 0.978494623655914 |
| Mortality | pink | GSE113883 | 93 | 91 | 0.978494623655914 |
| Mortality | red | GPL10379 | 117 | 91 | 0.7777777777777778 |
| Mortality | red | GPL15331 | 117 | 96 | 0.8205128205128205 |
| Mortality | red | GPL24324 | 117 | 100 | 0.8547008547008547 |
| Mortality | red | GSE113883 | 117 | 101 | 0.8632478632478633 |
| Mortality | salmon | GPL10379 | 12 | 9 | 0.75 |
| Mortality | salmon | GPL15331 | 12 | 10 | 0.8333333333333334 |
| Mortality | salmon | GPL24324 | 12 | 10 | 0.8333333333333334 |
| Mortality | salmon | GSE113883 | 12 | 10 | 0.8333333333333334 |
| Mortality | skyblue | GPL10379 | 20 | 18 | 0.9 |
| Mortality | skyblue | GPL15331 | 20 | 18 | 0.9 |
| Mortality | skyblue | GPL24324 | 20 | 18 | 0.9 |
| Mortality | skyblue | GSE113883 | 20 | 18 | 0.9 |
| Mortality | tan | GPL10379 | 61 | 51 | 0.8360655737704918 |
| Mortality | tan | GPL15331 | 61 | 52 | 0.8524590163934426 |
| Mortality | tan | GPL24324 | 61 | 59 | 0.9672131147540983 |
| Mortality | tan | GSE113883 | 61 | 58 | 0.9508196721311475 |
| Mortality | turquoise | GPL10379 | 332 | 275 | 0.8283132530120482 |
| Mortality | turquoise | GPL15331 | 332 | 276 | 0.8313253012048193 |
| Mortality | turquoise | GPL24324 | 332 | 293 | 0.8825301204819277 |
| Mortality | turquoise | GSE113883 | 332 | 292 | 0.8795180722891566 |
| Mortality | yellow | GPL10379 | 169 | 153 | 0.9053254437869822 |
| Mortality | yellow | GPL15331 | 169 | 142 | 0.8402366863905325 |
| Mortality | yellow | GPL24324 | 169 | 158 | 0.9349112426035503 |
| Mortality | yellow | GSE113883 | 169 | 157 | 0.9289940828402367 |

라이선스 제목은 MGB OPEN ACCESS LICENSE 1.0. 모듈 구성 유전자와 계수는 결과 파일과 이 보고서에 적지 않았다.

## 6. GSE113883 RNA-seq 개발

대상은 GSE113883 153샘플, 11명이다. 조건은 전부 control이다. GSE196728은 읽지 않았다. t = `clock_time_h`.

### 6.1 QC 재실행

샘플 중앙값은 mean TPM ≥ 1인 유전자의 `log2(TPM+1)`로 계산했다. MAD는 중앙값 절대편차이고 1.4826을 곱하지 않았다. 표시만 하고 샘플을 빼지 않았다.

| dataset | value_scale | gene_filter | n_genes | n_samples | median_of_sample_medians | mad_of_sample_medians | mad_constant | n_flagged_3mad | flagged_gsm |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE113883 | log2(TPM+1) | mean TPM >= 1 | 12390 | 153 | 3.69376571221778 | 0.18397853773122 | none | 6 | GSM3122552,GSM3122587,GSM3122617,GSM3122624,GSM3122625,GSM3122643 |

### 6.2 시계 값 (포함 수, 대체 수, 분포, official–lin 상관)

lin은 절편 + 계수 · `log2(TPM+1)`이다. 행렬에 없는 유전자는 빼고 대체하지 않는다. official은 `predict_age`에 `exprtype="FPKM"`, `idtype="SYMBOL"`, `tissue="blood"`, `stype="all"`을 넘겼다. 입력 값은 TPM이다. 대체 수는 서명 유전자 중 행렬에 없는 수다.

| clock | method | n | mean | sd | min | max | n_genes_used | n_genes_signature | n_genes_imputed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | lin | 153 | 62.7599299733162 | 13.2772548430179 | 32.7300518005297 | 95.849546253573 | 271 | 312 |  |
| RA_GTEx | official | 153 | 47.2878576885725 | 12.9942271536497 | 16.8577310425423 | 79.0150833774752 | 271 | 312 | 41 |
| RA_Peters | lin | 153 | 39.6065298952733 | 14.5136881716668 | 6.89986370272937 | 76.3563319503402 | 276 | 289 |  |
| RA_Peters | official | 153 | 32.6553543699098 | 14.4088509859448 | -0.26635732992775 | 69.1389729912341 | 276 | 289 | 13 |
| RA_default | lin | 153 | 37.8218395219532 | 9.41556090582936 | 14.1438176938915 | 63.5447753341867 | 242 | 279 |  |
| RA_default | official | 153 | 29.3950762596875 | 9.42280561208904 | 5.46839193952964 | 56.6544030928359 | 242 | 279 | 37 |

| clock | pearson_within | sd_within_official_minus_lin | n | use_lin_only |
| --- | --- | --- | --- | --- |
| RA_GTEx | 0.997677855578136 | 0.593570220352772 | 153 | TRUE |
| RA_Peters | 0.998982885935051 | 0.266674735766916 | 153 | TRUE |
| RA_default | 0.999014038243109 | 0.323284149297551 | 153 | TRUE |

세 시계의 개인 내 Pearson이 0.99 이상이라 이후 전체 분석은 lin이다. official은 점추정만 남겼다.

### 6.3 주 모형 표 (보고서 2 §5.2와 같은 열)

n_subjects 11, n_samples 153. lin 구간은 부트스트랩 2,000회, 시드 20261015, 2.5·97.5 백분위 type 7이다. 정점 괄호는 `peak_center`와 `halfwidth_95`다. 비수렴은 0이다. 점추정의 REML·ML·귀무 수렴은 TRUE다.

| clock | method | role | beta_cos | beta_sin | 2A | R2_within | 2A/SD_between | 2A/SD_resid | peak_h (center, half) | LRT p | ICC | boot_converged | reml_singular |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | lin | primary | -0.0956555930233394 | 2.16973955504679 | 4.34369415553687 [2.14402772768952, 10.7945031669817] | 0.0315491727978717 [0.00905836854543264, 0.14463204052845] | 0.405335937361681 [0.18429339925241, 1.25692476563122] | 0.488103052622962 [0.255318611205345, 1.10479018672001] | 6.1682879453552 (6.7371459498579, 4.82078466749953) | 0.104384139695722 | 0.554021910692076 | 2000/2000 | FALSE |
| RA_GTEx | official | appendix_point | -0.336248241223975 | 1.92722940101417 | 3.9126850340193 [미계산] | 0.0269579334362133 [미계산] | 0.374888182397042 [미계산] | 0.445661721431295 [미계산] | 6.65979393801801 (미계산) | 0.14519648763658 | 0.547197554132143 | 미계산 | FALSE |
| RA_Peters | lin | primary | -1.57614289141543 | 0.0676221643784579 | 3.15518568155639 [0.724064511975616, 6.78626343037065] | 0.03527416607544 [0.00202148111734601, 0.138905271852363] | 0.229721092919869 [0.0659678259118176, 0.660877698468735] | 0.516102685110757 [0.125429103071359, 1.03339812617933] | 11.8362208570621 (11.6741391650562, 4.62490864529145) | 0.077344312576993 | 0.818809345054716 | 2000/2000 | FALSE |
| RA_Peters | official | appendix_point | -1.68410808342564 | -0.0415097396670217 | 3.36923913971497 [미계산] | 0.0395884490849726 [미계산] | 0.247468227841035 [미계산] | 0.55064158884193 [미계산] | 12.0941290115822 (미계산) | 0.0561796678896831 | 0.815886288922892 | 미계산 | FALSE |
| RA_default | lin | primary | -3.38101122840144 | 0.165981334106488 | 6.77016594474564 [4.54731764000349, 9.17576797367988] | 0.134425171769798 [0.0636149012292298, 0.247718423922361] | 0.954019901666464 [0.635341349852203, 1.77170408171584] | 1.00862397715271 [0.708562523439538, 1.37792588545123] | 11.8126319723118 (11.8212982069323, 1.44667958367679) | 3.24675832907807e-05 | 0.486008334841039 | 2000/2000 | FALSE |
| RA_default | official | appendix_point | -3.51267648401214 | 0.0345667118031814 | 7.02569311567101 [미계산] | 0.137544944472577 [미계산] | 1.00723774274417 [미계산] | 1.02451533976697 [미계산] | 11.9624130340213 (미계산) | 2.49979371279294e-05 | 0.465115328178416 | 미계산 | FALSE |

### 6.4 무작위 서명, S2, S3

무작위 서명은 자료의 전체 샘플 `log2(TPM+1)` 평균으로 10분위를 나누고, 시계 유전자가 아닌 같은 분위에서 서명 안에서는 비복원으로 뽑았다. 계수는 섞고 절편은 0이다. 시드 20261015. 1,000개. 백분위는 100 × (무작위 R2_within ≤ 시계 R2_within)의 평균이다.

| clock | n_random | n_finite | n_not_converged | random_median | random_p05 | random_p95 | clock_r2_within | percentile_rank |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | 1000 | 1000 | 0 | 0.0819958688009636 | 0.00725006666787109 | 0.22889008286771 | 0.0315491727978717 | 20.3 |
| RA_Peters | 1000 | 1000 | 0 | 0.0832288408686809 | 0.00827492373235772 | 0.236624516362962 | 0.03527416607544 | 22.9 |
| RA_default | 1000 | 1000 | 0 | 0.0774110681374275 | 0.00625949107422514 | 0.25615075081854 | 0.134425171769798 | 72.9 |

S2 귀무는 `score ~ 1 + (1 + cos_t + sin_t | subject_id)`다. 기울기 항이 있어 `icc_null`은 비어 있다. S1, S4, S5는 하지 않았다.

| clock | n_subjects | n_samples | two_a | r2_within | peak_h | lrt_p | icc_null | two_a_over_sd_between | reml_converged | ml_converged | null_converged | reml_singular | null_formula |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | 11 | 153 | 4.3487585309192 | 0.0247387026502068 | 6.39413214510532 | 0.144144617714443 |  | 0.40580852435083 | TRUE | TRUE | TRUE | TRUE | score ~ 1 + (1 + cos_t + sin_t / subject_id) |
| RA_Peters | 11 | 153 | 3.22861812794821 | 0.00301353084523837 | 11.9308044538506 | 0.209886322370578 |  | 0.235067523698734 | TRUE | TRUE | TRUE | TRUE | score ~ 1 + (1 + cos_t + sin_t / subject_id) |
| RA_default | 11 | 153 | 6.71732945455258 | 0.00523418351006422 | 11.7929917360237 | 0.00407334556958206 |  | 0.946574432295436 | TRUE | TRUE | TRUE | TRUE | score ~ 1 + (1 + cos_t + sin_t / subject_id) |

S3은 참가자별 OLS다. 샘플 ≥ 5이고 `clock_time_h`의 max−min ≥ 16인 참가자만 썼다. 2A는 참가자별 2A의 중앙값, 괄호에 해당하는 열은 25·75 백분위다. R2_within 열은 참가자별 OLS R²의 중앙값이다. 정점은 참가자 정점의 원형 평균이다. 혼합모형 LRT는 맞추지 않았다. Rayleigh는 `circular::rayleigh.test`다.

| clock | n_subjects | n_samples | two_a | two_a_q25 | two_a_q75 | r2_within | peak_h | rayleigh_stat | rayleigh_p | clock_range_definition |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | 11 | 153 | 9.80688743022131 | 4.86201842318975 | 11.1197728611371 | 0.138013074484506 | 8.74321899382482 | 0.536617110488638 | 0.0382527207811252 | max(clock_time_h) - min(clock_time_h) |
| RA_Peters | 11 | 153 | 5.33913843038554 | 3.64466795792879 | 8.63352950534084 | 0.127018820698532 | 12.2667352011868 | 0.372335990733508 | 0.22151868254916 | max(clock_time_h) - min(clock_time_h) |
| RA_default | 11 | 153 | 7.92531712439221 | 4.88806800160775 | 11.2862449137707 | 0.224844881920846 | 12.1669063930325 | 0.724739069979986 | 0.00152800179344408 | max(clock_time_h) - min(clock_time_h) |

### 6.5 구성 점수와 하루 리듬, 투영 매개 비율, 시계–구성 상관

MCPcounter 1.2.0. 입력은 `log2(TPM+1)`, `featuresType="HUGO_symbols"`, 유전자 파일은 `data/processed/ref/MCPcounter_genes.txt`다. 분석한 혈액 집단은 8개다. Endothelial cells와 Fibroblasts는 점수만 있고 리듬 모형에는 넣지 않았다. |r| > 0.9인 쌍은 빼지 않았다.

| population | n_marker | n_marker_in_matrix |
| --- | --- | --- |
| T cells | 16 | 15 |
| CD8 T cells | 1 | 1 |
| Cytotoxic lymphocytes | 7 | 7 |
| B lineage | 9 | 9 |
| NK cells | 9 | 8 |
| Monocytic lineage | 7 | 7 |
| Myeloid dendritic cells | 6 | 5 |
| Neutrophils | 15 | 13 |
| Endothelial cells | 33 | 32 |
| Fibroblasts | 8 | 8 |

| population | n_subjects | n_samples | two_a | r2_within | peak_h | lrt_p | icc_null | reml_converged | ml_converged | null_converged |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| T cells | 11 | 153 | 0.501333019807266 | 0.285597123119417 | 1.74942881938859 | 4.99931417782781e-11 | 0.354974209283768 | TRUE | TRUE | TRUE |
| CD8 T cells | 11 | 153 | 0.419891128542506 | 0.226446703839558 | 0.619346944281585 | 1.15565194309793e-08 | 0.813778663913656 | TRUE | TRUE | TRUE |
| Cytotoxic lymphocytes | 11 | 153 | 0.243090121688901 | 0.0972174260922877 | 9.50166976042668 | 0.000732999198460932 | 0.474900414196802 | TRUE | TRUE | TRUE |
| NK cells | 11 | 153 | 0.485126215137497 | 0.306458361396565 | 11.3747771559795 | 5.23804105419162e-12 | 0.432139625912167 | TRUE | TRUE | TRUE |
| B lineage | 11 | 153 | 0.332784618573047 | 0.225011095585362 | 23.9510317757202 | 1.42327204249319e-08 | 0.615937120299878 | TRUE | TRUE | TRUE |
| Monocytic lineage | 11 | 153 | 0.123338435985521 | 0.0474975712537372 | 3.5632959277126 | 0.0315096138649207 | 0.739498660856843 | TRUE | TRUE | TRUE |
| Myeloid dendritic cells | 11 | 153 | 0.125222679215964 | 0.0578355659338845 | 0.75409377990084 | 0.014113491536299 | 0.606827459088031 | TRUE | TRUE | TRUE |
| Neutrophils | 11 | 153 | 0.28398112036914 | 0.0839649563291982 | 13.8007121554084 | 0.00187867425470755 | 0.504420106405589 | TRUE | TRUE | TRUE |

| population_a | population_b | pearson_r | abs_r_gt_0_9 |
| --- | --- | --- | --- |
| T cells | CD8 T cells | 0.317190619089936 | FALSE |
| T cells | Cytotoxic lymphocytes | 0.322017128370196 | FALSE |
| T cells | NK cells | 0.163744385411754 | FALSE |
| T cells | B lineage | 0.376145222159579 | FALSE |
| T cells | Monocytic lineage | 0.146144446365591 | FALSE |
| T cells | Myeloid dendritic cells | 0.0645149159924547 | FALSE |
| T cells | Neutrophils | -0.173207815354377 | FALSE |
| CD8 T cells | Cytotoxic lymphocytes | 0.16805375608272 | FALSE |
| CD8 T cells | NK cells | -0.184267790579194 | FALSE |
| CD8 T cells | B lineage | -0.364852743332035 | FALSE |
| CD8 T cells | Monocytic lineage | -0.584891072200442 | FALSE |
| CD8 T cells | Myeloid dendritic cells | -0.147964009862876 | FALSE |
| CD8 T cells | Neutrophils | -0.402062300321167 | FALSE |
| Cytotoxic lymphocytes | NK cells | 0.77296535386469 | FALSE |
| Cytotoxic lymphocytes | B lineage | 0.0614215992175894 | FALSE |
| Cytotoxic lymphocytes | Monocytic lineage | 0.333996068952386 | FALSE |
| Cytotoxic lymphocytes | Myeloid dendritic cells | 0.231904733979738 | FALSE |
| Cytotoxic lymphocytes | Neutrophils | -0.0347638759237555 | FALSE |
| NK cells | B lineage | 0.0485187033371649 | FALSE |
| NK cells | Monocytic lineage | 0.31674252202641 | FALSE |
| NK cells | Myeloid dendritic cells | -0.00863176666879634 | FALSE |
| NK cells | Neutrophils | 0.311434222537657 | FALSE |
| B lineage | Monocytic lineage | 0.693843644466909 | FALSE |
| B lineage | Myeloid dendritic cells | 0.518189045351404 | FALSE |
| B lineage | Neutrophils | 0.017634516440665 | FALSE |
| Monocytic lineage | Myeloid dendritic cells | 0.661693224105589 | FALSE |
| Monocytic lineage | Neutrophils | 0.23833440337456 | FALSE |
| Myeloid dendritic cells | Neutrophils | -0.0443538904508809 | FALSE |

투영은 비조정 REML의 (β_cos, β_sin)을 r_u, 8개 점수 조정 REML의 계수를 r_a로 둔다. `p_med` = 1 − (r_a · r_u) / |r_u|². `two_a_orth` = 2|r_a − proj|. 부트스트랩은 §6.3과 같은 2,000개 추출이다. 같은 추출에서 두 모형을 맞춘다. RNA-seq에는 2A 감소 비율을 계산하지 않았다.

| clock | n_boot | n_converged_both_reml | p_med | p_med_p025 | p_med_p975 | two_a_orth | two_a_orth_p025 | two_a_orth_p975 | lrt_p_unadj | ru_norm |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | 2000 | 2000 | -0.149683650058583 | -0.549633569275945 | 1.49242245851642 | 0.345308977625604 | 0.0611309367272339 | 5.19307959744199 | 0.104384139695722 | 2.17184707776843 |
| RA_Peters | 2000 | 2000 | 0.812219272485387 | -0.47918439806038 | 3.41052546250761 | 1.35899567483826 | 0.0712273990076862 | 4.23284571533323 | 0.077344312576993 | 1.57759284077819 |
| RA_default | 2000 | 2000 | 0.870626030579166 | 0.464188530693954 | 1.54833140363586 | 3.41508275044276 | 0.610507876317039 | 6.69625276967478 | 3.24675832907807e-05 | 3.38508297237282 |

| clock | population | n_samples | n_subjects | pearson_within |
| --- | --- | --- | --- | --- |
| RA_GTEx | T cells | 153 | 11 | 0.146042777905975 |
| RA_GTEx | CD8 T cells | 153 | 11 | 0.0830889599521202 |
| RA_GTEx | Cytotoxic lymphocytes | 153 | 11 | 0.119842167187801 |
| RA_GTEx | NK cells | 153 | 11 | 0.0823564010207206 |
| RA_GTEx | B lineage | 153 | 11 | 0.214293512754773 |
| RA_GTEx | Monocytic lineage | 153 | 11 | 0.152372073887921 |
| RA_GTEx | Myeloid dendritic cells | 153 | 11 | 0.188382558159267 |
| RA_GTEx | Neutrophils | 153 | 11 | -0.022125731803156 |
| RA_Peters | T cells | 153 | 11 | 0.0618202336372216 |
| RA_Peters | CD8 T cells | 153 | 11 | 0.0160246769751935 |
| RA_Peters | Cytotoxic lymphocytes | 153 | 11 | 0.339145631454113 |
| RA_Peters | NK cells | 153 | 11 | 0.259894222354301 |
| RA_Peters | B lineage | 153 | 11 | -0.0796431137270133 |
| RA_Peters | Monocytic lineage | 153 | 11 | -0.0151441096115187 |
| RA_Peters | Myeloid dendritic cells | 153 | 11 | 0.0525090285386409 |
| RA_Peters | Neutrophils | 153 | 11 | -0.0998062958140424 |
| RA_default | T cells | 153 | 11 | -0.497022431722507 |
| RA_default | CD8 T cells | 153 | 11 | -0.408634709083263 |
| RA_default | Cytotoxic lymphocytes | 153 | 11 | -0.120148271049008 |
| RA_default | NK cells | 153 | 11 | 0.157047603029415 |
| RA_default | B lineage | 153 | 11 | -0.314891654185659 |
| RA_default | Monocytic lineage | 153 | 11 | -0.310624717920735 |
| RA_default | Myeloid dendritic cells | 153 | 11 | -0.317494964712713 |
| RA_default | Neutrophils | 153 | 11 | 0.120906454884594 |

### 6.6 아침 − 저녁 시계 차이

창은 07:00–09:00과 19:00–21:00, 양 끝 포함이다. 차이는 아침 − 저녁이다. 두 창에 샘플이 있는 참가자만 썼다. 구간은 대응 t 검정의 95% 구간이다. 유전자 행은 단계 3에 있고 여기서는 lin 시계 행만 더했다.

| outcome | n_subjects | mean_diff | sd_diff | ci95_low | ci95_high | p |
| --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | 11 | 3.38675095929399 | 7.6092984988986 | -1.72524533014459 | 8.49874724873257 | 0.170680785168479 |
| RA_Peters | 11 | 0.956543945170376 | 7.27976563446149 | -3.93406915362397 | 5.84715704396472 | 0.672241284299487 |
| RA_default | 11 | 2.33773972761062 | 7.16235159340991 | -2.47399354510375 | 7.14947300032499 | 0.304432352848547 |

### 6.7 시각 축 비교 (ΔAIC 표)

T_clock과 T_wake는 153샘플, 11명 전부다. T_dlmo는 §4.1에서 하지 않으므로, 7명으로 다시 맞추지 않았다. 결과 변수는 lin 세 시계, PER1, PER3, BMAL1이다. BMAL1의 `matrix_symbol`은 BMAL1이다. ΔAIC는 그 모형의 ML AIC에서 T_clock의 ML AIC를 뺀 값이다.

| outcome | matrix_symbol | axis | time_col | n_subjects | n_samples | two_a | peak_h | lrt_p | r2_within | aic_ml | delta_aic_vs_T_clock | ml_converged | reml_converged |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | RA_GTEx | T_clock | clock_time_h | 11 | 153 | 4.34369415553687 | 6.1682879453552 | 0.104384139695722 | 0.0315491727978717 | 1140.52952888452 | 0 | TRUE | TRUE |
| RA_GTEx | RA_GTEx | T_wake | hours_since_wake | 11 | 153 | 4.3436941555369 | 0.168287945355171 | 0.104384139695722 | 0.0315491728067385 | 1140.52952888452 | 0 | TRUE | TRUE |
| RA_Peters | RA_Peters | T_clock | clock_time_h | 11 | 153 | 3.15518568155639 | 11.8362208570621 | 0.077344312576993 | 0.03527416607544 | 1038.79888821693 | 0 | TRUE | TRUE |
| RA_Peters | RA_Peters | T_wake | hours_since_wake | 11 | 153 | 3.15518568155705 | 5.83622085706282 | 0.0773443137299703 | 0.0352688623398265 | 1038.79888824674 | 2.98141458188184e-08 | TRUE | TRUE |
| RA_default | RA_default | T_clock | clock_time_h | 11 | 153 | 6.77016594474564 | 11.8126319723118 | 3.24675832907807e-05 | 0.134425171769798 | 1035.26720531803 | 0 | TRUE | TRUE |
| RA_default | RA_default | T_wake | hours_since_wake | 11 | 153 | 6.77016594474358 | 5.81263197231228 | 3.24675832907807e-05 | 0.134425171747783 | 1035.26720531803 | 0 | TRUE | TRUE |
| PER1 | PER1 | T_clock | clock_time_h | 11 | 153 | 0.620551996416151 | 6.88888783794359 | 2.75027621424795e-09 | 0.243712744626445 | 187.762080516524 | 0 | TRUE | TRUE |
| PER1 | PER1 | T_wake | hours_since_wake | 11 | 153 | 0.620551996416324 | 0.888887837941677 | 2.75027621424802e-09 | 0.243712749030786 | 187.762080516524 | 5.6843418860808e-14 | TRUE | TRUE |
| PER3 | PER3 | T_clock | clock_time_h | 11 | 153 | 0.449708658716999 | 4.3428731541087 | 7.58295067863991e-11 | 0.278992024827892 | 43.5253470814227 | 0 | TRUE | TRUE |
| PER3 | PER3 | T_wake | hours_since_wake | 11 | 153 | 0.449708658716996 | 22.3428731541088 | 7.5829506786398e-11 | 0.278992013629664 | 43.5253470814227 | -2.8421709430404e-14 | TRUE | TRUE |
| BMAL1 | BMAL1 | T_clock | clock_time_h | 11 | 153 | 0.270966203946098 | 2.29347352520855 | 0.00933679580233329 | 0.064122510682494 | 157.17446749686 | 0 | TRUE | TRUE |
| BMAL1 | BMAL1 | T_wake | hours_since_wake | 11 | 153 | 0.270966203946194 | 20.2934735252097 | 0.00933679580233315 | 0.0641225106819208 | 157.17446749686 | -2.8421709430404e-14 | TRUE | TRUE |

| axis | performed | reason | n_samples_other_axes |
| --- | --- | --- | --- |
| T_dlmo | FALSE | salivary_dlmo_circadian_time values are outside 18-26 | 153 |

### 6.8 그림 경로

`results/stage4/figures/GSE113883_RA_GTEx_lin.png`

`results/stage4/figures/GSE113883_RA_Peters_lin.png`

`results/stage4/figures/GSE113883_RA_default_lin.png`

각 그림은 control 샘플의 lin 점수와 `clock_time_h`, 참가자별 선, REML 코시너 곡선이다.

## 7. 모듈 시계 개발 분석

Part B에서 배치 파일 7개의 md5가 Zenodo와 같았고 채점 경로가 정해져 Part D를 수행했다. 대상은 GSE39445, GSE56931, GSE113883의 control이다. GSE196728은 쓰지 않았다.

### 7.1 입력 준비 (로그 밑, 오솔로그 대응 수, 래퍼 경로, 파이프라인별 특징 비율·대체 비율)

래퍼는 `scripts/27_partD_score.py`다. 마이크로어레이는 단계 1의 log2, RNA-seq은 `log2(TPM+1)`이다. 로그 밑 변환은 하지 않았다. 오솔로그는 `scmodage.preprocess.map_to_mouse_entrez(species="human")`다. 상대 발현은 대응 뒤 값에서 그 자료 control 샘플의 유전자별 평균을 뺀 것이다. 그다음 파이프라인 `predict`만 호출했다. 계수와 파이프라인은 고치지 않았다. 없는 특징은 파이프라인의 `SimpleImputer`가 채운다. `fraction_imputed` = 1 − `fraction_present`.

| dataset | value | log_base_conversion | ortholog | relative_expression | n_control_samples | n_symbols_input | n_symbols_finite | n_mouse_entrez |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | log2 microarray | none | scmodage.preprocess.map_to_mouse_entrez species=human | mapped value minus mean of control samples in this dataset | 221 | 19541 | 19541 | 14120 |
| GSE56931 | log2 microarray | none | scmodage.preprocess.map_to_mouse_entrez species=human | mapped value minus mean of control samples in this dataset | 130 | 20755 | 20755 | 15147 |
| GSE113883 | log2(TPM+1) | none | scmodage.preprocess.map_to_mouse_entrez species=human | mapped value minus mean of control samples in this dataset | 153 | 35581 | 35581 | 16063 |

| dataset | clock | outcome | module | n_features | n_features_present | fraction_present | fraction_imputed | imputer |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE113883 | TM_black_chronological | chronological | black | 42 | 38 | 0.9047619047619048 | 0.09523809523809523 | pipeline SimpleImputer |
| GSE113883 | TM_blue_chronological | chronological | blue | 357 | 318 | 0.8907563025210085 | 0.10924369747899154 | pipeline SimpleImputer |
| GSE113883 | TM_brown_chronological | chronological | brown | 255 | 245 | 0.9607843137254902 | 0.039215686274509776 | pipeline SimpleImputer |
| GSE113883 | TA_chrono | chronological | composite | 16517 | 14715 | 0.8909002845553067 | 0.10909971544469332 | pipeline SimpleImputer |
| GSE113883 | TM_cyan_chronological | chronological | cyan | 20 | 20 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE113883 | TM_darkgreen_chronological | chronological | darkgreen | 36 | 35 | 0.9722222222222222 | 0.02777777777777779 | pipeline SimpleImputer |
| GSE113883 | TM_darkgrey_chronological | chronological | darkgrey | 31 | 31 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE113883 | TM_darkorange_chronological | chronological | darkorange | 17 | 16 | 0.9411764705882353 | 0.05882352941176472 | pipeline SimpleImputer |
| GSE113883 | TM_darkred_chronological | chronological | darkred | 26 | 24 | 0.9230769230769231 | 0.07692307692307687 | pipeline SimpleImputer |
| GSE113883 | TM_green_chronological | chronological | green | 97 | 96 | 0.9896907216494846 | 0.010309278350515427 | pipeline SimpleImputer |
| GSE113883 | TM_greenyellow_chronological | chronological | greenyellow | 59 | 59 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE113883 | TM_grey60_chronological | chronological | grey60 | 52 | 48 | 0.9230769230769231 | 0.07692307692307687 | pipeline SimpleImputer |
| GSE113883 | TM_lightgreen_chronological | chronological | lightgreen | 59 | 40 | 0.6779661016949152 | 0.3220338983050848 | pipeline SimpleImputer |
| GSE113883 | TM_lightyellow_chronological | chronological | lightyellow | 16 | 13 | 0.8125 | 0.1875 | pipeline SimpleImputer |
| GSE113883 | TM_magenta_chronological | chronological | magenta | 54 | 50 | 0.9259259259259259 | 0.07407407407407407 | pipeline SimpleImputer |
| GSE113883 | TM_midnightblue_chronological | chronological | midnightblue | 52 | 42 | 0.8076923076923077 | 0.1923076923076923 | pipeline SimpleImputer |
| GSE113883 | TM_orange_chronological | chronological | orange | 30 | 28 | 0.9333333333333333 | 0.06666666666666665 | pipeline SimpleImputer |
| GSE113883 | TM_pink_chronological | chronological | pink | 93 | 91 | 0.978494623655914 | 0.021505376344086002 | pipeline SimpleImputer |
| GSE113883 | TM_red_chronological | chronological | red | 117 | 101 | 0.8632478632478633 | 0.1367521367521367 | pipeline SimpleImputer |
| GSE113883 | TM_salmon_chronological | chronological | salmon | 12 | 10 | 0.8333333333333334 | 0.16666666666666663 | pipeline SimpleImputer |
| GSE113883 | TM_skyblue_chronological | chronological | skyblue | 20 | 18 | 0.9 | 0.09999999999999998 | pipeline SimpleImputer |
| GSE113883 | TM_tan_chronological | chronological | tan | 61 | 58 | 0.9508196721311475 | 0.049180327868852514 | pipeline SimpleImputer |
| GSE113883 | TM_turquoise_chronological | chronological | turquoise | 332 | 292 | 0.8795180722891566 | 0.12048192771084343 | pipeline SimpleImputer |
| GSE113883 | TM_yellow_chronological | chronological | yellow | 169 | 157 | 0.9289940828402367 | 0.07100591715976334 | pipeline SimpleImputer |
| GSE113883 | TM_black_mortality | mortality | black | 42 | 38 | 0.9047619047619048 | 0.09523809523809523 | pipeline SimpleImputer |
| GSE113883 | TM_blue_mortality | mortality | blue | 357 | 318 | 0.8907563025210085 | 0.10924369747899154 | pipeline SimpleImputer |
| GSE113883 | TM_brown_mortality | mortality | brown | 255 | 245 | 0.9607843137254902 | 0.039215686274509776 | pipeline SimpleImputer |
| GSE113883 | TA_mort | mortality | composite | 16760 | 14822 | 0.8843675417661098 | 0.11563245823389023 | pipeline SimpleImputer |
| GSE113883 | TM_cyan_mortality | mortality | cyan | 20 | 20 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE113883 | TM_darkgreen_mortality | mortality | darkgreen | 36 | 35 | 0.9722222222222222 | 0.02777777777777779 | pipeline SimpleImputer |
| GSE113883 | TM_darkgrey_mortality | mortality | darkgrey | 31 | 31 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE113883 | TM_darkorange_mortality | mortality | darkorange | 17 | 16 | 0.9411764705882353 | 0.05882352941176472 | pipeline SimpleImputer |
| GSE113883 | TM_darkred_mortality | mortality | darkred | 26 | 24 | 0.9230769230769231 | 0.07692307692307687 | pipeline SimpleImputer |
| GSE113883 | TM_green_mortality | mortality | green | 97 | 96 | 0.9896907216494846 | 0.010309278350515427 | pipeline SimpleImputer |
| GSE113883 | TM_greenyellow_mortality | mortality | greenyellow | 59 | 59 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE113883 | TM_grey60_mortality | mortality | grey60 | 52 | 48 | 0.9230769230769231 | 0.07692307692307687 | pipeline SimpleImputer |
| GSE113883 | TM_lightgreen_mortality | mortality | lightgreen | 59 | 40 | 0.6779661016949152 | 0.3220338983050848 | pipeline SimpleImputer |
| GSE113883 | TM_lightyellow_mortality | mortality | lightyellow | 16 | 13 | 0.8125 | 0.1875 | pipeline SimpleImputer |
| GSE113883 | TM_magenta_mortality | mortality | magenta | 54 | 50 | 0.9259259259259259 | 0.07407407407407407 | pipeline SimpleImputer |
| GSE113883 | TM_midnightblue_mortality | mortality | midnightblue | 52 | 42 | 0.8076923076923077 | 0.1923076923076923 | pipeline SimpleImputer |
| GSE113883 | TM_orange_mortality | mortality | orange | 30 | 28 | 0.9333333333333333 | 0.06666666666666665 | pipeline SimpleImputer |
| GSE113883 | TM_pink_mortality | mortality | pink | 93 | 91 | 0.978494623655914 | 0.021505376344086002 | pipeline SimpleImputer |
| GSE113883 | TM_red_mortality | mortality | red | 117 | 101 | 0.8632478632478633 | 0.1367521367521367 | pipeline SimpleImputer |
| GSE113883 | TM_salmon_mortality | mortality | salmon | 12 | 10 | 0.8333333333333334 | 0.16666666666666663 | pipeline SimpleImputer |
| GSE113883 | TM_skyblue_mortality | mortality | skyblue | 20 | 18 | 0.9 | 0.09999999999999998 | pipeline SimpleImputer |
| GSE113883 | TM_tan_mortality | mortality | tan | 61 | 58 | 0.9508196721311475 | 0.049180327868852514 | pipeline SimpleImputer |
| GSE113883 | TM_turquoise_mortality | mortality | turquoise | 332 | 292 | 0.8795180722891566 | 0.12048192771084343 | pipeline SimpleImputer |
| GSE113883 | TM_yellow_mortality | mortality | yellow | 169 | 157 | 0.9289940828402367 | 0.07100591715976334 | pipeline SimpleImputer |
| GSE39445 | TM_black_chronological | chronological | black | 42 | 37 | 0.8809523809523809 | 0.11904761904761907 | pipeline SimpleImputer |
| GSE39445 | TM_blue_chronological | chronological | blue | 357 | 268 | 0.7507002801120448 | 0.24929971988795518 | pipeline SimpleImputer |
| GSE39445 | TM_brown_chronological | chronological | brown | 255 | 224 | 0.8784313725490196 | 0.1215686274509804 | pipeline SimpleImputer |
| GSE39445 | TA_chrono | chronological | composite | 16517 | 13280 | 0.8040201005025126 | 0.1959798994974874 | pipeline SimpleImputer |
| GSE39445 | TM_cyan_chronological | chronological | cyan | 20 | 20 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE39445 | TM_darkgreen_chronological | chronological | darkgreen | 36 | 33 | 0.9166666666666666 | 0.08333333333333337 | pipeline SimpleImputer |
| GSE39445 | TM_darkgrey_chronological | chronological | darkgrey | 31 | 31 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE39445 | TM_darkorange_chronological | chronological | darkorange | 17 | 15 | 0.8823529411764706 | 0.11764705882352944 | pipeline SimpleImputer |
| GSE39445 | TM_darkred_chronological | chronological | darkred | 26 | 22 | 0.8461538461538461 | 0.15384615384615385 | pipeline SimpleImputer |
| GSE39445 | TM_green_chronological | chronological | green | 97 | 93 | 0.9587628865979382 | 0.04123711340206182 | pipeline SimpleImputer |
| GSE39445 | TM_greenyellow_chronological | chronological | greenyellow | 59 | 58 | 0.9830508474576272 | 0.016949152542372836 | pipeline SimpleImputer |
| GSE39445 | TM_grey60_chronological | chronological | grey60 | 52 | 46 | 0.8846153846153846 | 0.11538461538461542 | pipeline SimpleImputer |
| GSE39445 | TM_lightgreen_chronological | chronological | lightgreen | 59 | 38 | 0.6440677966101694 | 0.35593220338983056 | pipeline SimpleImputer |
| GSE39445 | TM_lightyellow_chronological | chronological | lightyellow | 16 | 11 | 0.6875 | 0.3125 | pipeline SimpleImputer |
| GSE39445 | TM_magenta_chronological | chronological | magenta | 54 | 49 | 0.9074074074074074 | 0.09259259259259256 | pipeline SimpleImputer |
| GSE39445 | TM_midnightblue_chronological | chronological | midnightblue | 52 | 39 | 0.75 | 0.25 | pipeline SimpleImputer |
| GSE39445 | TM_orange_chronological | chronological | orange | 30 | 27 | 0.9 | 0.09999999999999998 | pipeline SimpleImputer |
| GSE39445 | TM_pink_chronological | chronological | pink | 93 | 86 | 0.9247311827956989 | 0.07526881720430112 | pipeline SimpleImputer |
| GSE39445 | TM_red_chronological | chronological | red | 117 | 96 | 0.8205128205128205 | 0.17948717948717952 | pipeline SimpleImputer |
| GSE39445 | TM_salmon_chronological | chronological | salmon | 12 | 10 | 0.8333333333333334 | 0.16666666666666663 | pipeline SimpleImputer |
| GSE39445 | TM_skyblue_chronological | chronological | skyblue | 20 | 18 | 0.9 | 0.09999999999999998 | pipeline SimpleImputer |
| GSE39445 | TM_tan_chronological | chronological | tan | 61 | 52 | 0.8524590163934426 | 0.14754098360655743 | pipeline SimpleImputer |
| GSE39445 | TM_turquoise_chronological | chronological | turquoise | 332 | 276 | 0.8313253012048193 | 0.1686746987951807 | pipeline SimpleImputer |
| GSE39445 | TM_yellow_chronological | chronological | yellow | 169 | 142 | 0.8402366863905325 | 0.1597633136094675 | pipeline SimpleImputer |
| GSE39445 | TM_black_mortality | mortality | black | 42 | 37 | 0.8809523809523809 | 0.11904761904761907 | pipeline SimpleImputer |
| GSE39445 | TM_blue_mortality | mortality | blue | 357 | 268 | 0.7507002801120448 | 0.24929971988795518 | pipeline SimpleImputer |
| GSE39445 | TM_brown_mortality | mortality | brown | 255 | 224 | 0.8784313725490196 | 0.1215686274509804 | pipeline SimpleImputer |
| GSE39445 | TA_mort | mortality | composite | 16760 | 13368 | 0.7976133651551313 | 0.20238663484486874 | pipeline SimpleImputer |
| GSE39445 | TM_cyan_mortality | mortality | cyan | 20 | 20 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE39445 | TM_darkgreen_mortality | mortality | darkgreen | 36 | 33 | 0.9166666666666666 | 0.08333333333333337 | pipeline SimpleImputer |
| GSE39445 | TM_darkgrey_mortality | mortality | darkgrey | 31 | 31 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE39445 | TM_darkorange_mortality | mortality | darkorange | 17 | 15 | 0.8823529411764706 | 0.11764705882352944 | pipeline SimpleImputer |
| GSE39445 | TM_darkred_mortality | mortality | darkred | 26 | 22 | 0.8461538461538461 | 0.15384615384615385 | pipeline SimpleImputer |
| GSE39445 | TM_green_mortality | mortality | green | 97 | 93 | 0.9587628865979382 | 0.04123711340206182 | pipeline SimpleImputer |
| GSE39445 | TM_greenyellow_mortality | mortality | greenyellow | 59 | 58 | 0.9830508474576272 | 0.016949152542372836 | pipeline SimpleImputer |
| GSE39445 | TM_grey60_mortality | mortality | grey60 | 52 | 46 | 0.8846153846153846 | 0.11538461538461542 | pipeline SimpleImputer |
| GSE39445 | TM_lightgreen_mortality | mortality | lightgreen | 59 | 38 | 0.6440677966101694 | 0.35593220338983056 | pipeline SimpleImputer |
| GSE39445 | TM_lightyellow_mortality | mortality | lightyellow | 16 | 11 | 0.6875 | 0.3125 | pipeline SimpleImputer |
| GSE39445 | TM_magenta_mortality | mortality | magenta | 54 | 49 | 0.9074074074074074 | 0.09259259259259256 | pipeline SimpleImputer |
| GSE39445 | TM_midnightblue_mortality | mortality | midnightblue | 52 | 39 | 0.75 | 0.25 | pipeline SimpleImputer |
| GSE39445 | TM_orange_mortality | mortality | orange | 30 | 27 | 0.9 | 0.09999999999999998 | pipeline SimpleImputer |
| GSE39445 | TM_pink_mortality | mortality | pink | 93 | 86 | 0.9247311827956989 | 0.07526881720430112 | pipeline SimpleImputer |
| GSE39445 | TM_red_mortality | mortality | red | 117 | 96 | 0.8205128205128205 | 0.17948717948717952 | pipeline SimpleImputer |
| GSE39445 | TM_salmon_mortality | mortality | salmon | 12 | 10 | 0.8333333333333334 | 0.16666666666666663 | pipeline SimpleImputer |
| GSE39445 | TM_skyblue_mortality | mortality | skyblue | 20 | 18 | 0.9 | 0.09999999999999998 | pipeline SimpleImputer |
| GSE39445 | TM_tan_mortality | mortality | tan | 61 | 52 | 0.8524590163934426 | 0.14754098360655743 | pipeline SimpleImputer |
| GSE39445 | TM_turquoise_mortality | mortality | turquoise | 332 | 276 | 0.8313253012048193 | 0.1686746987951807 | pipeline SimpleImputer |
| GSE39445 | TM_yellow_mortality | mortality | yellow | 169 | 142 | 0.8402366863905325 | 0.1597633136094675 | pipeline SimpleImputer |
| GSE56931 | TM_black_chronological | chronological | black | 42 | 38 | 0.9047619047619048 | 0.09523809523809523 | pipeline SimpleImputer |
| GSE56931 | TM_blue_chronological | chronological | blue | 357 | 275 | 0.7703081232492998 | 0.22969187675070024 | pipeline SimpleImputer |
| GSE56931 | TM_brown_chronological | chronological | brown | 255 | 239 | 0.9372549019607843 | 0.06274509803921569 | pipeline SimpleImputer |
| GSE56931 | TA_chrono | chronological | composite | 16517 | 13971 | 0.8458557849488406 | 0.15414421505115938 | pipeline SimpleImputer |
| GSE56931 | TM_cyan_chronological | chronological | cyan | 20 | 20 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE56931 | TM_darkgreen_chronological | chronological | darkgreen | 36 | 35 | 0.9722222222222222 | 0.02777777777777779 | pipeline SimpleImputer |
| GSE56931 | TM_darkgrey_chronological | chronological | darkgrey | 31 | 28 | 0.9032258064516129 | 0.09677419354838712 | pipeline SimpleImputer |
| GSE56931 | TM_darkorange_chronological | chronological | darkorange | 17 | 15 | 0.8823529411764706 | 0.11764705882352944 | pipeline SimpleImputer |
| GSE56931 | TM_darkred_chronological | chronological | darkred | 26 | 23 | 0.8846153846153846 | 0.11538461538461542 | pipeline SimpleImputer |
| GSE56931 | TM_green_chronological | chronological | green | 97 | 94 | 0.9690721649484536 | 0.030927835051546393 | pipeline SimpleImputer |
| GSE56931 | TM_greenyellow_chronological | chronological | greenyellow | 59 | 57 | 0.9661016949152542 | 0.03389830508474578 | pipeline SimpleImputer |
| GSE56931 | TM_grey60_chronological | chronological | grey60 | 52 | 45 | 0.8653846153846154 | 0.13461538461538458 | pipeline SimpleImputer |
| GSE56931 | TM_lightgreen_chronological | chronological | lightgreen | 59 | 40 | 0.6779661016949152 | 0.3220338983050848 | pipeline SimpleImputer |
| GSE56931 | TM_lightyellow_chronological | chronological | lightyellow | 16 | 12 | 0.75 | 0.25 | pipeline SimpleImputer |
| GSE56931 | TM_magenta_chronological | chronological | magenta | 54 | 49 | 0.9074074074074074 | 0.09259259259259256 | pipeline SimpleImputer |
| GSE56931 | TM_midnightblue_chronological | chronological | midnightblue | 52 | 42 | 0.8076923076923077 | 0.1923076923076923 | pipeline SimpleImputer |
| GSE56931 | TM_orange_chronological | chronological | orange | 30 | 28 | 0.9333333333333333 | 0.06666666666666665 | pipeline SimpleImputer |
| GSE56931 | TM_pink_chronological | chronological | pink | 93 | 88 | 0.946236559139785 | 0.053763440860215006 | pipeline SimpleImputer |
| GSE56931 | TM_red_chronological | chronological | red | 117 | 90 | 0.7692307692307693 | 0.23076923076923073 | pipeline SimpleImputer |
| GSE56931 | TM_salmon_chronological | chronological | salmon | 12 | 9 | 0.75 | 0.25 | pipeline SimpleImputer |
| GSE56931 | TM_skyblue_chronological | chronological | skyblue | 20 | 18 | 0.9 | 0.09999999999999998 | pipeline SimpleImputer |
| GSE56931 | TM_tan_chronological | chronological | tan | 61 | 51 | 0.8360655737704918 | 0.16393442622950816 | pipeline SimpleImputer |
| GSE56931 | TM_turquoise_chronological | chronological | turquoise | 332 | 275 | 0.8283132530120482 | 0.17168674698795183 | pipeline SimpleImputer |
| GSE56931 | TM_yellow_chronological | chronological | yellow | 169 | 152 | 0.8994082840236687 | 0.10059171597633132 | pipeline SimpleImputer |
| GSE56931 | TM_black_mortality | mortality | black | 42 | 38 | 0.9047619047619048 | 0.09523809523809523 | pipeline SimpleImputer |
| GSE56931 | TM_blue_mortality | mortality | blue | 357 | 275 | 0.7703081232492998 | 0.22969187675070024 | pipeline SimpleImputer |
| GSE56931 | TM_brown_mortality | mortality | brown | 255 | 239 | 0.9372549019607843 | 0.06274509803921569 | pipeline SimpleImputer |
| GSE56931 | TA_mort | mortality | composite | 16760 | 14061 | 0.8389618138424821 | 0.16103818615751786 | pipeline SimpleImputer |
| GSE56931 | TM_cyan_mortality | mortality | cyan | 20 | 20 | 1.0 | 0.0 | pipeline SimpleImputer |
| GSE56931 | TM_darkgreen_mortality | mortality | darkgreen | 36 | 35 | 0.9722222222222222 | 0.02777777777777779 | pipeline SimpleImputer |
| GSE56931 | TM_darkgrey_mortality | mortality | darkgrey | 31 | 28 | 0.9032258064516129 | 0.09677419354838712 | pipeline SimpleImputer |
| GSE56931 | TM_darkorange_mortality | mortality | darkorange | 17 | 15 | 0.8823529411764706 | 0.11764705882352944 | pipeline SimpleImputer |
| GSE56931 | TM_darkred_mortality | mortality | darkred | 26 | 23 | 0.8846153846153846 | 0.11538461538461542 | pipeline SimpleImputer |
| GSE56931 | TM_green_mortality | mortality | green | 97 | 94 | 0.9690721649484536 | 0.030927835051546393 | pipeline SimpleImputer |
| GSE56931 | TM_greenyellow_mortality | mortality | greenyellow | 59 | 57 | 0.9661016949152542 | 0.03389830508474578 | pipeline SimpleImputer |
| GSE56931 | TM_grey60_mortality | mortality | grey60 | 52 | 45 | 0.8653846153846154 | 0.13461538461538458 | pipeline SimpleImputer |
| GSE56931 | TM_lightgreen_mortality | mortality | lightgreen | 59 | 40 | 0.6779661016949152 | 0.3220338983050848 | pipeline SimpleImputer |
| GSE56931 | TM_lightyellow_mortality | mortality | lightyellow | 16 | 12 | 0.75 | 0.25 | pipeline SimpleImputer |
| GSE56931 | TM_magenta_mortality | mortality | magenta | 54 | 49 | 0.9074074074074074 | 0.09259259259259256 | pipeline SimpleImputer |
| GSE56931 | TM_midnightblue_mortality | mortality | midnightblue | 52 | 42 | 0.8076923076923077 | 0.1923076923076923 | pipeline SimpleImputer |
| GSE56931 | TM_orange_mortality | mortality | orange | 30 | 28 | 0.9333333333333333 | 0.06666666666666665 | pipeline SimpleImputer |
| GSE56931 | TM_pink_mortality | mortality | pink | 93 | 88 | 0.946236559139785 | 0.053763440860215006 | pipeline SimpleImputer |
| GSE56931 | TM_red_mortality | mortality | red | 117 | 90 | 0.7692307692307693 | 0.23076923076923073 | pipeline SimpleImputer |
| GSE56931 | TM_salmon_mortality | mortality | salmon | 12 | 9 | 0.75 | 0.25 | pipeline SimpleImputer |
| GSE56931 | TM_skyblue_mortality | mortality | skyblue | 20 | 18 | 0.9 | 0.09999999999999998 | pipeline SimpleImputer |
| GSE56931 | TM_tan_mortality | mortality | tan | 61 | 51 | 0.8360655737704918 | 0.16393442622950816 | pipeline SimpleImputer |
| GSE56931 | TM_turquoise_mortality | mortality | turquoise | 332 | 275 | 0.8283132530120482 | 0.17168674698795183 | pipeline SimpleImputer |
| GSE56931 | TM_yellow_mortality | mortality | yellow | 169 | 152 | 0.8994082840236687 | 0.10059171597633132 | pipeline SimpleImputer |

### 7.2 복합 시계 H1 표

부트스트랩 추출은 마이크로어레이가 `results/stage2/boot_draws_*.tsv`, GSE113883이 `results/stage4/rnaseq/boot_draws_GSE113883.tsv`다. 시드는 그 파일을 만들 때 20261015였다. 모듈 점수는 부트스트랩하지 않았다. 출력 단위는 §5.5와 같다.

| dataset | clock | output_unit | n_subjects | n_samples | 2A | R2_within | 2A/SD_between | 2A/SD_resid | peak_h (center, half) | LRT p | ICC | n_converged | n_not_converged |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | TA_chrono | fraction of maximum lifespan | 24 | 221 | 0.0752260169510138 [0.0333439302176083, 0.128101448514176] | 0.0818707783732363 [0.0181532226910966, 0.20123919586041] | 0.747395182534912 [0.378876871337183, 1.27347282361924] | 0.750240076435841 [0.353119847965313, 1.18271146036711] | 15.3360608213662 (15.2425147637526, 2.1292787060599) | 0.000235567727394529 | 0.463850671172076 | 2000 | 0 |
| GSE39445 | TA_mort | log10 hazard ratio | 24 | 221 | 0.470276403374039 [0.215647266233542, 0.771302409358931] | 0.0719801227245767 [0.0151416898877572, 0.176914923892424] | 0.720473745860167 [0.37214510427527, 1.17672572611328] | 0.705167502585696 [0.335383835193897, 1.10155636042881] | 15.5538826108025 (15.4706502081726, 1.95030763365726) | 0.000615890141275604 | 0.450506238667552 | 2000 | 0 |
| GSE56931 | TA_chrono | fraction of maximum lifespan | 14 | 130 | 0.0294906994760309 [0.0109107651095166, 0.0522992803151469] | 0.0589265110850511 [0.00816751804898923, 0.208099354293255] | 0.490736217939359 [0.183399595436075, 1.20225983716478] | 0.711643633227566 [0.253300640812028, 1.33849429320705] | 0.827909618359159 (0.605119785999099, 2.83506082730293) | 0.029920138819241 | 0.649531474992563 | 2000 | 0 |
| GSE56931 | TA_mort | log10 hazard ratio | 14 | 130 | 0.169978316903373 [0.0360959155157471, 0.359225336328813] | 0.032665947996622 [0.0014817349214344, 0.159495880841973] | 0.737461667077993 [0.184011096418156, 1.57682280114873] | 0.52558701889251 [0.108289006124987, 1.17230288643993] | 1.28654792915375 (1.10363376536451, 5.49568961854717) | 0.157968351876906 | 0.266356971585751 | 2000 | 0 |
| GSE113883 | TA_chrono | fraction of maximum lifespan | 11 | 153 | 0.0699539619458422 [0.0369519119878254, 0.101275123307544] | 0.0966791813140306 [0.0247800777387931, 0.235397109571688] | 0.7304551667482 [0.3827976679813, 1.3912026164744] | 0.904354140589051 [0.452195636186797, 1.39569854899423] | 15.667602590042 (15.6128080538023, 1.82336493839622) | 0.000751939790616891 | 0.569256878192063 | 2000 | 0 |
| GSE113883 | TA_mort | log10 hazard ratio | 11 | 153 | 0.0466007880250426 [0.0164682998603527, 0.180279015492726] | 0.00416979906927029 [0.000555824004425451, 0.0650123736562311] | 0.147617554907699 [0.0593244426425826, 0.697869644748395] | 0.179923991856036 [0.0663404129344585, 0.695219796595057] | 17.6460030950774 (17.9051883332584, 10.0827106730224) | 0.745144998027147 | 0.561185272755364 | 2000 | 0 |

### 7.3 모듈 시계 H1 표 (자료 × 목표 × 모듈: R2_within, LRT p, BH q, 정점, 2A)

BH q는 자료 × 목표 안의 모듈 23개 LRT p에 붙였다. 복합 시계에는 붙이지 않았다.

| dataset | outcome | module | n_subjects | n_samples | r2_within | lrt_p | q_bh | peak_h | two_a | two_a_over_sd_between |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE113883 | chronological | black | 11 | 153 | 0.0920649868165196 | 0.00104272227494776 | 0.00239826123237984 | 9.43426408390714 | 0.0668828044975735 | 0.552889459518017 |
| GSE113883 | chronological | blue | 11 | 153 | 0.0949335321810294 | 0.000889833155799814 | 0.00227401806482175 | 13.343483035576 | 0.198343438431693 | 0.864351511288424 |
| GSE113883 | chronological | brown | 11 | 153 | 0.00824164069931765 | 0.54229698366573 | 0.54229698366573 | 20.9352127938297 | 0.0387549991216469 | 0.408376940735253 |
| GSE113883 | chronological | cyan | 11 | 153 | 0.123362172471242 | 6.94656665720149e-05 | 0.000266285055192724 | 15.4118534698436 | 0.194305528738217 | 1.54259446993807 |
| GSE113883 | chronological | darkgreen | 11 | 153 | 0.0443816487496821 | 0.0394479954252844 | 0.0604869263187694 | 13.8172747297159 | 0.0574962711118861 | 0.368138816840889 |
| GSE113883 | chronological | darkgrey | 11 | 153 | 0.193471250411986 | 2.6468084551431e-07 | 1.52191486170728e-06 | 13.7807786019332 | 0.306212507749194 | 2.6800133791416 |
| GSE113883 | chronological | darkorange | 11 | 153 | 0.202292621573607 | 7.38058900210152e-08 | 5.65845156827783e-07 | 13.0407547772519 | 0.195738022457963 | 1.88060078653627 |
| GSE113883 | chronological | darkred | 11 | 153 | 0.0489837158397306 | 0.0281393678234745 | 0.0462289614242796 | 11.6969906951804 | 0.0791034816825276 | 0.455924233907432 |
| GSE113883 | chronological | green | 11 | 153 | 0.0283018515741664 | 0.140448938560477 | 0.161516279344548 | 12.7952611234869 | 0.0761758776937856 | 0.959380030734378 |
| GSE113883 | chronological | greenyellow | 11 | 153 | 0.0848094442133197 | 0.00184236749904285 | 0.00385222295254414 | 10.8979759962828 | 0.0917600732355018 | 0.560883354198333 |
| GSE113883 | chronological | grey60 | 11 | 153 | 0.114840502179783 | 0.000154386719262304 | 0.000507270649004714 | 13.6397119948119 | 0.191702775537332 | 1.38442202820411 |
| GSE113883 | chronological | lightgreen | 11 | 153 | 0.0217209816336884 | 0.211500951209 | 0.231643898943191 | 0.610438851335317 | 0.0821897005232743 | 0.524482995370835 |
| GSE113883 | chronological | lightyellow | 11 | 153 | 0.0360141107383472 | 0.0613543948049248 | 0.0855998284319558 | 13.4227516081739 | 0.0954029371223188 | 1.07838454316542 |
| GSE113883 | chronological | magenta | 11 | 153 | 0.126198905064618 | 6.18663055924296e-05 | 0.000266285055192724 | 3.24652801917783 | 0.0779897814266044 | 0.825408889106047 |
| GSE113883 | chronological | midnightblue | 11 | 153 | 0.0146870597153587 | 0.359631823335839 | 0.375978724396559 | 21.4272941407535 | 0.055560735068777 | 0.545974405473107 |
| GSE113883 | chronological | orange | 11 | 153 | 0.237566290178325 | 4.02441765520324e-09 | 9.25616060696745e-08 | 14.1583138939012 | 0.206804349325327 | 1.57559164550323 |
| GSE113883 | chronological | pink | 11 | 153 | 0.212319385973019 | 4.80387011478345e-08 | 5.52445063200097e-07 | 11.4550545518287 | 0.0878268285016077 | 0.645536077198604 |
| GSE113883 | chronological | red | 11 | 153 | 0.0763074511089803 | 0.00377567349522923 | 0.00723670753252269 | 1.70607891543789 | 0.136308924477375 | 0.663785796195427 |
| GSE113883 | chronological | salmon | 11 | 153 | 0.0377226599091116 | 0.0655159201986105 | 0.0855998284319558 | 17.9390325766023 | 0.0412458238189219 | 0.352801060980765 |
| GSE113883 | chronological | skyblue | 11 | 153 | 0.0282175443143632 | 0.122490396691927 | 0.14827784862707 | 11.5618933160481 | 0.0790365662960361 | 0.612951118618877 |
| GSE113883 | chronological | tan | 11 | 153 | 0.0372348524161215 | 0.0669911700771828 | 0.0855998284319558 | 12.752035911046 | 0.0855131638195602 | 0.515009921340629 |
| GSE113883 | chronological | turquoise | 11 | 153 | 0.0642818018652043 | 0.00909627741769173 | 0.0160934138928392 | 1.13969806986858 | 0.0604547155039366 | 0.41988553021994 |
| GSE113883 | chronological | yellow | 11 | 153 | 0.110856561888506 | 0.000240544783693973 | 0.000691566253120172 | 15.5766940777842 | 0.237181521625488 | 0.846807376873637 |
| GSE113883 | mortality | black | 11 | 153 | 0.180051879317388 | 7.51473389146692e-07 | 8.64194397518696e-06 | 10.9097923409538 | 0.477457555396201 | 0.693186186622412 |
| GSE113883 | mortality | blue | 11 | 153 | 0.0204800727249376 | 0.234320942124992 | 0.296158682354586 | 14.5296601202447 | 0.472926075533795 | 0.436327911455731 |
| GSE113883 | mortality | brown | 11 | 153 | 0.0601982394593815 | 0.0123276056507036 | 0.0330162364431858 | 13.752169203839 | 0.662702614302119 | 0.88447692064344 |
| GSE113883 | mortality | cyan | 11 | 153 | 0.12634953899778 | 7.2150772242239e-05 | 0.000331893552314299 | 15.9715619563109 | 0.631037722638229 | 1.55972847241303 |
| GSE113883 | mortality | darkgreen | 11 | 153 | 0.0113963005443716 | 0.448160024993602 | 0.490536100194936 | 9.01603504720726 | 0.122777369606789 | 0.173471314765207 |
| GSE113883 | mortality | darkgrey | 11 | 153 | 0.165501059230859 | 2.81529071490669e-06 | 2.15838954809513e-05 | 14.0262131726109 | 0.987738376435797 | 1.89034598325829 |
| GSE113883 | mortality | darkorange | 11 | 153 | 0.0974056889863783 | 0.000606000741228307 | 0.00232300284137518 | 12.5957894777605 | 0.383612786429356 | 1.0235247072312 |
| GSE113883 | mortality | darkred | 11 | 153 | 0.0593356600743192 | 0.0129193968690727 | 0.0330162364431858 | 11.7418291295997 | 0.490598058665928 | 0.54641074264881 |
| GSE113883 | mortality | green | 11 | 153 | 0.00727744333767699 | 0.614711725512297 | 0.614711725512297 | 10.1493538043249 | 0.114872720305081 | 0.27724380979572 |
| GSE113883 | mortality | greenyellow | 11 | 153 | 0.206825201478375 | 7.00536210296014e-08 | 1.61123328368083e-06 | 11.5406136291133 | 0.68757026208995 | 0.708410599635171 |
| GSE113883 | mortality | grey60 | 11 | 153 | 0.0186371248605184 | 0.274871285203076 | 0.316101977983537 | 1.253999661326 | 0.357887993137256 | 0.426347991670353 |
| GSE113883 | mortality | lightgreen | 11 | 153 | 0.0254412903465814 | 0.161792310366885 | 0.248081542562558 | 1.00127115461664 | 0.406057506559653 | 0.539685383744935 |
| GSE113883 | mortality | lightyellow | 11 | 153 | 0.0198339626171844 | 0.206866312809622 | 0.296158682354586 | 13.2076443735556 | 0.318950183598149 | 0.832362539580184 |
| GSE113883 | mortality | magenta | 11 | 153 | 0.0664444723016787 | 0.00675632441920848 | 0.0221993516631136 | 4.0424626850677 | 0.271854393994277 | 1.00221541578618 |
| GSE113883 | mortality | midnightblue | 11 | 153 | 0.0562925968431228 | 0.0163194315335967 | 0.0375346925272725 | 19.2183714558209 | 0.591693206272589 | 0.952308140474875 |
| GSE113883 | mortality | orange | 11 | 153 | 0.135688293041572 | 2.8459735767537e-05 | 0.000163643480663338 | 13.9279380615639 | 1.00072322314913 | 1.16893202365014 |
| GSE113883 | mortality | pink | 11 | 153 | 0.0399901537273766 | 0.0553988457722534 | 0.0910123894829878 | 7.39782198995527 | 0.296386857978733 | 0.379022947758606 |
| GSE113883 | mortality | red | 11 | 153 | 0.0107498435749711 | 0.469208443664722 | 0.490536100194936 | 2.85945718876094 | 0.302206507490001 | 0.0766163323795952 |
| GSE113883 | mortality | salmon | 11 | 153 | 0.0421977270777493 | 0.0471365433288297 | 0.0833954228125448 | 14.7114513603138 | 0.136072814106418 | 0.413533530038258 |
| GSE113883 | mortality | skyblue | 11 | 153 | 0.0204953298620161 | 0.230454562251633 | 0.296158682354586 | 5.24138228441545 | 0.162150102629153 | 0.264606932881049 |
| GSE113883 | mortality | tan | 11 | 153 | 0.050643987897681 | 0.0246558077612224 | 0.0515530525916469 | 11.0836100014759 | 0.55569516672134 | 0.557349096357486 |
| GSE113883 | mortality | turquoise | 11 | 153 | 0.01958278310867 | 0.244652824553788 | 0.296158682354586 | 8.39017465316885 | 0.12378609671399 | 0.189238295810394 |
| GSE113883 | mortality | yellow | 11 | 153 | 0.048984193053558 | 0.0291533549700461 | 0.0558772636925884 | 17.4615677991321 | 0.546002214568777 | 0.597853829916604 |
| GSE39445 | chronological | black | 24 | 221 | 0.00419662297051537 | 0.649316597852252 | 0.649316597852252 | 15.5430089372762 | 0.0253390071268415 | 0.139592254825397 |
| GSE39445 | chronological | blue | 24 | 221 | 0.237829014095836 | 3.15097496339781e-12 | 3.62362120790748e-11 | 12.8858696782194 | 0.306392863404453 | 1.20691697333535 |
| GSE39445 | chronological | brown | 24 | 221 | 0.120012185415971 | 3.1110076635881e-06 | 1.43106352525053e-05 | 2.41251767514426 | 0.178499819493741 | 0.849634547263129 |
| GSE39445 | chronological | cyan | 24 | 221 | 0.0784202942552863 | 0.000310206618827339 | 0.000891844029128599 | 17.0881605720823 | 0.0626366791287462 | 0.400904915236077 |
| GSE39445 | chronological | darkgreen | 24 | 221 | 0.0438650714444653 | 0.0141317219871215 | 0.0270858004753162 | 12.2297235510656 | 0.111673059359522 | 0.467261940074429 |
| GSE39445 | chronological | darkgrey | 24 | 221 | 0.0228754471660846 | 0.0995629987206397 | 0.14312181066092 | 13.4425874433925 | 0.100489270007581 | 0.483218333123922 |
| GSE39445 | chronological | darkorange | 24 | 221 | 0.252125912766848 | 6.43156705950182e-13 | 1.47926042368542e-11 | 13.1994196787556 | 0.203246818413915 | 1.34238501727207 |
| GSE39445 | chronological | darkred | 24 | 221 | 0.170268168064267 | 1.11038075040807e-08 | 6.38468931484642e-08 | 14.6035753971536 | 0.0938181297988642 | 0.637768917975413 |
| GSE39445 | chronological | green | 24 | 221 | 0.0142119820256309 | 0.192140166425451 | 0.220961191389269 | 23.1118641649033 | 0.0585566987935926 | 0.391268920705086 |
| GSE39445 | chronological | greenyellow | 24 | 221 | 0.0916708400447678 | 8.55966574357814e-05 | 0.000328120520170495 | 14.0126566688726 | 0.1837978926623 | 0.899391305042763 |
| GSE39445 | chronological | grey60 | 24 | 221 | 0.0452800114168385 | 0.0118589267499855 | 0.0247959377499697 | 1.61122714677413 | 0.0618588514015 | 0.304573317361092 |
| GSE39445 | chronological | lightgreen | 24 | 221 | 0.0113722292224852 | 0.337849225785314 | 0.353206008775556 | 19.5440607248531 | 0.0400270907317322 | 0.218856818236956 |
| GSE39445 | chronological | lightyellow | 24 | 221 | 0.016958944730711 | 0.188011532026613 | 0.220961191389269 | 16.3703421898219 | 0.0460476121930791 | 0.369636326696857 |
| GSE39445 | chronological | magenta | 24 | 221 | 0.031841709089661 | 0.0418042393016934 | 0.0640998335959299 | 4.83609712470271 | 0.0778257211224785 | 0.302928297593307 |
| GSE39445 | chronological | midnightblue | 24 | 221 | 0.0360812244223455 | 0.0263980738213452 | 0.0442020942617806 | 2.48971135004274 | 0.0875421274577231 | 0.706051123694089 |
| GSE39445 | chronological | orange | 24 | 221 | 0.0623826475609638 | 0.00194298914813481 | 0.00496541671190006 | 15.5211451674636 | 0.131986794103685 | 0.585104851126349 |
| GSE39445 | chronological | pink | 24 | 221 | 0.202568746239115 | 1.81628517628059e-10 | 1.39248530181512e-09 | 13.8622823951712 | 0.120770235943479 | 0.898521256156496 |
| GSE39445 | chronological | red | 24 | 221 | 0.0166988451905946 | 0.210605494933635 | 0.23066316111779 | 4.09149417114038 | 0.0631606030125016 | 0.289550951495543 |
| GSE39445 | chronological | salmon | 24 | 221 | 0.0207582463074085 | 0.12989270482601 | 0.175737188882249 | 14.2862780888814 | 0.0293646692822192 | 0.23333152794679 |
| GSE39445 | chronological | skyblue | 24 | 221 | 0.0367524785410869 | 0.0269056225941273 | 0.0442020942617806 | 2.12416317227273 | 0.0685859947947901 | 0.417156359041061 |
| GSE39445 | chronological | tan | 24 | 221 | 0.0197098625397636 | 0.139182162275866 | 0.177843874019162 | 16.4114414075102 | 0.0940545722906332 | 0.543800108714267 |
| GSE39445 | chronological | turquoise | 24 | 221 | 0.0482271919172995 | 0.00919790868805224 | 0.0211551899825202 | 11.6135068062167 | 0.0545605326288795 | 0.334577890279516 |
| GSE39445 | chronological | yellow | 24 | 221 | 0.0788763500807268 | 0.00030816196430124 | 0.000891844029128599 | 3.98006739489746 | 0.14377779990153 | 0.469842037649312 |
| GSE39445 | mortality | black | 24 | 221 | 0.0334487137985466 | 0.0347355763656406 | 0.0726289324008849 | 15.2074353819786 | 0.401688641198103 | 0.395269479689631 |
| GSE39445 | mortality | blue | 24 | 221 | 0.173384151988194 | 9.86053798875224e-09 | 1.13396186870651e-07 | 13.3350666192169 | 1.29787930191934 | 1.19478262617152 |
| GSE39445 | mortality | brown | 24 | 221 | 0.0180956475419376 | 0.179471889493929 | 0.275190230557358 | 12.7854394600048 | 0.284587699346108 | 0.284470197259363 |
| GSE39445 | mortality | cyan | 24 | 221 | 0.0387155662466691 | 0.0194919313586336 | 0.0448314421248573 | 14.963377925844 | 0.162727983866685 | 0.342989233983329 |
| GSE39445 | mortality | darkgreen | 24 | 221 | 0.0253755609023505 | 0.085967523581588 | 0.16477108686471 | 12.2070204131449 | 0.393949067041888 | 0.327762232124016 |
| GSE39445 | mortality | darkgrey | 24 | 221 | 0.0157714699716435 | 0.215049197736021 | 0.290948914584029 | 10.5754611115216 | 0.360362347622348 | 0.314085859039108 |
| GSE39445 | mortality | darkorange | 24 | 221 | 0.222732128205742 | 2.5440199674316e-11 | 5.85124592509269e-10 | 11.990507591373 | 0.619138288605142 | 1.20192715550285 |
| GSE39445 | mortality | darkred | 24 | 221 | 0.126281955122939 | 1.93892619480763e-06 | 8.91906049611511e-06 | 13.6277744098779 | 0.434254400139623 | 0.768350592970405 |
| GSE39445 | mortality | green | 24 | 221 | 0.0169112079772625 | 0.156286236015497 | 0.256755959168316 | 21.0399406078005 | 0.26324527402271 | 0.499942983293326 |
| GSE39445 | mortality | greenyellow | 24 | 221 | 0.121801675571493 | 3.18247726036681e-06 | 1.21994961647395e-05 | 14.1088658354686 | 1.05625203285525 | 1.17314852356444 |
| GSE39445 | mortality | grey60 | 24 | 221 | 0.0884199606988083 | 0.000127459660491946 | 0.000366446523914344 | 3.47713463219173 | 0.673735221223613 | 0.750987270589488 |
| GSE39445 | mortality | lightgreen | 24 | 221 | 0.00200577581082695 | 0.834398256894445 | 0.872325450389647 | 19.6939855837568 | 0.067653715841672 | 0.0933460204018298 |
| GSE39445 | mortality | lightyellow | 24 | 221 | 0.00352594888358004 | 0.703834950188363 | 0.770866850206302 | 16.0614428390137 | 0.0909032999156293 | 0.148973434275692 |
| GSE39445 | mortality | magenta | 24 | 221 | 0.0169441690686331 | 0.191846089160145 | 0.275778753167708 | 4.65555674736115 | 0.248220615866462 | 0.246404192346565 |
| GSE39445 | mortality | midnightblue | 24 | 221 | 0.0143688103662002 | 0.240032836220979 | 0.30670862406014 | 2.00118991393135 | 0.248041212443898 | 0.295007594614641 |
| GSE39445 | mortality | orange | 24 | 221 | 0.139130759332057 | 5.11987945961834e-07 | 2.94393068928055e-06 | 14.9948982783575 | 0.828371516049201 | 1.02927050874338 |
| GSE39445 | mortality | pink | 24 | 221 | 0.0863849327369744 | 0.00011647346563911 | 0.000366446523914344 | 14.6593137052705 | 0.337221237625616 | 0.587985659898256 |
| GSE39445 | mortality | red | 24 | 221 | 0.0214812002813103 | 0.118609431495671 | 0.20984745572311 | 3.15560654057404 | 0.438072352017299 | 0.116570470330096 |
| GSE39445 | mortality | salmon | 24 | 221 | 0.0535088229874735 | 0.00468818860828374 | 0.0119809264433918 | 14.3737687655268 | 0.133847496858067 | 0.367494202841448 |
| GSE39445 | mortality | skyblue | 24 | 221 | 0.00897552098664223 | 0.418185721734746 | 0.506224821047325 | 2.03710442275896 | 0.158796872671845 | 0.186820967388943 |
| GSE39445 | mortality | tan | 24 | 221 | 0.00106264679483215 | 0.912885926819834 | 0.912885926819834 | 1.85297786563406 | 0.105900290889614 | 0.172866594043816 |
| GSE39445 | mortality | turquoise | 24 | 221 | 0.144979647914123 | 2.59158198342643e-07 | 1.98687952062693e-06 | 13.079788251819 | 0.450458789316072 | 0.506352652960087 |
| GSE39445 | mortality | yellow | 24 | 221 | 0.00586286783638201 | 0.580147184298968 | 0.667169261943813 | 6.16149376718439 | 0.152796481707908 | 0.136368929997923 |
| GSE56931 | chronological | black | 14 | 130 | 0.00689307540415207 | 0.675530523241086 | 0.70623645611568 | 18.8430210798695 | 0.0176776988281424 | 0.159505747389247 |
| GSE56931 | chronological | blue | 14 | 130 | 0.0835709307842948 | 0.00551597823577023 | 0.0115334090384287 | 8.49374611650761 | 0.0930803773216562 | 0.513675643733922 |
| GSE56931 | chronological | brown | 14 | 130 | 0.103992356891917 | 0.00180315277649665 | 0.0082945027718846 | 23.3615266605128 | 0.0912619479334364 | 0.778558322290536 |
| GSE56931 | chronological | cyan | 14 | 130 | 0.0931102048235051 | 0.00367203908014258 | 0.0110523304077737 | 18.4887569220951 | 0.0485211119416131 | 0.444208237570394 |
| GSE56931 | chronological | darkgreen | 14 | 130 | 0.0467864370573649 | 0.0613364693244821 | 0.0940492529642058 | 5.47809913544504 | 0.0676059161557205 | 0.651786286915868 |
| GSE56931 | chronological | darkgrey | 14 | 130 | 0.000406752578562819 | 0.975456154614006 | 0.975456154614006 | 0.151978490634241 | 0.00714835219421273 | 0.0775549842817842 |
| GSE56931 | chronological | darkorange | 14 | 130 | 0.223188558409583 | 5.24850078796061e-07 | 1.20715518123094e-05 | 13.0257575359056 | 0.093590737918322 | 1.20157014592551 |
| GSE56931 | chronological | darkred | 14 | 130 | 0.0523315873682693 | 0.0517382344463975 | 0.0849985280190816 | 12.9319319521483 | 0.0279790359478605 | 0.756790859656383 |
| GSE56931 | chronological | green | 14 | 130 | 0.0225949280167521 | 0.253978923752743 | 0.307448170858583 | 8.61391421556577 | 0.0367077954645625 | 0.519563328394144 |
| GSE56931 | chronological | greenyellow | 14 | 130 | 0.0900155127962415 | 0.00454601902938541 | 0.0110523304077737 | 8.25986430634499 | 0.119851253410518 | 0.889462788607745 |
| GSE56931 | chronological | grey60 | 14 | 130 | 0.0320803366959598 | 0.176594826357309 | 0.238922412130477 | 14.7965669723146 | 0.0574069663463747 | 0.710885742491929 |
| GSE56931 | chronological | lightgreen | 14 | 130 | 0.0282868105573336 | 0.192796504568752 | 0.246351089171184 | 12.2750491530141 | 0.0311975009887524 | 0.290656222571811 |
| GSE56931 | chronological | lightyellow | 14 | 130 | 0.0919596079778938 | 0.00343578887956388 | 0.0110523304077737 | 5.21997300629944 | 0.0743556942209195 | 1.27516480698028 |
| GSE56931 | chronological | magenta | 14 | 130 | 0.0180737310695097 | 0.339025317633109 | 0.389879115278075 | 9.36134792638769 | 0.0317507616785345 | 0.460853426622173 |
| GSE56931 | chronological | midnightblue | 14 | 130 | 0.0712450857755046 | 0.0143823762637541 | 0.0275662211721953 | 2.31207348973244 | 0.108635191440427 | 1.20736821661093 |
| GSE56931 | chronological | orange | 14 | 130 | 0.206661212076165 | 1.59780340143946e-06 | 1.83747391165538e-05 | 15.2637127251035 | 0.15134038360775 | 2.3169798129481 |
| GSE56931 | chronological | pink | 14 | 130 | 0.126342340445824 | 0.00036335103065202 | 0.00208926842624911 | 11.1199417302676 | 0.040691216368205 | 0.459364143506234 |
| GSE56931 | chronological | red | 14 | 130 | 0.0868918422466177 | 0.00480536104685813 | 0.0110523304077737 | 22.3663160298919 | 0.0915921861699867 | 0.446978658782821 |
| GSE56931 | chronological | salmon | 14 | 130 | 0.0657790403304638 | 0.017234347310599 | 0.0304915375495213 | 19.4527771849152 | 0.029134566734471 | 0.320154726126958 |
| GSE56931 | chronological | skyblue | 14 | 130 | 0.166661129139752 | 2.50696964183316e-05 | 0.000192201005873875 | 4.23405108385446 | 0.118955304907104 | 1.47985879551545 |
| GSE56931 | chronological | tan | 14 | 130 | 0.00644165551802955 | 0.673349937764662 | 0.70623645611568 | 4.62056644992187 | 0.0264951903425614 | 0.227061695289138 |
| GSE56931 | chronological | turquoise | 14 | 130 | 0.091941341024543 | 0.00388621528189854 | 0.0110523304077737 | 3.51438834939511 | 0.0385955453481906 | 0.404394786271329 |
| GSE56931 | chronological | yellow | 14 | 130 | 0.0330910772045983 | 0.139519326354172 | 0.200559031634122 | 7.42135124811049 | 0.0460291589132994 | 0.410986324205495 |
| GSE56931 | mortality | black | 14 | 130 | 0.0172923543267151 | 0.364508273835368 | 0.381076831736976 | 15.0010467690942 | 0.157452465794652 | 0.266596244337249 |
| GSE56931 | mortality | blue | 14 | 130 | 0.0556424862973739 | 0.0333910809856194 | 0.063495321455112 | 6.38160062769946 | 0.310448500452738 | 0.35087304073573 |
| GSE56931 | mortality | brown | 14 | 130 | 0.0779639907500801 | 0.00936479530210405 | 0.0269237864935492 | 15.3192656248541 | 0.423567018181625 | 0.602972487850452 |
| GSE56931 | mortality | cyan | 14 | 130 | 0.0530021929012933 | 0.046891518458832 | 0.0770360660395097 | 17.4680181184238 | 0.193082800686707 | 0.534048410322706 |
| GSE56931 | mortality | darkgreen | 14 | 130 | 0.0453470709027217 | 0.0713540512810864 | 0.0965378340861757 | 4.53105837810403 | 0.280192892598758 | 0.417497038088834 |
| GSE56931 | mortality | darkgrey | 14 | 130 | 0.0107799719222351 | 0.522396032647298 | 0.522396032647298 | 3.67332251022731 | 0.198031870563805 | 0.51867367436756 |
| GSE56931 | mortality | darkorange | 14 | 130 | 0.0895951442676672 | 0.00477623855766983 | 0.0169232126065604 | 13.9271338124809 | 0.196176536281696 | 0.789178815472078 |
| GSE56931 | mortality | darkred | 14 | 130 | 0.108126184963575 | 0.00134898352061123 | 0.0103422069913528 | 11.7921623328549 | 0.222348810295337 | 0.885723452290527 |
| GSE56931 | mortality | green | 14 | 130 | 0.0183216635047599 | 0.32919155148734 | 0.360543127819467 | 7.43524324525559 | 0.146668004291273 | 0.425782432461275 |
| GSE56931 | mortality | greenyellow | 14 | 130 | 0.103721383650612 | 0.00215346916179761 | 0.0123824476803362 | 8.77025120655723 | 0.672416889831356 | 1.17341677815714 |
| GSE56931 | mortality | grey60 | 14 | 130 | 0.0576479338343247 | 0.0358886599528894 | 0.063495321455112 | 16.3703934608109 | 0.448728159479172 | 1.16538844766955 |
| GSE56931 | mortality | lightgreen | 14 | 130 | 0.046798175391677 | 0.0628541571264112 | 0.090352850869216 | 10.7480267161231 | 0.141404747357716 | 0.3180316894884 |
| GSE56931 | mortality | lightyellow | 14 | 130 | 0.119040291733443 | 0.000633691830898987 | 0.00728745605533835 | 5.18623631730539 | 0.36122281132082 | 1.39638556368888 |
| GSE56931 | mortality | magenta | 14 | 130 | 0.0462426563030093 | 0.0599008637003445 | 0.090352850869216 | 10.1520816986573 | 0.267613834635087 | 0.833135935539167 |
| GSE56931 | mortality | midnightblue | 14 | 130 | 0.0228218027858318 | 0.263429132501458 | 0.310861487998706 | 2.10206579880247 | 0.298046385016261 | 0.757845230141443 |
| GSE56931 | mortality | orange | 14 | 130 | 0.0356302204809017 | 0.123510868170419 | 0.157819442662202 | 15.0857874262409 | 0.295502970201542 | 0.729803548259059 |
| GSE56931 | mortality | pink | 14 | 130 | 0.17016687754802 | 1.82336978045962e-05 | 0.000419375049505713 | 11.6840808585176 | 0.258860108665723 | 0.587933080962635 |
| GSE56931 | mortality | red | 14 | 130 | 0.0695703562698264 | 0.0149398846306077 | 0.0381797051671085 | 0.162919136687949 | 0.356559726071841 | 0.124613950981799 |
| GSE56931 | mortality | salmon | 14 | 130 | 0.0213253750617666 | 0.270314337390179 | 0.310861487998706 | 19.0996621216644 | 0.0386024564685023 | 0.18026175501001 |
| GSE56931 | mortality | skyblue | 14 | 130 | 0.0866878748408164 | 0.0049047617602511 | 0.0169232126065604 | 4.64568014131419 | 0.346087229733568 | 0.939570100023015 |
| GSE56931 | mortality | tan | 14 | 130 | 0.0592320049884755 | 0.0298849040758074 | 0.0624866176130518 | 2.51201116672037 | 0.299610968731108 | 0.575341026861434 |
| GSE56931 | mortality | turquoise | 14 | 130 | 0.0864204056444807 | 0.00515054296721403 | 0.0169232126065604 | 8.61329057080909 | 0.161695606903725 | 0.264084645130497 |
| GSE56931 | mortality | yellow | 14 | 130 | 0.0671274627141432 | 0.0177643086151436 | 0.0408579098148302 | 6.89325179690463 | 0.326799139320041 | 0.508465040533105 |

### 7.4 복합 시계 H2

구성 점수는 마이크로어레이가 단계 2 MCPcounter 점수, GSE113883이 §6.5 점수다. 조정 모형은 8개 혈액 집단을 더한 REML이다. 부트스트랩 추출은 §7.2와 같다.

| dataset | clock | n_boot | n_converged_both_reml | p_med | p_med_p025 | p_med_p975 | two_a_orth | two_a_orth_p025 | two_a_orth_p975 | lrt_p_unadj |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | TA_chrono | 2000 | 2000 | 0.542336640271811 | -0.00294389882458876 | 1.24769032899823 | 0.0276903644889938 | 0.00261801410304637 | 0.0533707482509308 | 0.000235567727394529 |
| GSE39445 | TA_mort | 2000 | 2000 | 0.381125400689681 | -0.339161623894615 | 1.03190871771426 | 0.161409023182385 | 0.00912634683154319 | 0.354787891327057 | 0.000615890141275604 |
| GSE56931 | TA_chrono | 2000 | 2000 | -0.427352460718101 | -1.34929896396744 | 0.912007989005958 | 0.012640021206852 | 0.000714767991185486 | 0.0345436413608205 | 0.029920138819241 |
| GSE56931 | TA_mort | 2000 | 2000 | -0.231528920792587 | -1.65443398283297 | 2.15248411959897 | 0.0470034993375978 | 0.0028741154305547 | 0.171342128344774 | 0.157968351876906 |
| GSE113883 | TA_chrono | 2000 | 2000 | 0.705733971429448 | 0.361173919734201 | 1.14866299455836 | 0.0121381718104434 | 0.00059876077454297 | 0.0395192231484417 | 0.000751939790616891 |
| GSE113883 | TA_mort | 2000 | 2000 | 0.653179763882686 | -4.11325442494246 | 3.27695421305645 | 0.0632261677410873 | 0.00209492657032637 | 0.21774005276563 | 0.745144998027147 |

### 7.5 복합 시계와 RA_Peters의 개인 내 상관

RA_Peters lin은 절편 + Peters 계수 · log2 (RNA-seq은 `log2(TPM+1)`)이다. 상관은 control 샘플에서 참가자 평균을 뺀 값끼리의 Pearson이다.

| dataset | clock | other | n_samples | n_subjects | pearson_within |
| --- | --- | --- | --- | --- | --- |
| GSE39445 | TA_chrono | RA_Peters_lin | 221 | 24 | -0.286431165735805 |
| GSE39445 | TA_mort | RA_Peters_lin | 221 | 24 | -0.302321387372661 |
| GSE56931 | TA_chrono | RA_Peters_lin | 130 | 14 | -0.0509675037347897 |
| GSE56931 | TA_mort | RA_Peters_lin | 130 | 14 | -0.016011385011495 |
| GSE113883 | TA_chrono | RA_Peters_lin | 153 | 11 | -0.0238493224772275 |
| GSE113883 | TA_mort | RA_Peters_lin | 153 | 11 | -0.0136643443406104 |

## 8. H5 개발

재학습은 하지 않았다. 부트스트랩은 하지 않았다. limma 설계는 `~ subject_id + cos(2*pi*t/24) + sin(2*pi*t/24)`, t = `clock_time_h`, subject는 고정효과다. cos·sin 두 계수의 조절 F 검정 뒤 BH q다. 리듬 유전자는 q < 0.05이다. 마이크로어레이는 유한한 유전자 전부, GSE113883은 mean TPM ≥ 1이고 log2가 유한한 유전자다. GSE113883은 전체 153샘플이다.

### 8.1 리듬 유전자 수 (자료별)

| dataset | n_genes_tested | n_rhythmic_q_lt_0_05 | gene_filter |
| --- | --- | --- | --- |
| GSE39445 | 19541 | 10041 | all finite genes |
| GSE56931 | 20755 | 4087 | all finite genes |
| GSE113883 | 12390 | 8755 | mean TPM >= 1 and finite log2 |

거를 목록은 GSE39445가 R_56931 ∪ R_113883, GSE56931이 R_39445 ∪ R_113883, GSE113883이 R_39445 ∪ R_56931이다.

### 8.2 거른 시계의 유전자 수 (시계 × 평가 자료)

| dataset | clock | n_genes_original | n_genes_removed | n_genes_remaining |
| --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | 229 | 122 | 107 |
| GSE39445 | RA_Peters | 280 | 222 | 58 |
| GSE39445 | RA_default | 214 | 122 | 92 |
| GSE56931 | RA_GTEx | 255 | 163 | 92 |
| GSE56931 | RA_Peters | 278 | 241 | 37 |
| GSE56931 | RA_default | 235 | 157 | 78 |
| GSE113883 | RA_GTEx | 271 | 144 | 127 |
| GSE113883 | RA_Peters | 276 | 240 | 36 |
| GSE113883 | RA_default | 242 | 145 | 97 |

### 8.3 평가 표 (시계 × 평가 자료: 원래·거른 R2_within, LRT p, 2A, 2A/SD_between, 참가자 평균 상관)

참가자 평균 상관은 원래 점수와 거른 점수의 참가자 평균끼리 Pearson이다. 같은 시계의 original 행과 filtered 행에 같은 값이 들어 있다. 절편과 남은 계수는 원래 값 그대로다.

| dataset | clock | version | n_subjects | n_samples | r2_within | lrt_p | two_a | two_a_over_sd_between | pearson_subject_mean | reml_converged | ml_converged | null_converged |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | original | 24 | 221 | 0.00480935568643748 | 0.643663244864561 | 4.56221470455666 | 0.242731389057453 | 0.952675076502066 | TRUE | TRUE | TRUE |
| GSE39445 | RA_GTEx | filtered | 24 | 221 | 0.00839588294303228 | 0.452197903065829 | 6.15997040057254 | 0.294740639816505 | 0.952675076502066 | TRUE | TRUE | TRUE |
| GSE39445 | RA_Peters | original | 24 | 221 | 0.0237318052396572 | 0.0822398090507771 | 3.86882060477724 | 0.25915869900759 | 0.320866373280244 | TRUE | TRUE | TRUE |
| GSE39445 | RA_Peters | filtered | 24 | 221 | 0.0472106995689958 | 0.0093673944669514 | 2.58567618638333 | 0.339323291824199 | 0.320866373280244 | TRUE | TRUE | TRUE |
| GSE39445 | RA_default | original | 24 | 221 | 0.00680181474094321 | 0.571779550797631 | 1.98498899424706 | 0.211605758610427 | 0.700138202872276 | TRUE | TRUE | TRUE |
| GSE39445 | RA_default | filtered | 24 | 221 | 0.00817794646007386 | 0.469347076758069 | 1.70658886044549 | 0.186035983629465 | 0.700138202872276 | TRUE | TRUE | TRUE |
| GSE56931 | RA_GTEx | original | 14 | 130 | 0.0366339113390549 | 0.0977660921595171 | 8.4129561513483 | 0.69674420051673 | 0.702451762976866 | TRUE | TRUE | TRUE |
| GSE56931 | RA_GTEx | filtered | 14 | 130 | 0.0292048630098719 | 0.171004651604909 | 7.35272171637796 | 0.752032359637741 | 0.702451762976866 | TRUE | TRUE | TRUE |
| GSE56931 | RA_Peters | original | 14 | 130 | 0.145491804603583 | 0.00010511653009091 | 6.37997357324778 | 0.451478484596109 | 0.698616403013853 | TRUE | TRUE | TRUE |
| GSE56931 | RA_Peters | filtered | 14 | 130 | 0.0640713388376701 | 0.0202408622391958 | 3.61789049776709 | 0.39054572983295 | 0.698616403013853 | TRUE | TRUE | TRUE |
| GSE56931 | RA_default | original | 14 | 130 | 0.0658778147622046 | 0.0187399469939367 | 4.61694071235948 | 0.634553234472917 | 0.464065753956846 | TRUE | TRUE | TRUE |
| GSE56931 | RA_default | filtered | 14 | 130 | 0.0211941126866814 | 0.284301850441114 | 1.86363980430221 | 0.472768073002708 | 0.464065753956846 | TRUE | TRUE | TRUE |
| GSE113883 | RA_GTEx | original | 11 | 153 | 0.0315491727978717 | 0.104384139695722 | 4.34369415553687 | 0.405335937361681 | 0.358379331746022 | TRUE | TRUE | TRUE |
| GSE113883 | RA_GTEx | filtered | 11 | 153 | 0.003459203222695 | 0.802231496868176 | 1.3214850845559 | 0.259805977275101 | 0.358379331746022 | TRUE | TRUE | TRUE |
| GSE113883 | RA_Peters | original | 11 | 153 | 0.03527416607544 | 0.077344312576993 | 3.15518568155639 | 0.229721092919869 | 0.446910683302218 | TRUE | TRUE | TRUE |
| GSE113883 | RA_Peters | filtered | 11 | 153 | 0.0383896778926445 | 0.059349715619989 | 1.22716198363208 | 0.31122997811319 | 0.446910683302218 | TRUE | TRUE | TRUE |
| GSE113883 | RA_default | original | 11 | 153 | 0.134425171769798 | 3.24675832907807e-05 | 6.77016594474564 | 0.954019901666464 | 0.57106656133164 | TRUE | TRUE | TRUE |
| GSE113883 | RA_default | filtered | 11 | 153 | 0.026148980350796 | 0.141001085805726 | 2.36336937589158 | 0.747581988881185 | 0.57106656133164 | TRUE | TRUE | TRUE |

## 9. 버그·재실행 기록 (없으면 "없음")

`scripts/26_partC_rnaseq.R`의 첫 두 실행은 exit 1이었다. `data.table`의 `[i, j]` 안에서 `clock`이 반복 변수가 아니라 `clock` 열로 잡혀, official과 lin을 합칠 때 행이 1377개로 늘었다. `which()`를 `[` 밖에서 계산하도록 고친 뒤 세 번째 실행이 exit 0이었다. 실패한 실행은 부트스트랩 전에 멈췄고, 그 실행의 QC와 시계 점수 파일은 성공한 실행이 다시 썼다. Part D와 Part E의 첫 실행은 exit 0이었다. 혼합모형에서 `boundary (singular) fit` 경고가 나왔다. 수렴 판정은 optimizer 코드가 0인지이고, singular 경고와 별개다.

## 10. 판단이 필요한 사항

- Zenodo 메타데이터의 `version`이 빈칸이다. 공개일 순서로 22166800을 버전 2로 썼다. 그보다 새로운 기록은 없었다.
- 기록의 복합 pkl은 60개, 약 58.12 GB다. Part D가 쓰는 파일은 설치류 multitissue scaled-diff 복합 시계 2개와 설치류 모듈 압축이다. 이름·크기·md5가 맞은 EN multitissue scaleddiff 5개와, 이름을 맞춰 받은 모듈 압축 2개만 배치했다. 나머지 복합 pkl 55개(58,122,080,804바이트)와 학습 파일 8개(5,644,397,134바이트)는 받지 않았다.
- 로컬 `rodent_modules.zip`과 `multispecies_modules.zip`은 크기와 md5가 Zenodo와 같았으나 파일명이 달랐다. §2.2의 이름 조건을 통과시키지 않고, Zenodo 파일명으로 다시 받아 md5를 확인했다.
- 사람 채점 파일의 커밋터 날짜는 `clock-specificity`가 늦다. 그 스크립트의 목록은 다종 복합 시계와 다종 모듈 14개이고 `predict`가 아니다. Part D 입력은 지시서 3 §7.1이다. 로그 밑은 바꾸지 않았다. 오솔로그는 `map_to_mouse_entrez`다. 상대 발현은 control 평균을 뺀 값이다. 예측은 설치류 복합 2개와 모듈 46개(목표 2 × 23)의 `predict`다. `allmodulegenes`는 23개에 넣지 않았다.
- 복합 시계 특징 수는 16,517과 16,760이다. 논문에 적힌 10,487과 다르다. 모형을 고치지 않았다.
- QC 샘플 중앙값의 척도는 `log2(TPM+1)`이다. 유전자 거름은 mean TPM ≥ 1이다.
- `hours_since_wake`와 `clock_time_h`의 선형 차이는 6 또는 −18이고, 24로 나눈 나머지는 전부 6이다. T_clock과 T_wake는 153샘플로 맞췄다.

## 11. 지시서와 다르게 한 것

- 시계 파일이 로컬에 모두 있지 않았다. 시계 파일 전부를 받지 않고, §10의 7개만 배치했다. 학습 자료는 받지 않았다.
- `sc-module-aging/requirements.txt`에 적힌 패키지 중 시계 적재에 쓰이지 않는 것은 가상환경에 설치하지 않았다. 설치한 고정 버전은 §2에 적었다.
- `scripts/27_partD_score.py`는 `/usr/bin/time -v`로 재지 않았다. 경과 시간은 6.29초였고 RSS는 없다.
- 선택한 저장소의 `scripts/10_score_human.py`는 다종 시계를 가중치로 채점한다. Part D는 그 스크립트를 실행하지 않고, 지시서 3 §7.1의 입력과 설치류 파이프라인 `predict`를 썼다.

## 12. 산출물 목록, Part별 커밋 해시

| Part | 커밋 | 메시지 |
| --- | --- | --- |
| A | `c6a54b752f9d5f734ecfb0dedde92618f35adc16` | TOD4 Part A |
| B | `f1ca59128eff3b165bf86466def2e2fdc1c4b037` | TOD4 Part B |
| C | `96c5fbe2cfa45ec1228e21633c12b941440f8eed` | TOD4 Part C |
| D | `9ddbb62c940ae8ee3d4a906904626df3e01826c9` | TOD4 Part D |
| E | `963e75aee3e3b97870c0a537ae3e9d048f76b45c` | TOD4 Part E |

Part A: `config/decisions.tsv`, `results/stage4/partA_decisions.tsv`, `gate_rejudge.tsv`, `gse113883_dlmo_values.tsv`, `gse113883_dlmo_decision.tsv`, `gse113883_hours_since_wake.tsv`, `gse113883_hour_minus_wake_summary.tsv`, `scripts/23_partA_decisions.py`.

Part B: `config/paths.yaml`, `config/env_versions.txt`, `results/stage4/module_clocks/`의 Zenodo·md5·구조·포함률 표, `scripts/24_partB_clocks.py`, `scripts/25_partB_structure.py`. 시계 pkl과 zip은 TOD 밖에 두었고 커밋하지 않았다.

Part C: `results/stage4/rnaseq/`, `results/stage4/figures/`, `scripts/26_partC_rnaseq.R`.

Part D: `results/stage4/module_clocks/`의 점수·H1·H2 표, `scripts/27_partD_score.py`, `scripts/28_partD_analyze.R`.

Part E: `results/stage4/h5/`, `scripts/29_partE_h5.R`.
