# TOD 보고서 1 — 자료 기반 구축, 시각 변수 검증, 시계 준비, 개발 탐색

## 1. 요약 (수치만, 12줄 이내)

- GSE39445: 참가자 26, control 24명·221샘플, sleep_restriction 24명·217샘플. 참고 샘플 221; 217과 같다. 참고 참가자 24와 전체 참가자 26은 다르다.
- GSE56931: 참가자 14, control 130, sleep_deprivation 79, recovery 40. 참고 130; 119과 샘플 합이 같다.
- GSE48113: 참가자 22, control/D1 148, out_of_phase/D1 1, out_of_phase/D4 138. out_of_phase 합 139. 참고 147; 139. control 샘플은 참고보다 1개 많다.
- GSE122725: 참가자 11, 샘플 141. `clock_time_h`는 전원 NA.
- 잠금 파일 2개의 SHA-256을 기록했다. 헤더 추출 4건 모두 `!series_matrix_table_begin`에서 정지, 그 앞의 non-`!` 줄은 0. 잠금 발현 값 열람 0건.
- 양성 대조: GSE39445 review, GSE56931 review. 기준 (a) 미충족, (b) 충족. Part F 미실행.
- RNAAgeCalc blood 서명 유전자 수: GTExAge 312, Peters 289, DESeq2(기본) 279. official 대체 유전자: GSE39445 83, 9, 65. GSE56931 57, 11, 44. GPL24324 포함률은 기호 열이 없어 계산하지 않음.
- 연구자 모듈 시계: 미제공 (`module_clock_dir` 비어 있음).
- 개발 탐색 2A, 정점 시각, LRT: Part F 미실행.
- 마지막 분석 커밋: `d00573e95fce5987cb0d155685df51f7eefac437` (TOD Part F).

## 2. 환경과 자원

- OS: Ubuntu 22.04.5 LTS. 커널 Linux 6.8.0-138-generic.
- CPU: AMD Ryzen 5 3600 6-Core Processor, 스레드 12. `MemTotal` 65754920 kB.
- Python `/home/jkang/anaconda3/bin/python3` 3.12.7. pandas 2.3.3, numpy 1.26.4, scipy 1.13.1, statsmodels 0.14.6, matplotlib 3.9.2.
- R `/home/jkang/anaconda3/envs/r_env/bin/Rscript`, R 4.5.3 (2026-03-11). 라이브러리 `/ssd240/anaconda3/envs/r_env/lib/R/library`.
- RNAAgeCalc 1.22.0, GPL-2. impute 1.84.0, GPL-2. lme4 2.0-6, GPL (>= 2). data.table 1.18.6.1, MPL-2.0. circular 0.5-2, GPL-2.
- `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `NUMEXPR_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `BLIS_NUM_THREADS`는 설정되어 있지 않다. `ulimit -t`와 `ulimit -v`는 unlimited. GPU는 사용하지 않았다.
- RNAAgeCalc 실행을 `/usr/bin/time -v`로 감싼 기록 (`logs/partE_rnaage_time.log`): wall clock 0:44.56, user 46.54 s, system 4.19 s, CPU 113%, maximum resident set size 2194440 kB, exit status 0.

## 3. 자료 목록

series matrix는 네 자료 모두 1개였다. 개발 자료의 보충 파일은 series matrix에 발현 표가 있어 받지 않았다. 잠금 자료의 보충 파일은 받지 않고 Content-Length만 기록했다. GPL15331과 GPL24324의 `annot.gz`는 HTTP 404였다. family SOFT는 받지 않았다. 플랫폼 레코드는 `acc.cgi?targ=self&view=full&form=text`로 받았다.

| 경로 | 바이트 | SHA-256 |
|---|---:|---|
| data/raw/GSE39445_series_matrix.txt.gz | 99904009 | e7b4abce3b619f5a9df58611613d6ab14b636f1c148ed2032e4cf74d9c81ba42 |
| data/raw/GSE56931_series_matrix.txt.gz | 51487791 | cd20a18fbfb55eb9663ec1a2e8a847eb5b85808979d0fd84f89c6f952c0a5330 |
| data/locked/GSE48113_series_matrix.txt.gz | 72911994 | 6a1991afcd58f862446692fc1dca7ffd437910a8a725af8f2b59c100a1771500 |
| data/locked/GSE122725_series_matrix.txt.gz | 16328295 | af19a3ab11e62638ce6a735f9e4d7b8e9fe73a63a61b8b586d759014185604b9 |
| data/raw/gpl/GPL10379.annot.gz | 8197525 | c25fdff86787178767349e8a467fab42d780959cd035b0c3381eac2aefb4e602 |
| data/raw/gpl/GPL15331_platform.txt | 25434313 | 88295686b1a6b11cece7bbd49d88b8a02073de38e9e94a4eddf6494dec57c4cc |
| data/raw/gpl/GPL24324_platform.txt | 117045079 | 50f87c7ebf05ae1959fbc48d60c566344da13c67aeab50f6e80636855df4489d |

| 자료 | GPL | 역할 | 채널 수 | 행 수 필드 |
|---|---|---|---:|---:|
| GSE39445 | GPL15331 | 개발 | 1 | 43758 |
| GSE56931 | GPL10379 | 개발 | 1 | 52378 |
| GSE48113 | GPL15331 | 잠금 | 1 | 43758 |
| GSE122725 | GPL24324 | 잠금 | 1 | 27189 |

잠금 보충 파일 (미다운로드, Content-Length):

| 자료 | 이름 | 바이트 |
|---|---|---:|
| GSE48113 | GSE48113_RAW.tar | 2585128960 |
| GSE48113 | filelist.txt | 19132 |
| GSE122725 | GSE122725_RAW.tar | 157327360 |
| GSE122725 | filelist.txt | 8835 |

개발 보충 파일 (미다운로드): GSE39445_RAW.tar 3938027520, filelist.txt 28859, GSE56931_RAW.tar 1337272320, filelist.txt 23402.

## 4. 메타데이터

헤더 추출은 `!`로 시작하는 줄만 저장하고 `!series_matrix_table_begin`에서 멈췄다. 네 파일 모두 `stopped_at_table=True`, `non_bang_before_stop=0`.

### 4.1 열 정의와 유도 규칙

**GSE39445.** `subject_id`는 `Sample_characteristics_ch1`의 `subject`. 11개 샘플은 이 키가 없어 `Sample_title`의 `AF`+숫자에서 취했다. 예: `AF0033_SE_T2` → `AF0033`. 제목과 characteristic가 둘 다 있는 427개는 불일치 0. `condition_raw`=`sleepprotocol`. `Sleep Extension`→`control`, `Sleep Restriction`→`sleep_restriction`. `period`는 모든 샘플 `extended_wakefulness`. `sample_type`=`tissue` 값 `Blood`. `clock_time_h`는 `timesampletaken`의 `H:MM`를 시간 소수로 변환 (`given`). 범위 0–23.8333. `hours_awake`=`hoursawake` (`given`). `elapsed_h`는 같은 `subject_id`와 `condition` 안에서 `hours_awake` 최솟값을 0으로 둔 차이 (`derived`). `circadian_ref_h`=`circadianphase` 숫자. `age`, `sex`, `batch_fields`는 NA. 빈 토큰 `NA`는 결측으로 두었다.

Series summary 원문: "In a cross balanced design human subjects were given one week of sufficient sleep and insufficient sleep (6 h per day). Following each of these conditions they were then kept awake for a day, a night and the subsequent day. During each of these periods of extended wakefulness 10 blood samples were taken, RNA was extracted from leukocytes, labelled and hybridised to human whole-genome microarrays"

**GSE56931.** `subject_id`=`patient`. `condition_raw`=`biological group`. `baseline`→`control`/`baseline`, `sleep deprivation`→`sleep_deprivation`/`deprivation`, `recovery`→`recovery`/`recovery`. `sample_type`=`blood` (`Sample_source_name_ch1`의 "blood sample"). `clock_time_h`=`hour` (`given`), 관측값 0, 4, 8, 12, 16, 20. `elapsed_h`=`(day-1)*24+hour`에서 같은 `subject_id`와 `period`의 최솟값을 뺀 값 (`derived`). `hours_awake`, `circadian_ref_h`, `age`, `batch_fields`는 NA. `sex`=`Sex`.

`Series_overall_design` 원문 일부: "24-hour normal baseline (6 time points: day 1 8hr, 12hr, 16hr, 20hr, day 2 0hr, 4hr), 38 hours of continuous wakefulness (10 time points: day 2 8hr, 12hr, 16hr, 20hr, day 3 0hr, 4hr, 8hr, 12hr, 16hr, 20hr) and subsequent recovery sleep (3 time points: day 4 0hr, 4hr, 8hr)"

**GSE48113.** `subject_id`=`subject`. `condition_raw`=`sleep timing`. `In phase with respect to melatonin`→`control`, `Out of phase with respect to melatonin`→`out_of_phase`. `period`는 `day sample was taken in protocol` 1→`D1`, 4→`D4`. `sample_type`=`tissue` 값 `Blood`. `clock_time_h`=`time sample taken` (`given`). 값 `NA`인 1건은 결측. `elapsed_h`는 같은 `subject_id`와 `period`에서 `time point` 순으로 정렬한 뒤, 시계 시각의 전진 차이 `(t2-t1) mod 24`를 누적 (`derived`). `hours_awake`, `circadian_ref_h`, `age`, `sex`, `batch_fields`는 NA.

Series summary 원문 일부: "During D1 and D4 7 blood samples were taken, RNA was extracted from leukocytes, labelled and hybridised to human whole-genome microarrays"

**GSE122725.** `subject_id`=`subject`. `condition_raw`=`condition`. `DayOrientedSchedule`→`control`/`period1`, `NightOrientedSchedule`→`night_oriented`/`period2`. `light_group`: `BrightLight`→`bright`, `Control`→`dim`. Series summary 원문: "Participants in the control group were exposed to dim light throughout the waking periods; participants in the bright group were exposed to bright light (~6,000lux) for 8 hours during the waking period." `sample_type`=`PBMC`. `relclocktime` 관측 범위 0.2–42.3. `clock_time_h`는 NA, `time_source`=`NA`. `elapsed_h`는 같은 `subject_id`와 `period`에서 `relclocktime` 최솟값을 0으로 둔 차이 (`derived`). `hours_awake`, `circadian_ref_h`, `age`, `sex`, `batch_fields`는 NA. 샘플 설명 예: "PBMC sample taken at 0.8 h after habitual wake time during day-oriented schedule in Subject T10 in the BrightLight group"

### 4.2 자료별 요약 표

| 자료 | 참가자 | 샘플 | 참가자당 샘플 min / median / max | clock_time_h | 연속 elapsed 간격 (h) |
|---|---:|---:|---|---|---|
| GSE39445 | 26 | 438 | 6 / 18.5 / 20 | n=438, 0–23.8333 | n=390, min 3, q25 3, median 3, q75 3, max 12 |
| GSE56931 | 14 | 249 | 15 / 18 / 19 | n=249, 0–20 | n=207, min 4, q25 4, median 4, q75 4, max 8 |
| GSE48113 | 22 | 287 | 6 / 14 / 14 | n=286, 0–23.9167 | n=242, min 3.33333, q25 3.98333, median 4, q75 4.02917, max 8.65 |
| GSE122725 | 11 | 141 | 10 / 13 / 14 | NA | n=119, min 1.9, q25 3.9, median 4, q75 4.1, max 8.1 |

조건·구간·빛 그룹별 참가자 수와 샘플 수:

| 자료 | condition | period | light_group | 참가자 | 샘플 |
|---|---|---|---|---:|---:|
| GSE39445 | control | extended_wakefulness | NA | 24 | 221 |
| GSE39445 | sleep_restriction | extended_wakefulness | NA | 24 | 217 |
| GSE56931 | control | baseline | NA | 14 | 130 |
| GSE56931 | sleep_deprivation | deprivation | NA | 14 | 79 |
| GSE56931 | recovery | recovery | NA | 14 | 40 |
| GSE48113 | control | D1 | NA | 22 | 148 |
| GSE48113 | out_of_phase | D1 | NA | 1 | 1 |
| GSE48113 | out_of_phase | D4 | NA | 21 | 138 |
| GSE122725 | control | period1 | bright | 3 | 18 |
| GSE122725 | control | period1 | dim | 8 | 52 |
| GSE122725 | night_oriented | period2 | bright | 3 | 20 |
| GSE122725 | night_oriented | period2 | dim | 8 | 51 |

GSE48113 `clock_time_h` 결측: GSM1168837, subject BB0264, `time sample taken` 원문 `NA`, D1, control, time point 7. D1이면서 out_of_phase인 샘플: GSM1168735, subject BB0121, `time sample taken` `20:25`, `clock_time_h` 20.41666667.

### 4.3 참고값 대조

참고 샘플 수 (지시서 §3, Hughey 2017 Table 1; control; perturbation): GSE39445 221; 217. GSE48113 147; 139. GSE56931 130; 119.

- GSE39445 control 221, sleep_restriction 217. 참고와 같다. 전체 참가자 26. `Series_overall_design`: "A total of 438 samples comprising 26 human subjects". 참고 참가자 수는 24. control 구간의 참가자는 24.
- GSE56931 control 130, sleep_deprivation+recovery 79+40=119. 참고와 같다. 참가자 14. 참고 참가자 14와 같다.
- `biological group`=`baseline`이면서 day 2, hour 8/12/16/20인 샘플은 14+11+14+12=51개다. `Series_overall_design`은 이 시각을 38시간 각성 구간에 적는다. 이 보고서의 `condition`은 샘플 characteristic `biological group`을 따랐다.
- GSE48113 in phase 148, out of phase 139. day 1 샘플 149, day 4 샘플 138. 참고 control 147, perturbation 139. perturbation 합은 같다. control로 둔 in-phase 샘플은 148이다.
- GSE122725 참가자 11. 지시서의 참고 참가자 11과 같다. 지시서의 Hughey 표에는 이 자료의 샘플 수가 없다. `Series_overall_design`: "A total of 141 samples were collected from 11 subjects". 샘플 141.

### 4.4 잠금 자료

GSE48113, GSE122725의 series matrix는 `data/locked/`에 있다. 메타데이터는 §4.1의 헤더 추출로만 읽었다. 발현 표 행은 읽지 않았다. 보충 파일은 받지 않았다. 플랫폼 파일 GPL15331, GPL24324는 잠금 대상이 아니어서 플랫폼 표만 읽었다.

## 5. 개발 자료 처리와 QC

### 5.1 처리 문구와 채널

GSE39445 `Sample_channel_count` 1. `Sample_label_ch1` `Cy3`. 두 채널 비율 자료로 표시하지 않았다.

`Sample_data_processing` 두 문구. 427샘플: "Individual samples were filtered based on AgilentQC metrics of reproducibility statistics, minimum detection level estimates and feature flags. Samples with a median coefficient of variation of less than 10% in spike ins or non-control replicated probes (NCRPs) were retained and quantile-normalized using the R Bioconductor package limma. NCRPs along with their corresponding flags were averaged. Probes with more than 5 flagged samples within a subject in more than half of the total number of subjects were excluded. In addition, for time-series analyses, series with 2 consecutive or more than 2 missing time points were excluded."

11샘플은 같은 문장에서 `quantile-normalized`가 `quantile-normalised`이다.

GSE56931 `Sample_channel_count` 1. `Sample_label_ch1` `biotin`. `Sample_data_processing`: "Data were processed using R/affy software. In brief, the RMA method adjusts the background of perfect match (PM) probes, applies a quantile normalization of the corrected PM values, and calculates final expression measures with the Tukey median polish algorithm."

### 5.2 로그 판정

판정 규칙: `q99>100`이거나 (`q100-q0>50`이고 `q25>0`)이면 선형, 아니면 로그. 분위수는 숫자 값만 사용했다. 입력 NA는 두 자료 모두 0.

| 자료 | n_numeric | q0 | q25 | q99 | q100 | 판정 | 0 이하를 NA로 바꾼 수 |
|---|---:|---:|---:|---:|---:|---|---:|
| GSE39445 | 19166004 | 1.121515378 | 3.78166512 | 15.72362557 | 19.45222559 | log | 0 |
| GSE56931 | 13042122 | 1.789731 | 4.172186 | 12.5924099 | 14.571592 | log | 0 |

로그로 판정되어 추가 log2 변환은 하지 않았다. 추가 정규화는 하지 않았다.

### 5.3 프로브 → 유전자

기호가 없거나, `///`로 나눈 고유 기호가 1개가 아니면 제외했다. `/` 한 글자는 구분자로 쓰지 않았다 (예: `OK/SW-CL.36`). 한 유전자에 프로브가 여럿이면 전체 샘플 평균 log2가 가장 큰 프로브 하나. 같으면 probe ID 사전순 앞선 것.

| 자료 | 기호 열 | 행렬 프로브 | 기호 없음 | 복수 기호 | 미선택 | 선택 유전자 | 플랫폼에만 있는 행 |
|---|---|---:|---:|---:|---:|---:|---:|
| GSE39445 | GENE_SYMBOL | 43758 | 13253 | 0 | 10964 | 19541 | 1462 |
| GSE56931 | Gene symbol | 52378 | 15149 | 451 | 16023 | 20755 | 0 |

행렬 프로브 중 플랫폼에 없는 ID는 0이었다.

### 5.4 QC

샘플 중앙값들의 중앙값과 MAD(중앙값 절대편차, 계수 1.4826을 곱하지 않음). 샘플 중앙값이 그 중앙값에서 3 MAD를 넘으면 표시만 하고 빼지 않았다.

| 자료 | 중앙값들의 중앙값 | MAD | 3 MAD 표시 수 |
|---|---:|---:|---:|
| GSE39445 | 7.770000963 | 0.020362502 | 20 |
| GSE56931 | 6.822182 | 0.031288 | 24 |

샘플 중앙값 분포 (min, q25, median, q75, max): GSE39445 7.633000996, 7.750925845, 7.770000963, 7.792720936, 7.841708466. GSE56931 6.652634, 6.787204, 6.822182, 6.849853, 6.995107.

샘플 IQR 분포: GSE39445 5.062749477, 5.195219443, 5.254108216, 5.308512607, 5.381699061. GSE56931 4.5528685, 4.8710105, 4.9382855, 4.9880375, 5.1128625.

GSE39445 표시 GSM: GSM968833, GSM968844, GSM968861, GSM968872, GSM968884, GSM968886, GSM968978, GSM968979, GSM968980, GSM969017, GSM969035, GSM969070, GSM969094, GSM969142, GSM969151, GSM969152, GSM969168, GSM969197, GSM969202, GSM969203.

GSE56931 표시 GSM: GSM1371587, GSM1371596, GSM1371604, GSM1371609, GSM1371618, GSM1371629, GSM1371630, GSM1371644, GSM1371681, GSM1371700, GSM1371709, GSM1371725, GSM1371727, GSM1371732, GSM1371746, GSM1371759, GSM1371760, GSM1371787, GSM1371792, GSM1371798, GSM1371801, GSM1371822, GSM1371834, GSM1371835.

PCA: 결측이 없는 유전자 중 분산 상위 1000개. 유전자별로 평균을 빼고 스케일은 맞추지 않은 뒤 SVD. PC1–PC5에 대한 R². `batch_fields`는 모든 샘플이 NA라 R²를 계산하지 않았다.

GSE39445 (n=438):

| PC | subject_id | condition | period | cos·sin(clock) |
|---|---:|---:|---:|---:|
| PC1 | 0.8068949456 | 0.002939063199 | 0 | 0.01752077317 |
| PC2 | 0.5766034049 | 0.0001123102288 | 0 | 0.0006618830219 |
| PC3 | 0.7321092484 | 0.004157637825 | 0 | 0.004520282996 |
| PC4 | 0.9606407844 | 0.003272297815 | 0 | 0.002023888936 |
| PC5 | 0.9360354955 | 4.040146747e-05 | 0 | 0.000174312646 |

GSE56931 (n=249). `condition`과 `period`는 한 쌍으로 대응한다.

| PC | subject_id | condition | period | cos·sin(clock) |
|---|---:|---:|---:|---:|
| PC1 | 0.9617240692 | 0.003308193322 | 0.003308193322 | 0.002532593096 |
| PC2 | 0.2850642903 | 0.001809621583 | 0.001809621583 | 0.005847934994 |
| PC3 | 0.8541740857 | 0.01194587511 | 0.01194587511 | 0.008930864642 |
| PC4 | 0.9140286427 | 0.0006595007544 | 0.0006595007544 | 0.001904526739 |
| PC5 | 0.9388482253 | 5.362595488e-05 | 5.362595488e-05 | 0.002354867082 |

## 6. 양성 대조

자료: 개발 자료의 `condition=control`이고 `clock_time_h`가 있는 샘플. GSE39445 24명, 221샘플. GSE56931 14명, 130샘플. 유전자 10개는 두 자료 모두 있었다. 모형 `expr ~ cos(2*pi*t/24) + sin(2*pi*t/24) + (1|subject_id)`. 계수는 REML, LRT는 ML, df=2. 3 MAD 표시 샘플은 빼지 않았다.

| 자료 | 유전자 | β_cos | β_sin | A | 정점(h) | LRT p | REML singular |
|---|---|---:|---:|---:|---:|---:|---|
| GSE39445 | PER1 | -0.0514596622941169 | 0.21354531415074 | 0.219658139023263 | 6.90324597433964 | 9.73597605831301e-18 | FALSE |
| GSE39445 | PER2 | 0.0325659288766659 | 0.133904533913401 | 0.13780770633809 | 5.08872615108073 | 5.71478055338033e-17 | FALSE |
| GSE39445 | PER3 | -0.00351865762853524 | 0.156560302975653 | 0.156599838503541 | 6.08583286807854 | 1.177374716677e-06 | FALSE |
| GSE39445 | NR1D1 | 0.205288042638116 | 0.134238723922147 | 0.245281910157344 | 2.21205623059111 | 1.67152276874189e-24 | FALSE |
| GSE39445 | NR1D2 | 0.096616059297151 | 0.186293706172415 | 0.209857112992543 | 4.17251607760881 | 3.14695497771e-18 | FALSE |
| GSE39445 | ARNTL | -0.0416712806403648 | -0.0387617345474874 | 0.0569119292884886 | 14.8618875133315 | 0.0449426621692836 | FALSE |
| GSE39445 | CRY1 | -0.00311843685572919 | 0.0959261399979879 | 0.095976814821793 | 6.12413048116649 | 0.041580574897009 | FALSE |
| GSE39445 | CRY2 | 0.04219988857847 | 0.0417873639671845 | 0.059388672180486 | 2.98123863770304 | 4.168248177311e-06 | FALSE |
| GSE39445 | DBP | -0.390425266447806 | -0.0503719304365573 | 0.39366130119272 | 12.4901054717509 | 2.04983377976056e-07 | FALSE |
| GSE39445 | TEF | 0.145798757902324 | 0.0587803398964352 | 0.157201800766406 | 1.46382392673113 | 0.0194597641265044 | FALSE |
| GSE56931 | PER1 | 0.0195976933024504 | 0.358765878129688 | 0.35930074435345 | 5.79155389251083 | 1.58787990764694e-24 | FALSE |
| GSE56931 | PER2 | -0.0342918280011894 | 0.150559585206476 | 0.15441540779731 | 6.85539646901579 | 1.36583694162499e-10 | FALSE |
| GSE56931 | PER3 | 0.0662695246970921 | 0.222557801464182 | 0.232214609566558 | 4.8945603273957 | 3.48666322673004e-14 | FALSE |
| GSE56931 | NR1D1 | 0.130812916051635 | 0.0459216361233354 | 0.13863915633823 | 1.28956873301937 | 1.40348603486936e-09 | FALSE |
| GSE56931 | NR1D2 | 0.153020508544788 | 0.123553986482304 | 0.196674511848828 | 2.59457275718178 | 3.80424994659179e-18 | FALSE |
| GSE56931 | ARNTL | -0.0452411776051107 | -0.0399551472397526 | 0.0603587437083268 | 14.7633078773663 | 2.95392758328053e-05 | FALSE |
| GSE56931 | CRY1 | -0.00206435684460958 | 0.0334790231540025 | 0.0335426081354465 | 6.23523062419983 | 0.546337341639908 | FALSE |
| GSE56931 | CRY2 | -0.00497694401189 | 0.040087501615455 | 0.0403952689985671 | 6.471811516987 | 0.0447041058504479 | FALSE |
| GSE56931 | DBP | 0.102898854612372 | 0.0888378706467456 | 0.135942419948984 | 2.7203820522671 | 3.8762589576806e-08 | FALSE |
| GSE56931 | TEF | 0.0962703088005004 | 0.231485497334631 | 0.250706018740686 | 4.49456836942552 | 3.59117197797277e-22 | FALSE |

기준 (a): PER1, PER3, NR1D1, NR1D2 정점이 4 이상 12 이하. GSE39445는 NR1D1 2.21205623059111이 범위 밖. GSE56931은 NR1D1 1.28956873301937, NR1D2 2.59457275718178이 범위 밖. 기준 (b): PER3와 ARNTL의 원형 거리(0–12 h). GSE39445 8.77605464525298. GSE56931 9.86874754997058. 둘 다 6 이상.

| 자료 | (a) | (b) | 거리(h) | 판정 |
|---|---|---|---:|---|
| GSE39445 | FALSE | TRUE | 8.77605464525298 | review |
| GSE56931 | FALSE | TRUE | 9.86874754997058 | review |

`clock_time_h`는 바꾸지 않았다. REML 적합의 `isSingular`은 20개 모두 FALSE. ML 적합도 모두 FALSE. null 모형 `isSingular`은 GSE56931 PER1만 TRUE. 실행 중 stderr에 `boundary (singular) fit: see help('isSingular')`가 1회 출력되었다.

## 7. 시계 준비

### 7.1 RNAAgeCalc

버전 1.22.0. 설치 경로 `/ssd240/anaconda3/envs/r_env/lib/R/library/RNAAgeCalc`.

| 파일 | 바이트 | SHA-256 |
|---|---:|---|
| DESCRIPTION | 1695 | ec789113f0e224e22d512b9504ee129107f90cb15047f331d75dfa56c9bc23c3 |
| R/sysdata.rdb | 2782564 | 5193841f4e4ab035f4a55e708490f04682eaaed0363f1da8bca8539499b9c086 |
| R/sysdata.rdx | 188 | b855fd22b795d8aeb97641f000d006998e5d96afde4e142d83453712c369619e |
| data/fpkmExample.RData | 387664 | 7add4c089abff4a2103d97d222461135354a3e11e9c5feaaf8b8453a3e81a447 |
| data/countExample.RData | 164584 | 93687a8412516ff2694c1eca3683df754f565a98d83c1bac49b54a9ce94cdb1c |

`predict_age` 소스에서 확인한 내용:

- 계수는 `genelist_all[[tissue]][[signature]]`. `stype="all"`이면 `genelist_all`, `stype="caucasian"`이면 `genelist_cau`. 벡터의 첫 원소 이름이 `(Intercept)`이고 나머지가 유전자 계수다.
- `tissue`를 주고 그 조직이 `names(genelist_all)`에 있으며 `signature`를 생략하면 서명은 `DESeq2`다. `blood`는 그 이름 목록에 있다. 실행 메시지: "signature is not provided, using DESeq2 signature automatically."
- `exprtype="counts"`이면 먼저 `count2FPKM`. 그 다음, 결측 유전자가 있으면 `log2(expr+1)` 행렬에 없는 서명 유전자를 NA 행으로 붙여 `impute.knn`을 호출한다. 결측이 없으면 `log2(선택한 유전자 + 1)`에 계수를 곱한다. `maxp`를 지정하지 않으면 `maxp=nrow(log2 행렬)`.
- 예측값 `sum(x * coefs[-1]) + coefs[1]`.

`impute` 1.84.0 도움말: `impute.knn`은 유전자 공간에서 k=10 최근접 이웃의 유클리드 거리로 결측을 평균한다. `rowmax` 기본 0.5. 도움말은 `rowmax`를 넘는 행은 샘플별 전체 평균으로 채운다고 적는다. 실행 확인: GSE39445 log2 행렬에서 유전자 100개·샘플 30개 블록에 NA 행 1개를 붙여 `impute.knn(maxp=nrow)`를 호출하면 그 행은 모든 샘플에서 0이었다. 관측 유전자의 열 평균과의 최대 절대차는 7.32891805584894. 경고 문구: "1 rows with more than 50 % entries missing; mean imputation used for these rows". `predict_age`는 `impute.knn` 경고를 `suppressWarnings`로 감싼다.

서명 감사. 문서 유전자 수는 패키지 비네트 `RNAAge-vignette.Rmd`의 Peters "1,497 genes"만 해당한다. GTExAge와 DESeq2의 유전자 수는 그 비네트에 숫자로 적혀 있지 않다.

| clock | signature | 절편 | 유전자 수 | 문서 수 | 일치 | 다음 단계로 넘김 |
|---|---|---:|---:|---:|---|---|
| RA_GTEx | GTExAge | 64.2016047245696 | 312 | 없음 | NA | TRUE |
| RA_Peters | Peters | 59.6846198087602 | 289 | 1497 | FALSE | FALSE |
| RA_default | DESeq2 | 65.5365453354773 | 279 | 없음 | NA | TRUE |

플랫폼 포함. 기호는 `///` 토큰이 그 플랫폼 주석에 있으면 있는 것으로 셌다. GPL24324 플랫폼 표 열은 `ID, probeset_id, seqname, strand, start, stop, total_probes, category, SPOT_ID, SPOT_ID`이다. 유전자 기호 열이 없어 포함률을 계산하지 않았다.

| 서명 | GPL15331 | GPL10379 | GPL24324 |
|---|---|---|---|
| RA_GTEx 312 | 229 (0.733974359) | 261 (0.8365384615) | 계산 안 함 |
| RA_Peters 289 | 280 (0.9688581315) | 278 (0.9619377163) | 계산 안 함 |
| RA_default 279 | 214 (0.7670250896) | 236 (0.8458781362) | 계산 안 함 |

값 계산은 개발 자료만. 두 자료 모두 채널 수 1이라 `official`을 적용했다. `official` 입력은 유전자 log2 행렬의 `2^값`, `exprtype="FPKM"`, `tissue="blood"`, `stype="all"`, `idtype="SYMBOL"`. `lin`은 `log2(FPKM+1)` 자리에 마이크로어레이 log2 값을 그대로 넣고, 없는 유전자는 빼고, 대체는 하지 않았다.

`lin`에 쓴 유전자 수:

| 자료 | clock | 서명 유전자 | 사용 |
|---|---|---:|---:|
| GSE39445 | RA_GTEx | 312 | 229 |
| GSE39445 | RA_Peters | 289 | 280 |
| GSE39445 | RA_default | 279 | 214 |
| GSE56931 | RA_GTEx | 312 | 255 |
| GSE56931 | RA_Peters | 289 | 278 |
| GSE56931 | RA_default | 279 | 235 |

`official` 자동 대체 수 (`n_genes_imputed`, 서명에 있고 발현 행렬에 없는 유전자):

| 자료 | clock | 대체 수 | 메시지에 적힌 비율 | elapsed_sec |
|---|---|---:|---|---:|
| GSE39445 | RA_GTEx | 83 | 26.6026% | 5.17 |
| GSE39445 | RA_Peters | 9 | 3.1142% | 2.448 |
| GSE39445 | RA_default | 65 | 23.2975% | 2.139 |
| GSE56931 | RA_GTEx | 57 | 18.2692% | 1.705 |
| GSE56931 | RA_Peters | 11 | 3.8062% | 0.890000000000001 |
| GSE56931 | RA_default | 44 | 15.7706% | 0.957000000000001 |

예측값 분포:

| 자료 | clock | method | n | mean | sd | min | max |
|---|---|---|---:|---:|---:|---:|---:|
| GSE39445 | RA_GTEx | lin | 438 | 35.8936597173399 | 29.8357960667814 | -89.5649004175149 | 115.520719729785 |
| GSE39445 | RA_GTEx | official | 438 | 471.369215211439 | 27.2362568401215 | 361.799470821225 | 546.25311401459 |
| GSE39445 | RA_Peters | lin | 438 | -63.3236929793359 | 17.5286191013874 | -125.077471096823 | -18.7166934727708 |
| GSE39445 | RA_Peters | official | 438 | 14.9924944317794 | 16.9426377741756 | -45.4027188543937 | 59.157660734378 |
| GSE39445 | RA_default | lin | 438 | -71.5262223236156 | 12.3698832803123 | -102.459397791902 | -11.4602204158499 |
| GSE39445 | RA_default | official | 438 | -37.1492143898583 | 11.9873840099039 | -65.9200204929699 | 20.4051732293915 |
| GSE56931 | RA_GTEx | lin | 249 | 124.673321176093 | 17.3656992270821 | 81.123549514504 | 176.612205787625 |
| GSE56931 | RA_GTEx | official | 249 | 308.377362947738 | 16.6400205148744 | 267.767071416076 | 358.916789581689 |
| GSE56931 | RA_Peters | lin | 249 | 107.404709790173 | 15.1167681774903 | 69.6483335242973 | 143.053918431723 |
| GSE56931 | RA_Peters | official | 249 | 37.547759941717 | 14.9242541482509 | 0.106938341914706 | 72.5363904790724 |
| GSE56931 | RA_default | lin | 249 | -82.4905937571467 | 9.23312161833081 | -107.571775385825 | -54.4191731651428 |
| GSE56931 | RA_default | official | 249 | -40.7704396434804 | 9.10841109694217 | -65.5078187383084 | -13.2402972560872 |

출력: `results/stage1/clocks/GSE39445_rnaage.tsv`, `results/stage1/clocks/GSE56931_rnaage.tsv`.

### 7.2 연구자 모듈 시계

`config/paths.yaml`의 `module_clock_dir`는 비어 있다. 연구자 자료 미제공. 다른 경로에서 시계를 구하지 않았다.

### 7.3 세포 분해 참조

MCPcounter 1.2.0. DESCRIPTION License `GPL-3`. DESCRIPTION SHA-256 `a6c57cc11cce22b43d2d1fb2592bb86fa3348a0aa50ea452ce0cde0c190e8f8d`, 538바이트. 표지 파일 `Signatures/genes.txt`는 CC0 1.0. SHA-256 `408f6c5d02c8f9bd2f1599c598367881853520934f8ce398a4886a8b296922bb`, 5342바이트. HUGO 기호 111개, 세포집단 10개. 비율 추정은 하지 않았다.

| 집단 | n | GPL15331 | GPL10379 | GPL24324 |
|---|---:|---|---|---|
| B lineage | 9 | 8 (0.8888888889) | 9 (1) | 계산 안 함 |
| CD8 T cells | 1 | 1 (1) | 1 (1) | 계산 안 함 |
| Cytotoxic lymphocytes | 7 | 7 (1) | 7 (1) | 계산 안 함 |
| Endothelial cells | 33 | 30 (0.9090909091) | 30 (0.9090909091) | 계산 안 함 |
| Fibroblasts | 8 | 8 (1) | 8 (1) | 계산 안 함 |
| Monocytic lineage | 7 | 7 (1) | 7 (1) | 계산 안 함 |
| Myeloid dendritic cells | 6 | 5 (0.8333333333) | 6 (1) | 계산 안 함 |
| NK cells | 9 | 6 (0.6666666667) | 8 (0.8888888889) | 계산 안 함 |
| Neutrophils | 15 | 15 (1) | 15 (1) | 계산 안 함 |
| T cells | 16 | 13 (0.8125) | 14 (0.875) | 계산 안 함 |
| ALL | 111 | 100 (0.9009009009) | 105 (0.9459459459) | 계산 안 함 |

EPIC 1.1.7의 LICENSE 제목은 "SOFTWARE LICENSE AGREEMENT FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY"이고, 프로그램 재배포는 사전 서면 허가 없이 하지 않는다고 적혀 있다. EPIC는 내려받지 않았다. BRef, TRef 포함률은 없다.

## 8. 개발 탐색 (확인 검정 아님)

### 8.1 주 모형 표

실행하지 않음. Part D 판정이 두 개발 자료 모두 `review`이다.

### 8.2 민감도 S1–S3

실행하지 않음. GSE39445에는 `hours_awake`가 있다. GSE56931의 `hours_awake`는 NA다. 민감도 모형은 돌리지 않았다.

### 8.3 그림 경로

그림 없음.

## 9. 버그·재실행 기록

- 첫 R 패키지 소스 설치는 `x86_64-conda-linux-gnu-cc`가 PATH에 없어 컴파일이 실패했다. conda 환경의 `bin`을 PATH에 넣은 뒤에도 `zlib.h` 없음으로 `data.table` 컴파일이 실패했다. 이후 conda-forge/bioconda의 빌드된 패키지로 RNAAgeCalc 1.22.0, impute 1.84.0, lme4 2.0-6, data.table 1.18.6.1, circular 0.5-2를 설치했다.
- 중단된 BiocManager 프로세스가 `SummarizedExperiment` 디렉터리를 지운 상태가 있었다. `conda install --force-reinstall bioconductor-summarizedexperiment`로 복구했다. 그 다음 `RNAAgeCalc` 로드와 `predict_age` 실행은 exit status 0이었다.
- 양성 대조는 `isSingular` 열을 넣기 위해 같은 데이터로 다시 맞췄다. 계수와 LRT p는 첫 실행과 같은 값이 출력되었다.

## 10. 판단이 필요한 사항

- GSE56931의 `condition`은 샘플 필드 `biological group`을 따랐다. day 2 hour 8/12/16/20의 baseline 51샘플은 `Series_overall_design`의 각성 구간 시각과 필드 값이 다르다.
- GSE122725 `relclocktime`은 0.2–42.3이고 샘플 설명은 habitual wake 이후 시간이다. `clock_time_h`는 NA로 두었다. `hours_awake`도 NA로 두었다. `elapsed_h`만 `relclocktime` 차이로 계산했다.
- GSE122725 `experimentalgroup=Control`을 `light_group=dim`으로 두었다. 근거는 §4.1의 series summary 문장이다.
- GSE39445 `period`는 샘플 필드에 구간 이름이 없어 모든 행을 `extended_wakefulness`로 두었다.
- GSE48113 `elapsed_h`는 `time point` 순서의 시계 시각 전진 차이로 계산했다. GSM1168837은 시각이 `NA`라 `elapsed_h`도 NA다.
- GPL15331 family SOFT는 306306821바이트, GPL24324 family SOFT는 416883830바이트였다. family 파일은 받지 않았다. 플랫폼 표만 `targ=self`로 받았다.
- GPL24324에는 유전자 기호 열이 없다. 두 번째 `SPOT_ID` 자유 텍스트에서 기호를 뽑지 않았다. 그 플랫폼의 서명·MCPcounter 포함률은 비어 있다.
- Peters blood 계수 유전자 289와 비네트의 1,497은 다르다. `RA_Peters`는 다음 단계로 넘기지 않는다. Part E의 값은 계산해 두었다.
- EPIC는 비상업·재배포 제한 라이선스라 받지 않았다.
- `ulimit -v`는 unlimited이다. 측정된 최대 RSS는 2194440 kB이다.
- `impute.knn` 도움말의 "overall mean per sample"과 실행 결과(전부 0)가 다르다. 위 실행 결과를 기록했다.
- 양성 대조가 `review`라 Part F를 하지 않았다. 기준에 맞추려고 `clock_time_h`를 바꾸지 않았다.

## 11. 지시서와 다르게 한 것

- GPL15331, GPL24324는 `annot.gz`가 없어 플랫폼 SOFT family 대신 `acc.cgi`의 플랫폼 레코드를 받았다.
- 개발 자료 보충 RAW.tar는 series matrix에 발현 표가 있어 받지 않았다.
- EPIC BRef/TRef는 라이선스 때문에 받지 않아 포함률이 없다.
- GPL24324 포함률은 기호 열이 없어 비어 있다.
- 메모리 상한 48 GB를 `ulimit`으로 걸지 않았다. 측정 RSS는 §2에 있다.

## 12. 산출물 목록, Part별 커밋 해시

| Part | 해시 |
|---|---|
| A | fe821182063c122f3bd400adee6d8b020c8ef455 |
| B | ad3ca6fd6bdcf0fb3168b97ad61021b3d5ffabf9 |
| C | 6a159e44baf39bcae3ca40407c424adf104ba1e1 |
| D | 2d07594196e45a5c5493a405e51649b15f18f9e0 |
| E | 345ce1fa0b2ba0a5d3a9e6b64570a4b1c939409a |
| F | d00573e95fce5987cb0d155685df51f7eefac437 |

산출물: `config/paths.yaml`, `config/env_versions.txt`, `config/raw_checksums.tsv`, `scripts/01`–`07`, `data/processed/{GSE}_samples.tsv`, `data/processed/{GSE}_expr_gene_log2.tsv.gz`, `data/processed/{GSE}_probe_map.tsv`, `results/stage1/`의 요약 표와 `results/stage1/clocks/{GSE}_rnaage.tsv`. `data/`는 git에서 제외되어 있다.
