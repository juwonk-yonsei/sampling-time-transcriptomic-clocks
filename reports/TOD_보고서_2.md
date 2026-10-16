# TOD 보고서 2 — 판단 사항 처리, 시계 감사 보완, H1·H2 개발 탐색, RNA-seq 자료 탐색

## 1. 요약 (수치만, 12줄 이내)

- 양성 대조 정정 (a')(b): GSE39445 pass, GSE56931 pass. 재적합 없음. GSE48113 in-phase이면서 `clock_time_h`가 있는 샘플 147, 참고값 147, 일치 true. 제외 GSM1168837, GSM1168735.
- RA_Peters: 계수 유전자 289가 Peters 2015 `sign-REPL == y` 1497 기호의 부분집합 (289/289). pass. GPL24324 심볼이 있는 플랫폼 행 19483, 고유 심볼 20179.
- official–lin 개인 내 Pearson 최솟값 0.9977116 (GSE39445 RA_GTEx, n=221). 여섯 쌍 모두 ≥ 0.99. 전 분석은 lin. official은 점추정만.
- H1 lin, GSE39445: R2_within [구간] RA_GTEx 0.004809 [0.000539, 0.06557], RA_Peters 0.02373 [0.002608, 0.07939], RA_default 0.006802 [0.000577, 0.07044]. 2A/SD_between [구간] 0.2427 [0.08203, 0.9792], 0.2592 [0.07733, 0.5683], 0.2116 [0.07081, 0.7120]. 정점 시각 13.71, 6.845, 12.01. LRT p 0.6437, 0.08224, 0.5718. 무작위 R2 백분위 10.2, 28.4, 12.1.
- H1 lin, GSE56931: R2_within [구간] RA_GTEx 0.03663 [0.01200, 0.1101], RA_Peters 0.1455 [0.08426, 0.2474], RA_default 0.06588 [0.01266, 0.1682]. 2A/SD_between [구간] 0.6967 [0.3873, 1.520], 0.4515 [0.3236, 0.7274], 0.6346 [0.2553, 1.677]. 정점 시각 13.27, 4.542, 9.887. LRT p 0.09777, 0.0001051, 0.01874. 무작위 R2 백분위 52.9, 91, 68.4. lin 부트스트랩 수렴 실패 0/2000.
- H2 감소 비율 [구간], GSE39445: RA_GTEx 0.07214 [−5.117, 0.7870], RA_Peters 0.01077 [−1.878, 0.4447], RA_default −1.088 [−7.962, 0.7568]. GSE56931: RA_GTEx 0.2849 [−0.7247, 0.7906], RA_Peters 0.3118 [−0.3116, 0.7907], RA_default 0.5279 [−1.971, 0.8642]. 두 REML 모두 수렴 2000/2000.
- Part E 공개 후보 4 (GSE113883, GSE150337, GSE196728, GSE295758). 채혈 시각이 확인된 통제 접근 후보 0. GTEx v8 공개 표현형 열은 `SUBJID SEX AGE DTHHRDY`.
- 마지막 분석 커밋 `40f2746f9d7bf108e2f12c1daf097b4f2ccbead8`.

## 2. 환경과 자원 (바뀐 것, /usr/bin/time 기록)

단계 1 환경에 `clariomshumanhttranscriptcluster.db` 8.8.0 (Artistic-2.0)이 추가되어 있다. RNAAgeCalc 1.22.0, MCPcounter 1.2.0, lme4 2.0-6, circular 0.5-2. GPU를 쓰지 않았다. `ulimit`을 두지 않았다. OMP/MKL 스레드 상한을 두지 않았다.

| 실행 | elapsed | max RSS (kB) | exit |
| --- | --- | --- | --- |
| Part B GPL 대응 | 0:13.36 | 1039400 | 0 |
| Part C, lin 부트스트랩 재개 후 무작위 서명에서 중단 | 3:27.52 | 1428360 | 1 |
| Part C, 무작위 서명까지 마치고 S2에서 중단 | 5:17.42 | 1533504 | 1 |
| Part C, 기존 lin 부트스트랩·무작위 서명을 재사용해 민감도와 그림까지 | 0:21.55 | 1533444 | 0 |
| Part D, 저장된 결과 | 20:24.99 | 1647976 | 0 |

최대 RSS는 40 GB 아래였다.

## 3. 결정 기록 (config/decisions.tsv 전체)

| date | id | content | reason | timing |
| --- | --- | --- | --- | --- |
| 2026-10-16 | 1 | GSE56931 condition stays the sample field biological group. Sensitivity S4 uses only the design-text baseline times: day 1 hours 8, 12, 16, 20 and day 2 hours 0, 4. | Instruction 2 §1.1 item 1. | prior |
| 2026-10-16 | 2 | GSE122725 wake_ref_h = relclocktime, source given. clock_time_h and hours_awake stay NA. The time axis for this dataset is set at preregistration. | relclocktime 141/141, min 0.2, median 20.3, max 42.3. `results/stage2/GSE122725_wake_ref.tsv`. `data/processed`는 수정하지 않음. | prior |
| 2026-10-16 | 3 | GSE122725 light_group stays Control to dim. | Instruction 2 §1.1 item 3. | prior |
| 2026-10-16 | 4 | GSE39445 period stays extended_wakefulness. | Instruction 2 §1.1 item 4. | prior |
| 2026-10-16 | 5 | GSE48113 elapsed_h rule stays the forward clock difference within subject_id and period. | Instruction 2 §1.1 item 5. | prior |
| 2026-10-16 | 6 | GPL family SOFT stays not downloaded. | Instruction 2 §1.1 item 6. | prior |
| 2026-10-16 | 7 | GPL24324 symbols are mapped by instruction 2 §2.2. | Instruction 2 §1.1 item 7. | prior |
| 2026-10-16 | 8 | RA_Peters pass or fail follows the Part B audit in instruction 2 §2.1. | Instruction 2 §1.1 item 8. | prior |
| 2026-10-16 | 9 | EPIC is not used. | Instruction 2 §1.1 item 9. | prior |
| 2026-10-16 | 10 | ulimit is not set. Heavy steps record maximum RSS. Stop and report if RSS exceeds 40 GB. | Instruction 2 §1.1 item 10. | prior |
| 2026-10-16 | 11 | The impute.knn record from report 1 stays. Whether it affects within-person analysis is checked in Part C §3.1. | Instruction 2 §1.1 item 11. | prior |
| 2026-10-16 | 12 | Instruction 1 Part F, which was not run, is run in this stage as Part C after the §1.2 criterion change. | Instruction 2 §1.1 item 12. | posthoc |
| 2026-10-16 | 1.2 | Corrected positive-control criterion (a'): PER1 and PER3 peak hours are inside 04-12. Criterion (b) is unchanged. Both development datasets are pass. No refit. | 사후 정정. GSE39445 PER1 6.90324597433964, PER3 6.08583286807854, PER3–ARNTL 원형 거리 8.77605464525298. GSE56931 PER1 5.79155389251083, PER3 4.8945603273957, 거리 9.86874754997058. 잠금 자료 기준은 사전등록에 남김. | posthoc |
| 2026-10-16 | 1.3 | Confirmation analysis will drop GSM1168837 and GSM1168735. | 발현값을 읽기 전에 메타데이터로 정함. GSM1168837의 time sample taken은 문자 NA. GSM1168735는 D1, sleep timing은 Out of phase with respect to melatonin. in-phase이면서 시계 시각이 있는 샘플 147. 참고값 147. Match true. | prior |
| 2026-10-16 | 2.1 | RA_Peters passes to the next stage. | Ren and Kuan 2020이 후보 집합에 elastic net을 맞추고 계수를 저장. 패키지에 별도 후보 유전자 객체 없음. Peters 2015 Supplementary Data 1, sheet `Data.1`, `NEW-Gene-ID`, `sign-REPL == y` 1497. 계수 289/289가 부분집합. 저장 계수 중 0은 0개. | prior |

## 4. 보완

### 4.1 RA_Peters 감사

Ren and Kuan, PLoS ONE, PMC7402472, PMID 32750074. 2026-10-16에 Europe PMC 전문 XML에서 발췌했다.

> For each tissue, signature and sample type, we pre-trained the calculator using elastic net based on the GTEx samples. In the within-tissue age prediction section, we have demonstrated that elastic net model outperforms the other prediction algorithms. We saved the pre-trained model coefficients as internal data in the software.

> In all these models, chronological age was the response variable whereas the logarithm transformed FPKM were the predictors. The samples were first randomly split into 50%-50% training-test set, where the prediction algorithms were constructed on the training data and evaluated on the test data. The optimal parameters, namely alpha and lambda in elastic net, number of trees in GBM, cost and gamma in SVR, and bin size in ensemble LDA were selected by 10 fold cross validation in the training set. This 50%-50% training-test split and model evaluation were repeated 10 times.

> (iv) The 1,497 age-related genes of [18] (denoted by Peters).

> (vi) The 1,616 common age-related genes discussed in the meta-analysis section (denoted by GTExAge).

후보 집합 문단의 발췌: "differentially expressed genes ranked by the p-values from differential expression analysis." 이어서 "To reduce the computation cost, we considered the top 1,000 genes."

RNAAgeCalc 네임스페이스에 Peters, GTExAge, DESeq2의 별도 후보 유전자 객체는 없다 (`separate_candidate_object_in_package` False). 저장 계수 유전자 수(절편 제외)는 Peters 289, GTExAge 312, DESeq2 279. 세 혈액 서명에서 계수가 0인 유전자는 0개다.

후보 목록은 Peters et al. 2015 Nature Communications Supplementary Data 1, `ncomms9570-s2.xlsx`. 바이트 1757989. SHA256 `a1955742f15506cdf86153f50fd6d5ebafaa13d1a3c1e5abe862a8b09457cc6e`. 라이선스 CC BY 4.0. 계수로 쓰지 않음 (`used_as_coefficients` False). Sheet `Data.1`, 헤더는 Excel 3행, 심볼 열 `NEW-Gene-ID`, 필터 `sign-REPL == y`. 1497행, 고유 심볼 1497. 계수 유전자 289가 모두 들어 있다. `pass_to_next_stage` True.

GTExAge와 DESeq2는 패키지 안에 후보 목록 객체가 없어 같은 부분집합 검사를 하지 않았다. Peters 파일을 그 두 서명의 후보 목록으로 쓰지 않았다.

### 4.2 GPL24324

플랫폼 원문:

- `!Platform_title` = [Clariom_S_Human_HT] Affymetrix Clariom S Assay HT, Human (Includes Pico Assay)
- `!Platform_technology` = in situ oligonucleotide
- `!Platform_manufacturer` = Affymetrix
- 설명에 Clariom S Assay HT (902969), Clariom S Pico Assay HT (902964), `#%genome-version=hg38`, `#%netaffx-annotation-date=2016-06-03`, `#%netaffx-annotation-netaffx-build=36`

대응: `clariomshumanhttranscriptcluster.db` 8.8.0, `AnnotationDbi::select`, PROBEID → SYMBOL. 플랫폼 행 27189. DB 키 24351. 키가 DB에 있는 행 21448. 심볼이 있는 행 19483. 심볼이 있는 고유 probeset 19483. 대응 행 22903. 심볼이 있는 대응 행 20938. 고유 심볼 20179. `select()`는 1:many를 반환했다.

시계·MCPcounter 포함률 (심볼이 플랫폼 대응 심볼과 같은 유전자 수 / 참조 유전자 수):

| reference | n_reference | n_on_platform | fraction |
| --- | --- | --- | --- |
| RA_GTEx | 312 | 255 | 0.8173076923076923 |
| RA_Peters | 289 | 276 | 0.9550173010380623 |
| RA_default | 279 | 240 | 0.8602150537634409 |
| MCPcounter T cells | 16 | 14 | 0.875 |
| MCPcounter CD8 T cells | 1 | 1 | 1 |
| MCPcounter Cytotoxic lymphocytes | 7 | 7 | 1 |
| MCPcounter B lineage | 9 | 8 | 0.8888888888888888 |
| MCPcounter NK cells | 9 | 9 | 1 |
| MCPcounter Monocytic lineage | 7 | 7 | 1 |
| MCPcounter Myeloid dendritic cells | 6 | 5 | 0.8333333333333334 |
| MCPcounter Neutrophils | 15 | 14 | 0.9333333333333333 |
| MCPcounter Endothelial cells | 33 | 31 | 0.9393939393939394 |
| MCPcounter Fibroblasts | 8 | 8 | 1 |
| MCPcounter ALL | 111 | 104 | 0.9369369369369369 |

### 4.3 GSE39445 circadianphase

`circadianphase` 분포. 결측 0.

| subset | n | min | median | max | n_unique |
| --- | --- | --- | --- | --- | --- |
| all | 438 | −15 | 0 | 15 | 11 |
| control | 221 | −12 | 0 | 15 | 10 |

예시. 값은 `clock_time_h`, `hours_awake`, `circadian_ref_h` 순이다. AF0004는 샘플 10개, `hours_awake` 7.5–34.5, `circadian_ref_h` −12–15, 간격 3. AF0010은 샘플 9개이고 `hours_awake` 25.5, `circadian_ref_h` 6인 행이 없다. 단위와 0점 정의 문장은 이 단계 산출 파일에 없다.

### 4.4 연구자 모듈 시계

`config/paths.yaml`의 `module_clock_dir`이 비어 있다. 상태 `not_provided`. 모듈 시계는 계산하지 않았다.

## 5. H1 개발 탐색 (확인 검정 아님)

개발 자료의 control 샘플, `clock_time_h`가 있는 행. GSE39445 참가자 24, 샘플 221. GSE56931 참가자 14, 샘플 130. 시계는 RA_GTEx, RA_Peters, RA_default. 시드 20261015. 부트스트랩 추출은 자료마다 시드를 다시 20261015로 두고 참가자 2000회 복원추출. 같은 참가자가 여러 번 뽑히면 `subject#occurrence`를 붙였다. 추출 행 수는 GSE39445 48000, GSE56931 28000.

주 모형은 `score ~ cos(2*pi*t/24) + sin(2*pi*t/24) + (1|subject_id)`, t=`clock_time_h`. 추정은 REML. LRT(df=2)와 `R2_within`은 ML. `R2_within` = 1 − (주 모형 ML 잔차분산 / 귀무 모형 ML 잔차분산). 귀무 모형은 `score ~ 1 + (1|subject_id)`. `SD_between`은 참가자 평균 점수의 SD. `SD_resid`는 귀무 모형 잔차 SD. 2A = 2 × sqrt(β_cos² + β_sin²).

저장 lin 점수와 계수×발현으로 다시 만든 lin 점수의 최대 절대차와 상관:

| dataset | clock | n_genes_used | n_compared | max_abs_diff | cor |
| --- | --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | 229 | 438 | 4.547474e-13 | 1 |
| GSE39445 | RA_Peters | 280 | 438 | 4.831691e-13 | 1 |
| GSE39445 | RA_default | 214 | 438 | 4.831691e-13 | 1 |
| GSE56931 | RA_GTEx | 255 | 249 | 4.973799e-13 | 1 |
| GSE56931 | RA_Peters | 278 | 249 | 4.973799e-13 | 1 |
| GSE56931 | RA_default | 235 | 249 | 4.547474e-13 | 1 |

### 5.1 official–lin 일치

control, 참가자 평균을 뺀 뒤의 Pearson. `sd_within_official_minus_lin`은 (official−lin)을 참가자 안에서 평균을 뺀 값의 SD.

| dataset | clock | pearson_within | sd_within_official_minus_lin | n | use_lin_only |
| --- | --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | 0.9977116 | 2.6850734 | 221 | TRUE |
| GSE39445 | RA_Peters | 0.9987913 | 0.5000050 | 221 | TRUE |
| GSE39445 | RA_default | 0.9989172 | 0.5248264 | 221 | TRUE |
| GSE56931 | RA_GTEx | 0.9992513 | 0.8235418 | 130 | TRUE |
| GSE56931 | RA_Peters | 0.9994506 | 0.2155837 | 130 | TRUE |
| GSE56931 | RA_default | 0.9997851 | 0.1839905 | 130 | TRUE |

최솟값 0.9977116. 이후 부트스트랩, 무작위 서명, 민감도, 그림은 lin만 했다. official은 아래 점추정이다. 이 비교가 단계 1의 impute.knn 기록에 대한 §3.1 확인이다. impute.knn만 켜고 끈 별도 적합은 하지 않았다.

### 5.2 주 모형 표

lin의 구간은 부트스트랩 2.5·97.5 백분위 (type 7)다. 정점 구간은 부트스트랩 정점의 원형 평균 `peak_center`와, 그 평균까지 원형 거리(0–12시간)의 95백분위 `halfwidth_95`다. 수렴 실패는 REML, ML, 귀무 ML 중 하나라도 optimizer 코드가 0이 아닌 반복이다. lin은 여섯 시계 모두 0/2000. 점추정의 REML·ML·귀무 수렴은 모두 TRUE, REML singular는 FALSE.

GSE39445, n_subjects 24, n_samples 221.

| clock | method | β_cos | β_sin | 2A [구간] | R2_within [구간] | 2A/SD_between [구간] | 2A/SD_resid [구간] | peak_h (center, half) | LRT p | ICC |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | lin | −2.05534974435911 | −0.989438315904824 | 4.5622147045564 [1.61285125755799, 17.0458191350799] | 0.00480935569795549 [0.000539049792313514, 0.0655700363463853] | 0.242731389057439 [0.0820296542735993, 0.97918811506639] | 0.181791314709904 [0.0657195689371128, 0.681731950192155] | 13.7137336093391 (13.6618778047924, 9.4818705231949) | 0.643663244864414 | 0.283742283864298 |
| RA_GTEx | official | −2.15376105086438 | −0.78822420854626 | 4.58693102941774 [미계산] | 0.0058765595480359 [미계산] | 0.259811194937003 [미계산] | 0.201578188473544 [미계산] | 13.3400931979681 (미계산) | 0.587336245948952 | 0.303366536519521 |
| RA_Peters | lin | −0.424466616594337 | 1.8872655640872 | 3.86882060477993 [1.29322695284716, 7.38317071538597] | 0.0237318052532715 [0.00260828058054079, 0.0793917923749646] | 0.25915869900777 [0.077332978314821, 0.568255471420682] | 0.421217882293791 [0.146297295486235, 0.766168132741987] | 6.84503499428439 (6.70958325567532, 3.82180804350612) | 0.0822398090507864 | 0.70770638546879 |
| RA_Peters | official | −0.593114337131461 | 1.73507042482888 | 3.66729000545458 [미계산] | 0.0222427842780705 [미계산] | 0.25417386676187 [미계산] | 0.412097502834738 [미계산] | 7.25816113482645 (미계산) | 0.0955997779560659 | 0.706729126991898 |
| RA_default | lin | −0.992493083273389 | −0.00167525410351653 | 1.98498899424841 [0.654469343295716, 6.17298524100273] | 0.00680181476423569 [0.000576531365170474, 0.0704405878123816] | 0.211605758610572 [0.0708138863374765, 0.712026135834721] | 0.212875508690383 [0.0687469965019418, 0.736252745388078] | 12.0064473932836 (12.458465246635, 8.83488246830996) | 0.571779550797501 | 0.460120713994768 |
| RA_default | official | −1.11989490384088 | −0.00733481138992836 | 2.23983784690491 [미계산] | 0.00909814998077529 [미계산] | 0.242883169673744 [미계산] | 0.249736075218239 [미계산] | 12.0250170931648 (미계산) | 0.462314014136504 | 0.472362859761108 |

GSE56931, n_subjects 14, n_samples 130.

| clock | method | β_cos | β_sin | 2A [구간] | R2_within [구간] | 2A/SD_between [구간] | 2A/SD_resid [구간] | peak_h (center, half) | LRT p | ICC |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | lin | −3.97713729889457 | −1.36997692582768 | 8.41295616809185 [4.90563166315199, 13.602552044309] | 0.0366339113266437 [0.0119981991170903, 0.110143351985373] | 0.696744201903398 [0.387281124642329, 1.52027479452532] | 0.587062420973302 [0.335342768582716, 0.987805856017995] | 13.2671292133393 (13.3498596199379, 2.75243878714261) | 0.0977660921594727 | 0.352751240612617 |
| RA_GTEx | official | −3.71231792732727 | −1.18718435171451 | 7.79505255364227 [미계산] | 0.0345963072469041 [미계산] | 0.669765726077446 [미계산] | 0.571000567100919 [미계산] | 13.1822711571327 (미계산) | 0.110221959787183 | 0.358889520151384 |
| RA_Peters | lin | 1.1883334284853 | 2.96038500225555 | 6.37997357324785 [4.70090081289685, 8.77697299959059] | 0.145491804598675 [0.0842595575738163, 0.247437380369278] | 0.451478484596114 [0.323600524627966, 0.727421182415766] | 1.06466896629158 [0.811769161533672, 1.39129983014162] | 4.54192615314062 (4.53020538874562, 1.54219887476787) | 0.000105116530090887 | 0.834952575430877 |
| RA_Peters | official | 1.12012276972188 | 2.92609770591761 | 6.26632996380679 [미계산] | 0.146619946979854 [미계산] | 0.447911372278779 [미계산] | 1.06596730168653 [미계산] | 4.60352929557141 (미계산) | 9.72252317266578e-05 | 0.837517062704639 |
| RA_default | lin | −1.96418355918904 | 1.2128554452917 | 4.61694071236019 [1.97821187537065, 7.50896822888791] | 0.0658778147660831 [0.012661998886085, 0.168222912899573] | 0.634553234473014 [0.255304707663938, 1.67690719040041] | 0.673018542815358 [0.296945118590377, 1.07695973961528] | 9.88701727265479 (9.84702408091482, 2.49115920551054) | 0.0187399469939356 | 0.482693640222497 |
| RA_default | official | −1.96175182407929 | 1.17317371714593 | 4.57156725418337 [미계산] | 0.0669122502617618 [미계산] | 0.636172951438872 [미계산] | 0.679629190605548 [미계산] | 9.9413056001098 (미계산) | 0.0175416896990936 | 0.486736890127431 |

### 5.3 무작위 서명 비교

자료의 모든 유전자를 전체 샘플 평균 log2의 10분위로 나눴다. lin에 쓴 시계 유전자마다 같은 분위에서 시계 유전자가 아닌 유전자를 서명 안에서 중복 없이 뽑고, 계수를 섞어 붙였다. 절편은 0. 서명 1000개. 자료×시계마다 시드를 20261015로 다시 두었다. 맞춘 것은 주 모형의 `R2_within`뿐이다. 수렴하지 않은 서명은 0. 백분위 순위는 100 × (무작위 R2_within ≤ 시계 R2_within인 비율)이다. 5·95 백분위는 type 7.

| dataset | clock | median | p05 | p95 | clock R2_within | percentile_rank |
| --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | 0.0339428421885641 | 0.00266884265317641 | 0.160781517816631 | 0.00480935569795549 | 10.2 |
| GSE39445 | RA_Peters | 0.0502220062254742 | 0.00299734998119304 | 0.200087445435632 | 0.0237318052532715 | 28.4 |
| GSE39445 | RA_default | 0.0400305655760273 | 0.00236622094032111 | 0.181963130629645 | 0.00680181476423569 | 12.1 |
| GSE56931 | RA_GTEx | 0.0338608602551348 | 0.00208964864979564 | 0.145419854558692 | 0.0366339113266437 | 52.9 |
| GSE56931 | RA_Peters | 0.0409515902957628 | 0.0026331742693649 | 0.177083379503933 | 0.145491804598675 | 91 |
| GSE56931 | RA_default | 0.0402417339009208 | 0.00284461987505798 | 0.163829730915396 | 0.0658778147660831 | 68.4 |

### 5.4 민감도 S1–S5

부트스트랩은 하지 않았다. S1은 GSE39445만. 모형은 `score ~ cos_t + sin_t + hours_awake + (1|subject_id)`, 귀무는 `score ~ hours_awake + (1|subject_id)`. S2 귀무는 `score ~ 1 + (1 + cos_t + sin_t | subject_id)`. S2의 REML·ML·귀무는 모두 수렴 TRUE. S2의 `icc_null`은 비어 있다. S3은 참가자별 OLS다. 샘플 ≥ 5이고 `clock_time_h`의 max−min ≥ 16인 참가자만 썼고, 두 자료 모두 참가자 전원(24, 14)이 들어갔다. S3의 2A는 참가자별 2A의 중앙값, 괄호는 25·75 백분위. S3의 R2는 참가자별 OLS R²의 중앙값. S3의 정점은 참가자 정점의 원형 평균. S3 LRT는 맞추지 않았다. Rayleigh는 `circular::rayleigh.test`. S4는 GSE56931, day 1의 hour 8·12·16·20과 day 2의 hour 0·4. n_samples 79. S5는 단계 1의 3 MAD 표시 샘플을 control 안에서 뺀 것이다. GSE39445는 10개(221→211), GSE56931은 11개(130→119).

| dataset | clock | S | 2A | R2 | peak_h | LRT p 또는 Rayleigh p | n_samples |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | S1 | 2.56231423377487 | 0.00145101292931349 | 17.1109342722687 | 0.856738467893009 | 221 |
| GSE39445 | RA_Peters | S1 | 4.27852694464597 | 0.0275038489537423 | 7.60739182142201 | 0.0560156218723147 | 221 |
| GSE39445 | RA_default | S1 | 0.310501104206593 | 0.000218879567791874 | 16.0415515311994 | 0.983899575769617 | 221 |
| GSE39445 | RA_GTEx | S2 | 4.78509519871122 | 0.00323141549948947 | 14.0417004743891 | 0.666862993019604 | 221 |
| GSE39445 | RA_Peters | S2 | 3.76595125782325 | 0.00639982326133193 | 6.79457127152056 | 0.122275997221478 | 221 |
| GSE39445 | RA_default | S2 | 2.36899011751258 | 0.00707273101386641 | 11.3586752893535 | 0.540052558400389 | 221 |
| GSE56931 | RA_GTEx | S2 | 8.59655831481556 | 0.0337005309565047 | 13.2257701112999 | 0.0996715564148518 | 130 |
| GSE56931 | RA_Peters | S2 | 6.38276844320766 | 0.0560173122341118 | 4.48047059001304 | 0.00045899513695964 | 130 |
| GSE56931 | RA_default | S2 | 4.57295778630413 | 0.0337151457326444 | 9.84275414039334 | 0.032216706449245 | 130 |
| GSE39445 | RA_GTEx | S3 | 35.9481274122815 (17.9339875705045, 44.6496710682609) | 0.330448650292978 | 16.1762755006608 | Rayleigh stat 0.185630469844363, p 0.441858759090448 | 221 |
| GSE39445 | RA_Peters | S3 | 9.85664897967244 (6.36305874766871, 14.8120221426277) | 0.194487402242135 | 7.6596725556989 | Rayleigh stat 0.315734702121352, p 0.0904784901099544 | 221 |
| GSE39445 | RA_default | S3 | 11.8116924403024 (10.1320357110958, 13.6075347981965) | 0.274437679986311 | 11.33679242508 | Rayleigh stat 0.26842729345024, p 0.178321317225931 | 221 |
| GSE56931 | RA_GTEx | S3 | 14.3668282235998 (12.8548751407631, 17.4278252921441) | 0.17084147610703 | 13.9237054661072 | Rayleigh stat 0.485055786552104, p 0.0340847922157505 | 130 |
| GSE56931 | RA_Peters | S3 | 7.37981802972622 (5.52954857275055, 10.3474686335236) | 0.256810192164129 | 4.75791729818394 | Rayleigh stat 0.707202273146455, p 0.000370716511355789 | 130 |
| GSE56931 | RA_default | S3 | 9.70694408907529 (5.63209849917956, 10.5950000204829) | 0.225275498622277 | 9.86259080906412 | Rayleigh stat 0.507379775727391, p 0.0242048907430237 | 130 |
| GSE56931 | RA_GTEx | S4 | 9.45144315318197 | 0.0657576097470823 | 15.2079563224418 | 0.102054116182808 | 79 |
| GSE56931 | RA_Peters | S4 | 5.91454628100006 | 0.131459803409985 | 4.58487114312972 | 0.0100308192594568 | 79 |
| GSE56931 | RA_default | S4 | 3.71114318233623 | 0.0422030998705777 | 10.6271527650614 | 0.252718464718802 | 79 |
| GSE39445 | RA_GTEx | S5 | 4.50024814829927 | 0.00540194907635438 | 14.2207579983526 | 0.661882207283519 | 211 |
| GSE39445 | RA_Peters | S5 | 4.00783818201034 | 0.024870733263528 | 7.85667082971997 | 0.0832374554816746 | 211 |
| GSE39445 | RA_default | S5 | 1.67814938667061 | 0.00499023120235265 | 12.0388732864785 | 0.658727506723165 | 211 |
| GSE56931 | RA_GTEx | S5 | 9.23730638196741 | 0.047290471211853 | 13.3799729690707 | 0.0708000668072659 | 119 |
| GSE56931 | RA_Peters | S5 | 5.9008770683746 | 0.135121794172807 | 4.86497445910234 | 0.000434807958165243 | 119 |
| GSE56931 | RA_default | S5 | 5.67804320347146 | 0.0960484657022619 | 10.2477341205723 | 0.00483653414936823 | 119 |

### 5.5 그림 경로

lin만. x=`clock_time_h`, y=score, 참가자별 선, 주 모형 곡선.

- `results/stage2/figures/GSE39445_RA_GTEx_lin.png`
- `results/stage2/figures/GSE39445_RA_Peters_lin.png`
- `results/stage2/figures/GSE39445_RA_default_lin.png`
- `results/stage2/figures/GSE56931_RA_GTEx_lin.png`
- `results/stage2/figures/GSE56931_RA_Peters_lin.png`
- `results/stage2/figures/GSE56931_RA_default_lin.png`

## 6. H2 개발 탐색

MCPcounter 1.2.0. 입력은 개발 자료 전체 샘플의 유전자 log2. `featuresType = "HUGO_symbols"`. `genes`에는 `data/processed/ref/MCPcounter_genes.txt`를 넘겼다. 기본 인자의 GitHub `curl`은 호출하지 않았다. 모형에 넣은 8개: T cells, CD8 T cells, Cytotoxic lymphocytes, NK cells, B lineage, Monocytic lineage, Myeloid dendritic cells, Neutrophils. Endothelial cells와 Fibroblasts는 점수만 냈다.

행렬에 있는 표지 수 / 표지 수: GSE39445는 T 13/16, CD8 1/1, Cytotoxic 7/7, B 8/9, NK 6/9, Monocytic 7/7, Myeloid DC 5/6, Neutrophils 15/15, Endothelial 30/33, Fibroblasts 8/8. GSE56931은 T 14/16, CD8 1/1, Cytotoxic 7/7, B 9/9, NK 8/9, Monocytic 7/7, Myeloid DC 6/6, Neutrophils 15/15, Endothelial 30/33, Fibroblasts 8/8.

control이고 `clock_time_h`가 있는 샘플의 점수 평균과 SD. GSE39445 n_all 438, n_control 221. GSE56931 n_all 249, n_control 130.

| dataset | population | mean_control | sd_control |
| --- | --- | --- | --- |
| GSE39445 | T cells | 10.7703572465493 | 0.338250435495698 |
| GSE39445 | CD8 T cells | 12.5471287161086 | 0.514580217674011 |
| GSE39445 | Cytotoxic lymphocytes | 10.5558625229341 | 0.361890790582009 |
| GSE39445 | NK cells | 8.64335037507994 | 0.461844365063809 |
| GSE39445 | B lineage | 10.0579170526069 | 0.418672249404429 |
| GSE39445 | Monocytic lineage | 9.95051765697156 | 0.2743225482616 |
| GSE39445 | Myeloid dendritic cells | 6.39663518889683 | 0.347460088141755 |
| GSE39445 | Neutrophils | 11.3750405327816 | 0.410796476945196 |
| GSE39445 | Endothelial cells | 5.52629470262006 | 0.261349197851303 |
| GSE39445 | Fibroblasts | 4.55429039557975 | 0.213840007884758 |
| GSE56931 | T cells | 9.81700225879121 | 0.268972762513634 |
| GSE56931 | CD8 T cells | 10.1172517461538 | 0.475930174899725 |
| GSE56931 | Cytotoxic lymphocytes | 9.73619642417582 | 0.332319359492643 |
| GSE56931 | NK cells | 7.21695132115385 | 0.353289001732924 |
| GSE56931 | B lineage | 9.06274351025641 | 0.461830049513812 |
| GSE56931 | Monocytic lineage | 8.97298605934066 | 0.259345943602736 |
| GSE56931 | Myeloid dendritic cells | 6.22589409102564 | 0.188025887783432 |
| GSE56931 | Neutrophils | 11.3109632133333 | 0.254688693061023 |
| GSE56931 | Endothelial cells | 5.1349375074359 | 0.137947273068097 |
| GSE56931 | Fibroblasts | 5.36127001730769 | 0.144383608749363 |

표지와 시계 계수 유전자의 겹침이 0보다 큰 칸. 나머지 집단의 겹침은 0이다. 고유 심볼 겹침은 RA_GTEx 3 (`DCN`, `FLT3LG`, `KLRD1`), RA_Peters 4 (`ADAP2`, `CD8B`, `RASSF4`, `TRAT1`), RA_default 4 (`CD8A`, `COL6A2`, `FLT3LG`, `PBX4`).

| clock | population | n_marker | n_clock_coef | n_overlap |
| --- | --- | --- | --- | --- |
| RA_GTEx | T cells | 16 | 312 | 1 |
| RA_GTEx | Cytotoxic lymphocytes | 7 | 312 | 1 |
| RA_GTEx | Fibroblasts | 8 | 312 | 1 |
| RA_Peters | T cells | 16 | 289 | 1 |
| RA_Peters | CD8 T cells | 1 | 289 | 1 |
| RA_Peters | Monocytic lineage | 7 | 289 | 2 |
| RA_default | T cells | 16 | 279 | 2 |
| RA_default | Cytotoxic lymphocytes | 7 | 279 | 1 |
| RA_default | Fibroblasts | 8 | 279 | 1 |

8개 점수 쌍 56개(자료 2 × C(8,2)=28)에서 \|r\| > 0.9인 쌍은 0개다. 표본은 control이고 `clock_time_h`가 있는 행이다. 절대 상관이 가장 큰 쌍은 GSE39445에서 T cells와 NK cells, r = −0.540411711711575. GSE56931에서 T cells와 Neutrophils, r = −0.739963619808647. 모형에서 빼지 않았다.

### 6.2 구성 점수의 하루 리듬

control, 주 모형, 부트스트랩 없음. 8개 모두 REML·ML·귀무 수렴 TRUE.

| dataset | population | n | 2A | R2_within | peak_h | LRT p |
| --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | T cells | 221 | 0.370479518051608 | 0.364969539858517 | 2.01178453285395 | 4.67850977349069e-20 |
| GSE39445 | CD8 T cells | 221 | 0.456416154866187 | 0.327102559880078 | 1.46898840392452 | 8.86356306008075e-18 |
| GSE39445 | Cytotoxic lymphocytes | 221 | 0.231865236884549 | 0.12167128293316 | 13.4282218619447 | 1.92009443964674e-06 |
| GSE39445 | NK cells | 221 | 0.429616717109446 | 0.351259667157474 | 13.6363947099388 | 3.16805722470307e-19 |
| GSE39445 | B lineage | 221 | 0.0796877116596262 | 0.0233435735062996 | 0.54355485981716 | 0.0944215101801027 |
| GSE39445 | Monocytic lineage | 221 | 0.0943620945088407 | 0.0946686061628258 | 1.48106305219427 | 5.45963753264426e-05 |
| GSE39445 | Myeloid dendritic cells | 221 | 0.125350227546181 | 0.0338011414266206 | 23.7180924270176 | 0.0376237493997112 |
| GSE39445 | Neutrophils | 221 | 0.491398502744698 | 0.324558207564043 | 14.7133857365347 | 1.72783871707044e-17 |
| GSE56931 | T cells | 130 | 0.252352434185947 | 0.29961283868262 | 0.268450187543044 | 1.022722552105e-09 |
| GSE56931 | CD8 T cells | 130 | 0.365374367890902 | 0.302506681724577 | 0.821668717174426 | 7.44808106876266e-10 |
| GSE56931 | Cytotoxic lymphocytes | 130 | 0.274071086486271 | 0.25427993353199 | 10.8882969919055 | 3.92173321072183e-08 |
| GSE56931 | NK cells | 130 | 0.508119137506897 | 0.350078827602487 | 10.9420393769602 | 1.21575312140788e-11 |
| GSE56931 | B lineage | 130 | 0.196618375878881 | 0.228477826018004 | 22.5809101432575 | 2.98864519498162e-07 |
| GSE56931 | Monocytic lineage | 130 | 0.115108886651418 | 0.0907587888915912 | 13.2327615839502 | 0.00388916842751784 |
| GSE56931 | Myeloid dendritic cells | 130 | 0.0676595364044617 | 0.0464710274422268 | 17.8752932282879 | 0.059720272244732 |
| GSE56931 | Neutrophils | 130 | 0.184671515886958 | 0.227683731696244 | 14.2483972010014 | 3.37004394625701e-07 |

### 6.3 구성 보정 전후

보정 모형은 주 모형에 8개 점수를 고정효과로 더한 것이다. 보정 `R2_within`의 귀무는 8개 점수 + `(1|subject_id)`. 감소 비율 = 1 − 2A(보정) / 2A(보정 전). 구간은 Part C와 같은 추출이고, 반복마다 두 모형을 같은 표본에 맞췄다. 두 REML이 모두 수렴한 반복 2000/2000. 보정 전 2A는 §5.2 lin 2A와 같다.

| dataset | clock | 2A 전 | R2 전 | peak 전 | LRT p 전 | 2A 후 | R2 후 | peak 후 | LRT p 후 | 감소 비율 [구간] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | 4.5622147045564 | 0.00480935569795549 | 13.7137336093391 | 0.643663244864414 | 4.23307854134853 | 0.00395891797715864 | 22.0201400359367 | 0.741592663207429 | 0.0721439442293563 [−5.11654249519675, 0.787013921900512] |
| GSE39445 | RA_Peters | 3.86882060477993 | 0.0237318052532715 | 6.84503499428439 | 0.0822398090507864 | 3.82713657363215 | 0.029137558801865 | 8.65187574711682 | 0.05142521684644 | 0.0107743509989272 [−1.87801689952032, 0.444715202142757] |
| GSE39445 | RA_default | 1.98498899424841 | 0.00680181476423569 | 12.0064473932836 | 0.571779550797501 | 4.14551106787721 | 0.0252641486638288 | 0.897600707548435 | 0.164349430438606 | −1.08843025320997 [−7.9618711687648, 0.756788584700328] |
| GSE56931 | RA_GTEx | 8.41295616809185 | 0.0366339113266437 | 13.2671292133393 | 0.0977660921594727 | 6.01595169382857 | 0.0217031292741995 | 15.598217470026 | 0.490547103538976 | 0.284918217374589 [−0.724723456894724, 0.790564092846028] |
| GSE56931 | RA_Peters | 6.37997357324785 | 0.145491804598675 | 4.54192615314062 | 0.000105116530090887 | 4.3905565243361 | 0.0444599963186707 | 4.35827181241872 | 0.0393202367232618 | 0.311822145667446 [−0.31164556124922, 0.790687529744977] |
| GSE56931 | RA_default | 4.61694071236019 | 0.0658778147660831 | 9.88701727265479 | 0.0187399469939356 | 2.17987790567819 | 0.00598733704263055 | 4.22144685948071 | 0.462793070404024 | 0.527852307082402 [−1.97091941605562, 0.864180910821805] |

보정 모형의 REML·ML·귀무 수렴은 여섯 행 모두 TRUE.

## 7. RNA-seq 자료 탐색

검색일 2026-10-16. 발현 행렬과 보충 파일은 받지 않았다. series matrix는 `!series_matrix_table_begin` 앞에서 멈췄다. 이미 쓰는 GSE39445, GSE56931, GSE48113, GSE122725는 후보에 넣지 않았다. 순위를 매기지 않았다.

### 7.1 검색 기록

| database | query | hit_count | note |
| --- | --- | --- | --- |
| GEO gds | human AND (blood OR PBMC OR leukocyte) AND ("RNA-seq" OR "RNA sequencing") AND (circadian OR "time of day" OR "clock time") | 82 | 제목 82건을 저장 |
| GEO gds | "whole blood" AND (RNA-seq OR "high throughput sequencing") AND (circadian OR "time of day" OR "time-of-day") AND human | 160 | |
| GEO gds | PBMC AND "RNA-seq" AND (circadian OR "time of day") AND human | 1 | |
| GEO gds | TimeSignature AND human | 1 | |
| GEO gds | BodyTime AND human | 0 | |
| PubMed | TimeSignature blood RNA-seq | 0 | |
| PubMed | TimeSignature circadian blood | 2 | 33031153; 30201705 |
| PubMed | Wittenbrink TimeSignature | 0 | |
| PubMed | BodyTime blood transcriptome circadian | 1 | 29953415 |
| PubMed | BodyTime circadian transcriptome | 1 | 29953415 |
| PubMed | "whole blood" RNA-seq "time of day" circadian | 0 | |
| PubMed | blood RNA-seq circadian older adults | 10 | 저장된 id_list는 5개: 42589234; 41902670; 41285310; 41112764; 39097147 |
| BioStudies | blood AND RNA-seq AND circadian | 62 | |
| BioStudies | "whole blood" AND RNA-seq AND "time of day" | 4 | E-MTAB-12608, E-MTAB-14799, E-MTAB-15199, E-MTAB-14708 |
| BioStudies | PBMC AND RNA-seq AND circadian | 0 | |
| BioStudies | leukocyte AND RNA-seq AND circadian | 19 | |
| BioStudies | BodyTime AND blood | 1 | 연구 JSON에 GSE 없음 |
| BioStudies | TimeSignature AND blood | 1 | JSON에 GSE113883, GSE39445, GSE48113, GSE56931 |
| ArrayExpress | blood RNA-seq circadian | | HTTP 200, 본문은 HTML |
| ArrayExpress | whole blood time of day | | HTTP 200, 본문은 HTML |
| EGA | metadata/v2 studies query blood RNA-seq circadian | | HTTP 404 |
| EGA | EgaService wsdl | | HTTP 404 |
| PubMed elink gds | PMID 30201705 | 1 | GSE113883 |
| PubMed elink gds | PMID 29953415 | 0 | GEO 연결 없음 |
| PubMed elink gds | PMID 31473510 | 1 | GSE127530 |
| PubMed elink gds | PMID 41902670 | 1 | GSE309217. esummary 제목은 골격근. 헤더를 읽지 않음 |

82건 제목 목록에서 혈액으로 보이는 계열과 단서 논문에 연결된 계열의 헤더를 읽었다. GSE292324 series matrix URL은 HTTP 404라 샘플 characteristics를 읽지 않았다.

BioStudies 4건의 연구 JSON에서 본 Time 예: E-MTAB-12608 `0 day` (Sus scrofa), E-MTAB-14799 `day 64 after treatment initiation day`, E-MTAB-15199 `7 day`, E-MTAB-14708 `D60 day`. 이 네 값은 벽시계 시각이 아니다.

### 7.2 공개 후보 표

조건: 사람 혈액(전혈, PBMC, 백혈구, 혈액 세포 분획)의 bulk RNA-Seq, 샘플 시각 필드가 있고, 한 사람에게서 시각이 2개 이상. 접근은 GEO 공개.

| 접근 번호 | 기술 | 표본 | 참가자 | 참가자당 샘플 | 시각 필드와 예 | 나이 | 조건 (헤더 문장) | 처리 형식 | 논문 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE113883 | bulk RNA-Seq, GPL18573 | Whole blood | 11 (`subjectid`) | 13개 3명, 14개 6명, 15개 2명. 샘플 153 | `hour` 1,3,5,…,29 (15개 값). `local_time` 예 8:00(22), 10:00(20), 12:00(15). `salivary_dlmo_circadian_time` NA 58/153 | characteristics에 age 없음 | 중간 일주기 표현형 11명. 28시간에 2시간마다, 설계 문장상 165개 중 153개가 QA 통과 | TPM. Genome_build Gencode22 GRCh38 | 30201705 |
| GSE150337 | bulk RNA-Seq, GPL23227 | whole blood | 10 (제목의 MS01–MS09, HC01–HC03. 각 2샘플) | 2 | characteristic `sampling time`: daytime 10, nighttime 10. 제목 시각: 21_00×7, 14_00×3, dt_1×3, nt_1×3, 15_00×2, 14_30×1, 15_30×1. 설계 문장: 2 pm과 9 pm | age 없음 | RRMS 7명과 대조 3명. 하루 2회 | raw counts. GRCh38 | 34746708 |
| GSE196728 | bulk RNA-Seq, GPL24676 | Whole blood | 8 (`subject` A–H) | 18 | `time` 2, 6, 10, 14, 18, 22 각 24. `altitude` PU, RI, SL 각 48 | age 없음 | 4시간 간격, 24시간, 고도 3곳 (SL, PU 3800 m, RI 5100 m). 설계 문장의 지명은 Granoble | raw counts. Genome_build 문장은 mm10. `Sample_organism_ch1`은 Homo sapiens | 35977481 |
| GSE295758 | bulk RNA-Seq, GPL28038 | Whole peripheral blood | 58 (제목 `subjectN`) | 2가 56명, 4가 2명. 58/58이 timepoint D와 N을 모두 가짐. 샘플 120 | `timepoint` D 60, N 60. 설계 문장: 낮 7:00–9:00, 밤 19:00–21:00 | 20–75 (샘플 120, 고유 나이 36) | 만성 요통. 낮과 밤 채혈 | counts. Gencode GRCh38.p13. 보충 설명에 raw transcript counts와 batch 조정 VST | 41031896 |

헤더를 읽었으나 후보 표에 넣지 않은 자료:

| 접근 번호 | 이유 |
| --- | --- |
| GSE127530 | Whole blood RNA-Seq, 5명 × 9샘플. `time of blood draw` 값은 fasting, 3 hr postprandial, 6 hr postprandial |
| GSE346256 | source `Pheripheral Blood`, RNA-Seq. 설계 문장: 혈청 충격으로 맞춘 뒤 0h, 6h, 12h, 18h, 24h |
| GSE202737 | Plasma, miRNA-Seq, 6명 × 3 (morning, afternoon, evening) |
| GSE214977 | Plasma, miRNA-Seq, 6명, 샘플 17, morning/afternoon/evening |
| GSE254151 | Blood, ncRNA-Seq, 10명 × 6, time T1/T2/T3. 설계 문장은 29샘플·3명이고 헤더 샘플 수는 60 |
| GSE58122 | PBMC RNA-Seq, 48샘플, 시각 필드 없음. 제목 Serial RNA-seq of a single human |
| GSE323980 | Whole blood RNA-Seq, 40샘플, characteristics는 tissue뿐 |
| GSE320507 | blood/PBMC RNA-Seq, 9샘플, characteristics는 tissue뿐. 처리 문장에 Parse Biosciences splitpipe |
| GSE336230 | PBMC RNA-Seq, 9샘플, cell type·treatment만 |
| GSE290897 | blood RNA-Seq, 6샘플, 시각 필드 없음 |
| GSE96060 | blood RNA-Seq, 3명 × 1, 시각 필드 없음 |
| GSE180097 | 샘플 source가 CTC. RNA-Seq 306 |
| GSE254828 | FUCCI-hNSCs |
| GSE342887 | U2OS |
| GSE309217 | 골격근. 헤더를 읽지 않음 |
| GSE292324 | series matrix HTTP 404 |

### 7.3 통제 접근 후보 표

| 자료 | 접근 | 확인한 것 |
| --- | --- | --- |
| GTEx v8 공개 표현형 `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` | 공개 파일. 받은 것은 첫 줄 | 열 `SUBJID SEX AGE DTHHRDY`. 시계 시각 열은 없다 |
| dbGaP phs000424.v8.p2, 자료 `GTEx_Subject_Phenotypes`, pht002742.v8.p2 | 통제 접근. 개인 수준 표는 받지 않음 | 공개 설명: history at time of death (n=13 variables), death circumstances (n=26 variables). 그 변수 이름은 이 페이지에 없다 |
| DTHHRDY | 공개 변수 페이지는 phs000424.v4.p1, phv00169092.v4.p1 | Comment: Death Classification: 4-point Hardy Scale. 0 ventilator, 1 violent and fast, 2 fast natural, 3 intermediate, 4 slow. 시계 시각이 아니다 |

EGA 검색 URL 두 개는 HTTP 404였다. EGA에서 시각이 확인된 자료는 없다.

## 8. 버그·재실행 기록

- `use_lin_only`가 data.table 색인으로 길이 6 논리값을 돌려 `isTRUE`가 FALSE가 되었다. official 부트스트랩이 시작되었다. 그 프로세스를 멈추고 official 반복 파일을 지웠다. 이미 끝난 lin 반복 파일은 남겼다.
- 무작위 서명이 시계 유전자 행렬만으로 분위를 나눠 `decile pool smaller than clock-gene count`로 멈췄다. 전체 유전자의 평균 log2로 다시 계산했다.
- S2 귀무 모형의 무작위효과 SD가 3개여서 `icc_null` 대입이 멈췄다. SD 길이가 1이 아니면 `icc_null`을 비운다.
- `MCPcounter.estimate` 뒤에 검색 경로에서 `lmer`를 찾지 못했다. `lme4::lmer`로 다시 실행했다. 그 전 실행의 빈 계수 파일은 덮어썼다.
- GPL `select()`와 RNAAgeCalc 적재 때 `stack imbalance in '<-', 2 then 4` 경고가 났다. 대응 행 수는 §4.2와 같다.
- 첫 PMC 보충 파일 받기는 HTML이었다. SHA256 작업 증명 뒤 xlsx(1757989바이트)로 바꿨다. 계수가 아니다.

## 9. 판단이 필요한 사항

- S1의 LRT와 R2_within은 cos·sin만 빼서 비교한다. 귀무에도 `hours_awake`가 있다.
- S2의 귀무는 무작위 기울기 항을 유지한다. S2의 ICC는 비워 두었다.
- S3은 혼합모형 LRT를 맞추지 않는다. 시각 범위는 `max(clock_time_h) − min(clock_time_h)`이다. S3 표의 R2 열은 참가자별 OLS R²의 중앙값이다.
- 부트스트랩 시드는 자료마다, 무작위 서명 시드는 자료×시계마다 20261015로 다시 설정한다.
- 구성 점수 상관과 감소 비율의 표본은 control이고 `clock_time_h`가 결측이 아닌 행이다.
- 감소 비율 구간은 두 모형의 REML이 모두 수렴한 반복만 쓴다. 이번 파일은 2000/2000이다.
- GSE196728은 organism이 Homo sapiens이고 처리 문장의 Genome_build가 mm10이다. 후보 표에 둘 다 적고 제외하지 않았다.
- GSE127530의 채혈 시각 값과 GSE346256의 혈청 충격 후 시간은 후보의 시계 시각으로 세지 않았다.
- miRNA-Seq와 ncRNA-Seq는 후보 표에 넣지 않았다.
- `circadianphase`의 단위와 0점은 확인하지 못했다.
- GTEx의 death circumstances 26개 변수 이름은 공개 자료 설명에 없어 열지 않았다. DTHHRDY 설명은 v4.p1 변수 페이지에서 가져왔다.
- official 부트스트랩은 개인 내 상관이 모두 0.99 이상이라 하지 않았다.

## 10. 지시서와 다르게 한 것

없음. 위 판단은 9절에 있다.

## 11. 산출물 목록, Part별 커밋 해시

| Part | 커밋 | 내용 |
| --- | --- | --- |
| A | `541a64e5578477c38a8199de9b201b870a3f131f` | `config/decisions.tsv`, 양성 대조 재판정, GSE48113 제외 목록, GSE122725 `wake_ref_h` |
| B | `590714ebc55e4a5dc71230e21837875b47349ccf` | Peters 부분집합, GPL24324 심볼, circadianphase 분포, 모듈 시계 미제공 |
| C | `54b8055a5876c28f67739f697b6a8d9b5bff69bf` | H1 점추정, lin 부트스트랩, 무작위 서명, S1–S5, 그림 |
| D | `2c9a19620efa491fd289b5f2058e62cebbb7b9b1` | MCPcounter 점수, 구성 코사이너, 보정 모형, 감소 비율 |
| E | `40f2746f9d7bf108e2f12c1daf097b4f2ccbead8` | 검색 기록, 헤더, 후보 표, GTEx 공개 열과 dbGaP 페이지 |

단계 1 결과 `results/stage1/`와 `data/processed/`는 수정하지 않았다. 잠금 자료 GSE48113, GSE122725의 발현값은 읽지 않았다. H4, H5, 사전등록, 고른 RNA-seq의 정리는 하지 않았다.
