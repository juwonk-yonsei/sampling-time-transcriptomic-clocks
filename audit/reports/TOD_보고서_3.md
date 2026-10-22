# TOD 보고서 3 — 판단 사항 처리, 자료 역할, RNA-seq 개발, H1·H2·H4 개발, 관문 시험, 모듈 시계

## 1. 요약 (수치만, 16줄 이내)

- 자료 역할 8건 기록. GSE295758 잠금. 참가자 58, D와 N이 모두 있는 참가자 58, 나이 20–75, 샘플 120 (D 60, N 60).
- GSE39445 `circadianphase`: n=438, 최소 −15, 최대 15, 전부 3의 배수. `use_as_time_axis_in_section_5_1` TRUE. §5.1에서 T_circ를 맞췄다.
- GSE113883: TPM 행 58278, 선택 SYMBOL 35581, 다중 SYMBOL ENSEMBL 583 제외. 양성 대조 fail. PER1 정점 6.88888783794359, PER3 4.3428731541087, PER3–ARNTL 원형 거리 2.04939962890015. ARNTL의 행렬 기호는 BMAL1.
- GSE196728: 유전자 ID 12726, human-style 12559 (0.9868772592), ENSMUSG 0, mouse-style 160 (0.01257268584). `use_for_analysis` TRUE. SL 48샘플·8명. ARNTL `log2(FPKM+1)` 유한값 0. 양성 대조 fail.
- RNA-seq §3.2–§3.5 미실행. 시계 2A, R2_within, 정점, LRT, official–lin 상관, 무작위 서명, 구성 점수, 그림 없음.
- H2 `p_med` [구간]. GSE39445: RA_GTEx 1.52682501330863 [−2.78679825578026, 3.12129350554563], RA_Peters 0.119399199505249 [−0.622005757921436, 1.74440304777819], RA_default 3.03185040820847 [−3.61344425227156, 5.60947774217046]. GSE56931: RA_GTEx 0.413998323546127 [−0.489487133435497, 1.51854132384192], RA_Peters 0.31261747469437 [−0.220350765673326, 0.956070197020275], RA_default 0.95871462885305 [0.0195522302126154, 2.26398470535335]. 두 REML 수렴 2000/2000. RNA-seq 미계산.
- H4 수면 제한 평균 차이 (restriction − control) [구간]: RA_GTEx 3.21825858250506 [−3.6756284174559, 10.5636435279068], RA_Peters −0.887815897468002 [−3.55671400467413, 1.65945853998723], RA_default 2.4193949175103 [0.441165353256765, 4.75863844147715]. 두 조건을 모두 가진 참가자 22. 부트스트랩 유한값 2000/2000.
- GSE56931 0시 deprivation − baseline [Wald 95%]: RA_GTEx 9.55678469886526 [−1.33972523936263, 20.4532946370931], RA_Peters 4.25692107534861 [−0.550134480487721, 9.06397663118493], RA_default 1.48210725460707 [−3.76069025470135, 6.7249047639155]. 4시: 7.75465389906009 [−3.38152983393059, 18.8908376320508], −1.43969785540255 [−6.35165459335046, 3.47225888254535], −3.57972990558479 [−8.93747867408558, 1.77801886291601].
- 관문 G1: GSE39445 pass, GSE56931 pass, GSE113883 fail, GSE196728 fail. G2: fail, pass, pass, fail (참가자 0). G3 `hours_awake` 정점 PER1 22.7776132062867, PER3 21.9104219255107.
- 모듈 시계 경로 `not_found`. 체크섬 비교 없음. Part G 미수행. 복합 시계 점수 없음.
- 마지막 분석 커밋 `74328cfaa3484c4a87f08d5b1e8f39458f0eff25`.

## 2. 환경과 자원 (바뀐 것, 새 가상환경, /usr/bin/time 기록)

새 가상환경은 없다. `org.Hs.eg.db` 3.22.0 (Artistic-2.0)을 `config/env_versions.txt`에 붙였다. Part B의 R `sessionInfo()`도 같은 파일에 붙였다. RNAAgeCalc 1.22.0, MCPcounter 1.2.0, lme4 2.0-6. GPU를 쓰지 않았다. `ulimit`을 두지 않았다. OMP/MKL 스레드 상한을 두지 않았다.

| 실행 | elapsed | max RSS (kB) | exit |
| --- | --- | --- | --- |
| 16_partB_prepare.py | 0:22.47 | 187944 | 0 |
| 17_partB_expr.R | 0:17.58 | 1516952 | 0 |
| 18_partB_qc.py | 0:04.73 | 233708 | 0 |
| 19_partC_rnaseq.R, ARNTL 별칭 수정 전 | 0:15.55 | 1265916 | 1 |
| 19_partC_rnaseq.R | 0:16.55 | 1380360 | 0 |
| 20_partD_projection.R, data.table 필터 수정 전 | 0:15.39 | 1341828 | 1 |
| 20_partD_projection.R | 7:47.93 | 1505944 | 0 |
| 21_partE_h4.R | 3:54.91 | 1541912 | 0 |
| 22_partF_gates.R | 0:20.80 | 1444452 | 0 |

최대 RSS는 40 GB 아래였다.

## 3. 결정 기록 (config/decisions.tsv에 이번에 더한 행)

| date | id | content | reason | timing |
| --- | --- | --- | --- | --- |
| 2026-10-16 | 3.1.1 | S1 comparison method is approved. LRT and R2_within drop only cos and sin. The null keeps hours_awake. | Instruction 3 §1.1 item 1. Report 2 §9. | prior |
| 2026-10-16 | 3.1.2 | S2 null and ICC are approved. The null keeps the random slopes. ICC stays blank. | Instruction 3 §1.1 item 2. Report 2 §9. | prior |
| 2026-10-16 | 3.1.3 | S3 definition is approved. S3 is descriptive. No mixed-model LRT. The time range is max(clock_time_h) minus min(clock_time_h). The R2 column is the median of per-subject OLS R2. | Instruction 3 §1.1 item 3. Report 2 §9. | prior |
| 2026-10-16 | 3.1.4 | The seed reset method is approved. This stage uses the same method. Bootstrap seeds are reset to 20261015 per dataset. Random-signature seeds are reset to 20261015 per dataset by clock. | Instruction 3 §1.1 item 4. Report 2 §9. | prior |
| 2026-10-16 | 3.1.5 | The composition-analysis sample is approved: control samples with non-missing clock_time_h. | Instruction 3 §1.1 item 5. Report 2 §9. | prior |
| 2026-10-16 | 3.1.6 | The reduction-ratio interval rule is approved. The metric changes to the rhythm-vector projection ratio in instruction 3 §4. | Instruction 3 §1.1 item 6. Report 2 §9. | prior |
| 2026-10-16 | 3.1.7 | GSE196728 is designated a development dataset. Gene IDs in the processed file are checked in instruction 3 §2.3. It is used only if that check passes. | Instruction 3 §1.1 item 7 and §1.2. | prior |
| 2026-10-16 | 3.1.8 | GSE127530 and GSE346256 are approved as not candidates. | Instruction 3 §1.1 item 8. Report 2 §9. | prior |
| 2026-10-16 | 3.1.9 | miRNA-Seq and ncRNA-Seq are approved as not candidates. | Instruction 3 §1.1 item 9. Report 2 §9. | prior |
| 2026-10-16 | 3.1.10 | The unit of GSE39445 circadianphase is checked from the paper and the data in instruction 3 §2.1. It is not decided in this part. | Instruction 3 §1.1 item 10. Report 2 §9. | prior |
| 2026-10-16 | 3.1.11 | The GTEx judgment is approved. H3 uses GSE295758. | Instruction 3 §1.1 item 11. Report 2 §9. | prior |
| 2026-10-16 | 3.1.12 | Skipping the official bootstrap when all within-person correlations are at least 0.99 is approved. RNA-seq datasets follow instruction 3 §3.2. | Instruction 3 §1.1 item 12. Report 2 §9. | prior |
| 2026-10-16 | 3.2.GSE39445 | Role development. Technology microarray. Uses: H1, H2, H4 development. | Instruction 3 §1.2. | prior |
| 2026-10-16 | 3.2.GSE56931 | Role development. Technology microarray. Uses: H1, H2, H4 development. | Instruction 3 §1.2. | prior |
| 2026-10-16 | 3.2.GSE113883 | Role development, designated in this stage. Technology RNA-seq. Uses: year-scale H1, H2, H4 (DLMO) development, and the gate test. | Instruction 3 §1.2. | prior |
| 2026-10-16 | 3.2.GSE196728 | Role development, designated in this stage, conditional. Technology RNA-seq. Uses: sea-level condition H1 and H2 development, only if instruction 3 §2.3 passes. | Instruction 3 §1.2. | prior |
| 2026-10-16 | 3.2.GSE48113 | Role locked. Technology microarray. Uses: confirmation H1 and H4 (sleep phase misalignment). Expression values stay unread. | Instruction 3 §1.2. Instruction 1 §1.3. | prior |
| 2026-10-16 | 3.2.GSE122725 | Role locked. Technology microarray. Uses: confirmation H1 and H4 (schedule shift, light). Expression values stay unread. | Instruction 3 §1.2. Instruction 1 §1.3. | prior |
| 2026-10-16 | 3.2.GSE295758 | Role locked, designated in this stage. Technology RNA-seq. Uses: confirmation H1 (morning-evening) and H3 (age interaction). From this decision, instruction 1 §1.3 applies. Expression values and supplementary files are not opened. | Instruction 3 §1.2. | prior |
| 2026-10-16 | 3.2.GSE150337 | Role held. Technology RNA-seq. Not opened. | Instruction 3 §1.2. | prior |
| 2026-10-16 | 3.3 | RA_Peters is the provisional primary clock. The final choice is made at preregistration. This stage analyzes RA_GTEx, RA_Peters, and RA_default. | Instruction 3 §1.3. Coverages already computed: GPL15331 280/289 = 0.9688581315; GPL10379 278/289 = 0.9619377163; GPL24324 276/289 = 0.9550173010380623. No fourth platform coverage row is in those files. | prior |

## 4. 자료 준비

### 4.1 GSE39445 circadianphase (인용, offset 원형 SD, 판정)

출처: Möller-Levet et al. 2013 PNAS, PMID 23440187, SI (PMC3607048, `1217154110_pnas.201217154SI.pdf`).

영점: "We referenced each participant’s sampling times (t) to their melatonin peak (p) by following the transformation tp = 24t − 24p (both, t and p in fractional day), where tp is time away from the melatonin peak in hours."

구간: "Note that Circadian-Time-bin (3-h-wide circadian bins) reflects the melatonin phase-aligned sampling times sorted in 45° bins to create categorical values. Because of the variance between subjects in circadian phase, this resulted in the 10 time points to cover 11 circadian bins."

저장값 438개. 최소 −15, 최대 15. 고유값 −15,−12,−9,−6,−3,0,3,6,9,12,15. 전부 3의 배수. 인접 고유값 차이 3. `hours_awake − circadianphase`의 참가자 내 선형 SD 최댓값 1.133893419027682, 그 SD가 0보다 큰 참가자–조건 1 (control AF0182). 중앙값 0. `data_matches_hour_scale` TRUE. `use_as_time_axis_in_section_5_1` TRUE.

| condition | n | clock offset 원형 SD min, median, max | clock offset 원형 평균 min, median, max | awake offset 원형 SD max | awake offset 원형 평균 min, median, max |
| --- | --- | --- | --- | --- | --- |
| control | 24 | 0, 0.005952245801, 0.7134796627 | 2.261666478, 3.875, 7.170341323 | 1.042100947 | 19.5, 19.5, 22.09878254 |
| sleep_restriction | 24 | 0, 0.002500026697, 0.05248817133 | 1.258333301, 5.255832947, 8.166666669 | 8.049455132e-08 | 19.5, 22.5, 22.5 |

### 4.2 GSE113883 (파일·해시, 메타데이터 요약, DLMO 정의 인용, 프로토콜 인용, 유전자 대응, QC)

| 파일 | bytes | SHA-256 |
| --- | --- | --- |
| data/raw/GSE113883/GSE113883_series_matrix.txt.gz | 8995 | 541bfd45d6de012ba54b046355f82fe4a1532ed18356c51d316abd5050ce3ee9 |
| data/raw/GSE113883/GSE113883_TPM_processed_data.txt.gz | 7674591 | 0546e6bf916f3928668004b7c72be703473594155f67e4642590dd7aacda5d9b |

샘플 153. `condition`은 전부 `control`. `period`는 `time_course`. `clock_time_h`는 `local_time`. `salivary_dlmo_circadian_time`이 있는 샘플 95, 참가자 7. 나머지 4명·58샘플은 결측. 개인 나이는 특성에 없다. 논문 포함 기준 문장: "between the ages of 18–60 years".

DLMO 언급: "Chronotype was determined by the timing of preferred sleep and wake times, sleep logs, and actigraphy, and confirmed by the timing of salivary dim light melatonin onset (DLMO) from samples collected at home prior to the experiment."

본문: "The reference time was the clock time each sample was taken, as melatonin data are not universally available."

`salivary_dlmo_circadian_time`이 기상 후 시간인지, DLMO의 시계 시각인지, DLMO를 CT 14로 둔 시각인지는 가져온 본문과 SI에 없다. 상태: 확인 못 함. `circadian_ref_h`에는 그 필드의 숫자를 넣었다. DLMO가 있는 7명에서 `circadian_ref_h`의 참가자 내 고유값 수는 1이다 (§5.1 `n_unique_time_within_max`).

프로토콜: "During these 3 days, subjects remained at rest (either seated or semi-recumbent) during the 16 waking hours, and were provided with regularly spaced isocaloric snacks. On arrival to the CRU, light intensity was 100-300 lux during 16 waking hours and dark during the 8 sleep hours. Night 1 was an adaptation night to familiarize subjects with the laboratory environment. Beginning upon awakening on Day 2, light intensity was <20 lux during the waking hours. Blood samples were taken every 2 h starting 2 h after wake on on Day 2 for 29 h, comprising a total of 15 draws per subject."

GEO: "At each timepoint for each subject, 2.5 mL of whole blood was drawn into a PAXgene RNA tube and processed according to manufacturer protocol (BD Vacutainer)."

유전자 대응 (`org.Hs.eg.db` 3.22.0). ID 예: `TSPAN6__ENSG00000000003.14`. 행 58278, ENSEMBL 키 58233, 한 ENSEMBL에 SYMBOL이 둘 이상인 키 583은 제외, 미대응 22645, 대응 행 35633, 선택 SYMBOL 35581, 접두 기호와 DB SYMBOL이 같은 선택 행 29387. 한 SYMBOL에 ID가 여럿이면 평균 TPM이 가장 큰 행을 남겼다.

QC (`log2(TPM+1)`): 유전자 35581, 샘플 153, 샘플 중앙값의 중앙값 0, MAD 0, 3 MAD 표시 0. PCA는 분산 상위 1000유전자, n=153. PC1–PC5의 `subject_id` R²: 0.38081027, 0.7127455584, 0.8912413292, 0.9131099657, 0.9497519066. `cos_sin_clock_time_h` R²: 0.1083693546, 0.03438985255, 0.01166019137, 0.01854352578, 0.00266801604. `condition`과 `period` R²는 0. `batch_fields`는 비어 있다.

### 4.3 GSE196728 (종 확인 예시와 비율, 사용 여부, 쓴 경우 메타데이터·프로토콜·FPKM 방법)

| 파일 | bytes | SHA-256 |
| --- | --- | --- |
| data/raw/GSE196728/GSE196728_series_matrix.txt.gz | 12211 | d689f6b1a602692d4f55b63adabb2df57a806d396a38165ad5362a8338363082 |
| data/raw/GSE196728/GSE196728_Peru_raw.xlsx | 8544929 | e772a1eedc23b936206026b5958319fa7d82af589030e1c5295e2a0d6ad816fb |

유전자 ID 12726. human-style 12559 (0.9868772592). human ENSEMBL 0. mouse ENSEMBL 0. mouse-style 160 (0.01257268584). 그 밖 7. 예 10개: A1BG-AS1, A2M-AS1, A4GALT, AAAS, AACS, AAGAB, AAK1, AAMDC, AAMP, AAR2. 규칙: ENSMUSG이 0이고, 대문자 또는 ENSG 분율 ≥ 0.8이고, 첫 글자만 대문자인 분율 < 0.2이면 사용. `use_for_analysis` TRUE. 샘플 제목과 count 열 144/144 대응. 1열은 대소문자만 달랐다 (`sub6_PU2`와 `Sub6_PU2.raw`).

SL 샘플의 `condition`은 `control`. PU와 RI는 H1·H2에 넣지 않았다. SL 48샘플, 참가자 8. `clock_time_h`는 2, 6, 10, 14, 18, 22. `elapsed_h`는 비워 두었다.

프로토콜 (PMID 35977481): "Whole-blood samples were collected from every subject at 4 h intervals throughout 24 h in each altitude." "Lifestyle conditions were controlled and standardized for all participants, especially meal times (breakfast between 7:30 and 8:00, lunch between 13:00 and 13:30, dinner between 20:00 and 20:30), bedtime (22:30) and wake-up time (5:30)." "At each altitude level (SL, PU, and RI) and time around the clock (06:00, 10:00, 14:00, 18:00, 22:00 and 02:00), a whole blood sample (2.5 mL) was collected in a dedicated blood RNA tube (PAXgene, PreAnalytiX, Qiagen, Hombrechtikon, Swizerland) from a venipuncture of the arm for each participant (144 samples in total)." 조도(lux): 확인 못 함.

GEO 처리: "Genome_build: mm10". "Reads were mapped to genome (/shareDB/iGenomes/Mus_musculus/UCSC/mm10/Sequence/STAR_index) using STAR (default parameters)". Organism: Homo sapiens.

FPKM: `RNAAgeCalc::count2FPKM`, `idtype="SYMBOL"`, 길이 열 `bp_length`. 유전자 12726, 샘플 144, 모든 샘플이 NA인 유전자 858. 경고: "Can't find gene length for 6.7421% genes when converting raw count to FPKM." ARNTL count는 144샘플 모두 유한 (최소 8, 최대 157). ARNTL FPKM과 `log2(FPKM+1)`은 144샘플 모두 NA.

QC (`log2(FPKM+1)`): 유전자 12726, 샘플 144, 샘플 중앙값의 중앙값 1.425256153, MAD 0.204796836, 3 MAD 표시 2 (GSM5899766, GSM5899772). 표시만 하고 빼지 않았다. PCA n=144, 유전자 1000. PC1 `condition` R² 0.7891902071, `subject_id` R² 0.06102123954, `cos_sin_clock_time_h` R² 0.00512523754. PC2–PC5 `subject_id` R²: 0.610818112, 0.266376151, 0.7931764371, 0.5689087414.

### 4.4 GSE295758 (메타데이터 요약 표, 처리 문장 원문, 보충 파일 목록)

발현과 보충 파일은 열지 않았다. series matrix는 `!` 행만 읽었다.

| 파일 | bytes | SHA-256 |
| --- | --- | --- |
| data/locked/GSE295758/GSE295758_series_matrix.txt.gz | 9691 | 76f7eec9d779f82dc01a018ffee53e84bfe245d97150142e7025085d5d4fd7d9 |

| 항목 | 값 |
| --- | --- |
| 샘플 | 120 |
| 참가자 | 58 |
| 참가자당 샘플 min, median, max | 2, 2, 4 |
| D와 N이 모두 있는 참가자 | 58 |
| 나이 min, q1, median, q3, max | 20, 40, 52, 61, 75 |
| 성별 (참가자) | female 29, male 29 |
| timepoint | D 60, N 60 |
| opioids/group (참가자) | 0: 37, 1: 16, 2: 5 |

설계 원문: "RNA-sequencing of human peripheral blood cells. Blood samples were collected during the day (7:00-9:00) and night (19:00-21:00) from participants with chronic low back pain."

처리 원문: "Raw-sequencing reads were cleaned by BGI using SOAPnuke v1.5.2 with the parameters “-l 15 -q 0.2 -n 0.05”." "Reads were aligned to GENCODE’s primary human reference genome vGRCh38.p13 with Hisat2 v2.2.1 using parameters "--dta --sensitive --no-discordant --no-mixed"." "Samtools v1.10 was used to sort resulting BAM files, mark duplicates, index BAM files." "StringTie v2.1.5 assembled aligned reads into transcripts and quantified their abundance using GENCODE’s vGRCh38.p13 comprehensive gene annotation." "IsoformAnalyzeR v1.17.4 in R v4.2.1 assigned a known gene to “novel” transcripts based on their genomic coordinates when possible." "For network analysis, ComBat-seq from sva v3.36.0 adjusted raw and filtered transcript counts for batch effects. DESeq2 v1.26.0 applied the variance stabilizing transformation to these adjusted counts." "Assembly: Gencode GRCh38.p13".

보충 목록 (다운로드 안 함): `GSE295758_combatVstTranscriptCounts.csv.gz` 목록 크기 80M, 수정 2025-04-28 08:05. `GSE295758_rawTranscriptCounts.csv.gz` 목록 크기 29M, 수정 2025-04-28 08:06.

### 4.5 모듈 시계 (찾아본 위치, 정한 경로, 저장소 커밋·태그, setup 스크립트, 기대·실제 체크섬, 파일 목록, 라이선스, 읽기 환경, 구조 확인 표, 채점 코드 확인, 플랫폼 포함률)

`config/paths.yaml`의 `module_clock_dir`와 `module_clock_code_dir`는 비어 있다. 상태 `not_found`. Part G `skipped`. 다운로드하지 않았다. 기대 체크섬과 실제 체크섬은 없다. 파일 목록, 라이선스, 읽기 환경, 구조, 채점 코드, 플랫폼 포함률은 없다.

| 위치 | 확인 | 결과 |
| --- | --- | --- |
| /home/jkang | 이름 검색 깊이 4, git 검색 깊이 5 | v0.3.0 저장소 없음 |
| /HDD8T1/WORK/GeroSci/revision | 폴더 | git 저장소가 아님. `JAAA-D-26-02318.pdf`, `major_revision_1.txt`만 있음 |
| /HDD8T1/WORK/GeroSci/sc-module-aging | git, `scripts/00_download.sh`가 Zenodo 22166800과 tAge를 언급 | pyproject 버전 0.2.0. 태그 없음. v0.3.0 커밋 메시지 없음. 다운로드 스크립트에 체크섬 없음. 커밋 `a0e9d090cbabcdde441a599bb6443ed3f5ba0350` |
| sc-module-aging/data/external/tAge | 중첩 git | 태그 1.5.0 |
| /HDD8T1/WORK/GeroSci/clock-specificity | git | v0.3.0 태그 없음 |
| /SSD1T, /HDD6T, /HDD8T2, /ssd240 | git 깊이 5 | v0.3.0에 해당하는 저장소 없음 |

## 5. RNA-seq 개발 H1·H2

### 5.1 양성 대조 표와 판정

t = `clock_time_h`. GSE113883은 전체 153샘플·11명. GSE196728은 `condition=control` (SL) 48샘플·8명. (a')는 PER1과 PER3 정점이 4–12에 들어가는지이다. (b)는 PER3–ARNTL 원형 거리 ≥ 6이다. 둘 다 통과해야 §3.2–§3.5를 한다.

| dataset | gene | matrix_symbol | peak_h | lrt_p | two_a | n_samples | n_subjects | reml_converged |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE113883 | PER1 | PER1 | 6.88888783794359 | 2.75027621424795e-09 | 0.620551996416151 | 153 | 11 | TRUE |
| GSE113883 | PER2 | PER2 | 4.8732860858618 | 9.27215074557918e-09 | 0.320745822015945 | 153 | 11 | TRUE |
| GSE113883 | PER3 | PER3 | 4.3428731541087 | 7.58295067863991e-11 | 0.449708658716999 | 153 | 11 | TRUE |
| GSE113883 | NR1D1 | NR1D1 | 0.913464171710277 | 0.00052689207268246 | 0.280681915296824 | 153 | 11 | TRUE |
| GSE113883 | NR1D2 | NR1D2 | 3.21755304741653 | 3.13505667430625e-15 | 1.02092699550255 | 153 | 11 | TRUE |
| GSE113883 | ARNTL | BMAL1 | 2.29347352520855 | 0.00933679580233329 | 0.270966203946098 | 153 | 11 | TRUE |
| GSE113883 | CRY1 | CRY1 | 3.54047677289818 | 3.4799087393093e-05 | 0.458041061525944 | 153 | 11 | TRUE |
| GSE113883 | CRY2 | CRY2 | 5.54979841432233 | 0.0522701033256887 | 0.112952026576084 | 153 | 11 | TRUE |
| GSE113883 | DBP | DBP | 6.45559270353684 | 0.145623612973268 | 0.15735235377638 | 153 | 11 | TRUE |
| GSE113883 | TEF | TEF | 4.14074509211546 | 1.12522933891847e-16 | 0.475282772862979 | 153 | 11 | TRUE |
| GSE196728 | PER1 | PER1 | 6.41144900180597 | 0.187843896897575 | 0.244053056953198 | 48 | 8 | TRUE |
| GSE196728 | PER2 | PER2 | 8.66116903749671 | 0.779269043831126 | 0.0430049625371678 | 48 | 8 | TRUE |
| GSE196728 | PER3 | PER3 | 4.0865909290848 | 0.0573355762517646 | 0.205664314738903 | 48 | 8 | TRUE |
| GSE196728 | NR1D1 | NR1D1 | 3.0687819107217 | 0.177876684941138 | 0.269374842161372 | 48 | 8 | TRUE |
| GSE196728 | NR1D2 | NR1D2 | 5.60676333371992 | 0.173007965160013 | 0.17897050287304 | 48 | 8 | TRUE |
| GSE196728 | ARNTL | ARNTL |  |  |  | 0 | 0 | FALSE |
| GSE196728 | CRY1 | CRY1 | 2.07645926880775 | 0.238769501047531 | 0.180483693180549 | 48 | 8 | TRUE |
| GSE196728 | CRY2 | CRY2 | 2.82251648365978 | 0.28805771229822 | 0.169365592823513 | 48 | 8 | TRUE |
| GSE196728 | DBP | DBP | 1.90643315523274 | 0.277409026659748 | 0.232303927766801 | 48 | 8 | TRUE |
| GSE196728 | TEF | TEF | 4.6245806659025 | 0.120447853540006 | 0.274154718526839 | 48 | 8 | TRUE |

| dataset | (a') | (b) | PER3–ARNTL 원형 거리 | 판정 |
| --- | --- | --- | --- | --- |
| GSE113883 | TRUE | FALSE | 2.04939962890015 | fail |
| GSE196728 | TRUE | FALSE |  | fail |

### 5.2 시계 값 (포함 수, 대체 수, 분포, official–lin 상관)

미계산. 두 자료 모두 양성 대조 fail. §3.2를 실행하지 않았다.

### 5.3 주 모형 표 (보고서 2 §5.2와 같은 열)

미계산.

### 5.4 무작위 서명 비교, 민감도

미계산. S2와 S3을 실행하지 않았다.

### 5.5 구성 점수와 하루 리듬

미계산. MCPcounter를 RNA-seq 행렬에 실행하지 않았다.

### 5.6 그림 경로

없음.

## 6. H2 지표 교체

마이크로어레이 보정 전 2A는 단계 2 lin 점추정과 같다. 절대차 최댓값 1.67435256770432e-08 (GSE56931 RA_GTEx). 단계 2의 감소 비율 파일은 바꾸지 않았다. GSE113883과 GSE196728의 투영 비율과 개인 내 상관은 미계산이다. lin 점수와 MCPcounter 점수가 §3.2–§3.4에서 만들어지지 않았다.

부트스트랩은 단계 2의 `boot_draws_GSE39445.tsv`, `boot_draws_GSE56931.tsv`이다. 시드 20261015, 2000회. 반복마다 보정 전과 보정 후 REML을 같은 표본에 맞췄다. 구간은 백분위 95%, type 7. 두 REML이 수렴한 반복 2000/2000. `p_med`가 유한한 반복 2000/2000.

### 6.1 투영 매개 비율 표 (자료 × 시계: 보정 전 LRT p, |r_u|, p_med [구간], 2A_orth [구간])

| dataset | clock | LRT p | \|r_u\| | p_med | p_med 0.025 | p_med 0.975 | 2·A_orth | 2·A_orth 0.025 | 2·A_orth 0.975 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | 0.643663244864707 | 2.28110735226545 | 1.52682501330863 | −2.78679825578026 | 3.12129350554563 | 3.48456468800812 | 0.241478816531498 | 12.8600842877736 |
| GSE39445 | RA_Peters | 0.0822398090507771 | 1.93441030239009 | 0.119399199505249 | −0.622005757921436 | 1.74440304777819 | 1.74358784754434 | 0.0787539562310089 | 6.08656898539641 |
| GSE39445 | RA_default | 0.571779550797631 | 0.992494497123529 | 3.03185040820847 | −3.61344425227156 | 5.60947774217046 | 0.95841180193476 | 0.129140714408146 | 7.38862058830963 |
| GSE56931 | RA_GTEx | 0.0977660921595171 | 4.20647807567416 | 0.413998323546127 | −0.489487133435497 | 1.51854132384192 | 3.447711065702 | 0.208924349805865 | 9.09567352385791 |
| GSE56931 | RA_Peters | 0.000105116530090898 | 3.18998678662389 | 0.31261747469437 | −0.220350765673326 | 0.956070197020275 | 0.211018938462172 | 0.0309072888951493 | 3.38139640664492 |
| GSE56931 | RA_default | 0.0187399469939367 | 2.30847035606281 | 0.95871462885305 | 0.0195522302126154 | 2.26398470535335 | 2.17152819704037 | 0.0795287085460286 | 7.37246031430758 |

GSE39445 n=24/221. GSE56931 n=14/130. 시계는 lin. 보정 모형은 주 모형에 8개 구성 점수를 더한 REML이다.

### 6.2 시계–구성 개인 내 상관 표

control이고 `clock_time_h`가 있는 샘플. 참가자 평균을 뺀 값의 Pearson 상관.

| dataset | clock | T cells | CD8 T cells | Cytotoxic lymphocytes | NK cells | B lineage | Monocytic lineage | Myeloid dendritic cells | Neutrophils |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | −0.197505568551426 | −0.20993722988163 | −0.00652672795021659 | 0.153424147857917 | 0.00162958351418949 | −0.0751569093830547 | −0.127940487648064 | 0.0898933138743642 |
| GSE39445 | RA_Peters | 0.253419065593425 | 0.275264284282644 | 0.326561534806622 | −0.0590223528422122 | −0.395185396489041 | 0.0341233876811266 | 0.174823245199813 | −0.176828598312388 |
| GSE39445 | RA_default | −0.367031299924488 | −0.312152806466937 | −0.0465379780293033 | 0.164981241053637 | −0.0897983813311287 | −0.170386466714967 | −0.44024889823705 | 0.176782338997834 |
| GSE56931 | RA_GTEx | −0.167426190946125 | −0.103189848561146 | 0.121402272663986 | 0.182468589808423 | −0.117560502139278 | 0.0326210028360044 | −0.151479785815535 | 0.125995490260275 |
| GSE56931 | RA_Peters | 0.115665384915649 | 0.175872158322802 | 0.0569999551074129 | 0.0868243261043974 | −0.0268850275914042 | −0.103116034430659 | −0.21674582233898 | −0.2012098950184 |
| GSE56931 | RA_default | −0.210468177688468 | −0.142650068812006 | 0.268974030361793 | 0.317228017655641 | −0.175819804348665 | 0.219911846078777 | −0.180289193121932 | 0.158275743587801 |

## 7. H4 개발

### 7.1 시각 축 비교 (ΔAIC 표)

ΔAIC = AIC(해당 축) − AIC(T_clock). AIC는 ML. GSE39445 control, 세 시각이 모두 있는 샘플 221, 참가자 24. T_circ의 t는 `circadian_ref_h` (`circadianphase`). T_wake의 t는 `hours_awake`. GSE113883은 DLMO가 있는 샘플 95, 참가자 7. 시계(lin)는 미계산. T_dlmo의 참가자 내 고유 시각 수 최댓값은 1이다.

| dataset | outcome | axis | n | 2A | peak_h | LRT p | AIC_ML | ΔAIC |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | T_clock | 221 | 4.56221470453089 | 13.7137336093738 | 0.643663244864707 | 2097.46921733189 | 0 |
| GSE39445 | RA_GTEx | T_circ | 221 | 6.37145260149485 | 8.88818248511983 | 0.439035781907637 | 2096.70402781165 | −0.765189520243439 |
| GSE39445 | RA_GTEx | T_wake | 221 | 5.94282476146924 | 4.74571867363762 | 0.482338353793975 | 2096.89215766982 | −0.577059662066404 |
| GSE39445 | RA_Peters | T_clock | 221 | 3.86882060478018 | 6.84503499428447 | 0.0822398090507771 | 1687.69391330542 | 0 |
| GSE39445 | RA_Peters | T_circ | 221 | 4.28560487395079 | 2.24050372540086 | 0.0417687748552327 | 1686.33893244041 | −1.35498086501207 |
| GSE39445 | RA_Peters | T_wake | 221 | 4.22778858781387 | 21.9621992111359 | 0.0469387442921098 | 1686.57232122317 | −1.12159208225103 |
| GSE39445 | RA_default | T_clock | 221 | 1.98498899424706 | 12.0064473932843 | 0.571779550797631 | 1675.04549952352 | 0 |
| GSE39445 | RA_default | T_circ | 221 | 4.33642388183845 | 8.92866708954318 | 0.0613794960626038 | 1670.58214416881 | −4.46335535471258 |
| GSE39445 | RA_default | T_wake | 221 | 4.42485548061331 | 4.47604782501751 | 0.0542125472189409 | 1670.33381725105 | −4.71168227246403 |
| GSE39445 | PER1 | T_clock | 221 | 0.439316278046525 | 6.90324597433964 | 9.73597605831301e-18 | 49.3146296852052 | 0 |
| GSE39445 | PER1 | T_circ | 221 | 0.471935709121597 | 3.19855166600159 | 6.9006383979288e-21 | 34.8106911915784 | −14.5039384936268 |
| GSE39445 | PER1 | T_wake | 221 | 0.474871117706392 | 22.7776132062867 | 4.60083182863588e-21 | 33.9999375729718 | −15.3146921122334 |
| GSE39445 | PER3 | T_clock | 221 | 0.313199677007082 | 6.08583286807854 | 1.177374716677e-06 | 178.218994599916 | 0 |
| GSE39445 | PER3 | T_circ | 221 | 0.349396429394201 | 2.40991202668664 | 3.10631509855509e-08 | 170.948954276746 | −7.27004032316967 |
| GSE39445 | PER3 | T_wake | 221 | 0.343532957703934 | 21.9104219255107 | 5.96963454018461e-08 | 172.255451359859 | −5.9635432400562 |
| GSE39445 | ARNTL | T_clock | 221 | 0.113823858576977 | 14.8618875133315 | 0.0449426621692836 | 61.5515021516275 | 0 |
| GSE39445 | ARNTL | T_circ | 221 | 0.0846088865417176 | 12.118197907295 | 0.168602675502959 | 64.1958169737784 | 2.64431482215087 |
| GSE39445 | ARNTL | T_wake | 221 | 0.0836293206360964 | 7.37049546630091 | 0.177501452760856 | 64.2986847319189 | 2.7471825802914 |
| GSE113883 | PER1 | T_clock | 95 | 0.635894607850685 | 7.11236086937803 | 2.08593469333727e-05 | 139.445159334557 | 0 |
| GSE113883 | PER1 | T_dlmo | 95 | 3.18745438957537 | 1.8593424135732 | 0.427950004195629 | 159.303078360869 | 19.8579190263124 |
| GSE113883 | PER3 | T_clock | 95 | 0.417430954930491 | 4.15372302271565 | 1.28342626519682e-05 | 42.3758695718028 | 0 |
| GSE113883 | PER3 | T_dlmo | 95 | 5.82776738604978 | 1.37385729951439 | 0.00230210869050469 | 52.7547944506195 | 10.3789248788167 |
| GSE113883 | ARNTL (BMAL1) | T_clock | 95 | 0.228050037626824 | 2.46740778829001 | 0.132175815798831 | 102.230716954272 | 0 |
| GSE113883 | ARNTL (BMAL1) | T_dlmo | 95 | 8.14348339278813 | 1.36482673174668 | 0.0804424810245483 | 101.237535821041 | −0.993181133231374 |

REML 수렴은 위 25행 모두 TRUE. REML singular는 모두 FALSE. ML 수렴은 모두 TRUE.

### 7.2 수면 제한 효과 (LRT, 평균 차이 [구간], 조건별 2A·정점, 민감도)

GSE39445. `condition` 수준은 control, sleep_restriction. 샘플 438. control 참가자 24, sleep_restriction 참가자 24, 두 조건이 모두 있는 참가자 22. 평균 수준 차이는 `score ~ cos_t + sin_t + condition + (1 | subject_id)`의 `conditionsleep_restriction` 계수 (restriction − control). condition 주효과 LRT는 그 모형 대 `score ~ cos_t + sin_t + (1 | subject_id)`, df 1. 상호작용 LRT는 `condition:cos_t`와 `condition:sin_t`를 더한 모형 대 평균 수준 모형, df 2. 둘 다 ML. 95% 구간은 두 조건이 모두 있는 22명을 단위로 한 부트스트랩 2000회, 시드 20261015, 백분위 type 7. 유한한 반복 2000/2000. SD_between(control)은 control 샘플의 참가자 평균 점수의 SD이다.

| clock | 차이 | SD_between | 차이/SD_between | condition LRT p | 상호작용 LRT p |
| --- | --- | --- | --- | --- | --- |
| RA_GTEx | 3.21825858250506 | 18.7953223613647 | 0.171226570134304 | 0.201787186217236 | 0.794124215878378 |
| RA_Peters | −0.887815897468002 | 14.9283841120993 | −0.0594716675831267 | 0.341583324792914 | 0.526315716465803 |
| RA_default | 2.4193949175103 | 9.3806000710098 | 0.257914728183253 | 0.00875104390077384 | 0.0972621301754652 |

구간: RA_GTEx [−3.6756284174559, 10.5636435279068], RA_Peters [−3.55671400467413, 1.65945853998723], RA_default [0.441165353256765, 4.75863844147715]. LRT 수렴은 세 시계 모두 TRUE.

| clock | condition | n_subjects | n_samples | 2A | peak_h | LRT p |
| --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | control | 24 | 221 | 4.56221470453089 | 13.7137336093738 | 0.643663244864707 |
| RA_GTEx | sleep_restriction | 24 | 217 | 8.7988044863254 | 11.3725102001519 | 0.208577536917136 |
| RA_Peters | control | 24 | 221 | 3.86882060478018 | 6.84503499428447 | 0.0822398090507771 |
| RA_Peters | sleep_restriction | 24 | 217 | 2.6802593103669 | 3.04100854979083 | 0.303346676234885 |
| RA_default | control | 24 | 221 | 1.98498899424706 | 12.0064473932843 | 0.571779550797631 |
| RA_default | sleep_restriction | 24 | 217 | 7.32591167769245 | 12.6043511192968 | 0.0002877111080939 |

민감도: 같은 모형에 `hours_awake` 선형 항을 더했다. 샘플 438. 같은 22명 추출. 유한한 반복 2000/2000.

| clock | 차이 | 구간 | 차이/SD_between | condition LRT p | 상호작용 LRT p |
| --- | --- | --- | --- | --- | --- |
| RA_GTEx | 3.40021038660913 | [−3.43480614608607, 10.8122719564661] | 0.180907266246124 | 0.173002263840059 | 0.810818777423374 |
| RA_Peters | −0.903085622523648 | [−3.59046653798959, 1.68218554445904] | −0.0604945328136154 | 0.333232557562609 | 0.532580903447561 |
| RA_default | 2.51725004437555 | [0.560714909350125, 4.87660700300516] | 0.268346377131561 | 0.00530849676702354 | 0.0620187898990932 |

### 7.3 수면 상태 효과 (인용, 모형 1 LRT, 모형 2 시각별 차이 표)

GSE56931. baseline과 deprivation. recovery는 빼었다. 샘플 209, 참가자 14. baseline 130, deprivation 79.

인용 (PMID 25197809, PMC4173916): "parallel blood draws every 4 h during 24-h normal baseline (sleep and wakefulness), 38 h of continuous wakefulness, and during subsequent recovery sleep". "On Day 2 the sleep deprivation started at 08:00 and lasted until 22:00 on Day 3, a total of 38 h." "Time points for blood draws (every 4 h across wake [yellow] and sleep periods [black])". "Lights on and lights off periods were defined as times 8:00 – 22:00 and 22:00 – 8:00, respectively." GEO 설계: baseline은 day 1 8hr, 12hr, 16hr, 20hr와 day 2 0hr, 4hr. continuous wakefulness는 day 2 8hr, 12hr, 16hr, 20hr와 day 3 0hr, 4hr, 8hr, 12hr, 16hr, 20hr. baseline 0시·4시는 lights-off 22:00–08:00 안에 있고, 논문은 baseline 채혈이 sleep와 wakefulness에 걸쳐 있다고 적는다. deprivation 0시·4시는 Day 2 08:00에 시작해 Day 3 22:00에 끝나는 38시간 continuous wakefulness 안에 있다.

모형 1은 §7.2와 같은 구조이고 `condition` 자리에 `period`를 썼다. 계수는 deprivation − baseline. 부트스트랩은 하지 않았다. LRT 수렴은 TRUE.

| clock | 차이 | period LRT p (df 1) | 상호작용 LRT p (df 2) |
| --- | --- | --- | --- |
| RA_GTEx | 1.03664014830376 | 0.606102091094096 | 0.228378322567274 |
| RA_Peters | 0.361179241794387 | 0.685739682537868 | 0.148648500914668 |
| RA_default | −0.479127633605237 | 0.615941000310065 | 0.132933545372051 |

| clock | period | n_samples | 2A | peak_h | LRT p |
| --- | --- | --- | --- | --- | --- |
| RA_GTEx | baseline | 130 | 8.41295615134832 | 13.2671292169983 | 0.0977660921595171 |
| RA_GTEx | deprivation | 79 | 5.09288933643161 | 18.870925851358 | 0.526018825077103 |
| RA_Peters | baseline | 130 | 6.37997357324779 | 4.54192615314072 | 0.000105116530090898 |
| RA_Peters | deprivation | 79 | 3.63951824878276 | 1.02267711531618 | 0.270661255212625 |
| RA_default | baseline | 130 | 4.61694071212563 | 9.88701727310532 | 0.0187399469939367 |
| RA_default | deprivation | 79 | 5.9257656280342 | 13.5946406390903 | 0.0245803809626937 |

모형 2: `score ~ factor(clock_time_h) * period + (1 | subject_id)`. REML. Wald 구간은 추정 ± `qnorm(0.975)` × SE. REML 수렴 TRUE, singular FALSE. 시각별 n은 시계 행과 구성 행이 같다.

| outcome | 0시 차이 [구간] (n 13/12) | 4시 (n 11/13) | 8시 (n 28/14) | 12시 (n 25/14) | 16시 (n 27/13) | 20시 (n 26/13) |
| --- | --- | --- | --- | --- | --- | --- |
| RA_GTEx | 9.55678469886526 [−1.33972523936263, 20.4532946370931] | 7.75465389906009 [−3.38152983393059, 18.8908376320508] | −9.26595962836218 [−18.1414479643868, −0.390471292337528] | −0.501567458099604 [−9.56264944634051, 8.5595145301413] | 4.57693032260941 [−4.59526035370705, 13.7491209989259] | −1.36044794444658 [−10.5947240549448, 7.87382816605164] |
| RA_Peters | 4.25692107534861 [−0.550134480487721, 9.06397663118493] | −1.43969785540255 [−6.35165459335046, 3.47225888254535] | −2.33228387454639 [−6.24606982620346, 1.58150207711068] | 0.775171981312385 [−3.22093959457541, 4.77128355720017] | 2.05092261227367 [−1.99463290350574, 6.09647812805308] | 1.14490824052939 [−2.92832120580161, 5.21813768686039] |
| RA_default | 1.48210725460707 [−3.76069025470135, 6.7249047639155] | −3.57972990558479 [−8.93747867408558, 1.77801886291601] | −2.19640449967481 [−6.46605127274687, 2.07324227339725] | −1.18308376225736 [−5.54222747554174, 3.17605995102703] | 4.09517977109311 [−0.317612711742576, 8.5079722539288] | −1.07898662116354 [−5.52177610044393, 3.36380285811685] |
| T cells | −0.0419025641019562 [−0.155270238583353, 0.0714651103794401] | 0.0365241820802602 [−0.0793205737669574, 0.152368937927478] | −0.0268068852040762 [−0.119114330778416, 0.0655005603702635] | −0.0790589061388084 [−0.173306208498596, 0.015188396220979] | −0.0301424094090128 [−0.125554179759717, 0.0652693609416912] | −0.0963104249282952 [−0.192373769551135, −0.000247080305455685] |
| CD8 T cells | −0.027157806419721 [−0.185503107611987, 0.131187494772545] | 0.11966818337978 [−0.0421334380525407, 0.281469804812101] | −0.0401243928571446 [−0.169046941372883, 0.0887981556585939] | −0.0627830773581022 [−0.194416989672267, 0.0688508349560626] | −0.100250210134449 [−0.233512379068941, 0.0330119588000441] | −0.133962154421672 [−0.268135613527587, 0.000211304684243241] |
| Cytotoxic lymphocytes | −0.000965315697714677 [−0.12474234316202, 0.122811711766591] | 0.0992669307715681 [−0.0272136715409463, 0.225747533084083] | −0.210952178571431 [−0.311733330633637, −0.110171026509225] | −0.0225320797271655 [−0.125431714717351, 0.0803675552630204] | 0.0308325076763829 [−0.0733390013683564, 0.135004016721122] | −0.013453430772283 [−0.118336667648392, 0.0914298061038265] |
| NK cells | −0.0171262357225946 [−0.176616067650507, 0.142363596205318] | 0.167567726680763 [0.00458387060947152, 0.330551582752055] | −0.272310584821428 [−0.402190607218305, −0.14243056242455] | 0.0315501429016975 [−0.101054001798168, 0.164154287601564] | 0.148566200821871 [0.0143285393280798, 0.282803862315662] | 0.0887469084420288 [−0.0464042683598877, 0.223898085243945] |
| B lineage | −0.0554972237481346 [−0.16016281566064, 0.0491683681643703] | −0.0451425192522308 [−0.152090844694977, 0.0618058061905159] | −0.0165383253968319 [−0.101751837813468, 0.0686751870198046] | −0.0829256593699115 [−0.169932376229611, 0.00408105748978829] | 0.00398659160820949 [−0.0840973530666927, 0.0920705362831117] | −0.077978830409936 [−0.166665780679061, 0.010708119859189] |
| Monocytic lineage | 0.0799433173555889 [−0.0266492501517387, 0.186535884862917] | 0.0738231391963413 [−0.0350982668422249, 0.182744545234907] | −0.0323885663265216 [−0.119179096001373, 0.0544019633483292] | −0.0543053860479711 [−0.142919953258277, 0.0343091811623351] | 0.023831313352083 [−0.0658782362693031, 0.113540862973469] | −0.011619450093625 [−0.101941705462135, 0.078702805274885] |
| Myeloid dendritic cells | 0.0574883299429373 [−0.0410023515231731, 0.155979011409048] | 0.0500176081255331 [−0.0506289749703243, 0.150664191221391] | 0.040272119047611 [−0.0399298399000601, 0.120474077995282] | 0.00638828142437012 [−0.0754968833493785, 0.0882734461981187] | −0.0194323106505655 [−0.102327146014667, 0.0634625247135364] | −0.0152515935317092 [−0.0987111675729724, 0.0682079805095541] |
| Neutrophils | 0.0497811036362136 [−0.046249568626954, 0.145811775899381] | 0.0445779767241536 [−0.0535497498174975, 0.142705703265805] | 0.14941749999999 [0.0712288111032873, 0.227606188896692] | 0.0448831268081196 [−0.0349494125478725, 0.124715666164112] | 0.0474409595154467 [−0.0333785885253817, 0.128260507556275] | 0.0470600003421284 [−0.0343118953852535, 0.12843189606951] |

n은 baseline/deprivation이다.

## 8. 관문 시험

기준은 결과를 본 뒤에 바꾸지 않았다. G2 창은 07:00–09:00과 19:00–21:00, 양 끝 포함. 두 창에 샘플이 모두 있는 참가자만 썼다. 차이는 아침 평균 − 저녁 평균. 유전자 검정은 대응 t 검정 (양측). 두 유전자 모두 차이 > 0이고 p < 0.01이면 G2 통과. 시계 구간은 같은 대응 t 검정의 95% 구간이다.

### 8.1 G1 (자료별)

통과 조건: PER1과 PER3의 LRT p가 모두 < 0.01이고, PER3–ARNTL 원형 거리 ≥ 6.

| dataset | PER1 peak | PER1 LRT p | PER3 peak | PER3 LRT p | ARNTL peak | ARNTL 기호 | 거리 | 판정 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | 6.90324597433964 | 9.73597605831301e-18 | 6.08583286807854 | 1.177374716677e-06 | 14.8618875133315 | ARNTL | 8.77605464525298 | pass |
| GSE56931 | 5.79155389251083 | 1.58787990764694e-24 | 4.8945603273957 | 3.48666322673004e-14 | 14.7633078773663 | ARNTL | 9.86874754997058 | pass |
| GSE113883 | 6.88888783794359 | 2.75027621424795e-09 | 4.3428731541087 | 7.58295067863991e-11 | 2.29347352520855 | BMAL1 | 2.04939962890015 | fail |
| GSE196728 | 6.41144900180597 | 0.187843896897575 | 4.0865909290848 | 0.0573355762517646 |  | ARNTL |  | fail |

GSE196728 ARNTL의 유한한 `log2(FPKM+1)`은 0이다. REML 수렴 FALSE.

### 8.2 G2 (자료별 참가자 수, 유전자별 차이와 p, 시계 차이 [구간]과 SD)

| dataset | gene | n | 차이 | SD | p | 95% 구간 | 유전자 기준 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | PER1 | 15 | 0.457282198733333 | 0.390086823168289 | 0.000462323702321782 | [0.241259316455056, 0.673305081011611] | TRUE |
| GSE39445 | PER3 | 15 | 0.2552221425 | 0.421460938601263 | 0.0342660715272614 | [0.0218248542120165, 0.488619430787983] | FALSE |
| GSE56931 | PER1 | 14 | 0.666719392857143 | 0.274073921211973 | 5.29896048819085e-07 | [0.508473822048558, 0.824964963665728] | TRUE |
| GSE56931 | PER3 | 14 | 0.320380964285714 | 0.26166409536946 | 0.000514941233529122 | [0.169300613595117, 0.471461314976312] | TRUE |
| GSE113883 | PER1 | 11 | 0.663564636365152 | 0.438466332707098 | 0.000522275801438845 | [0.368998945957966, 0.958130326772337] | TRUE |
| GSE113883 | PER3 | 11 | 0.247747093162695 | 0.251436643678242 | 0.00845948309236336 | [0.0788296575149657, 0.416664528810423] | TRUE |
| GSE196728 | PER1 | 0 |  |  |  |  | FALSE |
| GSE196728 | PER3 | 0 |  |  |  |  | FALSE |

G2 판정: GSE39445 fail, GSE56931 pass, GSE113883 pass, GSE196728 fail. GSE196728 control의 `clock_time_h`는 2, 6, 10, 14, 18, 22라 두 창에 동시에 들어간 참가자가 0이다.

| dataset | clock | n | 차이 | SD | 95% 구간 |
| --- | --- | --- | --- | --- | --- |
| GSE39445 | RA_GTEx | 15 | −1.19848633720139 | 29.570255863457 | [−17.5739482137342, 15.1769755393314] |
| GSE39445 | RA_Peters | 15 | 2.94326940353098 | 13.8922894187166 | [−4.7500240466312, 10.6365628536932] |
| GSE39445 | RA_default | 15 | 0.0360875520089001 | 14.9192225221645 | [−8.22590249526463, 8.29807759928243] |
| GSE56931 | RA_GTEx | 14 | 0.641323562546064 | 14.2657338282151 | [−7.59546592354024, 8.87811304863237] |
| GSE56931 | RA_Peters | 14 | 3.21349555127013 | 5.39283490516242 | [0.0997653792416124, 6.32722572329864] |
| GSE56931 | RA_default | 14 | 2.23332903123105 | 5.88912814943408 | [−1.16695232881608, 5.63361039127819] |

GSE113883과 GSE196728의 시계 행은 미계산이다.

### 8.3 G3

GSE39445 control. t = `hours_awake`. n=24, 샘플 221. REML 수렴 TRUE.

| gene | peak_h | 2A | LRT p |
| --- | --- | --- | --- |
| PER1 | 22.7776132062867 | 0.474871117706392 | 4.60083182863588e-21 |
| PER3 | 21.9104219255107 | 0.343532957703934 | 5.96963454018461e-08 |

## 9. 모듈 시계 개발 분석

Part G를 하지 않았다. §2.5 경로를 찾지 못했고 체크섬을 비교하지 않았다. 모듈 유전자와 계수를 읽지 않았다.

### 9.1 입력 준비 (로그 밑, 오솔로그 대응 수, 래퍼 경로, 파이프라인별 특징 비율·대체 비율)

없음.

### 9.2 복합 시계 H1 표 (보고서 2 §5.2와 같은 열)

없음.

### 9.3 모듈 시계 H1 표 (자료 × 목표 × 모듈: R2_within, LRT p, BH q, 정점, 2A)

없음.

### 9.4 복합 시계 H2 (투영 매개 비율)

없음.

### 9.5 복합 시계와 RA_Peters의 개인 내 상관

없음.

## 10. 버그·재실행 기록 (없으면 "없음")

`19_partC_rnaseq.R` 첫 실행은 exit 1이다. `org.Hs.eg.db` 3.22.0에서 ENSG00000133794의 SYMBOL은 BMAL1이고, ARNTL은 그 SYMBOL의 ALIAS이다. GSE113883 행렬에 ARNTL 행이 없어 정점 슬롯이 비었고, `is.finite`가 길이 0을 반환했다. ARNTL 요청 유전자를 행렬의 BMAL1 행으로 읽고, 정점이 없으면 NA로 둔 뒤 다시 실행했다. 두 번째 실행은 exit 0이다.

`20_partD_projection.R` 첫 실행은 exit 1이다. data.table 부분집합 안에서 `dataset`이 열 이름으로 평가되어 단계 2 점추정 행을 찾지 못했다. `which()`로 걸러 낸 뒤 다시 실행했다. 두 번째 실행은 exit 0이다. 보정 전 2A는 단계 2 lin과 절대차 1e-6 미만으로 맞았다.

GSE196728 양성 대조에서 singular-fit 경고가 기록되었다. 유한값이 있는 유전자의 저장 REML 수렴은 TRUE이다. ARNTL은 유한값이 0이라 모형을 맞추지 않았다. Part F의 GSE56931 적합 중 singular-fit 경고가 1회 기록되었다. 저장 REML 수렴은 TRUE이다.

## 11. 판단이 필요한 사항

- 지시서 §1.3은 네 플랫폼 포함률 ≥ 0.95를 적는다. 계산된 포함률 행은 세 플랫폼이다. 세 수치를 결정 3.3에 적고, 넷째 행은 만들지 않았다. 잠정 주 시계는 지시서대로 RA_Peters이다.
- GSE113883에 condition 칸이 없다. 전 샘플을 `condition=control`, `period=time_course`로 두었다.
- GSE196728 SL만 `condition=control`로 두었다. PU와 RI는 H1·H2에 넣지 않았다.
- 종 판정 문턱은 human-style ≥ 0.8, mouse-style < 0.2, ENSMUSG = 0이다. `C#orf#` 형태가 첫 글자만 대문자인 규칙에 걸려, mouse-style이 0이어야 한다는 규칙이면 사람 기호 파일이 탈락한다.
- count 열 1개는 제목과 대소문자만 달랐고 대응이 하나라 그 열을 썼다.
- 한 ENSEMBL이 SYMBOL을 둘 이상 가지면 그 키를 제외했다. 583키.
- `circadianphase`를 §5.1의 시간 축으로 썼다. SI의 tp는 시간이고, 저장값은 3의 배수이며 45° 문장은 3시간 구간을 가리킨다. 값은 이산 구간 라벨이다. control AF0182의 awake−phase 오프셋은 상수가 아니다.
- `salivary_dlmo_circadian_time`의 영점: 확인 못 함. `circadian_ref_h`에는 그 숫자를 그대로 넣었다.
- 발현 행렬은 `results/stage3/expr/`에 두고 git에서 제외했다.
- 모듈 시계 저장소 v0.3.0을 찾지 못했다. 경로는 비웠다. `sc-module-aging`의 시계 파일을 쓰지 않았다. Part G를 하지 않았다.
- GSE113883 샘플 중앙값의 MAD가 0이라 3 MAD 규칙은 아무도 표시하지 않는다.
- 양성 대조가 실패한 자료에는 §3.2–§3.5를 하지 않았다. §4, §5.1, §6.2가 그 자료의 lin 시계를 적고 있어도, 시계 점수와 RNA-seq 구성 점수는 만들지 않았다. 유전자 모형(PER1, PER3, ARNTL)은 맞췄다.
- GSE113883에서 지시서의 ARNTL은 행렬 기호 BMAL1 행이다. ENSG00000133794.

## 12. 지시서와 다르게 한 것

§4, §5.1, §6.2는 GSE113883의 lin 시계를, §2.3을 통과한 GSE196728의 lin 시계를 결과 변수로 적는다. 그 점수는 §3.2에서만 만들도록 되어 있고, 양성 대조 실패로 §3.2를 실행하지 않았다. 그래서 그 시계의 투영 비율, 개인 내 상관, 시각 축, 아침−저녁 차이는 비어 있다. §9의 보수적 선택으로 적는다.

## 13. 산출물 목록, Part별 커밋 해시

| Part | 커밋 | 메시지 |
| --- | --- | --- |
| A | c56a1bb43a27daa119024c25953a5919dad41c96 | TOD3 Part A |
| B | 0859337d58ba64dbbd574c108d4cfda8035f9e0c | TOD3 Part B |
| C | dbdedebb906818a0544401e31dac2220081f5b26 | TOD3 Part C |
| D | 69c028b1ee891702bc3cc895926c2374c172104d | TOD3 Part D |
| E | 8618fd97829916cfbae88f48a0dc25f6fd0574ec | TOD3 Part E |
| F | 74328cfaa3484c4a87f08d5b1e8f39458f0eff25 | TOD3 Part F |
| G | 없음 | 경로를 찾지 못해 수행하지 않음 |

단계 3에서 커밋한 스크립트: `scripts/15_partA_decisions.py`, `scripts/16_partB_prepare.py`, `scripts/17_partB_expr.R`, `scripts/18_partB_qc.py`, `scripts/19_partC_rnaseq.R`, `scripts/20_partD_projection.R`, `scripts/21_partE_h4.R`, `scripts/22_partF_gates.R`.

표는 `results/stage3/`의 `partA_decisions.tsv`, circadianphase 판정·offset, GSE113883·GSE196728 샘플·대응·QC, GSE295758 메타데이터와 보충 목록, `module_clock_status.tsv`, `rnaseq/positive_control.tsv`, `h2/projection_ci.tsv`, `h2/clock_composition_within.tsv`, `h4/time_axis_comparison.tsv`, `h4/sleep_restriction_bootstrap.tsv`, `h4/sleep_state_model2_wald.tsv`, `gates/g1_decision.tsv`, `gates/g2_genes.tsv`, `gates/g2_clocks.tsv`, `gates/g3_peaks.tsv`에 있다. 발현 행렬은 `results/stage3/expr/`에 있고 커밋하지 않았다.
