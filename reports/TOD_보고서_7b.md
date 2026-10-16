# TOD 보고서 7b — 이탈 처리, 확인 검정, 사후 탐색

## 1. 요약 (수치만, 16줄 이내)
- 이탈 계획 커밋 `7ffde8edc6f8ffbb395482fcf430ad91457d16e1` (2026-10-16 17:29:12 +0900). 수정 커밋 `8e38cffef0cc475f2267b176e8aaac9705a98445` (2026-10-16 18:34:56 +0900). 확인 출력 커밋은 없다.
- 실행 4 관문: L1 `dataset_pass` TRUE, L2 FALSE, L3 TRUE. `lrt_p`, `peak_h`, L3 `mean_delta`와 `p`는 보고서 7과 같다.
- 가족 P, S1, H1 보조 점검, 재현 표 파일은 H2 전에 써졌다. 실행 4는 종료 코드 1이다. 끝까지 끝난 실행이 없으므로 사전등록 결과가 아니다. 그 표의 p와 Holm q는 여기에 옮기지 않는다.
- H2, H3, H4, H5, S4, D1 완전 사례 민감도, D2, D6, 독립 검사는 실행 4에서 끝나지 않았다.
- D1 개발 검증: 유한한 log2 값 4318561개 중 863712개를 가렸다. 유전자 19541개 중 탈락 0. 가림 방식별 모형 M은 §4.
- 실행 4는 `H2 L1 TAms_chrono` 다음 줄에서 멈췄다. Day 1 샘플 `GSM1168771`의 `CD8B`가 NA이고, MCPcounter `CD8 T cells`가 비유한이다.
- 사후 탐색은 하지 않았다. L2 원형 평균 차이와 L2 H1 LRT p는 없다.

## 2. 환경과 자원 (/usr/bin/time 기록)

R은 `/home/jkang/anaconda3/envs/r_env/bin/Rscript`이다. 점수화 Python은 `/HDD8T1/WORK/GeroSci/tod4_clock_env/bin/python`이다. D1 개발 검증과 실행 4의 최대 RSS는 40GB 아래이다. `impute` 1.84.0의 `impute.knn` 기본 인자는 `k=10`, `rowmax=0.5`, `colmax=0.8`, `maxp=1500`, `rng.seed=362436069`이다. 시드 20261015는 가림에만 썼고, `impute.knn`의 `rng.seed`는 바꾸지 않았다.

| 작업 | 로그 | 경과 | 최대 RSS (kB) | 종료 |
|---|---|---|---|---|
| D1 개발 검증 | `logs/d1_development_check.log` | 58:26.91 | 1882620 | 0 |
| run 4 | `logs/confirm_run4.log` | 12:19.31 | 2297912 | 1 |

GENCODE 38 comprehensive GTF는 `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz`에서 받았다. `gzip -t`는 통과했다.

## 3. 이탈 계획 (D1–D7 표, 실행 3 부분 산출물 manifest)

기록 시각은 모두 2026-10-16, 시점 `결과 보기 전`이다. 같은 문장이 `config/decisions.tsv`와 `prereg/deviations.tsv`에 있다.

| 번호 | 결정 |
|---|---|
| D1 | L1 RA 주 분석: 285샘플 유전자 log2에서 비유한 값이 50%를 넘는 유전자를 빼고, 남은 비유한 값을 `impute.knn` k=10, 그 밖은 패키지 기본값으로 채운 뒤 점수를 낸다. 민감도는 비유한 값이 하나라도 있는 유전자를 모두 뺀 완전 사례. |
| D2 | L1 TAms와 모듈의 주 분석은 동결 코드. 민감도는 D1로 채운 행렬로 다시 계산. |
| D3 | L3 RA: 모든 샘플에서 FPKM이 비유한인 행을 뺀다. |
| D4 | L2는 가족 P, S1, S2, S3, S4에서 빠지고 gate fail. 이후 L2 분석은 사후 탐색. |
| D5 | L3 batch-2 민감도는 수행 불가. |
| D6 | L3 주 대응은 GENCODE 38 basic. 민감도는 comprehensive 대응으로 주 시계 3개의 L3 가족 P를 다시 계산. Holm 가족 밖. |
| D7 | 이탈 7.3.1–7.3.5는 버그 수정. 분석 정의를 바꾸지 않는다. |
| D8 | L2가 빠지면 가족 P는 L1과 L3의 6개 검정. 재현은 두 자료 이상 유의이므로 L1과 L3가 모두 유의해야 한다. |

실행 3 부분 산출물은 `results/confirm_run3_partial/`에 있다. 경로, 바이트, SHA-256은 `results/stage7b/run3_partial_manifest.tsv`의 30행이다. 그 안의 TAms 점추정과 부트스트랩 파일은 열지 않았다.

## 4. 수정과 D1 개발 검증 (바뀐 파일과 해시, 가림 방식별 비교 표)

`results/stage7b/changed_files.tsv`의 SHA-256은 다음과 같다. manifest 26개 경로 중 동결 해시와 다른 파일은 `scripts/confirm/run_pipeline.R` 하나이고, 바이트 13578, SHA-256 `d0e9903f19e797644bdcd5b2d306f1e2787dc9d0883bd0938ac9bf986ff8ea4c`이다.

| 경로 | 바이트 | SHA-256 |
|---|---|---|
| `scripts/confirm/run_pipeline.R` | 13578 | `d0e9903f19e797644bdcd5b2d306f1e2787dc9d0883bd0938ac9bf986ff8ea4c` |
| `scripts/confirm/run_locked_analysis.R` | 45245 | `2fe718c2875209f211e2e2de0e4482308d616f6ab664be6831897d4140cccbf7` |
| `scripts/47_part7b_d1_check.R` | 4299 | `94c7428d43c4cc3f759ccf90df38caa24e85c539315d5e9ded7dc1d7766f56ca` |
| `prereg/deviations.tsv` | 6680 | `0ca6488025a3b52baa60239a1280495333d647617f4e4340baf90f95e201f58e` |
| `prereg/ref/gencode.v38.annotation.gtf.gz` | 46556621 | `22020df0d3356e965868f4b193e89fa13e838b950a574349f7fcd461ac01c050` |
| `results/stage7b/gencode_v38_comprehensive_transcripts.tsv` | 18101039 | `b0b433334aef1b8b5ebd2594d922f1f1078dfbf77cd976cc71e148208cc9dd7a` |

comprehensive 표는 전사체 237012행이다. 헤더를 더하면 237013행이다. 열은 `transcript_id`, `gene_id`, `gene_name`, `transcript_id_base`, `gene_id_base`이다. 동결 basic 표는 바꾸지 않았다. 구현 diff는 `prereg/deviations.tsv`의 7b.B1–7b.B4이다.

개발 자료는 GSE39445 control, 참가자 24, 샘플 221, 유전자 19541이다. 잠금 파일은 열지 않았다. 두 가림은 각각 `set.seed(20261015)`이다. 낮은 발현 가림의 가중치는 평균 발현 순위의 역순이다. D1 절차는 가림 결과에 따라 바꾸지 않았다.

| scheme | n_masked | n_finite_before_mask | n_genes_dropped | n_genes_kept | n_values_imputed |
|---|---|---|---|---|---|
| random | 863712 | 4318561 | 0 | 19541 | 863712 |
| low_expression | 863712 | 4318561 | 0 | 19541 | 863712 |

모형 M. `n_subjects` 24, `n_samples` 221.

| scheme | clock | two_a | lrt_stat | lrt_p |
|---|---|---|---|---|
| unmasked | RA_default | 1.98498899424706 | 1.11800352518344 | 0.571779550797631 |
| unmasked | RA_Peters | 3.86882060477724 | 4.99623159826274 | 0.0822398090507771 |
| unmasked | RA_GTEx | 4.56221470455666 | 0.88115920268956 | 0.643663244864561 |
| random | RA_default | 3.59732072810026 | 4.05301218414456 | 0.131795197933585 |
| random | RA_Peters | 2.71773343430116 | 1.85034526385562 | 0.396462971000013 |
| random | RA_GTEx | 4.38646857517034 | 0.898765464439748 | 0.638021860434852 |
| low_expression | RA_default | 4.3409417895527 | 6.92911326014269 | 0.0312868740107463 |
| low_expression | RA_Peters | 3.55436088898673 | 3.47797892472477 | 0.175697860237701 |
| low_expression | RA_GTEx | 7.00189641952909 | 2.99852054281814 | 0.223295276970678 |

개인 내 Pearson 상관. 시계와 가림마다 참가자 24명이고, 각 참가자의 샘플 수는 6 이상 10 이하이다. 아래 네 값은 오름차순 1번째, 12번째, 13번째, 24번째이다.

| scheme | clock | 1st | 12th | 13th | 24th |
|---|---|---|---|---|---|
| random | RA_default | 0.626305914179446 | 0.827598231601562 | 0.854806842318687 | 0.986296888827285 |
| random | RA_Peters | 0.618986442570281 | 0.835427722056753 | 0.846627617620954 | 0.988561023189761 |
| random | RA_GTEx | 0.267436003969967 | 0.922440825952411 | 0.922551679440376 | 0.997764062619681 |
| low_expression | RA_default | 0.0848787029081507 | 0.815186061182987 | 0.815973021268296 | 0.970329635842812 |
| low_expression | RA_Peters | -0.296219801746613 | 0.769080045586621 | 0.773151983849318 | 0.969502956425735 |
| low_expression | RA_GTEx | 0.550274450188178 | 0.845913017770249 | 0.846211678198251 | 0.947078062107116 |

0보다 작은 상관은 `low_expression`의 `RA_Peters` 1개이다. 나머지 143개는 0보다 크다.

실행 4가 H2 전에 쓴 D1·D3 감사는 다음과 같다. L1 주 RA 점수 855/855의 value가 비어 있지 않다. L3 주 RA 점수 348/348의 value가 비어 있지 않다.

| 표 | n_genes_in | n_genes_dropped | n_genes_kept | n_values_imputed | k |
|---|---|---|---|---|---|
| L1 `d1_impute_audit.tsv` | 19541 | 2983 | 16558 | 281938 | 10 |
| L1 `d1_complete_audit.tsv` |  | 18387 | 1154 |  |  |
| L3 `d3_ra_audit.tsv` |  | 27614 | 31679 |  |  |

## 5. 관문

숫자는 실행 4의 `gate.tsv`이다. 보고서 7의 실행 3 관문과 같다.

L1 Day 1. 참가자 22, 샘플 147. 전처리: scale `log`, probes 43758, 샘플 285, 유전자 19541, `n_dropped_multi` 0, `n_dropped_no_symbol` 13253.

| gene | lrt_p | peak_h | pass | dataset_pass |
|---|---|---|---|---|
| PER1 | 1.35631359279823e-11 | 7.39244745511937 | TRUE | TRUE |
| PER3 | 2.79164389620025e-06 | 4.93116485780503 | TRUE | TRUE |

L2 period 1. 참가자 11, 샘플 70. 전처리: scale `log`, probes 27189, `n_multi` 911, 유전자 20179, 샘플 141. 빛 조건 참가자 수는 bright 3, dim 8이다. `stopped.tsv`는 `gate fail`이다.

| gene | lrt_p | peak_h | pass | dataset_pass |
|---|---|---|---|---|
| PER1 | 1.33395937684947e-05 | 8.62707238405378 | FALSE | FALSE |
| PER3 | 2.56672604887414e-06 | 6.37573752707371 | FALSE | FALSE |

L3. 참가자 58, 샘플 116. 분기 `transcript_to_gene`, 입력 행 336781, 대응 행 110352, 유전자 59293, batch-2 샘플 0. PER1은 `log2(CPM+1)`의 morning − evening이다. 로그에 `count2FPKM` 경고가 있다. 유전자 길이 없음 46.5721%.

| gene | mean_delta | p | pass | dataset_pass |
|---|---|---|---|---|
| PER1 | 0.565581638383806 | 0.00672475968762667 | TRUE | TRUE |

## 6. 가족 P

실행 4는 `family_P.tsv`, `family_S1.tsv`, `h1_secondary.tsv`, `replication.tsv`를 쓴 뒤 H2에서 멈췄다. 이 파일들은 `results/confirm_run4_partial/`에 있다. 종료 코드가 1이므로 이탈을 포함한 사전등록 결과가 아니다. p, Holm q, 2A, 정점, mean Δ, 재현 판정, 위상·부호는 이 보고서에 적지 않는다.

### 6.1 L1 모형 M

표 없음. 실행 4 로그는 L1 가족 P의 `TAms_chrono`, `TAms_mort`, `RA_default`와 가족 S1의 `RA_Peters`, `RA_GTEx`에 대해 부트스트랩 2000회 완료 줄을 찍었다.

### 6.2 L3 대응 t

표 없음. `paired_deltas.tsv`는 부분 산출물에 있다.

### 6.3 재현 판정

표 없음.

### 6.4 보조 점검 (L1 정점 ±3시간, L3 부호)

표 없음. `h1_secondary.tsv`는 부분 산출물에 있다.

## 7. 가족 S1

`family_S1.tsv`는 부분 산출물에 있다. Holm q와 판정은 여기에 적지 않는다. L2 S1은 관문 실패로 수행되지 않았다.

## 8. H2

수행되지 않았다. 로그의 마지막 R 오류는 다음과 같다.

`Error in h2_after_holm("L1", day1, l1_scores, l1_comp, primary_clocks,  : composition rows do not match the score samples`

이 함수는 가족 P 판정이 `significant`인 시계에서만 `H2` 줄을 찍는다. 찍힌 줄은 `H2 L1 TAms_chrono` 하나이고, 그 다음이 위 오류이다.

Day 1 147샘플의 MCPcounter 8개 열 중 비유한 값은 1개이다. 샘플 `GSM1168771`, 참가자 `BB0148`, `analysis_set` day1, 열 `CD8 T cells`. 점수화 표본 285개에서는 비유한 값이 3개이고, 나머지는 Day 4의 `GSM1168705`와 `GSM1168708`, 참가자 `BB0100`, 같은 열이다. MCPcounter 유전자 표에서 `CD8 T cells`에 해당하는 유전자는 `CD8B` 하나이다. 세 샘플의 `CD8B` log2는 NA이다. L3 `log2(FPKM+1)` 116샘플의 MCPcounter 값 중 비유한 값은 0개이다.

## 9. H3

수행되지 않았다.

## 10. H4 (L1)

수행되지 않았다.

## 11. H5

수행되지 않았다.

## 12. 모듈 수준 S4

수행되지 않았다. 확인 수 / 8은 없다.

## 13. 민감도 (사전등록 §8 항목 + D1·D2·D6)

수행되지 않았다. 3-MAD 제외, 임의 기울기, L3 combat, batch-2, D1 완전 사례 모형 M, D2 대치 행렬 재점수, D6 comprehensive 대응 t의 출력 파일은 없다. L1 완전 사례 점수 파일 `ra_scores_complete_case.tsv`와 D1 대치 행렬 `log2_genes_d1_imputed.tsv`는 H2 전에 써졌다. 그 파일의 모형 M은 없다.

## 14. 독립 검사

수행되지 않았다. 가족 P를 끝까지 마친 실행이 없다.

## 15. 사후 탐색 (관문 실패 자료)

수행되지 않았다. 지시서 7b는 확인 출력을 커밋한 뒤에만 Part D를 한다. 확인 출력 커밋은 없다.

### 15.1 L2 시각 기준 진단 (인용, 유전자별 정점과 원형 차이)

없음.

### 15.2 L2 H1 탐색 (LRT, BH q)

없음.

### 15.3 L2 H4 탐색

없음.

## 16. 버그·재실행·이탈 기록 (어느 실행의 출력이 사전등록 결과인지)

사전등록 결과로 쓸 끝까지 성공한 실행은 없다. 실행 4의 `run_id.tsv`는 `run4`, 시작 2026-10-16 18:35:33 +0900이다. 경과 12:19.31, 최대 RSS 2297912 kB, 종료 코드 1. 부분 산출물 43개 파일, 204173477바이트는 `results/confirm_run4_partial/`로 옮겼다. 경로와 SHA-256은 `results/stage7b/run4_partial_manifest.tsv`이다. `results/confirm/`은 비어 있는 상태가 아니라 경로 자체가 없다.

실행 4는 D1 주 RA와 D3 주 RA를 지나 가족 P 부트스트랩을 쓴 뒤 H2에서 멈췄다. 멈춘 이유는 §8이다. 표본을 빼거나 `CD8 T cells`를 채우는 규칙은 결과 전에 정해져 있지 않으므로 코드를 고치지 않았다.

## 17. 판단이 필요한 사항

L1 H2는 사전등록에서 조성 점수 8개를 모형 M의 고정 효과로 넣는다. Day 1 샘플 147개 중 1개, `GSM1168771`, 의 `CD8 T cells`가 비유한이다. 이 값은 표식 유전자 `CD8B` 하나의 결측에서 나온다. 이 샘플을 빼거나 값을 채우면 표본 또는 조성이 바뀐다. 그 규칙이 생기기 전에는 실행 5를 시작하지 않는다.

## 18. 지시서와 다르게 한 것

없음. 수정 커밋의 diff 표시는 지워진 `results/confirm/L3/log2_fpkm.tsv`와 `gene_counts.tsv`를 새 GENCODE 파일과 비슷한 내용으로 짝지었다. 커밋에 들어간 파일의 SHA-256은 §4의 표와 같다. 실행 3 부분 산출물 경로는 그 커밋 뒤에도 `results/confirm_run3_partial/`에 있다.

## 19. 산출물 목록, 커밋 해시

| 산출 | 경로 |
|---|---|
| 이탈 계획 | `config/decisions.tsv`, `prereg/deviations.tsv`, `results/stage7b/run3_partial_manifest.tsv` |
| 수정과 개발 검증 | `scripts/confirm/run_locked_analysis.R`, `scripts/47_part7b_d1_check.R`, `results/stage7b/changed_files.tsv`, `results/stage7b/d1_development_audit.tsv`, `results/stage7b/d1_development_model_m.tsv`, `results/stage7b/d1_development_within_pearson.tsv` |
| GENCODE 38 comprehensive | `prereg/ref/gencode.v38.annotation.gtf.gz`, `results/stage7b/gencode_v38_comprehensive_transcripts.tsv` |
| 실행 4 부분 산출 | `results/confirm_run4_partial/`, `results/stage7b/run4_partial_manifest.tsv` |
| 로그 (커밋하지 않음) | `logs/d1_development_check.log`, `logs/confirm_run4.log` |

| 커밋 | 해시 | 시각 |
|---|---|---|
| TOD7b deviation plan (before results) | `7ffde8edc6f8ffbb395482fcf430ad91457d16e1` | 2026-10-16 17:29:12 +0900 |
| TOD7b code fixes | `8e38cffef0cc475f2267b176e8aaac9705a98445` | 2026-10-16 18:34:56 +0900 |
| TOD7b confirm output run4 | 없음 |  |
| TOD7b posthoc L2 | 없음 |  |
| TOD7b report | 이 보고서를 넣는 커밋 |  |
