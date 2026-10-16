# TOD 보고서 7 — 잠금 해제와 확인 검정

## 1. 요약 (수치만, 14줄 이내)
- OSF 상태: not uploaded. 동결 해시 `f19f79c1988d75e7085a3f570b80a96a2e3350b368335ac80b0ae4fb162120ba`. manifest 26개 파일 일치, 불일치 0.
- 잠금 해제 전 커밋 `082ed49ed6f4733c28dcdcce9b502f7fb02c0862` (2026-10-16 16:46:36 +0900). 잠금 해제 커밋 `dc2bed1a64adb374884128eb91720ffbed3ad3e2` (2026-10-16 16:49:28 +0900).
- 관문 L1: PER1과 PER3 모두 pass TRUE. L2: PER1과 PER3 모두 pass FALSE, dataset_pass FALSE. L3: PER1 pass TRUE.
- 가족 P 표, Holm q, 재현 판정, H2, H3, H4, H5, S4는 만들어지지 않았다.
- L1 `RA_default` 점수가 비어 있어 실행 3이 `prepare_score`에서 멈췄다. L1 RA 점수 파일의 비어 있지 않은 value는 0/855이다. L3 RA 점수 파일은 0/348이다.
- 실행 3회(종료 코드 1, 1, 1). 이탈 5건. 사전등록 결과로 쓸 끝낸 실행은 없다.
- 확인 출력 커밋 `a868eee221933a57196936fcacbada4798485c82` (2026-10-16 17:21:12 +0900). 이 커밋은 멈춘 실행의 산출물이다.

## 2. 환경과 자원 (/usr/bin/time 기록)

R은 `/home/jkang/anaconda3/envs/r_env/bin/Rscript`이다. 점수화 Python은 `/HDD8T1/WORK/GeroSci/tod4_clock_env/bin/python`이다. 세 실행 모두 최대 RSS가 40GB 아래이다.

| 실행 | 로그 | 경과 | 최대 RSS (kB) | 종료 |
|---|---|---|---|---|
| L3 보충 파일 다운로드 | `logs/l3_download.log` | 0:30.88 | 79760 | 0 |
| run 1 | `logs/confirm_run1.log` | 0:01.75 | 235852 | 1 |
| run 2 | `logs/confirm_run2.log` | 1:31.22 | 2099556 | 1 |
| run 3 | `logs/confirm_run3.log` | 5:38.77 | 2097204 | 1 |

## 3. 잠금 해제 전 기록 (OSF 상태, 판단 사항 결정, 동결 확인 표)

`config/osf_status.txt`는 `not uploaded`이다. `config/decisions.tsv`에 2026-10-16, 시점 `사전`으로 7.1.1–7.1.13을 적었다. 7.1.5는 가족 P와 S1–S3의 Holm을 관문 통과 뒤 실제로 수행한 검정에만 적용하고, m은 그 검정 수이며, 관문 실패 자료는 가족에서 빠지고 gate fail로 보고하며, 두 자료 이상이라는 재현 규칙은 그대로라고 적는다. 7.1.13은 리허설 L2 H4 `RA_default` 귀무 제1종 오류 0.078을 알려진 성질로 적고 바꾸지 않는다고 적는다.

동결 확인은 잠금 해제 파일이 생기기 전에 했다. `sha256sum prereg/manifest.tsv`는 `f19f79c1988d75e7085a3f570b80a96a2e3350b368335ac80b0ae4fb162120ba`였다. manifest 26개 경로의 바이트와 SHA-256이 모두 일치했다. `git diff --stat 8a203ffd -- scripts/confirm prereg`는 비어 있었다. 표는 `results/stage7/preunlock_freeze_check.tsv`에 있다.

## 4. 잠금 해제 (시각, L3 파일 바이트와 SHA-256, 형식 분기와 열 이름)

`config/unlock_approved.txt` 첫 줄은 동결 해시이다. 파일을 쓴 시각은 2026-10-16 16:47:07 +0900, 2026-10-16 07:47:07 UTC이다.

| 파일 | 바이트 | SHA-256 |
|---|---|---|
| `data/locked/GSE295758/GSE295758_rawTranscriptCounts.csv.gz` | 30024282 | `0a2d1ebce4802a5a228059b197386466fdab0f2acdbbdf9bf23cafa80b5cd7f7` |
| `data/locked/GSE295758/GSE295758_combatVstTranscriptCounts.csv.gz` | 83783329 | `57fcd34f4199b413ac930f0084fe83f5955d07205bc125f1d0a5c191ce853d0a` |

실행 3이 기록한 L3 분기는 `transcript_to_gene`, 구분자는 comma, 헤더 열 수는 117이다. 첫 필드는 빈 문자열이다. 나머지 116개 열 이름은 다음과 같다.

`1_D,2_D,3_D,4_D,5_D,6_D,7_D,8_D_batch1,9_D,10_D,11_D,12_D,13_D,14_D,15_D,16_D_batch1,17_D,18_D,19_D,20_D,21_D,22_D,23_D,24_D,25_D,26_D,27_D,28_D,29_D,30_D,31_D,32_D,1_N,2_N,3_N,4_N,5_N,6_N,7_N,8_N_batch1,9_N,10_N,11_N,12_N,13_N,14_N,15_N,16_N_batch1,17_N,18_N,19_N,20_N,21_N,22_N,23_N,24_N,25_N,26_N,27_N,28_N,29_N,30_N,31_N,32_N,33_D,34_D,35_D,36_D,37_D,38_D,39_D,40_D,41_D,42_D,43_D,44_D,45_D,46_D,47_D,48_D,49_D,50_D,51_D,52_D,53_D,54_D,55_D,56_D,57_D,58_D,33_N,34_N,35_N,36_N,37_N,38_N,39_N,40_N,41_N,42_N,43_N,44_N,45_N,46_N,47_N,48_N,49_N,50_N,51_N,52_N,53_N,54_N,55_N,56_N,57_N,58_N`

`batch`가 들어간 열은 `8_D_batch1`, `16_D_batch1`, `8_N_batch1`, `16_N_batch1` 네 개이다. `batch2` 열은 0개이다. GENCODE 38 전사체 식별자로 합친 뒤 유전자 수는 59293, 입력 행 336781, 대응된 행 110352, 1차 표본 116, 피험자 58이다.

## 5. 관문 (자료별 유전자 결과와 판정)

아래 숫자는 실행 3의 `gate.tsv`이다. 실행 2의 같은 파일과 값이 같다. 가족 P 검정은 이 관문 뒤에 끝나지 못했다.

L1 Day 1. 표본은 참가자 22, 샘플 147. 전처리: scale `log`, probes 43758, 샘플 285, 유전자 19541, `n_dropped_multi` 0, `n_dropped_no_symbol` 13253.

| gene | lrt_p | peak_h | pass | dataset_pass |
|---|---|---|---|---|
| PER1 | 1.35631359279823e-11 | 7.39244745511937 | TRUE | TRUE |
| PER3 | 2.79164389620025e-06 | 4.93116485780503 | TRUE | TRUE |

L2 period 1. 표본은 참가자 11, 샘플 70. 시간축은 `wake_ref_h %% 24`이고, 기록된 `clock_time_h`의 최소는 0.199999999999999, 최대는 22.3이다. 전처리: scale `log`, probes 27189, `n_multi` 911, 유전자 20179, 점수화에 넣은 샘플 141. 빛 조건 참가자 수는 bright 3, dim 8이다. 관문 창은 peak_mod가 18 이상이거나 2 이하인 구간이다.

| gene | lrt_p | peak_h | pass | dataset_pass |
|---|---|---|---|---|
| PER1 | 1.33395937684947e-05 | 8.62707238405378 | FALSE | FALSE |
| PER3 | 2.56672604887414e-06 | 6.37573752707371 | FALSE | FALSE |

L3. 표본은 참가자 58, 샘플 116. PER1은 `log2(CPM+1)`의 morning − evening이다.

| gene | mean_delta | p | pass | dataset_pass |
|---|---|---|---|---|
| PER1 | 0.565581638383806 | 0.00672475968762667 | TRUE | TRUE |

## 6. 가족 P

가족 P 표는 없다. 실행 3은 L1 `TAms_chrono`와 `TAms_mort`의 부트스트랩 2000회를 마친 뒤 `RA_default`에서 `missing score RA_default`로 멈췄다. `results/confirm/family_P.tsv`는 없다. Holm q, 재현 판정, 위상·부호 보조 점검 표도 없다.

L1 `ra_scores.tsv`는 `RA_default`, `RA_Peters`, `RA_GTEx` 각 285행이고 value가 비어 있는 행이 855/855이다. L3 `ra_scores.tsv`는 세 시계 각 116행이고 비어 있지 않은 value가 0/348이다.

L1 로그 행렬에서 RNAAgeCalc blood 시그니처와 행 이름이 겹치는 유전자 중, 값이 유한하지 않은 샘플이 있는 유전자 수와 모든 샘플이 비유한인 유전자 수는 다음과 같다. 샘플 285개 모두가 시그니처 안에 비유한 값을 하나 이상 가진다.

| 시그니처 | 행이 있는 유전자 | 일부 샘플이 비유한 | 모든 샘플이 비유한 | 비유한 값이 있는 샘플 | 한 샘플의 최대 비유한 유전자 수 |
|---|---|---|---|---|---|
| DESeq2 | 214 | 201 | 0 | 285 | 43 |
| Peters | 280 | 260 | 0 | 285 | 34 |
| GTExAge | 229 | 212 | 0 | 285 | 58 |

L1 원자료 행렬은 43758행, 열 288개(식별자 열 포함)이고 수치로 읽은 뒤 비유한 칸의 비율은 0.2402063이다. 유전자 로그 행렬의 비유한 칸 비율은 0.1799938이다.

L3 `log2(FPKM+1)`의 비유한 칸 비율은 0.4657211이다. 시그니처와 겹친 유전자 중 모든 샘플이 비유한인 유전자 수는 DESeq2 20/261, Peters 7/283, GTExAge 27/296이다. 이 유전자들은 일부 샘플만 비유한 경우가 0개이다. 샘플 116개 모두가 그 행 때문에 시그니처 안에 비유한 값을 가진다. `count2FPKM` 경고는 유전자 길이 결측 46.5721%이다.

`ra_linear_scores`는 행 이름이 있는 유전자의 계수와 발현을 `crossprod`로 곱한다. 그 곱에 비유한 값이 있으면 그 샘플의 점수는 빈 값이 된다.

## 7. 가족 S1 (보조 시계)

표가 없다. L1과 L3의 `RA_Peters`, `RA_GTEx` value도 비어 있다.

## 8. H2 (자료 × 시계: p_med [구간], 판정)

표가 없다.

## 9. H3 (시계별: 나이 계수 [구간], p, Holm q, 방향)

표가 없다.

## 10. H4 (자료 × 시계: LRT p, Holm q, 조건별 2A와 정점, 평균 차이)

표가 없다. L2는 관문 실패로 이후 검정이 호출되지 않았고 `stopped.tsv`의 status는 `gate fail`이다.

## 11. H5 (시계별: c, mean(Δ) − c [구간], mean(Δ)/c, 판정)

표가 없다.

## 12. 모듈 수준 S4 (L1 확인 표와 확인 수, L2 LRT 표, L3 부호 일치 표, 탐색 BH 표 요약)

표가 없다. 확인 수 / 8은 계산되지 않았다. L1과 L3의 tAge 점수 파일은 시계 16개 행을 가지고 있다. L1 성분 행은 13680, L3 성분 행은 5568이다. 이 점수에 대한 모형 M과 대응 t는 실행되지 않았다.

## 13. 민감도

표가 없다. L3 헤더의 batch-2 열 개수는 0이다. combat 파일은 받아 두었고, 민감도 점수화는 실행 3이 그 단계에 도달하기 전에 멈췄다.

## 14. 독립 검사 (값별 차이)

가족 P 검정이 수행되지 않아 독립 검사를 하지 않았다.

## 15. 버그·재실행·이탈 기록 (어느 실행의 출력이 사전등록 결과인지)

사전등록 결과로 지정할 끝낸 실행은 없다. 실행 3의 관문 표가 마지막으로 계산된 관문이다. 가족 P 이후는 없다.

| 실행 | 멈춘 곳 | 남긴 산출 |
|---|---|---|
| run 1 | `data.table::fread`가 connection을 읽지 못함 | `logs/confirm_run1.log` |
| run 2 | L3 `ra_linear_scores`의 `as.matrix` 인자에서 subscript out of bounds. `count2FPKM`의 `data.frame`이 `1_D`를 `X1_D`로 바꿈 | `results/confirm_run2_partial/` |
| run 3 | L1 `prepare_score`가 `RA_default` 점수를 찾지 못함. 점수 행은 있고 value가 비어 있음 | `results/confirm/` |

`prereg/deviations.tsv`의 7.3.1–7.3.5는 읽기 함수 교체, L3 쉼표 헤더, 3–12단계 연결, batch-2 열 부재, FPKM 열 이름 복원이다. `scripts/confirm/lib_lock.R`의 디스크 해시는 실행 3의 manifest 검사에서 동결 해시와 같았다. 바뀐 manifest 파일은 `scripts/confirm/run_pipeline.R` 하나이고, 실행 3이 기록한 그 파일의 SHA-256은 `d0e9903f19e797644bdcd5b2d306f1e2787dc9d0883bd0938ac9bf986ff8ea4c`이다.

## 16. 판단이 필요한 사항

가족 P를 끝내려면 L1에서 일부 샘플만 비유한인 시그니처 유전자를 선형 점수에 어떻게 넣을지 규칙이 필요하다. 현재 `crossprod`는 그 샘플의 점수를 빈 값으로 만든다. 모든 샘플이 그런 유전자를 하나 이상 가지므로 `RA_default`, `RA_Peters`, `RA_GTEx`의 L1 점수가 전부 빈다. 이 규칙은 사전등록 문장과 현재 함수에 없다. 분석 정의를 바꿔야 하므로 실행을 멈췄다.

L3에서는 비유한 시그니처 유전자가 모든 샘플에서 비유한이다. 유전자 길이가 없어 FPKM이 비는 행이다. 그 행을 없는 유전자로 빼면 나머지 시그니처 유전자는 116개 샘플 모두에서 유한하다. 그 처리도 현재 함수는 하지 않는다. L1의 일부 결측과는 행의 패턴이 다르다.

## 17. 지시서와 다르게 한 것

끝낸 확인 실행이 없어 `TOD7 confirm output run1`이라는 커밋은 만들지 않았다. 멈춘 실행의 산출물은 `a868eee221933a57196936fcacbada4798485c82`에 있다. 독립 검사와 가족 P 이후 표는 없다.

## 18. 산출물 목록, 커밋 해시

| 커밋 | 해시 | 시각 |
|---|---|---|
| TOD7 pre-unlock record | `082ed49ed6f4733c28dcdcce9b502f7fb02c0862` | 2026-10-16 16:46:36 +0900 |
| TOD7 unlock | `dc2bed1a64adb374884128eb91720ffbed3ad3e2` | 2026-10-16 16:49:28 +0900 |
| TOD7 confirm stopped before family P | `a868eee221933a57196936fcacbada4798485c82` | 2026-10-16 17:21:12 +0900 |

실행 3 관문과 빈 RA 점수, L3 열 이름, tAge 점수 파일은 `results/confirm/`에 있다. 실행 2의 부분 산출은 `results/confirm_run2_partial/`에 있다. 로그는 `logs/confirm_run1.log`, `logs/confirm_run2.log`, `logs/confirm_run3.log`이다.
