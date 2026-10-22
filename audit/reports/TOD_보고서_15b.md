# TOD 보고서 15b — 공개 저장소 준비 계속

| 항목 | 결과 |
|---|---|
| 잠금 자료 발현 파일 수 | 70 (유전자×샘플 발현 행렬 19, 표본별 값 파일 51) |
| 그중 가장 이른 첫 커밋 시각 (KST) | 2026-10-16 17:21:12 (`a868eee221933a57196936fcacbada4798485c82`) |
| 잠금 해제 커밋 시각 (KST) | 2026-10-16 16:49:28 (`dc2bed1a64adb374884128eb91720ffbed3ad3e2`) |
| 판정 (모두 해제 뒤 / 앞선 파일 있음) | 모두 해제 뒤 |

## 1. 점검 결과 (§2: 큰 파일, 기록 전체, 공개 금지 내용, 멈춘 경우 그 이유)

멈추지 않았다.

지시서 15b §1에 따라 보고서 15의 유전자×샘플 행렬은 추적에서 빼지 않았고 기록도 고치지 않았다.

잠금 순서 확인은 `scripts/check_lock_order.py`로 했다. `/usr/bin/time -v` 기록은 경과 3:18.88, 최대 RSS 2949516 kB, 종료 상태 0이다.

- 잠금 자료 샘플 ID는 저장된 확인 행렬 머리에서 읽었다. GSE48113 285개, GSE122725 141개, GSE295758 116개.
- 대상: `git rev-list --all`의 커밋 65개, 그 트리의 blob 933개. GSE295758 ID는 앞에 `X`가 붙은 형태(`X1_D`)도 셌다.
- 한 자료의 ID가 10개 이상 나온 경로는 96개다.

| 분류 | 파일 수 | 가장 이른 첫 커밋 (`--diff-filter=A`) | 모두 해제 뒤 |
|---|---|---|---|
| 발현 행렬 | 19 | 2026-10-16 17:21:12 `a868eee` | True |
| 표본별 값 (점수, 성분, 플래그, 그림 원자료) | 51 | 2026-10-16 17:21:12 `a868eee` | True |
| 합성 리허설 (R2) | 2 | 2026-10-16 16:20:54 `967723e` | False |
| 표본 주석 | 23 | 2026-10-16 08:51:22 `541a64e` | False |
| 보고서 | 1 | 2026-10-16 17:22:26 `126351a` | True |

| 첫 커밋 | 시각 (KST) | 메시지 | 발현 행렬 | 표본별 값 파일 |
|---|---|---|---|---|
| `a868eee` | 2026-10-16 17:21:12 | TOD7 confirm stopped before family P | 8 | 17 |
| `7ffde8e` | 2026-10-16 17:29:12 | TOD7b deviation plan (before results) | 4 | 9 |
| `eac6837` | 2026-10-16 18:54:17 | TOD7b report | 5 | 10 |
| `e56a1c4` | 2026-10-17 11:12:06 | TOD7c confirm output run5 | 2 | 10 |
| `f6d8879` | 2026-10-17 11:20:50 | TOD7c posthoc L2 | 0 | 2 |
| `122d987` | 2026-10-19 18:33:36 | TOD10 figures | 0 | 3 |

잠금 해제 커밋까지(포함)의 커밋은 42개, blob은 493개다. ID가 1개 이상 나온 blob은 4개다.

| blob | 경로 | ID 수 |
|---|---|---|
| `debbd8162721` | `results/stage2/GSE122725_wake_ref.tsv` | GSE122725 141 |
| `82fa91d0d472` | `results/stage2/rnaseq_headers/GSE295758_header.txt`, `results/stage3/GSE295758_bang_lines.txt` | GSE295758 116 |
| `055091f1a5b5` | `results/stage6/rehearsal/R2/gene_log2.tsv.gz` | GSE122725 141 |
| `1eb19a2c411a` | `results/stage6/rehearsal/R2/samples_input.tsv` | GSE122725 141 |

발현 행렬과 표본별 값 파일 70개의 blob 중 해제 전 트리에 있는 것은 0개다.

`GSE295758_header.txt`의 비어 있지 않은 줄 74개 중 `!`로 시작하지 않는 줄은 0개다.

R2. `scripts/43_partD_pipeline.R`의 `run_r2`는 `set.seed(20261015)`과 `stats::rnorm`으로 행렬을 만든다. 같은 코드와 시드로 `/tmp/r2check/`에 다시 만든 표와 커밋된 `results/stage6/rehearsal/R2/gene_log2.tsv.gz`의 압축 해제 SHA-256은 둘 다 `cbcf0c2db76c092224ee0177d729d8e5a58258c87c12d9df1021721101976887`이다. 셸 명령 경과 7213 ms. 다시 만든 파일은 커밋하지 않았다.

`git log --all --diff-filter=A --follow`는 발현 행렬 18개(R2 포함)의 첫 커밋을 `553a2bf` 2026-10-16 15:37:46 (TOD6 Part B)로 냈다. `553a2bf`의 name-status에서 추가된 파일은 `prereg/` 아래 10개와 `scripts/35_partB_tables.py`이고, 행렬 경로는 없다. `results/confirm/L2/log2_genes.tsv`의 `--follow --name-status`는 `a868eee`에서 `C068 prereg/ref/gencode.v38.basic.annotation.gtf.gz → results/confirm_run2_partial/L2/log2_genes.tsv`를 보이고, 그다음 GTF의 `e0d37b6`(M)과 `553a2bf`(A)로 이어진다. `results/confirm/L1/log2_genes.tsv.gz`, `results/confirm/L2/log2_genes.tsv`, `results/confirm/L3/gene_counts.tsv`는 `553a2bf`, `082ed49`, `dc2bed1`의 트리에 없다. 두 열은 `audit/locked_expression_first_commits.tsv`에 같이 있다.

GTF. `553a2bf`의 `prereg/ref/gencode.v38.basic.annotation.gtf.gz` blob `f7599230268637b058e8210145df7d4c28f027b9`(24158208바이트)는 `gzip -t`에서 "unexpected end of file"이다. `e0d37b6`(2026-10-16 16:02:26)에서 blob `ab6059d515d869c10b9c2619b1978def40c3da26`로 바뀌었고, 현재 파일은 `gzip -t`를 통과한다. 스크립트는 끊긴 스트림을 풀 수 있는 데까지 풀어 검사한다. 그 blob에서 잠금 ID는 0개다.

공개될 유전자×샘플 표는 26개, 796343221바이트다. 보고서 15의 24개에 `results/stage6/rehearsal/R3/inputs/counts_gene_column.tsv.gz`와 `counts_transcript_column.tsv.gz`가 더해졌다. 두 파일의 첫 열은 `gene_name`, `transcript_id`다.

## 2. .gitignore와 추적에서 뺀 파일

`git rm --cached`를 쓰지 않았다. 추적에서 뺀 파일은 0개다.

더한 규칙:

- 원 GEO 자료: `*_series_matrix.txt`, `*_series_matrix.txt.gz`, `*_RAW.tar`, `*.CEL`, `*.CEL.gz`, `*.cel`, `*.cel.gz`, `*.fastq`, `*.fastq.gz`, `*.fq`, `*.fq.gz`, `*_rawTranscriptCounts.csv.gz`, `*_combatVstTranscriptCounts.csv.gz`, `*_TPM_processed_data.txt.gz`, `GPL*_platform.txt`, `GPL*.annot.gz`, `GPL*_family.soft.gz`
- tAge 모델 파일: `*.pkl`, `*module clocks.zip`, `*module clocks/`
- 캐시와 중간 파일: `.cache/`, `logs/`, `results/stage7c/l1_build/`, `results/stage2/partE_search_counts_partial.tsv`, `__pycache__/`, `*.pyc`, `.Rhistory`, `.RData`, `.Rproj.user/`, `.ipynb_checkpoints/`
- 개인 설정: `.vscode/`, `.idea/`, `.cursor/`, `.cursorignore`, `.env`, `.env.*`, `*.swp`, `*~`, `.DS_Store`
- 원고와 투고 파일: `/manuscript/`, `/submission/`

기존 두 줄 `data/`, `results/stage3/expr/`는 그대로 두었다.

`git ls-files -ci --exclude-standard`는 0개다. 처음 쓴 `manuscript/`, `submission/`(앞 `/` 없음)은 `figures/manuscript/`와 `figures/submission/`의 추적 파일 104개에 걸렸다. 앞에 `/`를 붙인 뒤 0개다.

무시되는 추적 밖 파일: `data/` 30, `logs/` 35, `manuscript/` 1, `results/stage2/` 1, `results/stage3/` 5, `results/stage7c/` 2, `scripts/__pycache__/` 6, `scripts/confirm/__pycache__/` 1.

## 3. 새 파일 목록 (README, LICENSE, CITATION.cff, .zenodo.json, verify_manifest.py, audit/)

| 파일 | 바이트 |
|---|---|
| `README.md` | 32585 |
| `LICENSE` | 1066 |
| `CITATION.cff` | 676 |
| `.zenodo.json` | 2809 |
| `scripts/verify_manifest.py` | 5014 |
| `scripts/00_download_data.sh` | 3361 |
| `scripts/check_lock_order.py` | 9042 |
| `environment/R_packages.tsv` | 9560 |
| `environment/python_analysis_requirements.txt` | 7162 |
| `environment/python_clock_env_requirements.txt` | 228 |
| `audit/locked_expression_first_commits.tsv` | 21165 |
| `audit/instructions/` (19개) | 210850 |
| `audit/reports/` (보고서 1–15, 18개) | 770196 |
| `audit/reports/TOD_보고서_15b.md`, `audit/INDEX.md` | 이 보고서 뒤에 만듦 |

- `README.md`: 지시서 15 §4.1의 1–8 순서. 15b §2의 "Processed data included" 절(26행 표, 중단 실행 폴더 설명, R1–R3 출처), "Lock order" 절(문장, 커밋별 표, 70행 표)을 넣었다.
- `LICENSE`: MIT, 저작권자 `[Authors]`.
- `CITATION.cff`: cff-version 1.2.0, `version: 1.0.0`, `license: MIT`, 키워드 6개(원고 키워드 줄).
- `.zenodo.json`: `upload_type` software, 원고 제목, README 요약 문단, MIT, 키워드 6개. `related_identifiers`: `10.5281/zenodo.22166800` requires, GEO 6개 references, 논문 `[Article DOI]` isSupplementTo.
- `scripts/verify_manifest.py`: manifest 26개 파일의 SHA-256을 다시 계산한다. 출력은 `match`, `mismatch`, `missing (not redistributed)`. 불일치 파일은 manifest를 넣은 커밋의 사본 해시와 그 뒤 바꾼 커밋을 같이 적는다. `--tage-dir`는 tAge 파일 7개(md5)와 모듈 압축 안 파일 78개(SHA-256)를 확인한다.
- `scripts/00_download_data.sh`: `config/raw_checksums.tsv`의 15개 파일. 자동 14개, 수동 1개(Peters 2015 표, PMC).
- `scripts/check_lock_order.py`: §1의 잠금 순서 확인.
- `environment/`: 2026-10-22에 만든 목록. R 패키지 195개(`r_env`), 분석 Python 392줄, 시계 Python 14줄. `results/stage4/module_clocks/venv_freeze.txt`의 `installed_pins` 5개(joblib 1.6.0, numpy 1.26.4, pandas 2.2.3, scikit-learn 1.3.2, scipy 1.13.1)는 시계 Python 목록과 같다. 시계 Python 목록의 `clockspec==0.1.0`, `scmodage==0.2.0`은 원격 저장소가 없는 로컬 git 설치다.
- `audit/`: 지시서 1–15b 19개와 보고서 1–15 18개를 이름을 바꾸지 않고 복사했다. 저장소에 없는 지시서와 보고서는 0개다.

## 4. 기록 대조표 (README에 넣은 표 그대로)

| Event | Commit | Date and time (KST) | Message |
|---|---|---|---|
| Rehearsal of the pipeline on development and synthetic data (R1–R3) | `967723e47b76d4094634bd52d9b056e56d9d210a` | 2026-10-16 16:20:54 | TOD6 Part D |
| Analysis plan freeze (manifest SHA-256 f19f79c1…) | `8a203ffdc5a5f1e933cd4f64154c0d5972480e32` | 2026-10-16 16:21:17 | TOD6 freeze |
| Pre-unlock record | `082ed49ed6f4733c28dcdcce9b502f7fb02c0862` | 2026-10-16 16:46:36 | TOD7 pre-unlock record |
| Unlock | `dc2bed1a64adb374884128eb91720ffbed3ad3e2` | 2026-10-16 16:49:28 | TOD7 unlock |
| Runs 1–3 stopped; bug fixes 7.3.1–7.3.5 recorded (later confirmed as D7) | `a868eee221933a57196936fcacbada4798485c82` | 2026-10-16 17:21:12 | TOD7 confirm stopped before family P |
| Deviations D1–D8 recorded before family-P results | `7ffde8edc6f8ffbb395482fcf430ad91457d16e1` | 2026-10-16 17:29:12 | TOD7b deviation plan (before results) |
| Code for D1–D8 | `8e38cffef0cc475f2267b176e8aaac9705a98445` | 2026-10-16 18:34:56 | TOD7b code fixes |
| Run 4 partial output | `eac683709ee4a63ce302c01b05fd3d629b02801f` | 2026-10-16 18:54:17 | TOD7b report |
| Deviations D9–D12 recorded before H2 and later results | `6904715f53f7d694ab5f710ec53b9e04ce07da34` | 2026-10-17 09:49:01 | TOD7c deviation plan (before H2+ results) |
| Code for D9–D11 | `bd728bacd95ee487679d0c0380cc99f00d780b88` | 2026-10-17 10:52:31 | TOD7c code fixes |
| Confirmation run 5 (prespecified output) | `e56a1c489a74453ca33a832835fdaa312e3e8371` | 2026-10-17 11:12:06 | TOD7c confirm output run5 |
| Post hoc, instruction 7c: GSE122725 time axis | `f6d8879c19bc710a7353dd9835dedee5efce5010` | 2026-10-17 11:20:50 | TOD7c posthoc L2 |
| Instruction 7c report and independent recomputation | `248b37ec68635bc6527a8da092f29daf70f680b2` | 2026-10-17 11:23:42 | TOD7c report |
| Post hoc, instruction 8: Part A | `67d4c2d7b3570f2b72aba06242706057cc4e635d` | 2026-10-19 11:09:06 | TOD8 Part A |
| Post hoc, instruction 8: Part B | `9317da205567e399907d98a2ae1ef21e50f4aab8` | 2026-10-19 12:03:27 | TOD8 Part B |
| Post hoc, instruction 8: report | `b97d159b3f1e86ce55c6020c8b91b29629f6c7da` | 2026-10-19 12:06:57 | TOD8 report |
| Post hoc, instruction 9: simulation and mitigation | `50e5f00580c580d77fa67215b9dd097528ca8268` | 2026-10-19 16:33:32 | TOD9 Part A–C |
| Post hoc, instruction 9: figures | `c31ddf74158a59a282e8f3546780b50081d74d55` | 2026-10-19 16:33:39 | TOD9 figures |
| Post hoc, instruction 9: report | `226949bc5cbe49c74dfd444f8d33bb680bcfd09b` | 2026-10-19 16:35:21 | TOD9 report |
| Instruction 10: figures | `122d9876599b2651709a8ebc718e75240e81b26b` | 2026-10-19 18:33:36 | TOD10 figures |
| Instruction 11: figures | `653eb14b0b751b45868b7132466076995ee08ac5` | 2026-10-21 16:28:52 | TOD11 figures |
| Instruction 12: figures | `64b482b816479fb03633fffdce9575f1527c6404` | 2026-10-21 18:00:58 | TOD12 figures |
| Instruction 13: figures | `694678a0913bcccb66a7979c97c68259ba7b0dbf` | 2026-10-21 18:33:01 | TOD13 figures |
| Instruction 14: figures | `927a5f35f6dd1ce4db9387fc2a93ce451fc338e6` | 2026-10-21 20:02:14 | TOD14 figures |
| Instruction 14b: figures and supplementary files | `d2dbad34eb5c21e03a2b6581da668fe99cfbf99c` | 2026-10-21 21:12:11 | TOD14b figures and ESM |
| Instruction 15: release check | `85486cce0e97cf0eadd7df88254b6e086d11d8bf` | 2026-10-22 08:49:51 | TOD15 release prep |

## 5. 최종 점검 (§5의 1–4)

1. 공개 금지 내용 검사 (공개될 파일과 기록의 blob 933개)
   - tAge 모델 파일: `zenodo_v2_files.tsv`의 md5와 `module_archive_members.tsv`의 SHA-256에 맞는 파일은 공개될 파일 0개, 기록 blob 0개. 이름이 `.pkl`이거나 `module clocks`를 포함하는 경로는 0개.
   - 유전자 열과 `coef`, `coefficient`, `weight`, `beta` 열이 같이 있는 표: 0개.
   - 비밀 정보 패턴(`ghp_`, `github_pat_`, `AKIA`, `PRIVATE KEY`, `api_key`, `password=`, `secret=`, `zenodo_token`, `ACCESS_TOKEN`, `Bearer `): 걸린 파일은 `reports/TOD_보고서_15.md`와 그 사본 `audit/reports/TOD_보고서_15.md`이고, 걸린 문자열은 보고서 15의 §1에 적힌 검색 패턴 목록이다. 기록 blob에서도 같은 보고서뿐이다. 이름이 `.env`, `id_rsa`, `.pem`, `credentials`인 경로는 0개.
   - 원 GEO 파일 이름(`series_matrix`, `fastq`, `.fq`, `.cel`): 공개될 파일 0개, 기록 0개.
   - 크기: 100 MB 넘는 파일 0개, 기록 blob 0개. 50 MB 넘는 파일 2개(79808508바이트, 같은 blob).
2. 공개될 파일(`git ls-files -co --exclude-standard`): 948개, 1067791734바이트. 이 보고서, 그 `audit/` 사본, `audit/INDEX.md` 3개는 이 수에 들어 있지 않다.
3. `python scripts/verify_manifest.py`
   - `manifest_sha256` `f19f79c1988d75e7085a3f570b80a96a2e3350b368335ac80b0ae4fb162120ba`, `config/unlock_approved.txt`와 같음.
   - match 25, mismatch 1, missing (not redistributed) 0. 종료 코드 1.
   - mismatch: `scripts/confirm/run_pipeline.R`. 동결 커밋 `8a203ff`의 사본은 manifest와 같다. 동결 뒤 바꾼 커밋은 `a868eee` 하나다.
   - `--tage-dir /HDD8T1/WORK/GeroSci/zenodo_22166800_v2`: match 85, mismatch 0, missing 0. 빈 폴더: missing 85.
   - `bash scripts/00_download_data.sh`(기존 파일): match 15. 내려받기 URL 14개의 `curl -I -L` 응답은 모두 200.
4. 자리표시

| 파일 | 줄 | 내용 |
|---|---|---|
| `README.md` | 5 | `[Authors]`, `[Year]`, `[DOI]` (논문) |
| `README.md` | 6 | `[DOI]` (Zenodo) |
| `README.md` | 194 | `[clock-specificity repository URL]`, `[sc-module-aging repository URL]` |
| `README.md` | 279 | `[Confirm before submission.]` (원고 2.11 그대로) |
| `LICENSE` | 3 | `[Authors]` |
| `CITATION.cff` | 6 | `[YYYY-MM-DD]` |
| `CITATION.cff` | 7 | `[https://github.com/OWNER/REPOSITORY]` |
| `CITATION.cff` | 16–19 | `[Family name]`, `[Given name]`, `[Affiliation]`, `[https://orcid.org/0000-0000-0000-0000]` |
| `.zenodo.json` | 6–8 | `[Family name, Given name]`, `[Affiliation]`, `[0000-0000-0000-0000]` |
| `.zenodo.json` | 66 | `[Article DOI]` |

## 6. 지시서와 다르게 한 것

- 잠금 순서 판정은 경로 기준 `git log --all --diff-filter=A`의 첫 커밋과, 해제 전 트리 42개의 blob 검사로 했다. `--follow`는 18개 행렬에서 GTF 기록으로 이어졌다(§1). `--follow` 결과는 TSV에 같이 두었다. 멈추지 않았다.
- "발현 값이 들어간 파일"에 유전자×샘플 행렬 19개와 함께 표본별 값 파일 51개(시계 점수, 모듈 성분, `flags_3mad`의 `median_log`, 그림 원자료)를 넣었다. 표본 주석 23개와 보고서 1개는 표에서 뺐다.
- R2 합성 행렬은 해제 전 커밋에 있다. 같은 코드와 시드로 다시 만들어 내용 해시를 비교했고(§1), `check_lock_order.py`에서 `results/stage6/rehearsal/` 경로를 합성 리허설로 분류했다.
- README의 가공 행렬 표는 26개다. 보고서 15의 24개에 R3 입력 2개를 더했다.
- 보고서는 `reports/`에 그대로 두고 `audit/reports/`에 같은 파일을 복사했다.
- 잠금 파일(`renv.lock`, `requirements.txt`)이 없어 `environment/`에 지금 환경의 목록을 만들었다. 분석 Python 목록에는 지시서 14b에서 설치한 pypdf 6.21.1이 들어 있다.
- 지시서에 이름이 없는 `scripts/00_download_data.sh`와 `scripts/check_lock_order.py`를 만들었다. 분석 코드는 고치지 않았다.
- README의 AI 사용 절은 원고 2.11 문단을 `[Confirm before submission.]`까지 그대로 옮겼다.
- 원고 초안 폴더 `/manuscript/`와 `/submission/`을 무시 목록에 넣었다.
- `.zenodo.json`의 논문 항목은 자리표시 `[Article DOI]`로 두었다.
- README 기록 대조표에 7c 사후(`f6d8879`), 7c 보고서(`248b37e`), 지시서 15(`85486cc`) 행을 더했다.
- 연구자 요청(대화)에 따라 이 커밋 뒤 저장소 밖에 GitHub용 경로와 투고 파일 보관 경로를 만든다. 그 결과는 이 보고서의 수치에 들어 있지 않다.
- push, 태그, 릴리스는 하지 않았다.

## 7. 커밋 해시

이 보고서를 넣는 커밋
