# TOD 보고서 15 — 공개 저장소 준비

## 1. 점검 결과 (§2: 큰 파일, 기록 전체, 공개 금지 내용, 멈춘 경우 그 이유)

§3으로 가지 않았다.

추적 파일은 899개, 바이트 합은 1066720917이다.

50 MB(52428800바이트)를 넘는 추적 파일은 2개다. 둘 다 79808508바이트이고 blob은 `8191b3e4d06a5f85689190afd1a11b814a496fc7` 하나다.

| 바이트 | 경로 |
| --- | --- |
| 79808508 | results/confirm/L1/log2_genes_d1_imputed.tsv |
| 79808508 | results/confirm_run4_partial/L1/log2_genes_d1_imputed.tsv |

100 MB를 넘는 추적 파일은 0개다. `git rev-list --objects --all`의 blob 중 100 MB를 넘는 것은 0개다. 50 MB를 넘는 blob은 위 1개다.

tAge 모델 파일. `results/stage5/code_vs_zenodo_md5.tsv`의 파일 이름 64개에 해당하는 `*.pkl` blob은 현재 트리와 기록에 0개다. 그 파일 이름과 md5 문자열은 표 안에 있다.

모듈 구성 유전자와 계수. 유전자 ID 열과 모듈·가중치·계수 값이 한 표에 같이 있는 파일은 0개다. `results/stage2/rnaage_coef_genes.tsv`의 열은 signature, gene이다. 계수 값은 없다.

비밀 정보. 현재 트리에서 `ghp_`, `github_pat_`, `AKIA`, `PRIVATE KEY`, `api_key`, `password=`, `secret=` 검색은 0건이다. 파일 이름 `.env`, `id_rsa`, `.pem`, `credentials`는 현재 트리와 기록에 0개다.

원 GEO 자료. 이름에 `series_matrix`, `fastq`, `.cel`이 있는 파일은 현재 트리와 기록에 0개다. 유전자×샘플 수치 행렬은 24개이고 기록에 있다. 이 때문에 멈췄다. 기록을 지우지 않았다.

| 처음 넣은 커밋 | 시각 (KST) | 메시지 | 파일 수 |
| --- | --- | --- | --- |
| 967723e47b76d4094634bd52d9b056e56d9d210a | 2026-10-16 16:20:54 +0900 | TOD6 Part D | 5 |
| a868eee221933a57196936fcacbada4798485c82 | 2026-10-16 17:21:12 +0900 | TOD7 confirm stopped before family P | 8 |
| 7ffde8edc6f8ffbb395482fcf430ad91457d16e1 | 2026-10-16 17:29:12 +0900 | TOD7b deviation plan (before results) | 4 |
| eac683709ee4a63ce302c01b05fd3d629b02801f | 2026-10-16 18:54:17 +0900 | TOD7b report | 5 |
| e56a1c489a74453ca33a832835fdaa312e3e8371 | 2026-10-17 11:12:06 +0900 | TOD7c confirm output run5 | 2 |

967723e47b76d4094634bd52d9b056e56d9d210a: `results/stage6/rehearsal/R1/gene_log2.tsv.gz`, `results/stage6/rehearsal/R1/scores_input.tsv.gz`, `results/stage6/rehearsal/R2/gene_log2.tsv.gz`, `results/stage6/rehearsal/R3/gene_column/gene_counts.tsv.gz`, `results/stage6/rehearsal/R3/transcript_column/gene_counts.tsv.gz`.

a868eee221933a57196936fcacbada4798485c82: `results/confirm/L1/log2_genes.tsv.gz`, `results/confirm/L2/log2_genes.tsv`, `results/confirm/L3/gene_counts.tsv`, `results/confirm/L3/log2_fpkm.tsv`, `results/confirm_run2_partial/L1/log2_genes.tsv.gz`, `results/confirm_run2_partial/L2/log2_genes.tsv`, `results/confirm_run2_partial/L3/gene_counts.tsv`, `results/confirm_run2_partial/L3/log2_fpkm.tsv`.

7ffde8edc6f8ffbb395482fcf430ad91457d16e1: `results/confirm_run3_partial/L1/log2_genes.tsv.gz`, `results/confirm_run3_partial/L2/log2_genes.tsv`, `results/confirm_run3_partial/L3/gene_counts.tsv`, `results/confirm_run3_partial/L3/log2_fpkm.tsv`.

eac683709ee4a63ce302c01b05fd3d629b02801f: `results/confirm_run4_partial/L1/log2_genes.tsv.gz`, `results/confirm_run4_partial/L1/log2_genes_d1_imputed.tsv`, `results/confirm_run4_partial/L2/log2_genes.tsv`, `results/confirm_run4_partial/L3/gene_counts.tsv`, `results/confirm_run4_partial/L3/log2_fpkm.tsv`.

e56a1c489a74453ca33a832835fdaa312e3e8371: `results/confirm/L1/log2_genes_d1_imputed.tsv`, `results/confirm/L3/combat_gene_counts.tsv`.

## 2. .gitignore와 추적에서 뺀 파일

하지 않았다. `git rm --cached`를 쓰지 않았다.

## 3. 새 파일 목록 (README, LICENSE, CITATION.cff, .zenodo.json, verify_manifest.py, audit/)

만들지 않았다.

## 4. 기록 대조표 (README에 넣은 표 그대로)

README를 만들지 않았다. 표 없음.

## 5. 최종 점검 (§5의 1–4)

하지 않았다. `scripts/verify_manifest.py`는 없다.

## 6. 지시서와 다르게 한 것

- §2에서 멈춰 §3–§6의 파일 작업을 하지 않았다.
- 유전자×샘플 표를 수치 행렬로 보았다. `series_matrix`, FASTQ, CEL 이름의 파일은 0개다.
- 커밋 해시 칸은 커밋 뒤에 알 수 있으므로 아래 문장으로 두었다.

## 7. 커밋 해시

이 보고서를 넣는 커밋
