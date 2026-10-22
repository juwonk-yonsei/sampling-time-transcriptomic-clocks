# TOD 지시서 14b — 그림 용어 통일과 보충 자료 묶음 (프로젝트 TOD, 14b단계)

- 대상: Cursor agent
- 선행 조건: TOD 지시서 14 완료. 14가 아직 시작 전이면 함께 해도 된다.
- 산출 보고서: `reports/TOD_보고서_14b.md` (템플릿은 §5)
- TOD 지시서 1의 §1을 적용한다.
- **새 계산은 하지 않는다.** 이미 있는 표와 그림만 쓴다. 필요한 값이 저장된 표에 없으면 계산하지 않는다. 멈추고 보고서에 적는다.
- 모듈 구성 유전자와 계수는 어떤 파일에도 넣지 않는다(라이선스).
- 커밋: "TOD14b figures and ESM"

---

## 1. 이유 (이해용)

1. 원고에서 "preregistered/preregistration"을 "prespecified"로 바꿨다. 분석 계획을 공개 등록소에 올리지 않았기 때문이다. 그림 속 글자도 같은 말로 맞춘다.
2. GeroScience는 보충 자료를 "Online Resource"라고 부른다. 파일 이름은 ESM_1, ESM_2처럼 매긴다. 파일마다 논문 제목, 학술지 이름, 저자 이름, 교신저자의 소속과 이메일을 넣어야 한다.

---

## 2. 그림 글자

| 그림 | 지금 | 바꿀 말 |
|---|---|---|
| Fig. 1 | Preregistration and locked confirmation | Prespecified plan and locked confirmation |
| Fig. 5a 범례, Fig. 5b 열 묶음 머리글 | Preregistered scoring | Prespecified scoring |
| Fig. 6b 범례 | Preregistered scoring | Prespecified scoring |

- 그림 스크립트 전체에서 "reregist"를 찾는다. 남은 곳이 없게 고치고, 찾은 곳을 모두 보고서에 적는다.
- 그 밖에는 아무것도 바꾸지 않는다.
- 지시서 12 §3과 지시서 13 §2의 점검을 다시 한다. 투고용 사본(`figures/submission/`)과 검토 묶음도 다시 만든다.

---

## 3. 보충 자료 묶음 (`figures/submission/esm/`)

### 3.1 공통 표지

두 파일 모두 맨 앞(PDF는 1쪽, XLSX는 첫 시트 "Cover")에 아래 내용을 둔다. 대괄호는 자리표시로 그대로 둔다.

- Title: Time of blood sampling biases transcriptomic aging clocks: development, prespecified confirmation, and practical guidance
- Journal: GeroScience
- Authors: [Author names]
- Corresponding author: [Affiliation], [e-mail]

### 3.2 ESM_1.pdf (Online Resource 1, 보충 그림)

- 표지 다음에 Fig. S1, S2, S3를 한 쪽에 하나씩 둔다. 각 그림 아래에 설명문을 둔다.
- 설명문은 아래 문장을 그대로 쓴다. 굵게는 "Fig. S1"처럼 그림 번호만이다.
  - **Fig. S1** Effect sizes and composition mediation for the secondary clocks (RNAAgeCalc-Peters, RNAAgeCalc-GTEx), as in Fig. 3a, b and Fig. 4a
  - **Fig. S2** Bias simulation for the secondary clocks, as in Fig. 7; in d, the laboratory lines of the two clocks coincide because their biases at s = 4 h have equal size (0.19 of the total SD) and opposite sign
  - **Fig. S3** Mitigation for the secondary clocks, as in Fig. 8; the time correction c was defined for the primary clocks only
- 그림은 지시서 14(와 이번 §2)를 반영한 최종 PDF를 그대로 넣는다(벡터 유지, 글꼴 포함).

### 3.3 ESM_2.xlsx (Online Resource 2, 보충 표)

- 시트: Cover, S1, S2, S3, S4, S5, S6. S3이 너무 크면 S3a–S3f로 나눠도 된다(아래 표).
- 각 시트의 1행에는 아래 표 제목을 그대로 쓴다. 표는 3행부터 시작한다.

| 시트 | 표 제목 (1행) | 내용과 출처 |
|---|---|---|
| S1 | Table S1. Deviations D1–D12: content, timing relative to results, and reason. | 지시서 8 Table S1 (`results/stage8/tables/`) |
| S2 | Table S2. Sensitivity analyses: prespecified (flagged-sample exclusion, random slopes, batch-corrected counts, comprehensive annotation) and deviation-defined (complete-case, kNN and participant-mean versions; imputed composition). | 지시서 8 Table S2 |
| S3 (또는 S3a–S3f) | Table S3. Post hoc analyses: GSE122725 time-axis correction and exploratory H1 and H4, module outcomes under participant-mean scoring, masking experiment, simulation parameters and results, mitigation. | 지시서 8 Table S3 (a: GSE122725, b: 모듈), 마스킹 실험 결과(c), 지시서 9 시뮬레이션 매개변수(d)와 결과(e: 설계 1–3 요약), 지시서 9 완화 표(f) |
| S4 | Table S4. Frozen development parameters: peak times, amplitudes, variance components, expected morning−evening effects and module outcomes. | `prereg/frozen/`의 F1, F2, F3 (모듈 유전자·계수 제외) |
| S5 | Table S5. Software versions, model-file checksums and the frozen manifest. | 동결 환경 표, tAge 모델 파일 체크섬, manifest 파일 목록과 SHA-256 |
| S6 | Table S6. Prospective power of the confirmatory tests. | 사전 검정력 시뮬레이션 결과 |

- 나눈 시트(S3a–S3f)의 1행에는 "Table S3a. GSE122725 time-axis correction and exploratory H1 and H4"처럼 짧은 부제목을 쓴다.
- 열 이름과 값에서 내부 표기를 원고 표기로 바꾼다.

| 내부 표기 | 원고 표기 |
|---|---|
| TAms_chrono / TAms_mort | tAge-chronological / tAge-mortality |
| RA_default / RA_Peters / RA_GTEx | RNAAgeCalc-DESeq2 / RNAAgeCalc-Peters / RNAAgeCalc-GTEx |
| L1 / L2 / L3 | GSE48113 / GSE122725 / GSE295758 |
| CMP_TAms_chrono_orange 등 | chrono-orange 등 (Fig. 5와 같은 짧은 이름) |
| TM_turquoise_mortality | TM-turquoise |
| frozen / D11 | prespecified scoring / participant-mean imputation |

- 숫자는 저장값을 그대로 둔다(반올림하지 않는다).
- 열 이름 아래에 단위가 필요하면 괄호로 적는다(예: "peak (h)").

---

## 4. 점검

1. `pdffonts esm/ESM_1.pdf`: 글꼴이 모두 포함(emb yes)되어야 한다.
2. ESM_1.pdf 쪽 수와 각 쪽의 그림 번호.
3. ESM_2.xlsx: 시트 이름과 각 시트의 행·열 수, 1행 제목.
4. 두 파일 모두에서 "reregist", "TAms_", "RA_", "CMP_", "L1", "L2", "L3"(열 이름과 값) 검색 결과 수. 모두 0이어야 한다.
5. 모듈 구성 유전자 이름이나 계수가 들어간 열이 없음을 확인한 방법을 적는다.

---

## 5. 보고서 템플릿 (`reports/TOD_보고서_14b.md`)

~~~markdown
# TOD 보고서 14b — 그림 용어 통일과 보충 자료 묶음

## 1. 바꾼 그림과 "reregist" 검색 결과
## 2. 그림 점검 (지시서 12 §3, 13 §2 표)
## 3. ESM 파일 (경로, 크기, 쪽·시트 구성, 표별 출처 파일)
## 4. ESM 점검 (§4의 1–5)
## 5. 지시서와 다르게 한 것 (빠짐없이)
## 6. 커밋 해시
~~~
