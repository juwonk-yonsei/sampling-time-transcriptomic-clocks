# TOD 지시서 11 — 원고 그림 최종 수정 (프로젝트 TOD, 11단계)

- 대상: Cursor agent
- 선행 조건: TOD 지시서 10 완료 (커밋 "TOD10 figures")
- 산출 보고서: `reports/TOD_보고서_11.md` (템플릿은 §6)
- TOD 지시서 1의 §1을 적용한다.
- **새 계산은 하지 않는다.** 저장된 표만 읽어서 그린다. 필요한 값이 저장된 표에 없으면 계산하지 않는다. 멈추고 보고서에 적는다. `results/`는 고치지 않는다.
- 출력: `figures/manuscript/`를 덮어쓴다(이전 판은 git에 남는다). 투고용 사본은 `figures/submission/`에 둔다(§2.12).
- 원고 범례는 이 지시서대로 그린 그림을 기준으로 미리 고쳐 두었다. **지시서와 다르게 그린 곳은 하나도 빠짐없이 보고서 §6에 적는다.**
- 순서: §2(공통 규칙) → §3(본문 그림) → §4(보충 그림) → §5(점검) → 보고서 → 커밋 "TOD11 figures"

---

## 1. 이번 수정의 이유 (이해용)

1. 투고할 학술지(GeroScience, Springer)의 그림 지침에 맞춘다.
   - 패널 기호는 소문자 a, b, c를 쓴다.
   - 글꼴은 Helvetica나 Arial 계열, 크기는 8–12 pt다.
   - 벡터 파일은 EPS를 권장하며 글꼴을 포함한다. 비트맵 선 그림은 1200 dpi다.
   - 선 굵기는 0.3 pt 이상이다. 그림 안에 제목과 설명을 두지 않는다.
2. Fig. 2의 p 값과 정점이 저장값 전체 자릿수로 찍혀 있다(예: `p = 2.96936394371969e-10`).
3. 글자가 겹치거나 잘린 곳이 있다.
   - Fig. 3b: 왼쪽 시계 이름이 잘렸다.
   - Fig. 3c: 원판 눈금 6과 18이 겹친다.
   - Fig. 5a: 범례가 × 표식을 덮고, "chrono-orange"가 축 밖으로 나갔다.
   - Fig. 6: a와 b 사이가 겹친다.
   - Fig. 8a: 비율 0.67과 0.50이 겹친다. Fig. S3a: 0.63과 0.62가 겹친다.
   - Fig. S3b: "composition"이 잘렸다.
4. 한 그림 안에서 같은 색이 두 가지 뜻을 가진다.
   - Fig. 6의 분홍은 a에서 GSE122725, b에서 D11, c에서 RNAAgeCalc-Peters다.
   - Fig. 8b와 S3b는 시계 색을 쓴다. 그런데 범례는 같은 색을 자료(GSE56931, GSE113883, GSE39445)로 적는다.
5. 독자가 모르는 내부 용어가 그림에 보인다: D2, D11, frozen, GEO axis, relclocktime − 8, n_entered, before/after, p_med, 범례 제목 "post hoc". 원판 범례의 "r = 1.0"은 상관계수로 읽힌다.
6. Fig. 3c 원판은 자료를 색으로만 구분하고, 구간 호가 서로 겹친다. 정점이 모두 9–15시에 있으므로 0–24시 직선 축으로 그리면 감싸기 문제 없이 읽힌다.
7. Fig. 4b의 12시→12시 축은 NK 세포, 세포독성 림프구, 호중구의 정점(9.4–14.7시)을 축 양 끝으로 가른다. 지시서 10에서 정한 이 축을 18시→18시로 바꾼다.
8. Fig. 4a에 보조 시계(RNAAgeCalc-Peters) 행이 섞여 있다. 본문 그림에는 주 시계만 둔다.

---

## 2. 공통 규칙 (지시서 10 §2를 이것으로 바꾼다)

1. **패널 기호**: 굵은 소문자 a, b, c…를 10 pt로 패널 왼쪽 위에 둔다. 패널이 하나뿐인 Fig. 1에는 두지 않는다.
2. **글꼴**: Arial을 쓴다. 없으면 Liberation Sans(Arial과 글자 폭이 같은 무료 글꼴)를 쓴다. 둘 다 없으면 DejaVu Sans를 그대로 쓰고 보고서에 적는다.
   - Liberation Sans 설치: `apt-get install fonts-liberation`이 되면 그것을 쓴다. 안 되면 `liberation-fonts` GitHub 릴리스의 TTF를 **저장소 밖** 사용자 글꼴 폴더에 둔다. 그리고 `matplotlib.font_manager.fontManager.addfont()`로 등록한다.
   - 수식 글자도 같은 글꼴로 맞춘다: `mathtext.fontset = "custom"`, `mathtext.rm`, `mathtext.it`, `mathtext.bf`를 같은 글꼴의 보통, 기울임, 굵게로 둔다.
   - `pdf.fonttype = 42`, `ps.fonttype = 42`.
3. **글자 크기**(최종 크기 기준)
   - 눈금, 범례, 그림 안 글자: 8 pt
   - 축 이름: 9 pt
   - 패널 기호: 10 pt
   - 8 pt보다 작은 글자는 없게 한다. 선 굵기 0.5 pt 이상은 그대로 지킨다.
4. **투명도 금지**: EPS는 투명도를 지원하지 않는다. `alpha`를 쓰지 않는다. 연한 색은 미리 섞은 단색으로 쓴다(예: 검정 25% → `#BFBFBF`, 띠 → `#E6E6E6`).
5. **숫자 형식**
   - p는 유효숫자 2자리로 쓴다. p < 0.001이면 `p = 3.0 × 10⁻¹⁰` 형식(mathtext 위첨자)으로 쓴다.
   - 정점은 소수 1자리에 "h"를 붙인다. 비율은 소수 2자리다.
6. **범례 제목 "post hoc"를 지운다.** 그림 전체가 사후인 경우는 원고 범례에 적혀 있다. 확인 결과와 사후 결과가 한 그림에 섞일 때만 범례 항목 끝에 "(post hoc)"를 붙인다.
7. **그림에 쓰는 말** (내부 용어 대신 쓴다)

| 지금 | 바꿀 말 |
|---|---|
| frozen | Preregistered scoring |
| D2 kNN | kNN imputation |
| D11, D11 participant mean | Participant-mean imputation |
| GEO axis | GEO metadata definition |
| relclocktime − 8 | Publication definition |
| Features entered, n_entered | Features used |
| before / after (Fig. 8a, S3a) | Unadjusted / Composition-adjusted |
| before / composition / time (Fig. 8b, S3b) | Unadjusted / Composition-adjusted / Time-corrected |
| before / kNN / mean (Fig. 6c) | Unmasked / kNN / Participant mean |
| lab / free-living / free-living upper | Laboratory / Free-living / Free-living, upper bound |
| p_med (축 이름) | Proportion mediated by composition |
| b(s) / total SD | Bias (fraction of total SD) |
| d (축 이름) | Effect size d (morning − evening) |
| 2A | 2A (score) 또는 2A (years) |
| Mean Δ | Mean Δ, morning − evening (score) 또는 (years) |
| score (centered) | tAge: Centered score / RNAAgeCalc: Centered estimate (years) |

8. **색 하나에 뜻 하나**
   - 자료 색(지시서 10 §2 표)은 자료를 나타낼 때만 쓴다.
   - 시계 색은 시계를 선으로 구분하는 Fig. 7과 Fig. S2에서만 쓴다.
   - 그 밖의 구분(채점 판, 마스킹 그림의 시계)은 회색조와 표식 모양으로 한다.
9. **채움 규칙**: 코사이너 결과를 보이는 모든 패널에서 같은 규칙을 쓴다. 그 모형의 LRT p < 0.05이면 채운 표식, p ≥ 0.05이면 빈 표식이다. 범례 항목은 "Open symbols: LRT p ≥ 0.05" 하나로 쓴다.
10. **겹침과 잘림 없음**: 글자와 글자, 글자와 표식, 범례와 자료가 겹치지 않게 한다. 그림 밖으로 나간 글자도 없어야 한다. §5 점검으로 확인한다.
11. **파일**
    - 그림마다 PDF(벡터), EPS(벡터, 글꼴 포함), TIFF(1200 dpi, LZW), PNG(300 dpi, 검토용)를 낸다.
    - 크기는 너비 174 mm, 높이 234 mm 이하다.
    - 원자료 표 규칙(`figures/manuscript/source_data/`)은 지시서 10과 같다. 바뀐 패널의 원자료 표는 새로 쓴다. 파일 이름은 그대로 둬도 된다.
12. **투고용 사본**: `figures/submission/`에 둔다.
    - 본문 그림: `Fig1.eps`…`Fig8.eps`, `Fig1.tif`…`Fig8.tif`
    - 보충 그림: `FigS1.pdf`, `FigS2.pdf`, `FigS3.pdf`

---

## 3. 본문 그림

### Fig. 1 — 연구 설계
- 왼쪽 위의 "A"를 지운다.
- 64자 해시를 지운다. Freeze 상자 안에 `manifest SHA-256 f19f79c1…`로 쓴다.
- 첫 자료 상자의 첫 줄에 "(participants / samples)"를 적는다.
- "5 clocks"를 "5 clocks (3 primary, 2 secondary)"로 바꾼다.
- "Frozen predictions" 상자는 "Frozen predictions: peak φ̄, amplitude, composition mediation, morning−evening d, 8 module outcomes"로 쓴다.
- 줄 사이를 잇는다.
  - "Frozen predictions"에서 아래 줄 "Freeze" 상자로 화살표를 둔다.
  - 확인 줄에서 Post hoc 줄로 점선 화살표를 두고, 옆에 "after confirmation"을 적는다.
- 잠금 표시
  - 점선 테두리는 그대로 둔다.
  - 상자 왼쪽 위 모서리에 자물쇠 그림을 그린다(사각형과 반원 호, `matplotlib.patches`).
  - 상자 첫 줄은 "Locked until freeze"로 한다.
- 관문 상자는 쉼표 뒤를 띄운다("PER1, PER3").
- 빈 공간을 줄여 높이를 약 120 mm로 한다.

### Fig. 2 — 개인 내 궤적
- 패널 머리글 첫 줄은 "GSE113883 · tAge-chronological" 형식으로 쓴다. 아래 행은 "GSE48113 Day 1 · …"다.
- 머리글 둘째 줄(축 밖)에 `p = 3.0 × 10⁻¹⁰; peak 12.2 h`를 쓴다(§2.5). 축 안의 p와 정점 글자는 지운다.
- y축 이름을 모든 패널에 둔다: tAge는 "Centered score", RNAAgeCalc는 "Centered estimate (years)".
- 참가자 선과 점은 `alpha` 대신 `#BFBFBF` 단색으로 그린다. 선은 0.5 pt다.

### Fig. 3 — 효과 크기와 정점
- **a**
  - §2.9 채움 규칙을 적용한다.
  - 축 안의 범례("GSE122725 (post hoc)")를 지운다.
  - 그림 맨 아래에 표식 범례를 하나 둔다(검정과 회색): ● Development, ■ Confirmation, ◆ Post hoc (GSE122725), "Open symbols: LRT p ≥ 0.05".
- **b**
  - 왼쪽 시계 이름이 잘리지 않게 한다(왼쪽 여백을 늘리거나 "RNAAgeCalc-" 다음에 줄을 바꾼다).
  - 패널 안 빈 곳(오른쪽 아래)에 범례를 둔다.
    - ■ Observed, GSE295758 (95% CI)
    - ◇ Predicted (median of development)
    - 회색 선: Range of development predictions
  - x 이름은 "Effect size d (morning − evening)"다.
- **c** — 원판을 지우고 새로 그린다.
  - 형식은 a와 같은 가로 구간 그림이다.
  - 시계 묶음 3개를 위에서부터 tAge-chronological, tAge-mortality, RNAAgeCalc-DESeq2 순서로 둔다.
  - 묶음마다 자료 4줄을 둔다: GSE39445, GSE56931, GSE113883, GSE48113. y 눈금 이름은 자료 이름이다. 묶음 이름은 오른쪽에 둔다(a와 같게).
  - x는 "Clock time (h)"이다. 범위 0–24, 눈금 0, 6, 12, 18, 24.
  - 점은 정점 점추정이다. 막대는 `peak_center ± peak_halfwidth_95`다. 지시서 10의 Fig. 3C 원자료와 같은 값이며, Fig. 4b와 같은 방식이다.
  - 0이나 24를 넘는 막대는 축 끝에서 화살표로 자른다.
  - 묶음마다 예측 정점 φ̄에 검은 세로 점선을 둔다. φ̄ ± 3 h에는 회색 띠(`#E6E6E6`, 자료 뒤)를 둔다. 이 띠는 사전등록 보조 점검 창이다.
  - §2.9 채움 규칙을 적용한다. 자료 색과 표식은 a와 같다.
  - GSE122725는 넣지 않는다. 시각 축이 기상 기준이라 벽시계 축에 둘 수 없다.
  - 원자료 열: 시계, 자료, peak_h, peak_center, peak_halfwidth_95, LRT p, φ̄

### Fig. 4 — 세포 구성
- **a**
  - GSE56931 RNAAgeCalc-Peters 행을 빼서 Fig. S1c로 옮긴다(§4).
  - 행을 시계별로 묶는다. 묶음 순서와 오른쪽 묶음 이름은 Fig. 3a와 같다. y 눈금 이름은 자료 이름만 쓴다.
  - x 이름은 "Proportion mediated by composition"이다. 0과 1의 세로선, −0.5–1.5 범위, 화살표 자르기는 그대로 둔다.
- **b**
  - x축을 18시에서 다음날 18시까지로 바꾼다. 눈금은 18, 0, 6, 12, 18이다.
  - 정점과 구간 끝을 이 축으로 옮긴다(t < 18이면 t + 24). 축을 넘는 구간은 화살표로 자른다.
  - 범례의 "LRT p ≥ 0.05" 항목은 "Open symbols: LRT p ≥ 0.05"로 쓴다. GSE48113 항목은 "GSE48113 Day 1 (post hoc)"로 쓴다.
  - 각 집단 줄 안에서 자료의 위아래 순서를 범례 순서와 같게 둔다.

### Fig. 5 — 모듈
- **a**
  - 범례를 축 밖(아래)으로 옮긴다. 항목은 다음과 같다.
    - 채운 원: Participant-mean imputation (post hoc)
    - Open symbols: LRT p ≥ 0.05
    - 회색 ×: Preregistered scoring
    - 크기 기준 원 3개: p = 0.05, 10⁻⁴, 10⁻⁸
  - 축 이름: x는 "Predicted peak φ̄ (h)", y는 "Observed peak, GSE48113 (h)"다.
  - 같은 모듈의 ×와 원을 회색 가는 점선으로 잇는다(x가 같다).
  - 모듈 이름은 모두 축 안에 두고, 짧은 지시선으로 점과 잇는다.
    - 이름끼리, 이름과 표식이 겹치지 않게 한다.
    - `adjustText`를 쓰거나 좌표를 손으로 정한다. 손으로 정한 좌표표는 원자료 폴더에 둔다.
  - × 옆 특징 수는 "n = 14" 형식으로 쓴다.
- **b** (근거 표)
  - 열은 다음과 같다.

| 열 묶음 | 칸 |
|---|---|
| Development −log10 q | GSE39445, GSE56931, GSE113883 |
| Preregistered scoring | Features, Confirmed |
| Participant-mean imputation (post hoc) | Features, Confirmed |
| GSE295758 | Predicted, Observed, Agree |

  - q 색칸 안에 −log10 q 값을 소수 1자리로 적는다. 짙은 칸은 흰 글자로 쓴다. 값을 적었으므로 색막대는 지운다.
  - Confirmed 칸은 "Yes"나 "No"로 쓴다. 특징 수가 0이면 "—"로 쓴다(확인 판정에서는 "아님"으로 셌다. 원고 범례에 적었다).
  - 부호는 "+"와 "−"(U+2212)로 쓴다. Agree는 "Yes"나 "No"로 쓰고 빈칸을 두지 않는다.
  - ✓, ✗ 기호는 쓰지 않는다(Arial에 없다).
  - 오른쪽 빈 공간을 없앤다.

### Fig. 6 — 공개 자료의 함정
- **a**
  - 축 이름: x는 "Peak, GSE39445 (h since wake)", y는 "Peak, GSE122725 (h since wake)"다.
  - **x에 쓴 GSE39445 시각 변수 이름을 보고서에 적는다.** 기상 기준이 아니면 축 이름을 그 변수에 맞게 쓴다. 그리고 보고서 §6에 적는다.
  - 같은 유전자의 두 판(GEO metadata definition, Publication definition)을 회색 가는 선으로 잇는다.
  - 범례: 회색 빈 원 "GEO metadata definition", 채운 분홍 마름모 "Publication definition (post hoc)", "Open symbols: LRT p ≥ 0.05".
  - a와 b 사이를 넓힌다. b의 패널 기호와 축이 a의 눈금이나 축 이름과 겹치지 않게 한다.
- **b**
  - 패널 머리글은 "GSE48113"이다.
  - 채점 판은 회색조와 모양으로 구분한다. 자료 색은 쓰지 않는다.
    - Preregistered scoring: ● 검정
    - kNN imputation: ■ `#808080`
    - Participant-mean imputation: ◆ 흰 바탕에 검은 테두리
  - 기준선
    - 그림 전체를 가로지르는 세로선 2개를 지운다.
    - tAge 두 행에만, 그 행 높이의 짧은 세로 막대로 0이 아닌 가중치 수(1,839와 10,487)를 표시한다.
    - 범례는 "Nonzero-weight features of the clock"이다.
  - 특징 수 0은 지금처럼 축 왼쪽 끝에 "0"으로 적는다. 범례는 "0 = no feature available"이다.
  - x 이름은 "Features used (log scale)"다.
- **c**
  - 패널 머리글은 "GSE39445, 20% of values masked"다.
  - x 위치: Unmasked │ Random masking: kNN, Participant mean │ Low-expression masking: kNN, Participant mean
    - 두 가림 방식의 Unmasked 값이 같으면 한 번만 둔다.
    - 다르면 두 번 두고, 그 이유를 보고서에 적는다.
  - 시계는 모양으로 구분한다: RNAAgeCalc-DESeq2 ●, RNAAgeCalc-Peters ▲, RNAAgeCalc-GTEx ▼. 모두 짙은 회색(`#404040`)이다. 시계 색은 쓰지 않는다(a의 분홍과 겹친다).
  - p = 0.05 파선의 오른쪽 끝에 "p = 0.05"라고 적는다.
- b와 c의 범례는 각 패널 옆에 따로 둔다. 두 범례를 섞지 않는다.

### Fig. 7 — 편향 시뮬레이션
- 범례 제목 "post hoc"를 지운다. 시나리오 이름은 §2.7을 따른다.
- **a**
  - y 이름은 "Bias (fraction of total SD)"다. y = 0에 가는 회색 선을 둔다.
  - RNAAgeCalc-DESeq2의 Laboratory 선과 Free-living 선에서 s = 3인 점에 작은 빈 원을 둔다.
  - 년 환산은 선과 겹치지 않는 빈 곳에 글 상자로 둔다: "RNAAgeCalc-DESeq2, s = 3 h: 1.1 years (laboratory), 0.6 years (free-living)"
  - 선 위의 "1.1 years"와 "0.6 years"는 지운다.
- **b–e**: 패널 안 왼쪽 위에 고정 조건을 적는다.
  - b: "Design 1, n = 100"
  - c: "Design 1, s = 3 h"
  - d: "Design 2, s = 4 h"
  - e: "Design 3, laboratory, n = 50"
- **b, c, d**: 0.05(점선)와 0.10(파선) 가로선을 모두 둔다(지금은 b에만 있다). 선 오른쪽 끝에 "0.05"와 "0.10"을 적는다.

### Fig. 8 — 완화
- **a**
  - x 범주는 "Unadjusted", "Composition-adjusted"다.
  - y 이름은 패널마다 둔다: tAge는 "2A (score)", RNAAgeCalc는 "2A (years)"다.
  - **채움 규칙을 바꾼다.** 왼쪽 점은 보정 전 LRT p로, 오른쪽 점은 보정 후 LRT p로 채움을 정한다(§2.9). 값은 `results/stage9/composition_cosinor.tsv`의 `lrt_p_before`와 `lrt_p_after`다.
  - 비율 글자는 오른쪽 점 옆에 y 순서대로 둔다. 글자 사이에 최소 간격(8 pt 글자 높이)을 둔다. 겹치면 위나 아래로 밀고, 가는 지시선으로 점과 잇는다.
- **b**
  - GSE295758의 자료 색 `#E69F00`과 사각형(확인 자료 표식)으로 그린다. 시계 색은 쓰지 않는다.
  - x 범주는 "Unadjusted", "Composition-adjusted", "Time-corrected"다.
  - y 이름: tAge는 "Mean Δ, morning − evening (score)", RNAAgeCalc는 "Mean Δ, morning − evening (years)"다.
  - x 양쪽에 여백을 두어 범주 글자가 잘리지 않게 한다.
- 범례
  - 제목을 지운다.
  - 자료 항목에 "GSE295758"(■ `#E69F00`)을 더한다.
  - 빈 표식 항목은 "Open symbols: LRT p ≥ 0.05"로 쓴다.

---

## 4. 보충 그림

| 그림 | 수정 |
|---|---|
| Fig. S1 | a와 b에 Fig. 3a, 3b의 수정을 그대로 한다. b의 빈 공간을 줄인다(높이 축소). **c를 새로 둔다**: Fig. 4a에서 옮긴 보조 시계의 p_med [95% 구간]이다. 개발 LRT p < 0.05인 보조 시계 × 자료만 넣는다(지금은 GSE56931 RNAAgeCalc-Peters 한 줄). 형식은 Fig. 4a와 같다. |
| Fig. S2 | Fig. 7의 수정을 그대로 한다. d에서 RNAAgeCalc-Peters와 RNAAgeCalc-GTEx의 Laboratory 선이 겹친다. 두 시계의 s = 4에서 \|b\|/총 SD 값을 보고서에 적는다. |
| Fig. S3 | Fig. 8의 수정을 그대로 한다(a 채움 규칙, b의 GSE295758 색과 사각형, 범주 이름, 여백). b에 Time-corrected가 없는 것은 그대로 둔다(시각 보정 c는 주 시계에만 정의했다). |

---

## 5. 점검 (그림을 다 그린 뒤)

모든 그림에서 `fig.canvas.draw()` 뒤에 아래를 계산해 표로 낸다.

1. **겹침**: 아래 네 가지를 센다. 목표는 모두 0이다. 남는 것은 그림과 글자별로 보고한다. 같은 축의 눈금 글자끼리는 빼도 된다.
   - 보이는 글자(`Text`)의 `get_window_extent()` 상자끼리 겹치는 쌍
   - 글자 상자와 표식 위치가 겹치는 쌍
   - 범례 상자 안에 들어간 자료 점
   - 그림 경계 밖으로 나간 글자
2. **글자 크기**: 가장 작은 `get_fontsize()`를 적는다. 8 미만이면 고친다.
3. **투명도**: `alpha`가 None도 1도 아닌 artist 수를 센다. 0이어야 한다.
4. **글꼴**: PDF는 `pdffonts`로 본다. EPS는 `ps2pdf -dEPSCrop`로 PDF로 바꾼 뒤 `pdffonts`로 본다. 모두 emb가 yes여야 한다. 글꼴 이름을 적는다.
5. **크기**: 그림별 크기(mm, 너비 × 높이)와 TIFF 파일 크기(MB)를 적는다.
6. 검토 묶음 `figures/manuscript/all_figures_review.pdf`를 다시 만든다.

---

## 6. 보고서 템플릿 (`reports/TOD_보고서_11.md`)

**보고서를 읽는 사람은 파일을 볼 수 없다.** 표는 본문에 넣는다. 해석과 결론은 쓰지 않는다. 그림 묶음 PDF는 연구자가 따로 넘긴다.

~~~markdown
# TOD 보고서 11 — 원고 그림 최종 수정

## 1. 파일 목록 (그림별 PDF·EPS·TIFF·PNG, 크기 mm, TIFF MB, submission 사본)
## 2. 글꼴 (사용 글꼴, 설치 방법, pdffonts 결과, 가장 작은 글자 크기)
## 3. 점검 결과 (§5의 1–5 표)
## 4. 확인용 값 (원고 수치 대조용, 저장값 그대로)
### 4.1 Fig. 3c: 시계 × 자료별 peak_h, peak_center, peak_halfwidth_95, LRT p, φ̄
### 4.2 Fig. 6a: x에 쓴 GSE39445 시각 변수 이름과 정의
### 4.3 Fig. 6c: 시계 × 가림 방식 × 처리별 LRT p
### 4.4 Fig. 8a와 S3a: 채움에 쓴 파일과 열 이름
### 4.5 Fig. S2d: RNAAgeCalc-Peters와 RNAAgeCalc-GTEx Laboratory의 s = 4 |b|/총 SD
## 5. 판단이 필요한 사항
## 6. 지시서와 다르게 한 것 (빠짐없이)
## 7. 커밋 해시
~~~
