# TOD 보고서 6 — 사전등록 배치, 확인 코드, 리허설, 동결

## 1. 요약 (수치만, 12줄 이내)

- 받은 사전등록 문서 SHA-256 `613cb42f1997f62aa1f28a1e37003cc830302c4cc63393125bbfaac307256b6b`. 표를 붙인 문서 SHA-256 `6be4e76d4345a1d99a32a18629a1a3afced3137f1cf8e1c5907e546c5b1b1e71`.
- 동결 해시(manifest SHA-256) `f19f79c1988d75e7085a3f570b80a96a2e3350b368335ac80b0ae4fb162120ba`. 동결 커밋 `8a203ffdc5a5f1e933cd4f64154c0d5972480e32`, 시각 2026-10-16 16:21:17 +0900. ots 없음. OSF zip SHA-256 `e9a1c31ce9d50f54587467a5dc7140f6eaba0430a4fec3416183cc85d650cfe8`.
- null에서 p < 0.05 비율 범위 0.038–0.078 (18행, H1·H3·H4). planted 비율 범위 0.125–1 (18행). H2는 `n_finite_p` 0이라 비율이 비어 있다. H5 `planted_1c` 0.03–0.075, `planted_0.5c` 0.485–0.775.
- 파이프라인 종료 코드: R1 0, R2 0, R3 0.
- 독립 검사 최대 차이: LRT 상대 차이 `1.2674368553051373e-09`, t 통계량 절대 차이 `8.881784197001252e-15`.
- 잠금 시험: L1·L2 모두 종료 코드 1, `openat` 횟수 0, 메시지 `config/unlock_approved.txt is missing; locked expression stays closed`.
- 판단 사항 14건 (6.2.1–6.2.10, 6.3.1–6.3.4).

## 2. 환경과 자원 (바뀐 것, /usr/bin/time 기록)

버전은 `prereg/frozen/environment.tsv`와 같다. R 4.5.3, lme4 2.0.6, RNAAgeCalc 1.22.0, MCPcounter 1.2.0, impute 1.84.0, limma 3.66.0, clariomshumanhttranscriptcluster.db 8.8.0, org.Hs.eg.db 3.22.0, Python 3.12.7, numpy 1.26.4, pandas 2.2.3, scikit-learn 1.3.2, scipy 1.13.1, joblib 1.6.0. 독립 검사의 statsmodels 0.14.6은 `/home/jkang/anaconda3/bin/python3`에 있다.

| 실행 | 종료 | 경과 | 최대 RSS (kB) |
|---|---|---|---|
| 점수 리허설 | 0 | 0:08:37.07 | 237344 |
| R1 | 0 | 0:12:17.70 | 1663016 |
| R2 | 0 | 0:06.96 | 677968 |
| R3 | 0 | 0:04:22.69 | 1696684 |

R1이 끝나는 구간과 R3 성공 실행이 겹쳤다. 각 작업의 최대 RSS는 40 GB 아래다. `ulimit`은 바꾸지 않았다. `ots` 명령은 없다.

## 3. 결정 기록 (config/decisions.tsv에 이번에 더한 행)

승인 5건. 6.1.1 사람 자료 모듈은 다종 mortality 14개이고 설치류 재채점은 `TMrd_`. 6.1.2 분해는 설치류 23개 모듈 유전자 집합. 6.1.3 σ는 모형 M의 REML 잔차 SD, σ_u는 영모형 임의절편의 REML SD. 6.1.4 3시간 조건은 원형 거리 3 이하. 6.1.5 `TMrd_` 이름.

구현 기록 6.2.1–6.2.10과 리허설 기록 6.3.1–6.3.4는 §9에 적는다.

## 4. 사전등록 문서와 고정 매개변수

### 4.1 배치와 해시

받은 파일은 `prereg/TOD_preregistration_v1.md` 앞에 그대로 있고, 그 뒤에 표 F1–F6만 붙었다. 받은 파일 SHA-256과 표를 붙인 파일 SHA-256은 §1과 같다. 확인 코드가 읽는 표는 `prereg/frozen/F1.tsv` … `F6.tsv`다.

### 4.2 표 F1–F6 (본문에 전부, 원천 경로와 SHA-256 포함)

F1. 원천 `results/stage5/clock_summary_long.tsv` SHA-256 `bb2be685a97ab6a2141ff49c7412902b3d936bfbd803a3c3ee48be00301b9dad`, `results/stage5/clock_summary_phase.tsv` SHA-256 `38077898da876ea43de2166396e270ded605aa11c393727872b3f3dadd3d7d7b`.

```
clock	role	phi_bar	max_circ_dist	dataset	two_a	sigma	sigma_u	peak_h	lrt_p
TAms_chrono	primary	12.616768998085	1.61677078040933	GSE39445	0.0345210353298852	0.0171269287358502	0.0276429862743754	14.2335397784943	5.49745798912395e-20
TAms_chrono	primary	12.616768998085	1.61677078040933	GSE56931	0.00943066206760042	0.0197534498990921	0.0331124542924689	11.4579688514277	0.181165188487945
TAms_chrono	primary	12.616768998085	1.61677078040933	GSE113883	0.0418702301338288	0.0253117943835499	0.0416315613478055	12.188050742382	2.96936394371969e-10
TAms_mort	primary	11.8599407935094	1.79984172390618	GSE39445	0.148755555495313	0.0860831456543238	0.161562476919705	13.4288594848278	3.45851149483147e-15
TAms_mort	primary	11.8599407935094	1.79984172390618	GSE56931	0.0629468359700225	0.0740143176007128	0.177501208368076	12.0688456403715	0.00673707816031792
TAms_mort	primary	11.8599407935094	1.79984172390618	GSE113883	0.102562103238955	0.0897840656082817	0.184575482929547	10.0600990696032	4.44466652817908e-06
RA_default	primary	11.2423863477089	1.35536907505411	GSE39445	1.98498899424841	9.34019861854301	8.81869510215541	12.0064473932836	0.571779550797501
RA_default	primary	11.2423863477089	1.35536907505411	GSE56931	4.61694071236019	6.68803391836383	6.90600678789726	9.88701727265479	0.0187399469939356
RA_default	primary	11.2423863477089	1.35536907505411	GSE113883	6.770165944749	6.28931435261543	6.86875908783177	11.812631972311	3.24675832907807e-05
RA_Peters	secondary	7.56604107260516	4.27017978445731	GSE39445	3.86882060477993	9.12159876450608	14.613080804435	6.84503499428439	0.0822398090507864
RA_Peters	secondary	7.56604107260516	4.27017978445731	GSE56931	6.37997357324785	5.58777765048951	13.9973116863226	4.54192615314062	0.000105116530090887
RA_Peters	secondary	7.56604107260516	4.27017978445731	GSE113883	3.1551856815567	6.04743769494992	13.6399297750127	11.8362208570625	0.0773443137299791
RA_GTEx	secondary	11.518486825228	5.35019887986941	GSE39445	4.5622147045564	25.1603887545092	16.2536695399875	13.7137336093391	0.643663244864414
RA_GTEx	secondary	11.518486825228	5.35019887986941	GSE56931	8.41295616809185	14.1878999046968	11.0606886507247	13.2671292133393	0.0977660921594727
RA_GTEx	secondary	11.518486825228	5.35019887986941	GSE113883	4.34369415553422	8.81987355742966	10.4317015658613	6.16828794535857	0.104384139695734
```

F2. 원천 `results/stage5/gse295758_expected_d.tsv` SHA-256 `b6c1ef3c4db6202a467ae3636d0de3aa7e55b47bf766bcd67ece61ca538c7e98`.

```
clock	role	d_GSE39445	d_GSE56931	d_GSE113883	d_exp_median	predicted_sign
TAms_chrono	primary	-0.087085971883512	0.208443813199012	0.534281736707176	0.208443813199012	+
TAms_mort	primary	0.182025471981377	0.291250817026961	0.693084284817919	0.291250817026961	+
RA_default	primary	0.0749177677231186	0.429771444888433	0.412448894924069	0.412448894924069	+
RA_Peters	secondary	0.286305297017999	0.498487573681539	0.197988187071947	0.286305297017999	+
RA_GTEx	secondary	0.00960009963837758	0.0799546478017167	0.308963171982023	0.0799546478017167	+
```

F3. 원천 `results/stage5/reproducible_modules.tsv` SHA-256 `6b6c9a1b27cc34b060358fc656d51e3ba329c3a3acdac3c3ee536c6c2ab4dadd`, `results/stage5/gse295758_expected_d.tsv` SHA-256 `b6c1ef3c4db6202a467ae3636d0de3aa7e55b47bf766bcd67ece61ca538c7e98`, `results/stage5/h1_point_estimates.tsv` SHA-256 `60f24af5ff25abce59570d1f79d391adbcd9093adad2f8faf3f90c34e2c5c2f1`.

```
clock	family	module	parent_clock	phi_bar	q_GSE39445	q_GSE56931	q_GSE113883	peak_GSE39445	peak_GSE56931	peak_GSE113883	d_exp_median	predicted_sign
TM_turquoise_mortality	TM	turquoise		11.5724011994242	1.7030021051141e-07	0.00657277863416067	3.693684980855e-06	12.6215119266462	11.1568084467969	10.9484328837053	0.467000081226947	+
CMP_TAms_chrono_turquoise	CMP	turquoise	TAms_chrono	13.2837328415609	2.42711785115231e-08	0.00626983142167074	0.0107579084724153	13.1432263174321	12.1422987334514	14.5730749724508	0.20086732450169	+
CMP_TAms_chrono_yellow	CMP	yellow	TAms_chrono	12.6128269162576	0.00155470474860942	4.69738653218747e-08	9.00796184650676e-09	13.045845902978	11.2503142784199	13.523129022035	0.1496592336641	+
CMP_TAms_chrono_orange	CMP	orange	TAms_chrono	1.07332762742524	4.17596855270051e-09	0.00616665634539307	9.86133709312738e-06	2.44005820482784	0.0303979642477727	0.765259718398529	-0.280140424913859	-
CMP_TAms_chrono_brown	CMP	brown	TAms_chrono	1.32505139804729	2.43965028951066e-08	0.0359875447558034	2.67476690758748e-05	0.91282893412168	0.84038477850989	2.22823394895771	-0.168830899708354	-
CMP_TAms_chrono_red	CMP	red	TAms_chrono	0.781691026210673	6.00047155091236e-06	0.00206822511320336	0.00113591207810166	2.08110231751181	23.2474694413108	1.00061811727099	-0.169417336437209	-
CMP_TAms_mort_turquoise	CMP	turquoise	TAms_mort	13.003811149926	6.77796115376162e-09	0.000125036076009668	3.88973233356303e-07	13.2360241903645	11.9113647277402	13.8564175468774	0.1888870282562	+
CMP_TAms_mort_darkred	CMP	darkred	TAms_mort	14.9666596497573	0.00299942175721486	0.0275912785277717	9.48016174502148e-08	17.0651005937284	13.327427235357	14.5608059570561	-0.167328334678284	-
```

F4. TAms 원천 `results/stage5/h1_point_estimates.tsv` SHA-256 `60f24af5ff25abce59570d1f79d391adbcd9093adad2f8faf3f90c34e2c5c2f1`. RA_default 원천 `results/stage4/rnaseq/h1_point_estimates.tsv` SHA-256 `f0ad35b89eb5b0a97979b809448bde321cd97c4ec658af6cae87cb835e6ff85e`.

```
clock	dataset	beta_cos	beta_sin	formula	c
TAms_chrono	GSE113883	-0.0209097494754014	-0.00103025227051267	c = f(8) - f(20); f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24); c = -beta_cos + beta_sin*sqrt(3)	0.019125300198260260499475555777287663859727086340733
TAms_mort	GSE113883	-0.0448086053054778	0.0249386275840794	c = f(8) - f(20); f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24); c = -beta_cos + beta_sin*sqrt(3)	0.088003575352142004086395806850464311947738779118951
RA_default	GSE113883	-3.38101122840144	0.165981334106488	c = f(8) - f(20); f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24); c = -beta_cos + beta_sin*sqrt(3)	3.6684993321819421764583897146097149703131816543407
```

F5. 소프트웨어 행은 R 4.5.3, lme4 2.0.6, RNAAgeCalc 1.22.0, MCPcounter 1.2.0, impute 1.84.0, limma 3.66.0, clariomshumanhttranscriptcluster.db 8.8.0, org.Hs.eg.db 3.22.0, python 3.12.7, numpy 1.26.4, pandas 2.2.3, scikit-learn 1.3.2, scipy 1.13.1, joblib 1.6.0, clock_python `/HDD8T1/WORK/GeroSci/tod4_clock_env`, clock_specificity_commit `9d109499f0bf47a7869cd7b4bcffe311b8be9b75`, tAge_tag `1.5.0`, md5_TAms_chrono `7dc8a470f5db20a75d7052829efbd37a`, md5_TAms_mort `87307d7f5572d9932382912d3fd9571c`. 나머지 행은 위 원천 파일의 경로와 SHA-256이고, 추가로 `results/stage4/rnaseq/g2_clocks.tsv` `9984f8ad52ecab1499a5d12e36e2923c03d0acf675a193f5416dd9b3cf4181ff`, `results/stage5/g2_tams.tsv` `2fd4c6b63eed1cd60c3adf4b265803e7ef660e670769823d4df39e8aa64e257a`, `results/stage5/locked_design_power.tsv` `e508d01b1b3f1ab9893c0e0332be9df90140ac49defb3faeee071c3525bdfb37`.

F6. 원천 `results/stage5/gse295758_expected_d.tsv` SHA-256 `b6c1ef3c4db6202a467ae3636d0de3aa7e55b47bf766bcd67ece61ca538c7e98`, `results/stage5/locked_design_power.tsv` SHA-256 `e508d01b1b3f1ab9893c0e0332be9df90140ac49defb3faeee071c3525bdfb37`.

```
block	clock	design	param_dataset	d_exp_median	power_alpha_0.05	power_alpha_0.01	frac_p_lt_0.05	frac_p_lt_0.01	n_rep	n_converged
L3	TAms_chrono	GSE295758	median_of_three	0.208443813199012	0.345077590733327	0.150545763753431				
L3	TAms_mort	GSE295758	median_of_three	0.291250817026961	0.587339593367421	0.336502748601378				
L3	RA_default	GSE295758	median_of_three	0.412448894924069	0.870300205086949	0.682114842365414				
L3	RA_Peters	GSE295758	median_of_three	0.286305297017999	0.572873321105643	0.323267563613161				
L3	RA_GTEx	GSE295758	median_of_three	0.0799546478017167	0.0919710023400413	0.0243768273749048				
L1L2	RA_GTEx	GSE48113	GSE39445				0.0935	0.027	2000	2000
L1L2	RA_GTEx	GSE48113	GSE56931				0.579	0.3245	2000	2000
L1L2	RA_GTEx	GSE48113	GSE113883				0.468	0.2485	2000	2000
L1L2	RA_Peters	GSE48113	GSE39445				0.361	0.1545	2000	2000
L1L2	RA_Peters	GSE48113	GSE56931				0.998	0.968	2000	2000
L1L2	RA_Peters	GSE48113	GSE113883				0.477	0.244	2000	2000
L1L2	RA_default	GSE48113	GSE39445				0.112	0.0375	2000	2000
L1L2	RA_default	GSE48113	GSE56931				0.773	0.5335	2000	2000
L1L2	RA_default	GSE48113	GSE113883				0.9795	0.931	2000	2000
L1L2	TAms_chrono	GSE48113	GSE39445				1	1	2000	2000
L1L2	TAms_chrono	GSE48113	GSE56931				0.431	0.201	2000	2000
L1L2	TAms_chrono	GSE48113	GSE113883				1	1	2000	2000
L1L2	TAms_mort	GSE48113	GSE39445				1	1	2000	2000
L1L2	TAms_mort	GSE48113	GSE56931				0.8865	0.7175	2000	2000
L1L2	TAms_mort	GSE48113	GSE113883				0.9945	0.9775	2000	2000
L1L2	RA_GTEx	GSE122725	GSE39445				0.0755	0.0225	2000	2000
L1L2	RA_GTEx	GSE122725	GSE56931				0.3395	0.1505	2000	2000
L1L2	RA_GTEx	GSE122725	GSE113883				0.1995	0.067	2000	2000
L1L2	RA_Peters	GSE122725	GSE39445				0.1635	0.059	2000	2000
L1L2	RA_Peters	GSE122725	GSE56931				0.83	0.6185	2000	2000
L1L2	RA_Peters	GSE122725	GSE113883				0.2995	0.1155	2000	2000
L1L2	RA_default	GSE122725	GSE39445				0.0975	0.03	2000	2000
L1L2	RA_default	GSE122725	GSE56931				0.456	0.22	2000	2000
L1L2	RA_default	GSE122725	GSE113883				0.842	0.644	2000	2000
L1L2	TAms_chrono	GSE122725	GSE39445				1	0.9985	2000	2000
L1L2	TAms_chrono	GSE122725	GSE56931				0.237	0.092	2000	2000
L1L2	TAms_chrono	GSE122725	GSE113883				0.997	0.976	2000	2000
L1L2	TAms_mort	GSE122725	GSE39445				0.9975	0.988	2000	2000
L1L2	TAms_mort	GSE122725	GSE56931				0.6175	0.3815	2000	2000
L1L2	TAms_mort	GSE122725	GSE113883				0.8635	0.6735	2000	2000
```

### 4.3 문서와 사실이 다른 점

문서 문장은 고치지 않았다. 표기 차이는 §9에 있다. 분석 결과를 바꾸는 차이로 동결을 멈춘 항목은 없다. 6.2.3은 PER1 프로브셋을 지우게 되어 6.2.10에서 철회했다.

## 5. 확인 분석 코드

### 5.1 파일 목록과 실행 순서

진입점 `scripts/confirm/run_confirm.sh`는 `scripts/confirm/run_pipeline.R --mode locked`를 호출한다. 순서는 manifest 대조(`check_manifest`), 잠금 해제 확인(`load_locked_expression`), 그 다음 전처리·관문·H1–H5·S4·민감도다. 이번 단계에서는 이 진입점을 잠금 자료에 대해 실행하지 않았다. 리허설은 `scripts/43_partD_pipeline.R`가 같은 함수를 잠금 경로가 아닌 파일에 호출한다.

확인 코드 파일: `preprocess_array.py`, `preprocess_l2.R`, `preprocess_l3.R`, `score_clocks.py`, `lib_tests.R`, `lib_lock.R`, `run_pipeline.R`, `download_l3.R`, `build_gencode_table.py`, `run_confirm.sh`.

### 5.2 implementation_map 표

`prereg/implementation_map.tsv`는 헤더를 빼고 46행이다. 열은 문서 절, 항목, 파일, 함수. SHA-256 `8329d444099a2f271edbeb84afbce708199f37e033ccaac08e21fc570db359c0`.

### 5.3 잠금 장치 구현

`load_locked_expression()`만 잠금 발현 경로를 연다. `config/unlock_approved.txt`가 없거나, 그 첫 줄과 환경변수 `TOD_UNLOCK`이 `prereg/manifest.tsv`의 SHA-256과 같지 않으면 경로를 열기 전에 오류를 낸다. `config/unlock_approved.txt`는 만들지 않았다. `download_l3_supplements()`는 `scripts/confirm/download_l3.R`에 있고 호출하지 않았다.

### 5.4 L3 형식 분기, GENCODE 38 표

`l3_layout_from_header()`는 첫 줄에 유전자 열이 있으면 `gene_column`, 없으면 `transcript_to_gene`을 반환한다. 유전자 열이 Ensembl이면 GENCODE로 기호에 연결한다. GENCODE release 38 basic GTF는 `prereg/ref/gencode.v38.basic.annotation.gtf.gz`, 바이트 28152190, SHA-256 `c99260412682fdce39dfb729ba8b436b43084bfec7857c5d71c9884f7cc343e6`. transcript 행 111832. 표 `prereg/frozen/gencode_v38_basic_transcripts.tsv`의 열은 `transcript_id`, `gene_id`, `gene_name`, `transcript_id_base`, `gene_id_base`. SHA-256 `d057a8d2f7a4e9dbef5884fe5ccd7f6bd887f7d3d193bebcdd0a76d6a3c2b8a1`.

## 6. 리허설

### 6.1 점수 수준 시뮬레이션 (시나리오 × 설계 × 시계: 1종 오류 [구간], 검정력, 편향)

설계 수. L1 Day 1 147샘플·22명, Day 4 138샘플·21명. L2 period 1 70샘플·11명, period 2 71샘플·11명. L3 58명. 시드는 시나리오 시작에서 `set.seed(20261015 + scenario_id)`, scenario_id 1–16. 부트스트랩은 하지 않았다. null 500회, planted 200회.

null의 p < 0.05 비율은 0.038 (H1 L1, TAms_mort, 95% 구간 0.0230307708062476–0.0587072904632285)에서 0.078 (H4 L2, RA_default, 구간 0.0560511737755604–0.105087726430066)이다. planted 비율은 0.125 (H4 L2, RA_default)에서 1 (H1 L1 TAms_chrono, H1 L2 TAms_chrono, H1 L3 TAms_mort)이다. H2의 `frac_p_lt_0.05`는 비어 있고 `bias_p_med`는 null에서 0.00105043169314353, −0.00154597763975779, 0.000720508167687195, planted에서 0.00684606961010625, −0.0374925221003727, −0.0183749811453217 (시계 순서 TAms_chrono, TAms_mort, RA_default)이다.

H1 planted의 `bias_2a`는 L1에서 0.00073074898826607, 0.00183270875756049, 0.136660183042706이고, `bias_peak`는 0.0423099572992793, −0.0090934022162017, −0.0785926124014953이다. H1 L3 planted의 `bias_mean_delta_minus_c`는 8.97808094107671e-06, 0.00016941668488401, −0.0440727484077078이다. H3 planted의 `bias_age`는 −0.000232035725788012, −0.000474429444318525, −0.0432238681914999이다. H5 `planted_1c`의 `bias_mean_delta_minus_c`는 −0.000151554961360311, 0.000845965166444997, −0.00793508729628999이고, `planted_0.5c`는 −0.00963640470408257, −0.0441596310774048, −1.79829584262428이다. 48행 전체는 `results/stage6/rehearsal/scores/score_rehearsal.tsv`에 있다.

### 6.2 파이프라인 리허설 R1–R3 (종료 코드, 시간, RSS, 산출물, R1 관문·가족 P 표)

시간과 RSS는 §2와 같다.

R1 관문. PER1 LRT p `9.73597605831301e-18`, peak `6.90324597433964`, pass TRUE. PER3 LRT p `1.177374716677e-06`, peak `6.08583286807854`, pass TRUE. `dataset_pass` TRUE. 전처리 audit: scale `log`, n_probes 43758, n_samples 438, n_genes 19541, n_dropped_multi 0, n_dropped_no_symbol 13253.

R1 가족 P. Holm은 이 실행의 3검정을 대상으로 한다. 부트스트랩 2000회는 세 시계 모두 `n_converged` 2000.

| clock | lrt_stat | lrt_p | q_holm | two_a | peak_h |
|---|---|---|---|---|---|
| TAms_chrono | 88.6948321165556 | 5.49745798912489e-20 | 1.64923739673747e-19 | 0.0345210353298852 | 14.2335397784943 |
| TAms_mort | 66.5958762044224 | 3.45851149483138e-15 | 6.91702298966275e-15 | 0.148755555495314 | 13.4288594848275 |
| RA_default | 1.11800352518344 | 0.571779550797631 | 0.571779550797631 | 1.98498899424706 | 12.0064473932843 |

H2는 TAms_chrono와 TAms_mort에 기록됐다. TAms_chrono `p_med` 0.604554545853913, 구간 0.360490943298219–0.8814079662504. TAms_mort `p_med` 0.648619163581135, 구간 0.3284627436383–1.00587179922408. H4 LRT p는 0.0252018157335846, 0.303286988044451, 0.0972621301754652이고 q_holm은 0.0756054472007539, 0.303286988044451, 0.19452426035093이다. `stopped.tsv`는 `complete`.

R2. `keys()` 24351, scale `log`, multi-symbol probeset 911, 유전자 20179. 빛 그룹 참가자 bright 3, dim 8. 관문: PER1 p 0.301053683543012, peak 2.60721425135562, pass FALSE. PER3 p 0.977450944381598, peak 12.0091723385739, pass FALSE. `dataset_pass` FALSE. `stopped.tsv`는 `gate fail`. 이후 검정 파일은 없다.

R3. 피험자 11, 유전자 31637, 샘플 22. `gene_column` 분기는 열 `gene_name`. `transcript_to_gene` 분기는 열 `transcript_id`. 두 분기의 PER1 관문은 mean_delta `0.710676514064405`, p `0.00787278002669761`, pass TRUE. 두 분기의 가족 P·S1 표는 같은 6행이다.

| clock | family | n | mean_delta | t_stat | p | q_holm |
|---|---|---|---|---|---|---|
| TAms_chrono | P | 11 | 0.0148068238558497 | 1.04370691660465 | 0.321195068142229 | 0.321195068142229 |
| TAms_mort | P | 11 | 0.115439828157331 | 3.3329176159176 | 0.00758047833344554 | 0.0227414350003366 |
| RA_default | P | 11 | 1.69809697431048 | 1.73013215149272 | 0.114290611511807 | 0.228581223023614 |
| RA_Peters | S1 | 11 | -2.63785951827066 | -1.47709132755082 | 0.1704372325043 | 0.3408744650086 |
| RA_GTEx | S1 | 11 | 0.0517913785462715 | 0.038238616846036 | 0.970250042078384 | 0.970250042078384 |

`count2FPKM`은 유전자 길이 결측 1.081%를 경고로 남겼다. batch-2 쌍 수는 0. combat-vst 파일은 없고 내려받지 않았다. 두 분기 모두 `stopped.tsv`는 `complete`.

### 6.3 독립 검사 (값별 차이)

statsmodels `MixedLM`, `reml=False`, optimizer `bfgs`. LBFGS는 RA_default에서 로그우도가 유한하지 않았다. 기준은 LRT 상대 차이 < 1e-2, t 절대 차이 < 1e-8이다.

| 자료 | clock | R | Python | 차이 |
|---|---|---|---|---|
| R1 LRT | TAms_chrono | 88.6948321165556 | 88.6948321167265 | 1.9268285956101325e-12 |
| R1 LRT | TAms_mort | 66.5958762044224 | 66.59587612001633 | 1.2674368553051373e-09 |
| R1 LRT | RA_default | 1.11800352518344 | 1.118003525874883 | 6.18462383528543e-10 |
| R3 t, 두 분기 | TAms_chrono | 1.04370691660465 | 1.043706916604641 | 8.881784197001252e-15 |
| R3 t, 두 분기 | TAms_mort | 3.3329176159176 | 3.3329176159176064 | 6.217248937900877e-15 |
| R3 t, 두 분기 | RA_default | 1.73013215149272 | 1.7301321514927264 | 6.439293542825908e-15 |

### 6.4 잠금 장치 시험

`strace -f -e openat`으로 `load_locked_expression("L1")`과 `load_locked_expression("L2")`를 호출했다. 종료 코드 1, 계열 행렬 경로의 `openat` 0회. stderr는 `config/unlock_approved.txt is missing; locked expression stays closed`.

## 7. 동결 (manifest 표, 동결 해시, 커밋, ots, OSF 꾸러미)

`prereg/manifest.tsv`는 26개 파일의 경로, 바이트, SHA-256이다. 사전등록 문서, `implementation_map.tsv`, `prereg/frozen/*.tsv` 14개, GENCODE GTF와 transcript 표, `scripts/confirm/` 10개 파일이 들어 있다. manifest 자체 SHA-256이 동결 해시 `f19f79c1988d75e7085a3f570b80a96a2e3350b368335ac80b0ae4fb162120ba`다. 커밋 `8a203ffdc5a5f1e933cd4f64154c0d5972480e32`, 시각 2026-10-16 16:21:17 +0900. ots는 없음. `prereg/osf_package/TOD_osf_package.zip`에는 사전등록 문서, manifest, `implementation_map.tsv`, `confirm_code.tar.gz`가 있다. zip SHA-256 `e9a1c31ce9d50f54587467a5dc7140f6eaba0430a4fec3416183cc85d650cfe8`, 바이트 25391. 업로드하지 않았다.

## 8. 버그·재실행 기록 (없으면 "없음")

- Part B 커밋 `553a2bf21762aa37d9742d1d2c86af2b2906413e`에 GENCODE gzip이 받기가 끝나기 전에 들어갔다. Part C에서 완전한 파일로 바꿨다. 완전한 파일의 `gzip -t`는 통과했고 SHA-256은 §5.4와 같다.
- 점수 리허설 첫 실행은 1.41초에 종료 코드 1이었다. `set.seed`에 문자 id가 들어갔고, 설계 수 표가 `stats::cov`를 피험자 표로 읽었다. 고친 뒤의 실행이 §2의 0:08:37.07이다.
- R2 첫 실행은 다중 기호 프로브셋을 빼서 PER1이 `gene_absent`였다. 경과 0:06.95, 최대 RSS 676552 kB, 종료 코드 0, `gate fail`. 6.2.10 이후 재실행이 §6.2의 R2다.
- R3 첫 실행은 종료 코드 1, 경과 0:03:46.43, 최대 RSS 1217356 kB. `count2FPKM` 결과가 행렬이 아니어서 `crossprod`가 멈췄다. 둘째 실행은 종료 코드 1, 경과 0:03:41.94, 최대 RSS 1261536 kB. `q_holm` 대입 길이가 달랐다. `set()`으로 고친 셋째 실행이 §2의 R3다.

## 9. 판단이 필요한 사항

1. 6.2.1. L2 관문 18–26시간을 24로 나눈 나머지는 peak가 18 이상이거나 2 이하일 때 통과로 읽었다.
2. 6.2.2. GENCODE 식별자는 첫 점 앞을 결합키로 쓴다. 버전 있는 식별자도 표에 남긴다.
3. 6.2.3. 다중 SYMBOL 프로브셋을 빼는 읽기. 6.2.10이 철회했다.
4. 6.2.4. Holm은 그 실행이 실제로 수행한 검정에 적용한다. 리허설은 자료 하나다. 관문 실패는 그 자료의 검정을 빼는 기록이다.
5. 6.2.5. L2 시각은 `wake_ref_h`를 24로 나눈 나머지다. F6은 다시 계산하지 않았다.
6. 6.2.6. 확인 코드의 성분은 가중 점수의 부분합이다. F3 숫자는 5단계 값을 유지한다.
7. 6.2.7. 점수 기준 표본은 그 자료의 분석 표본 전체다.
8. 6.2.8. RNA-seq에서 MCPcounter 입력은 log2(FPKM+1)이다. 관문 유전자는 log2(CPM+1)이다.
9. 6.2.9. `confirm_root`는 `scripts/lib_cosinor.R`가 있는 디렉터리까지 부모를 찾는다.
10. 6.2.10. 6.2.3을 철회했다. clariom 데이터베이스 PROBEID는 24351개이고 PER1 프로브셋은 1개이며 그 프로브셋의 SYMBOL이 둘 이상이다. L2는 기호마다 평균이 가장 큰 프로브셋을 고르고, 다중 기호 프로브셋도 각 기호의 후보로 남긴다. 다중 기호 프로브셋 수는 audit에 남긴다.
11. 6.3.1. H2의 p < 0.05 비율은 비우고 `p_med` 편향을 보고한다. 심은 구성은 `alpha = 0.6 * A`, `C = cos + N(0, 0.5)`이다. null은 alpha 0이다. H4 null은 F1 GSE113883 진폭을 양일에 두고, planted는 Day 4 또는 period 2 진폭을 절반으로 한다. H3 잔차 SD는 개발 자료의 아침−저녁 차이 SD다.
12. 6.3.2. R2는 clariom 데이터베이스의 PROBEID 키 24351개 전부다. 프로브 행렬은 저장하지 않는다. 시드 20261015.
13. 6.3.3. 독립 MixedLM은 BFGS다. LBFGS는 RA_default에서 유한하지 않은 로그우도를 냈다.
14. 6.3.4. statsmodels는 시계용 환경에 없다. 독립 검사는 `/home/jkang/anaconda3/bin/python3`, statsmodels 0.14.6을 썼다. 시계용 환경은 바꾸지 않았다.

## 10. 지시서와 다르게 한 것

독립 검사의 Python은 시계용 환경이 아니라 §9의 14번 인터프리터다. MixedLM 최적화는 BFGS다. R2 프로브 행렬은 저장하지 않았다. H2 리허설의 p < 0.05 칸은 비어 있다. 확인 코드의 FPKM 행렬 변환과 Holm `set()` 수정은 Part D 커밋에 있다.

## 11. 산출물 목록, Part별 커밋 해시

- Part A `7ec509225e4926819c825c80d43a17fc26f3bd73`
- Part B `553a2bf21762aa37d9742d1d2c86af2b2906413e`
- Part C `e0d37b6d15193d0e787f11340318967e7d5407e2`
- Part D `967723e47b76d4094634bd52d9b056e56d9d210a`
- 동결 `8a203ffdc5a5f1e933cd4f64154c0d5972480e32`

산출. `prereg/TOD_preregistration_v1.md`, `prereg/frozen/F1.tsv`–`F6.tsv`, `prereg/frozen/gencode_v38_basic_transcripts.tsv`, `prereg/frozen/exclusions.tsv`, `prereg/frozen/clock_list.tsv`, `prereg/frozen/environment.tsv`, `prereg/implementation_map.tsv`, `prereg/manifest.tsv`, `prereg/osf_package/TOD_osf_package.zip`, `scripts/confirm/`, `results/stage6/rehearsal/scores/score_rehearsal.tsv`, `results/stage6/rehearsal/R1/`, `results/stage6/rehearsal/R2/`, `results/stage6/rehearsal/R3/`, `results/stage6/rehearsal/independent/independent_check.tsv`, `results/stage6/rehearsal/lock/lock_test.tsv`.

잠금은 풀지 않았다. `config/unlock_approved.txt`는 없다. L3 보충 파일은 받지 않았다.
