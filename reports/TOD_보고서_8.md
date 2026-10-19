# TOD 보고서 8 — 확인 결과 확정, 사후 분석, 원고 그림·표

## 1. 요약 (수치만, 12줄 이내)

- Part A `67d4c2d7b3570f2b72aba06242706057cc4e635d` (2026-10-19 11:09:06 +0900). 사전등록 확인 결과는 실행 5 `e56a1c489a74453ca33a832835fdaa312e3e8371`의 `results/confirm/`. 이탈은 D1–D12. 8.A2는 L2 기호 BMAL1을 ARNTL 별칭으로 기록. 8.A3는 원형 SD 식을 기록.
- Part B `9317da205567e399907d98a2ae1ef21e50f4aab8` (2026-10-19 12:03:27 +0900). `scripts/51_part8_posthoc.R`의 `/usr/bin/time -v`: 경과 `46:00.17`, 최대 RSS `1560928` kB, 종료 상태 0.
- 습관적 수면 시작 시각의 평균과 범위는 확인 못 함. `time_transform`은 `not performed`.
- BMAL1 별칭: n 70/11, peak_h `21.3692664519028`, lrt_p `0.00925089859203897`. 유전자 10개의 circ_mean_h `9.69765168583892`, circ_sd_h `5.49841393915787`.
- L3 검정력, 쌍 58, 양측 α=0.05: `TAms_chrono` `0.176777740279187`, `TAms_mort` `0.291920013796665`, `RA_default` `0.686538271492044`, `RA_Peters` `0.511452502493332`, `RA_GTEx` `0.108269826724681`.
- D11 모형 M `posthoc_confirmed` TRUE 6/8. 동결 S4 L1 확인 수는 `1`/`8`.
- 그림 PDF·PNG 300 dpi: `figures/Fig2_trajectories`, `figures/Fig3A_effect_forest`, `figures/Fig3B_peak_clock`, `figures/Fig4_composition`, `figures/Fig5_modules`. 표 TSV·docx: `results/stage8/tables/Table1`, `Table2`, `Table3`, `TableS1`, `TableS2`, `TableS3`.

## 2. 확정 기록

`config/decisions.tsv`의 2026-10-19 행이다.

| id | content |
|---|---|
| 8.A1 | The preregistered confirmation result is run 5, results/confirm/, commit e56a1c489a74453ca33a832835fdaa312e3e8371. The deviations that apply are D1 through D12. |
| 8.A2 | In the L2 post-hoc clock-gene table, the L2 matrix symbol BMAL1 is used as the alias of ARNTL and the circular summary is recomputed. The recomputed table is post hoc only. results/confirm/ is not changed. |
| 8.A3 | The circular SD formula in the stage-7c post-hoc table is approved: sqrt(-2*log(R))*24/(2*pi), where R is the mean resultant length of the signed differences. |

BMAL1 적합 (`results/stage8/l2_bmal1_as_arntl.tsv`):

| gene | n_samples | n_subjects | peak_h | lrt_p | status | alias_of | title |
| --- | --- | --- | --- | --- | --- | --- | --- |
| BMAL1 | 70 | 11 | 21.3692664519028 | 0.00925089859203897 | fit | ARNTL | L2 BMAL1 used as ARNTL post hoc, exploratory |

원형 요약 (`results/stage8/gene_peak_circular_summary_arntl_alias.tsv`):

| n_genes | circ_mean_h | circ_sd_h | signed_diff_definition | circ_sd_definition | arntl_source | title |
| --- | --- | --- | --- | --- | --- | --- |
| 10 | 9.69765168583892 | 5.49841393915787 | (peak_l2 - peak_gse39445) mod 24, then subtract 24 when the result is greater than 12; range (-12, 12] | sqrt(-2*log(R)) * 24 / (2*pi), R = mean resultant length of the signed differences | L2 symbol BMAL1; GSE39445 symbol ARNTL | circular mean and SD after ARNTL alias post hoc, exploratory |

## 3. 사후 분석

### 3.1 L2 시각 기준 보정 (인용, 변환, 사후 관문, 정점 거리)

출처는 Kervezee et al., J Biol Rhythms, PMID 30621487, DOI 10.1177/0748730418821776. 습관적 수면 시작 시각의 평균과 범위는 확인 못 함. `time_transform`은 `not performed`. 사후 관문과 보정 축의 정점 거리는 계산하지 않았다. Fig 3B에는 L2를 넣지 않았다.

원문:

- Briefly, participants maintained a regular, selfselected sleep/wake schedule with a nocturnal sleep period of 8 h for at least one week before admission to the laboratory, which was verified by daily sleep logs and actigraphy.
- Upon admission to the laboratory, participants were involved in a 24-h baseline measurement period that consisted of an 8-h sleep period during their habitual sleep times and a 16-h constant posture procedure, involving a constant semi-recumbent position, minimal activity, and hourly isocaloric snacks.
- relative clock time refers to the time of day relative to the start of each participant's habitual sleep period.

검색 표:

| source | citation | sentence | habitual_sleep_start_mean_range | reason | time_transform | title |
| --- | --- | --- | --- | --- | --- | --- |
| Kervezee et al. 2019 methods | Kervezee et al., J Biol Rhythms, PMID 30621487, DOI 10.1177/0748730418821776 | Briefly, participants maintained a regular, selfselected sleep/wake schedule with a nocturnal sleep period of 8 h for at least one week before admission to the laboratory, which was verified by daily sleep logs and actigraphy. | 확인 못 함 | The paper states an 8 h habitual sleep period and defines relative clock time from its start. It does not state a clock time mean or range for habitual sleep onset. Protocol details are referred to Cuesta et al., 2017. | not performed | L2 time-axis source post hoc, exploratory |
| Kervezee et al. 2019 methods | Kervezee et al., J Biol Rhythms, PMID 30621487, DOI 10.1177/0748730418821776 | Upon admission to the laboratory, participants were involved in a 24-h baseline measurement period that consisted of an 8-h sleep period during their habitual sleep times and a 16-h constant posture procedure, involving a constant semi-recumbent position, minimal activity, and hourly isocaloric snacks. | 확인 못 함 | The paper states an 8 h habitual sleep period and defines relative clock time from its start. It does not state a clock time mean or range for habitual sleep onset. Protocol details are referred to Cuesta et al., 2017. | not performed | L2 time-axis source post hoc, exploratory |
| Kervezee et al. 2019 figure 1 legend | Kervezee et al., J Biol Rhythms, PMID 30621487, DOI 10.1177/0748730418821776 | relative clock time refers to the time of day relative to the start of each participant's habitual sleep period. | 확인 못 함 | The paper states an 8 h habitual sleep period and defines relative clock time from its start. It does not state a clock time mean or range for habitual sleep onset. Protocol details are referred to Cuesta et al., 2017. | not performed | L2 time-axis source post hoc, exploratory |

### 3.2 자료 전체 효과 요약 표

개발 자료 3개와 L2의 코사이너 구간은 `scripts/51_part8_posthoc.R`가 확인 코드 `boot_effects`로 계산한 값이다. 2A/σ는 `two_a/sigma_reml`이다. 부트스트랩은 참가자 재추출 2000회, seed `20261015`, 백분위 type 7이다. 개발 15개 시계의 `two_a`와 `prereg/frozen/F1.tsv`의 `two_a` 상대 차이의 최댓값은 `1.9902099137501764e-09`이다. 네 자료 모두 `n_converged`는 2000이다. L1과 L3의 점추정과 L1 구간, L3의 mean Δ 구간과 p는 `results/confirm/family_P.tsv`와 `family_S1.tsv`의 문자열이다. L3의 d 구간은 `results/stage8/l3_d_vs_prediction.tsv`이다. 표를 합치지 않았다.

코사이너. `two_a_over_sigma`의 구간은 `ratio_p025`, `ratio_p975`. `r2_within`의 구간은 `r2_p025`, `r2_p975`. 정점 구간은 `peak_center`와 `peak_halfwidth_95`이다.

| dataset | clock | role | n_subjects | n_samples | two_a_over_sigma | ratio_p025 | ratio_p975 | r2_within | r2_low | r2_high | peak_h | peak_center | peak_halfwidth_95 | p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | TAms_chrono | development | 24 | 221 | 2.01559986979017 | 1.50682441272597 | 2.68093603216236 | 0.362193176953973 | 0.23947737962447 | 0.501480783091199 | 14.2335397784943 | 14.2186064321564 | 0.820399743143368 | 5.49745798912395e-20 |
| GSE39445 | TAms_mort | development | 24 | 221 | 1.72804507043291 | 1.37105685104705 | 2.23430522926863 | 0.28699925739554 | 0.1997711867053 | 0.40680716136472 | 13.4288594848278 | 13.4014342140422 | 1.16438541795076 | 3.45851149483147e-15 |
| GSE39445 | RA_default | development | 24 | 221 | 0.2125210689102 | 0.0684349337069235 | 0.758860045770011 | 0.00680181474094321 | 0.000576531349617138 | 0.0704405878265855 | 12.0064473932843 | 12.4584652559485 | 8.83488247674881 | 0.571779550797631 |
| GSE39445 | RA_Peters | development | 24 | 221 | 0.424138432819158 | 0.145714846186837 | 0.79419196216918 | 0.0237318052396572 | 0.00260828051439757 | 0.0793917895796013 | 6.84503499428336 | 6.70958325579563 | 3.82180804362507 | 0.0822398090507771 |
| GSE39445 | RA_GTEx | development | 24 | 221 | 0.181325286706445 | 0.0653971731430674 | 0.700835756775039 | 0.00480935568643748 | 0.000539049327682534 | 0.0655700363537038 | 13.7137336093387 | 13.6618778041425 | 9.48187052265864 | 0.643663244864561 |
| GSE56931 | TAms_chrono | development | 14 | 130 | 0.477418482127208 | 0.127260837991756 | 1.29387535750537 | 0.0286364463384637 | 0.00208177552906286 | 0.173287788864715 | 11.4579688514277 | 11.2197113952556 | 7.28242031911783 | 0.181165188487945 |
| GSE56931 | TAms_mort | development | 14 | 130 | 0.850468368966173 | 0.368231277608618 | 1.65873262933755 | 0.0821416660164724 | 0.0173816008067029 | 0.2609873466887 | 12.0688456403715 | 12.1544296381706 | 3.06206316312207 | 0.00673707816031792 |
| GSE56931 | RA_default | development | 14 | 130 | 0.690328543295956 | 0.301519143354035 | 1.18386540924029 | 0.0658778147622046 | 0.0123684442626431 | 0.168237969754778 | 9.88701727265618 | 9.77493370481417 | 2.411697769131 | 0.0187399469939367 |
| GSE56931 | RA_Peters | development | 14 | 130 | 1.14177298602567 | 0.823758107326477 | 1.57754107188006 | 0.145491804603583 | 0.0798925712241155 | 0.252084121730364 | 4.54192615314065 | 4.54819791251166 | 1.53720197882062 | 0.00010511653009091 |
| GSE56931 | RA_GTEx | development | 14 | 130 | 0.592966981425032 | 0.317722565787296 | 1.04422829814623 | 0.0366339113390549 | 0.0116063389326546 | 0.111998247627982 | 13.2671292169983 | 13.2821882675504 | 2.73432939449719 | 0.0977660921595171 |
| GSE113883 | TAms_chrono | development | 11 | 153 | 1.65417866072112 | 1.24631116203654 | 2.26147658965948 | 0.265399481367344 | 0.165601193416533 | 0.404691023457529 | 12.188050742382 | 12.2120885405226 | 1.30484190182366 | 2.96936394371969e-10 |
| GSE113883 | TAms_mort | development | 11 | 153 | 1.14231965933045 | 0.731497214201511 | 1.77961485125071 | 0.158916853730887 | 0.0696695591092995 | 0.312364575276871 | 10.0600990696032 | 10.0668531335917 | 2.22504509606368 | 4.44466652817908e-06 |
| GSE113883 | RA_default | development | 11 | 153 | 1.07645532806507 | 0.726629518999648 | 1.57921635089029 | 0.134425171735901 | 0.063614901229032 | 0.247718427378464 | 11.812631972311 | 11.8212982069481 | 1.44667958366095 | 3.24675832907807e-05 |
| GSE113883 | RA_Peters | development | 11 | 153 | 0.521739262264997 | 0.124662755016569 | 1.10322302796303 | 0.0352688623430064 | 0.00202148111120034 | 0.13890527185579 | 11.8362208570625 | 11.6741391650501 | 4.62490864528075 | 0.0773443137299791 |
| GSE113883 | RA_GTEx | development | 11 | 153 | 0.492489390834316 | 0.254651161294862 | 1.18952772966499 | 0.0315491728014815 | 0.00905836853320457 | 0.144632046247262 | 6.16828794535857 | 6.73714594981626 | 4.82078466754032 | 0.104384139695734 |
| L2 | TAms_chrono | post hoc | 11 | 70 | 0.955955620895999 | 0.609633260585316 | 1.54933916967279 | 0.126783501953014 | 0.0516852449557172 | 0.272480085279994 | 10.800630536763 | 10.6691771391267 | 2.6741144433617 | 0.0177344777316928 |
| L2 | TAms_mort | post hoc | 11 | 70 | 0.954078993760188 | 0.389231414563289 | 1.96642479136579 | 0.125546783188962 | 0.0228342911254596 | 0.378907332292164 | 10.9821548000952 | 10.9995786263699 | 3.47941809422713 | 0.0174608076774516 |
| L2 | RA_default | post hoc | 11 | 70 | 1.07219562404825 | 0.346623245997025 | 2.09518632559289 | 0.136499241214436 | 0.0166364765606075 | 0.370305583281508 | 8.29862480973111 | 8.51106573293619 | 2.0806503386807 | 0.0116625636160948 |
| L2 | RA_Peters | post hoc | 11 | 70 | 1.34606228166385 | 0.718236299251354 | 2.29714268449238 | 0.206353947286669 | 0.068226216886336 | 0.4367403473941 | 9.718280100379 | 9.72709169400464 | 2.24431173078793 | 0.000688235638502205 |
| L2 | RA_GTEx | post hoc | 11 | 70 | 1.39064762099237 | 0.78612637238873 | 2.16396450187175 | 0.211532946980103 | 0.085730680048668 | 0.390401196331784 | 7.24447680202368 | 7.30471326906638 | 1.48670379924488 | 0.00107043094705134 |
| L1 | TAms_chrono | confirmation | 22 | 147 | 1.47644030704959 | 1.03384190558836 | 2.03769515562535 | 0.239604998243655 | 0.133334366044785 | 0.377239951848059 | 10.879869555482 | 10.855130663467 | 1.31149483591321 | 2.70474257565244e-08 |
| L1 | TAms_mort | confirmation | 22 | 147 | 1.1381830590621 | 0.700981308705566 | 1.76922760462403 | 0.156244931376745 | 0.0679705423411425 | 0.305423320613559 | 11.4858304321163 | 11.4467106578545 | 2.03534016526948 | 2.67642593880552e-05 |
| L1 | RA_default | confirmation | 22 | 147 | 0.716426862045517 | 0.315800747882794 | 1.31329139603625 | 0.066357183458594 | 0.0136587510552551 | 0.20321009677178 | 11.5918483633937 | 11.5611538419167 | 2.79211694956184 | 0.0125285842066838 |
| L1 | RA_Peters | confirmation | 22 | 147 | 0.72926928693365 | 0.282409354356446 | 1.30560038410114 | 0.0671837397484756 | 0.0111381513630943 | 0.185619898185969 | 3.52591635753771 | 3.66899071507773 | 2.95933916983926 | 0.0131608178088439 |
| L1 | RA_GTEx | confirmation | 22 | 147 | 0.422296335991481 | 0.135450842918859 | 0.961108073220587 | 0.0245786114049672 | 0.0025435807168417 | 0.115501532769887 | 16.2142152061473 | 16.2659536981038 | 5.36702827578855 | 0.214719525083053 |

L3. d의 구간은 `d_p025`, `d_p975`. mean Δ의 구간은 `mean_low`, `mean_high`.

| dataset | clock | role | n_subjects | n_samples | d | d_p025 | d_p975 | mean_delta | mean_low | mean_high | p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| L3 | TAms_chrono | confirmation | 58 | 116 | 0.137127542568716 | -0.134951777137601 | 0.423111690763127 | 0.0120385849794561 | -0.0110449218090888 | 0.0351220917680009 | 0.300740570646651 |
| L3 | TAms_mort | confirmation | 58 | 116 | 0.188455112293294 | -0.072646855125881 | 0.472594447464266 | 0.0827171714412529 | -0.0326915922512471 | 0.198125935133753 | 0.156685969793716 |
| L3 | RA_default | confirmation | 58 | 116 | 0.326704860352488 | 0.0837673566577421 | 0.60081434729321 | 2.28347265426549 | 0.445702236643247 | 4.12124307188773 | 0.0157896851797883 |
| L3 | RA_Peters | confirmation | 58 | 116 | 0.265595107143801 | 0.0142213143420831 | 0.523158041348489 | 2.86233409356227 | 0.0286509268443877 | 5.69601726028016 | 0.047803206296958 |
| L3 | RA_GTEx | confirmation | 58 | 116 | 0.0938927716485769 | -0.171819370151669 | 0.350465488181871 | 0.713702823968326 | -1.2849450811811 | 2.71235072911775 | 0.477486960820646 |

### 3.3 L3 관측 효과와 예측

관측 d는 참가자 58명의 morning−evening에 대한 부트스트랩 2000회, seed `20261015`, quantile type 7이다. 검정력은 paired t, df = n−1, ncp = d*sqrt(n), 양측 α=0.05이다. `d_exp_median`과 `d_GSE113883`은 `prereg/frozen/F2.tsv`이다.

| clock | n_pairs | d | d_p025 | d_p975 | d_exp_median | d_GSE113883 | predicted_sign | power_alpha_0.05 | power_definition | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TAms_chrono | 58 | 0.137127542568716 | -0.134951777137601 | 0.423111690763127 | 0.208443813199012 | 0.534281736707176 | + | 0.176777740279187 | paired t, df=n-1, ncp=d*sqrt(n), two-sided alpha=0.05 | L3 observed d and development prediction post hoc, exploratory |
| TAms_mort | 58 | 0.188455112293294 | -0.072646855125881 | 0.472594447464266 | 0.291250817026961 | 0.693084284817919 | + | 0.291920013796665 | paired t, df=n-1, ncp=d*sqrt(n), two-sided alpha=0.05 | L3 observed d and development prediction post hoc, exploratory |
| RA_default | 58 | 0.326704860352488 | 0.0837673566577421 | 0.60081434729321 | 0.412448894924069 | 0.412448894924069 | + | 0.686538271492044 | paired t, df=n-1, ncp=d*sqrt(n), two-sided alpha=0.05 | L3 observed d and development prediction post hoc, exploratory |
| RA_Peters | 58 | 0.265595107143801 | 0.0142213143420831 | 0.523158041348489 | 0.286305297017999 | 0.197988187071947 | + | 0.511452502493332 | paired t, df=n-1, ncp=d*sqrt(n), two-sided alpha=0.05 | L3 observed d and development prediction post hoc, exploratory |
| RA_GTEx | 58 | 0.0938927716485769 | -0.171819370151669 | 0.350465488181871 | 0.0799546478017167 | 0.308963171982023 | + | 0.108269826724681 | paired t, df=n-1, ncp=d*sqrt(n), two-sided alpha=0.05 | L3 observed d and development prediction post hoc, exploratory |

`RA_default`의 d는 이 표에서 `0.326704860352488`이고 `results/confirm/family_P.tsv`에서 `0.326704860352489`이다.

### 3.4 L1 모듈 결과의 D11 판

D11 점수는 `results/confirm/L1/tage_scores_d11.tsv`와 `tage_components_d11.tsv`의 period D1, 참가자 22, 샘플 147이다. 기준은 LRT p < 0.05이고 φ̄까지 원형 거리가 3시간 이내이다.

| clock | n_subjects | n_samples | two_a | peak_h | lrt_stat | p | phi_bar | circ_dist | posthoc_confirmed | criterion | status | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TM_turquoise_mortality | 22 | 147 | 0.255393380643979 | 13.0160343833073 | 28.4363179329579 | 6.68547113479287e-07 | 11.5724011994242 | 1.4436331838831 | TRUE | LRT p < 0.05 and circular distance to phi_bar <= 3 h | fit | L1 F3 outcomes on D11 scores post hoc, exploratory |
| CMP_TAms_chrono_turquoise | 22 | 147 | 0.0021030946520009 | 9.41647658790674 | 1.52390038187855 | 0.466755276928446 | 13.2837328415609 | 3.86725625365416 | FALSE | LRT p < 0.05 and circular distance to phi_bar <= 3 h | fit | L1 F3 outcomes on D11 scores post hoc, exploratory |
| CMP_TAms_chrono_yellow | 22 | 147 | 0.000782866342744818 | 10.7776829888493 | 8.30580812044627 | 0.015718702080819 | 12.6128269162576 | 1.8351439274083 | TRUE | LRT p < 0.05 and circular distance to phi_bar <= 3 h | fit | L1 F3 outcomes on D11 scores post hoc, exploratory |
| CMP_TAms_chrono_orange | 22 | 147 | 0.00382671362552297 | 2.62378198893345 | 43.0815069562091 | 4.41539568022842e-10 | 1.07332762742524 | 1.55045436150821 | TRUE | LRT p < 0.05 and circular distance to phi_bar <= 3 h | fit | L1 F3 outcomes on D11 scores post hoc, exploratory |
| CMP_TAms_chrono_brown | 22 | 147 | 0.00208036672250036 | 0.0401997235020163 | 25.5874382289953 | 2.77816723822928e-06 | 1.32505139804729 | 1.28485167454527 | TRUE | LRT p < 0.05 and circular distance to phi_bar <= 3 h | fit | L1 F3 outcomes on D11 scores post hoc, exploratory |
| CMP_TAms_chrono_red | 22 | 147 | 0.00038041525453051 | 2.49280571765594 | 17.2648672748608 | 0.000178230368847034 | 0.781691026210673 | 1.71111469144526 | TRUE | LRT p < 0.05 and circular distance to phi_bar <= 3 h | fit | L1 F3 outcomes on D11 scores post hoc, exploratory |
| CMP_TAms_mort_turquoise | 22 | 147 | 0.00983352590115993 | 13.0566980605596 | 3.33111022264541 | 0.18908566524961 | 13.003811149926 | 0.0528869106335694 | FALSE | LRT p < 0.05 and circular distance to phi_bar <= 3 h | fit | L1 F3 outcomes on D11 scores post hoc, exploratory |
| CMP_TAms_mort_darkred | 22 | 147 | 0.003662200552302 | 15.4115823911225 | 8.43470140444924 | 0.0147376372199395 | 14.9666596497573 | 0.444922741365206 | TRUE | LRT p < 0.05 and circular distance to phi_bar <= 3 h | fit | L1 F3 outcomes on D11 scores post hoc, exploratory |

`posthoc_confirmed` TRUE는 6행, FALSE는 2행이다.

특징 수 (`results/stage8/f3_feature_counts.tsv`). 유전자 이름은 적지 않았다.

| version | clock | n_entered | title |
| --- | --- | --- | --- |
| frozen | TM_turquoise_mortality | 14 | post hoc, exploratory |
| frozen | CMP_TAms_chrono_turquoise | 4 | post hoc, exploratory |
| frozen | CMP_TAms_chrono_yellow | 1 | post hoc, exploratory |
| frozen | CMP_TAms_chrono_orange | 1 | post hoc, exploratory |
| frozen | CMP_TAms_chrono_brown | 0 | post hoc, exploratory |
| frozen | CMP_TAms_chrono_red | 0 | post hoc, exploratory |
| frozen | CMP_TAms_mort_turquoise | 15 | post hoc, exploratory |
| frozen | CMP_TAms_mort_darkred | 1 | post hoc, exploratory |
| D2_knn | TM_turquoise_mortality | 200 | post hoc, exploratory |
| D2_knn | CMP_TAms_chrono_turquoise | 58 | post hoc, exploratory |
| D2_knn | CMP_TAms_chrono_yellow | 5 | post hoc, exploratory |
| D2_knn | CMP_TAms_chrono_orange | 4 | post hoc, exploratory |
| D2_knn | CMP_TAms_chrono_brown | 5 | post hoc, exploratory |
| D2_knn | CMP_TAms_chrono_red | 1 | post hoc, exploratory |
| D2_knn | CMP_TAms_mort_turquoise | 201 | post hoc, exploratory |
| D2_knn | CMP_TAms_mort_darkred | 15 | post hoc, exploratory |
| D11_participant_mean | TM_turquoise_mortality | 202 | post hoc, exploratory |
| D11_participant_mean | CMP_TAms_chrono_turquoise | 58 | post hoc, exploratory |
| D11_participant_mean | CMP_TAms_chrono_yellow | 7 | post hoc, exploratory |
| D11_participant_mean | CMP_TAms_chrono_orange | 4 | post hoc, exploratory |
| D11_participant_mean | CMP_TAms_chrono_brown | 6 | post hoc, exploratory |
| D11_participant_mean | CMP_TAms_chrono_red | 1 | post hoc, exploratory |
| D11_participant_mean | CMP_TAms_mort_turquoise | 204 | post hoc, exploratory |
| D11_participant_mean | CMP_TAms_mort_darkred | 17 | post hoc, exploratory |

D11 행렬을 같은 규칙으로 다시 센 감사 (`results/stage8/d11_matrix_audit_recount.tsv`):

| n_missing | n_filled_participant | n_filled_gene_mean | n_nonfinite_after |
| --- | --- | --- | --- |
| 1002419 | 645943 | 356476 | 8265 |

## 4. 그림 (파일 목록, 그림별 원자료와 n)

PDF와 PNG는 300 dpi이다. 제목에 post hoc, exploratory가 있다. 개발은 녹색 원, 확인은 주황 사각형, 사후는 보라 다이아몬드이다. 주 시계 이름은 tAge chronological, tAge mortality, RNAAgeCalc blood이다.

| 파일 | 내용 |
|---|---|
| figures/Fig2_trajectories.pdf, figures/Fig2_trajectories.png | GSE113883과 L1 Day 1의 참가자 궤적과 집단 코사이너. x는 벽시계 시각. |
| figures/Fig3A_effect_forest.pdf, figures/Fig3A_effect_forest.png | §3.2. 코사이너 자료는 2A/σ, L3는 d. |
| figures/Fig3B_peak_clock.pdf, figures/Fig3B_peak_clock.png | 주 시계 3개의 정점과 개발 φ̄. L2는 없다. |
| figures/Fig4_composition.pdf, figures/Fig4_composition.png | p_med와 MCPcounter 8개 집단의 정점 시각. |
| figures/Fig5_modules.pdf, figures/Fig5_modules.png | F3 8개 결과의 개발 φ̄. 바깥 원 L1 동결 확인, 가운데 사각형 L1 D11 사후 확인, 안쪽 다이아몬드 L3 부호 일치. 회색은 해당 상태가 아님. |
| figures/figure_sources.tsv | 그림별 원자료와 n. |

| figure | panel | n_subjects | n_samples | role | source |
| --- | --- | --- | --- | --- | --- |
| Fig2 | GSE113883 TAms_chrono | 11 | 153 | development | results/stage5/clock_scores.tsv; results/stage4/rnaseq/clock_scores.tsv; results/stage3/GSE113883_samples.tsv |
| Fig2 | GSE113883 TAms_mort | 11 | 153 | development | results/stage5/clock_scores.tsv; results/stage4/rnaseq/clock_scores.tsv; results/stage3/GSE113883_samples.tsv |
| Fig2 | GSE113883 RA_default | 11 | 153 | development | results/stage5/clock_scores.tsv; results/stage4/rnaseq/clock_scores.tsv; results/stage3/GSE113883_samples.tsv |
| Fig2 | L1 Day 1 TAms_chrono | 22 | 147 | confirmation | results/confirm/L1/tage_scores.tsv; results/confirm/L1/ra_scores.tsv; results/confirm/L1/samples_scoring.tsv; results/confirm/family_P.tsv |
| Fig2 | L1 Day 1 TAms_mort | 22 | 147 | confirmation | results/confirm/L1/tage_scores.tsv; results/confirm/L1/ra_scores.tsv; results/confirm/L1/samples_scoring.tsv; results/confirm/family_P.tsv |
| Fig2 | L1 Day 1 RA_default | 22 | 147 | confirmation | results/confirm/L1/tage_scores.tsv; results/confirm/L1/ra_scores.tsv; results/confirm/L1/samples_scoring.tsv; results/confirm/family_P.tsv |
| Fig3A | 2A/sigma and d | GSE39445:24/221; GSE56931:14/130; GSE113883:11/153; L2:11/70; L1:22/147; L3:58/116 | subjects/samples in n_subjects | mixed | results/stage8/effect_summary.tsv |
| Fig3B | primary peaks, L2 omitted | GSE39445:24/221; GSE56931:14/130; GSE113883:11/153; L1:22/147; L3:58/116 | subjects/samples in n_subjects | mixed | prereg/frozen/F1.tsv; results/stage8/effect_summary.tsv |
| Fig4 | p_med and cell peaks | GSE39445 24/221; GSE56931 14/130; GSE113883 11/153; L1 H2 146 samples | subjects/samples in n_subjects | development and confirmation | results/stage5/projection_ci.tsv; results/stage3/h2/projection_ci.tsv; results/stage4/rnaseq/projection_ci.tsv; results/confirm/h2.tsv; results/stage2/h2_population_cosinor.tsv; results/stage4/rnaseq/h2_population_cosinor.tsv |
| Fig5 | F3 | L1 22; L3 58 | L1 147; L3 116 | confirmation and post hoc | prereg/frozen/F3.tsv; results/confirm/s4.tsv; results/stage8/l1_f3_d11_model_m.tsv |

Fig 2의 집단 곡선은 저장된 `beta_cos`, `beta_sin`으로 계산한 파동을 점수에서 뺀 값의 평균을 세로 위치로 둔 곡선이다.

## 5. 표 (파일 목록, Table 1–3 본문 수록)

| 파일 | 내용 |
|---|---|
| results/stage8/tables/Table1.tsv, Table1.docx | 자료 8개 |
| results/stage8/tables/Table2.tsv, Table2.docx | 관문, 가족 P, 재현, H1, S1 |
| results/stage8/tables/Table3.tsv, Table3.docx | H2, H3, H4, H5, S4 |
| results/stage8/tables/TableS1.tsv, TableS1.docx | D1–D12 |
| results/stage8/tables/TableS2.tsv, TableS2.docx | 사전등록 §8과 D1, D2, D6, D9, D11 민감도 |
| results/stage8/tables/TableS3.tsv, TableS3.docx | §3.1–§3.4와 L2 H1·H4 탐색 |

### Table 1

| dataset | role | description | tissue | n_participants | n_samples | design | time_axis |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | development | Agilent microarray, leukocytes, extended wakefulness | leukocytes | 24 | 221 | repeated measures | clock_time_h |
| GSE56931 | development | Affymetrix microarray, whole blood, baseline control | whole blood | 14 | 130 | repeated measures | clock_time_h |
| GSE113883 | development | RNA-seq TPM, whole blood, sampling every 2 h | whole blood | 11 | 153 | repeated measures | clock_time_h |
| GSE48113 | confirmation L1 | Microarray GPL15331, forced desynchrony, Day 1 | leukocytes | 22 | 147 | repeated measures | clock_time_h |
| GSE122725 | confirmation L2, gate fail | Microarray GPL24324, period 1 day-oriented schedule | PBMC | 11 | 70 | repeated measures | wake_ref_h %% 24 |
| GSE295758 | confirmation L3 | RNA-seq raw counts, morning and evening | whole blood | 58 | 116 | paired | morning 8.0, evening 20.0 |
| GSE196728 | excluded | RNA-seq. Sea-level control subset. Decision 4.3 excludes it from later analyses. | whole blood | 8 | 48 | repeated measures | clock_time_h: 2, 6, 10, 14, 18, 22 |
| GSE150337 | held | RNA-seq. Decision 3.2.GSE150337: not opened. | 확인 못 함 | 확인 못 함 | 확인 못 함 | 확인 못 함 | 확인 못 함 |

### Table 2

관문:

| dataset | gene | lrt_p | peak_h | pass | reason | dataset_pass | mean_delta | p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| L1 | PER1 | 1.35631359279823e-11 | 7.39244745511937 | TRUE |  | TRUE | nan | nan |
| L1 | PER3 | 2.79164389620025e-06 | 4.93116485780503 | TRUE |  | TRUE | nan | nan |
| L2 | PER1 | 1.33395937684947e-05 | 8.62707238405378 | FALSE |  | FALSE | nan | nan |
| L2 | PER3 | 2.56672604887414e-06 | 6.37573752707371 | FALSE |  | FALSE | nan | nan |
| L3 | PER1 | nan | nan | TRUE | nan | TRUE | 0.565581638383806 | 0.00672475968762667 |

가족 P:

| dataset | clock | n_subjects | n_samples | two_a_over_sigma | ratio_p025 | ratio_p975 | r2_within | r2_p025 | r2_p975 | peak_h | peak_center | peak_halfwidth_95 | p | q_holm | decision | d | mean_delta | ci_low | ci_high |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| L1 | TAms_chrono | 22 | 147 | 1.47644030704959 | 1.03384190558836 | 2.03769515562535 | 0.239604998243655 | 0.133334366044785 | 0.377239951848059 | 10.879869555482 | 10.855130663467 | 1.31149483591321 | 2.70474257565244e-08 | 1.62284554539147e-07 | significant |  |  |  |  |
| L1 | TAms_mort | 22 | 147 | 1.1381830590621 | 0.700981308705566 | 1.76922760462403 | 0.156244931376745 | 0.0679705423411425 | 0.305423320613559 | 11.4858304321163 | 11.4467106578545 | 2.03534016526948 | 2.67642593880552e-05 | 0.000133821296940276 | significant |  |  |  |  |
| L1 | RA_default | 22 | 147 | 0.716426862045517 | 0.315800747882794 | 1.31329139603625 | 0.066357183458594 | 0.0136587510552551 | 0.20321009677178 | 11.5918483633937 | 11.5611538419167 | 2.79211694956184 | 0.0125285842066838 | 0.0501143368267352 | not significant |  |  |  |  |
| L2 | TAms_chrono |  |  |  |  |  |  |  |  |  |  |  |  |  | gate fail |  |  |  |  |
| L2 | TAms_mort |  |  |  |  |  |  |  |  |  |  |  |  |  | gate fail |  |  |  |  |
| L2 | RA_default |  |  |  |  |  |  |  |  |  |  |  |  |  | gate fail |  |  |  |  |
| L3 | TAms_chrono | 58 | 116 |  |  |  |  |  |  |  |  |  | 0.300740570646651 | 0.313371939587432 | not significant | 0.137127542568716 | 0.0120385849794561 | -0.0110449218090888 | 0.0351220917680009 |
| L3 | TAms_mort | 58 | 116 |  |  |  |  |  |  |  |  |  | 0.156685969793716 | 0.313371939587432 | not significant | 0.188455112293294 | 0.0827171714412529 | -0.0326915922512471 | 0.198125935133753 |
| L3 | RA_default | 58 | 116 |  |  |  |  |  |  |  |  |  | 0.0157896851797883 | 0.0501143368267352 | not significant | 0.326704860352489 | 2.28347265426549 | 0.445702236643247 | 4.12124307188773 |

재현:

| clock | n_significant | n_gate_fail | decision |
| --- | --- | --- | --- |
| TAms_chrono | 1 | 1 | not replicated |
| TAms_mort | 1 | 1 | not replicated |
| RA_default | 0 | 1 | not replicated |

H1:

| dataset | clock | phi_bar | peak_h | circ_dist | within_3h | predicted_sign | observed_sign | sign_agrees |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| L1 | TAms_chrono | 12.616768998085 | 10.879869555482 | 1.736899442603 | TRUE |  |  |  |
| L3 | TAms_chrono |  |  |  |  | + | + | TRUE |
| L1 | TAms_mort | 11.8599407935094 | 11.4858304321163 | 0.374110361393063 | TRUE |  |  |  |
| L3 | TAms_mort |  |  |  |  | + | + | TRUE |
| L1 | RA_default | 11.2423863477089 | 11.5918483633937 | 0.349462015684843 | TRUE |  |  |  |
| L3 | RA_default |  |  |  |  | + | + | TRUE |

가족 S1:

| dataset | clock | n_subjects | n_samples | two_a_over_sigma | ratio_p025 | ratio_p975 | r2_within | r2_p025 | r2_p975 | peak_h | peak_center | peak_halfwidth_95 | p | q_holm | decision | d | mean_delta | ci_low | ci_high |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| L1 | RA_Peters | 22 | 147 | 0.72926928693365 | 0.282409354356446 | 1.30560038410114 | 0.0671837397484756 | 0.0111381513630943 | 0.185619898185969 | 3.52591635753771 | 3.66899071507773 | 2.95933916983926 | 0.0131608178088439 | 0.0526432712353756 | not significant |  |  |  |  |
| L1 | RA_GTEx | 22 | 147 | 0.422296335991481 | 0.135450842918859 | 0.961108073220587 | 0.0245786114049672 | 0.0025435807168417 | 0.115501532769887 | 16.2142152061473 | 16.2659536981038 | 5.36702827578855 | 0.214719525083053 | 0.429439050166105 | not significant |  |  |  |  |
| L2 | RA_Peters |  |  |  |  |  |  |  |  |  |  |  |  |  | gate fail |  |  |  |  |
| L2 | RA_GTEx |  |  |  |  |  |  |  |  |  |  |  |  |  | gate fail |  |  |  |  |
| L3 | RA_Peters | 58 | 116 |  |  |  |  |  |  |  |  |  | 0.047803206296958 | 0.143409618890874 | not significant | 0.265595107143803 | 2.86233409356227 | 0.0286509268443877 | 5.69601726028016 |
| L3 | RA_GTEx | 58 | 116 |  |  |  |  |  |  |  |  |  | 0.477486960820646 | 0.477486960820646 | not significant | 0.0938927716485764 | 0.713702823968326 | -1.2849450811811 | 2.71235072911775 |

### Table 3

H2:

| p_med | lrt_p | two_a | p_med_p025 | p_med_p975 | dataset | clock | n_dropped | n_samples | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.457228366383518 | 2.41897966602237e-08 | 0.0120084580513382 | 0.172340120664341 | 0.712303720087498 | L1 | TAms_chrono | 1 | 146 | supported |
| 1.04559654136252 | 3.24731259384485e-05 | 0.0361347645463202 | 0.646579898044605 | 1.40568236076452 | L1 | TAms_mort | 1 | 146 | supported |

H3:

| age_coef | age_p | n | clock | ci_low | ci_high | direction | p | q_holm | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -0.00637449895696949 | 0.450198895779863 | 58 | TAms_chrono | -0.0231818913377684 | 0.0104328934238295 | - | 0.450198895779863 | 1 | not significant |
| -0.0178521222339094 | 0.673095298812265 | 58 | TAms_mort | -0.102250942349652 | 0.0665466978818336 | - | 0.673095298812265 | 1 | not significant |
| -0.467466723286264 | 0.48209655938365 | 58 | RA_default | -1.7919373596678 | 0.857003913095277 | - | 0.48209655938365 | 1 | not significant |

H4:

| lrt_stat | lrt_p | converged | reml_converged | dataset | clock | p | mean_difference | q_holm | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 7.72981589070241 | 0.0209648522603376 | TRUE | TRUE | L1 | TAms_chrono | 0.0209648522603376 | 0.000307361838231626 | 0.0628945567810128 | not significant |
| 5.68242972879443 | 0.0583547298006714 | TRUE | TRUE | L1 | TAms_mort | 0.0583547298006714 | 0.00530602444996918 | 0.116709459601343 | not significant |
| 0.126705701224637 | 0.938612225353535 | TRUE | TRUE | L1 | RA_default | 0.938612225353535 | 1.30859604147651 | 0.938612225353535 | not significant |
|  |  |  |  | L2 | TAms_chrono |  |  |  | gate fail |
|  |  |  |  | L2 | TAms_mort |  |  |  | gate fail |
|  |  |  |  | L2 | RA_default |  |  |  | gate fail |

H5:

| clock | c | c_numeric | mean_delta | mean_delta_minus_c | ci_low | ci_high | ratio | l3_significant | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TAms_chrono | 0.019125300198260260499475555777287663859727086340733 | 0.0191253001982603 | 0.0120385849794561 | -0.00708671521880421 | -0.0301702220073491 | 0.0159967915697406 | 0.629458615271887 | FALSE | not calibrated |
| TAms_mort | 0.088003575352142004086395806850464311947738779118951 | 0.088003575352142 | 0.0827171714412529 | -0.00528640391088912 | -0.120695167603389 | 0.110122359781611 | 0.939929668882931 | FALSE | not calibrated |
| RA_default | 3.6684993321819421764583897146097149703131816543407 | 3.66849933218194 | 2.28347265426549 | -1.38502667791645 | -3.2227970955387 | 0.452743739705789 | 0.622454155636259 | FALSE | not calibrated |

S4, `in_f3` TRUE:

| dataset | clock | n_subjects | n_samples | two_a | peak_h | p | confirmed | circ_dist | phi_bar | mean_delta | observed_sign | predicted_sign | sign_agrees |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| L1 | TM_turquoise_mortality | 22 | 147 | 0.0390746756069444 | 21.7889728924429 | 0.00164866962432351 | FALSE | 10.2165716930187 | 11.5724011994242 |  |  |  |  |
| L1 | CMP_TAms_chrono_brown | 22 | 147 | 0 | 0 |  | FALSE | 1.32505139804729 | 1.32505139804729 |  |  |  |  |
| L1 | CMP_TAms_chrono_orange | 22 | 147 | 0.00297698520908832 | 2.77928506086391 | 9.888135840689e-09 | TRUE | 1.70595743343867 | 1.07332762742524 |  |  |  |  |
| L1 | CMP_TAms_chrono_red | 22 | 147 | 0 | 0 |  | FALSE | 0.781691026210673 | 0.781691026210673 |  |  |  |  |
| L1 | CMP_TAms_chrono_turquoise | 22 | 147 | 0.000401168491997726 | 20.5994829668058 | 0.355877614252857 | FALSE | 7.31575012524486 | 13.2837328415609 |  |  |  |  |
| L1 | CMP_TAms_chrono_yellow | 22 | 147 | 0.000158800838647645 | 3.03741437413107 | 8.37244837103072e-07 | FALSE | 9.57541254212653 | 12.6128269162576 |  |  |  |  |
| L1 | CMP_TAms_mort_darkred | 22 | 147 | 0.000591230920396307 | 3.04307785612592 | 0.0400291314216308 | FALSE | 11.9235817936314 | 14.9666596497573 |  |  |  |  |
| L1 | CMP_TAms_mort_turquoise | 22 | 147 | 0.00608539245989998 | 0.523167710801014 | 4.63877250086027e-07 | FALSE | 11.519356560875 | 13.003811149926 |  |  |  |  |
| L3 | TM_turquoise_mortality | 58 |  |  |  | 0.0367459632068223 |  |  |  | 0.18561858764589 | + | + | TRUE |
| L3 | CMP_TAms_chrono_brown | 58 |  |  |  | 0.543751743051306 |  |  |  | -0.0003247132009814 | - | - |  |
| L3 | CMP_TAms_chrono_orange | 58 |  |  |  | 0.0383882761729438 |  |  |  | -0.000825802934541766 | - | - | TRUE |
| L3 | CMP_TAms_chrono_red | 58 |  |  |  | 0.967330930103839 |  |  |  | 5.55046917861069e-06 | + | - |  |
| L3 | CMP_TAms_chrono_turquoise | 58 |  |  |  | 0.0686526251496304 |  |  |  | 0.00514397164720193 | + | + | TRUE |
| L3 | CMP_TAms_chrono_yellow | 58 |  |  |  | 0.953453644653339 |  |  |  | -3.30195431830906e-05 | - | + |  |
| L3 | CMP_TAms_mort_darkred | 58 |  |  |  | 0.762045749116664 |  |  |  | 0.00107588556555147 | + | - |  |
| L3 | CMP_TAms_mort_turquoise | 58 |  |  |  | 0.647755299503364 |  |  |  | 0.0040473650870528 | + | + |  |

S4 L1 확인 수:

| n_confirmed | n_f3 |
| --- | --- |
| 1 | 8 |

## 6. 판단이 필요한 사항

없음.

## 7. 지시서와 다르게 한 것

- docx는 python-docx, officer, flextable이 설치되어 있지 않아 OOXML zip으로 썼다. 패키지는 설치하지 않았다.
- Fig 4의 MCPcounter 정점은 시각 0–24 가로축에 집단 이름을 세로축으로 두었다.
- L2 `RA_default` LRT p는 이번 적합에서 `0.0116625636160948`이고 `results/posthoc_stage7c/h1_exploratory.tsv`에서 `0.0116625636160945`이다.
- `scripts/52_part8_feature_counts.py`의 `/usr/bin/time -v` 기록은 없다.

## 8. 산출물 목록, 커밋 해시

- Part A: `67d4c2d7b3570f2b72aba06242706057cc4e635d`, 2026-10-19 11:09:06 +0900, `TOD8 Part A`.
- Part B: `9317da205567e399907d98a2ae1ef21e50f4aab8`, 2026-10-19 12:03:27 +0900, `TOD8 Part B`.
- 이 보고서를 넣는 커밋: 보고서 작성 후 기록.
- 사후 표: `results/stage8/l2_sleep_start_search.tsv`, `l2_bmal1_as_arntl.tsv`, `gene_peak_differences_arntl_alias.tsv`, `gene_peak_circular_summary_arntl_alias.tsv`, `cosinor_effects_boot.tsv`, `effect_summary.tsv`, `l3_d_vs_prediction.tsv`, `l1_f3_d11_model_m.tsv`, `f3_feature_counts.tsv`, `d11_matrix_audit_recount.tsv`.
- 스크립트: `scripts/51_part8_posthoc.R`, `scripts/52_part8_feature_counts.py`, `scripts/53_part8_figures.py`, `scripts/54_part8_tables.py`, `scripts/55_part8_report.py`.
- `results/confirm/`은 변경하지 않았다.
