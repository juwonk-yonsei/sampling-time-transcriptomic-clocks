# TOD 보고서 9 — L2 사후 관문, 편향 시뮬레이션, 완화

## 1. 요약 (수치만, 14줄 이내)

- PER1 peak_wake_corr `0.627072384054571`, in_prereg_window `TRUE`, lrt_p_relclock `1.33395937684947e-05`, lrt_p_wake_corr `1.33395937684951e-05`, lrt_p_difference `3.79470760369927e-19`, prereg_gate_pass `TRUE`. PER3 peak_wake_corr `22.3757375270708`, in_prereg_window `TRUE`, lrt_p_relclock `2.56672604887414e-06`, lrt_p_wake_corr `2.56672604887414e-06`, lrt_p_difference `0`, prereg_gate_pass `TRUE`.
- 유전자 10개 circ_mean_h `1.69765168037113`, circ_sd_h `5.49841394966151`.
- 설계 1, s = 3, lab, n = 100: TAms_chrono b `0.00161130715539329`, b/total_sd `0.110609413254012`, fpr `0.287597392213601`; TAms_mort b `0.00693560844477734`, b/total_sd `0.112080414032841`, fpr `0.328504596537995`; RA_default b `1.09380254924844`, b/total_sd `0.0936705280285436`, fpr `0.169984918325999`; RA_Peters b `-2.40846522220964`, b/total_sd `-0.14465637064403`, fpr `0.497739451786296`; RA_GTEx b `4.28408351756162`, b/total_sd `0.141640575768605`, fpr `0.199274647655761`.
- 설계 1, s = 3, free-living, n = 100: TAms_chrono b `0.00947216425426115`, b/total_sd `0.112256607571981`, fpr `0.187568364077656`; TAms_mort b `0.0345051688024083`, b/total_sd `0.0890048454721773`, fpr `0.121801022888275`; RA_default b `0.582482849833315`, b/total_sd `0.0740603616754689`, fpr `0.130935658522455`; RA_Peters b `-0.534648480913311`, b/total_sd `-0.0310113426725348`, fpr `0.0780800617596416`; RA_GTEx b `0.227642526504202`, b/total_sd `0.0279017210378773`, fpr `0.0601364956539231`.
- n = 100에서 위양성률이 threshold `0.1`을 넘는 가장 작은 s: lab TAms_chrono `1` (exceeds); lab TAms_mort `1` (exceeds); lab RA_default `2` (exceeds); lab RA_Peters `2` (exceeds); lab RA_GTEx `2` (exceeds); free-living TAms_chrono `2` (exceeds); free-living TAms_mort `3` (exceeds); free-living RA_default `2` (exceeds); free-living RA_Peters `4` (exceeds); free-living RA_GTEx `12` (exceeds); free-living upper TAms_chrono `1` (exceeds); free-living upper TAms_mort `1` (exceeds); free-living upper RA_default `1` (exceeds); free-living upper RA_Peters `3` (exceeds); free-living upper RA_GTEx `2` (exceeds).
- 설계 2, s = 4, 군당 n = 100: lab TAms_chrono `0.100718309824148`; lab TAms_mort `0.12335153033894`; lab RA_default `0.102724735193889`; lab RA_Peters `0.261543396224426`; lab RA_GTEx `0.262672673073346`; free-living TAms_chrono `0.14593686612444`; free-living TAms_mort `0.101363966459141`; free-living RA_default `0.0783884169876189`; free-living RA_Peters `0.0641311541987358`; free-living RA_GTEx `0.0544567318102879`; free-living upper TAms_chrono `0.743969698644312`; free-living upper TAms_mort `0.353393392127528`; free-living upper RA_default `0.144830603170745`; free-living upper RA_Peters `0.107089928918182`; free-living upper RA_GTEx `0.116258479110633`.
- 설계 3, W = 0: TAms_chrono power `0.56`, bias `4.48391679889163e-05`; TAms_mort power `0.546`, bias `-5.84019623913225e-05`; RA_default power `0.563`, bias `0.0386793122491103`; RA_Peters power `0.5485`, bias `0.0175604272730201`; RA_GTEx power `0.569`, bias `0.253167915116043`.
- 설계 3, W = 2: TAms_chrono power `0.536`, bias `1.16925958745752e-05`; TAms_mort power `0.538`, bias `-7.32373131568809e-05`; RA_default power `0.5455`, bias `-0.0194796656702789`; RA_Peters power `0.546`, bias `-0.0469179667765909`; RA_GTEx power `0.552`, bias `0.151759439941641`.
- 설계 3, W = 4: TAms_chrono power `0.5205`, bias `-4.73959810330495e-05`; TAms_mort power `0.5325`, bias `-0.000219091549425118`; RA_default power `0.536`, bias `-0.0344874395344781`; RA_Peters power `0.5345`, bias `-0.0331392348678246`; RA_GTEx power `0.5445`, bias `0.123551863870455`.
- 설계 3, W = 8: TAms_chrono power `0.524`, bias `-1.72512069822784e-05`; TAms_mort power `0.5445`, bias `3.47417744971188e-05`; RA_default power `0.553`, bias `0.0311306280871482`; RA_Peters power `0.517`, bias `-0.0199246300612885`; RA_GTEx power `0.539`, bias `-0.170163487908811`.
- L1 two_a_ratio: TAms_chrono `0.49552049306934`; TAms_mort `0.335858870871784`; RA_default `0.185701803638067`; RA_Peters `0.754695041061693`; RA_GTEx `1.00312677926618`.
- L3 mean Δ: TAms_chrono before `0.0120385849794561`, composition `0.0166461557079569`, time `-0.00708671521880421`; TAms_mort before `0.0827171714412529`, composition `0.0828733287438787`, time `-0.00528640391088912`; RA_default before `2.28347265426549`, composition `1.81652100332312`, time `-1.38502667791645`; RA_Peters before `2.86233409356227`, composition `2.02172266582033`, time ``; RA_GTEx before `0.713702823968326`, composition `0.159401851182051`, time ``.
- L3 나이 Pearson r: TAms_chrono morning `0.271730437294783` / `0.25049954468322`, evening `0.344473699826594` / `0.379420830334003`, person mean `0.357659887896734`; TAms_mort morning `0.459242940208403` / `0.454368246593447`, evening `0.462942410027015` / `0.504173016355484`, person mean `0.542417587627288`; RA_default morning `0.446171503104569` / `0.484858641020211`, evening `0.517349079713589` / `0.499271914644175`, person mean `0.536867692830739`; RA_Peters morning `0.716354841847154` / `0.731704823753873`, evening `0.709903447516876` / `0.70297198458081`, person mean `0.740682571652942`; RA_GTEx morning `0.295924783031412` / `0.306451158604095`, evening `0.482587756628502` / `0.495235502983008`, person mean `0.438835111423925`.

## 2. 환경과 자원

- 분석: `scripts/56_part9_analysis.R`. R version 4.5.3 (2026-03-11), conda 환경 `r_env`.
- 명령: `/usr/bin/time -v Rscript scripts/56_part9_analysis.R`.
- 사용자 시간 `108.79` s, 시스템 시간 `2.18` s, CPU `102%`, 경과 `1:48.72`, 최대 RSS `2106748` kB, 종료 상태 0.
- 결정 10의 중지 기준은 최대 RSS 40 GB다. 관측 최대 RSS는 `2106748` kB다.
- 몬테카를로 시드 `20261019`. 설계 1 몬테카를로 시작에 `set.seed(20261019)`, 설계 3 시작에 다시 `set.seed(20261019)`.
- 로그 `logs/stage9_analysis.log`의 Part A 구간에 `boundary (singular) fit` 5줄, `Warning: stack imbalance` 4줄이 있다. 설계 1·2 출력 뒤에 `There were 25 warnings`가 있다.
- 그림: `/home/jkang/anaconda3/bin/python3`, `MPLBACKEND=Agg`, `scripts/57_part9_figures.py`. `savefig` dpi 300. PNG info dpi는 Fig 6과 Fig 7 모두 `(299.9994, 299.9994)`.
- 시계 표기: `TAms_chrono` tAge chronological, `TAms_mort` tAge mortality, `RA_default` RNAAgeCalc blood, `RA_Peters`, `RA_GTEx`.
- 표와 그림 제목은 `post hoc`다. `results/confirm/`은 수정하지 않았다.

## 3. L2 사후 관문 (표)

`wake_corr = (relclocktime - 8) mod 24`. period 1. LRT p는 `relclocktime` 축과 `wake_corr` 축에서 각각 계산했다.

사후 관문 (`results/stage9/l2_wake_corr_gate.tsv`):

| gene | n_samples | n_subjects | peak_relclock | peak_wake_corr | in_prereg_window | window_rule | lrt_p_relclock | lrt_p_wake_corr | lrt_stat_relclock | lrt_stat_wake_corr | lrt_p_difference | lrt_p_below_0.01 | prereg_gate_pass | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PER1 | 70 | 11 | 8.62707238405378 | 0.627072384054571 | TRUE | peak >= 18 or peak <= 2, the preregistered 18-26 h window modulo 24 | 1.33395937684947e-05 | 1.33395937684951e-05 | 22.4495479401546 | 22.4495479401545 | 3.79470760369927e-19 | TRUE | TRUE | post hoc |
| PER3 | 70 | 11 | 6.37573752707371 | 22.3757375270708 | TRUE | peak >= 18 or peak <= 2, the preregistered 18-26 h window modulo 24 | 2.56672604887414e-06 | 2.56672604887414e-06 | 25.7457587642308 | 25.7457587642308 | 0 | TRUE | TRUE | post hoc |

유전자 10개 (`results/stage9/l2_wake_corr_genes.tsv`):

| gene | symbol_used | n_samples | n_subjects | peak_wake_corr | lrt_p_wake_corr | status | peak_gse39445_hours_awake | signed_diff_h | title | circ_dist_0_12 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PER1 | PER1 | 70 | 11 | 0.627072384054571 | 1.33395937684951e-05 | fit | 22.7776132062867 | 1.84945917776787 | post hoc | 1.84945917776787 |
| PER2 | PER2 | 70 | 11 | 22.7975702798119 | 5.47315115180425e-05 | fit | 20.7330572724929 | 2.064513007319 | post hoc | 2.064513007319 |
| PER3 | PER3 | 70 | 11 | 22.3757375270708 | 2.56672604887414e-06 | fit | 21.9104219255107 | 0.4653156015601 | post hoc | 0.4653156015601 |
| NR1D1 | NR1D1 | 70 | 11 | 18.3535171570868 | 5.56394793297238e-07 | fit | 17.7523074513839 | 0.601209705702875 | post hoc | 0.601209705702875 |
| NR1D2 | NR1D2 | 70 | 11 | 20.2703663646416 | 7.14824571359563e-05 | fit | 19.5452527557332 | 0.725113608908448 | post hoc | 0.725113608908448 |
| ARNTL | BMAL1 | 70 | 11 | 13.3692664519889 | 0.00925089859203762 | fit | 7.37049546630091 | 5.99877098568803 | post hoc | 5.99877098568803 |
| CRY1 | CRY1 | 70 | 11 | 8.96706876018854 | 0.399939321666829 | fit | 20.3669852307156 | -11.3999164705271 | post hoc | 11.3999164705271 |
| CRY2 | CRY2 | 70 | 11 | 3.67074663411381 | 0.181311703326108 | fit | 18.5308484498846 | 9.13989818422921 | post hoc | 9.13989818422921 |
| DBP | DBP | 70 | 11 | 19.7658069493448 | 0.00200569749239259 | fit | 4.06024592076617 | -8.29443897142142 | post hoc | 8.29443897142142 |
| TEF | TEF | 70 | 11 | 14.8844839167013 | 0.252474488903234 | fit | 17.4419627275353 | -2.55747881083399 | post hoc | 2.55747881083399 |

원형 요약 (`results/stage9/l2_wake_corr_gene_circular.tsv`):

| n_genes | circ_mean_h | circ_sd_h | signed_diff_definition | circ_sd_definition | arntl_source | title |
| --- | --- | --- | --- | --- | --- | --- |
| 10 | 1.69765168037113 | 5.49841394966151 | (peak_wake_corr - peak_gse39445_hours_awake) mod 24, then subtract 24 when the result is greater than 12; range (-12, 12] | sqrt(-2*log(R)) * 24 / (2*pi), R = mean resultant length of the signed differences | L2 symbol BMAL1; GSE39445 symbol ARNTL, hours_awake axis, stage 8 table | post hoc |

시계 정점 (`results/stage9/l2_wake_corr_clocks.tsv`):

| clock | n_subjects | n_samples | peak_relclock | peak_wake_corr | lrt_p_relclock | lrt_p_wake_corr | lrt_p_difference | role | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TAms_chrono | 11 | 70 | 10.800630536763 | 2.80063053676038 | 0.0177344777316928 | 0.0177344777316921 | -7.56339435525888e-16 | reference | post hoc |
| TAms_mort | 11 | 70 | 10.9821548000952 | 2.98215480009343 | 0.0174608076774516 | 0.0174608076774515 | -6.24500451351651e-17 | reference | post hoc |
| RA_default | 11 | 70 | 8.29862480973111 | 0.298624823104874 | 0.0116625636160948 | 0.0116625636160945 | -3.3133218391157e-16 | reference | post hoc |
| RA_Peters | 11 | 70 | 9.718280100379 | 1.71828010038026 | 0.000688235638502205 | 0.000688235638502225 | 1.95156391047391e-17 | reference | post hoc |
| RA_GTEx | 11 | 70 | 7.24447680202368 | 23.2444768020377 | 0.00107043094705134 | 0.00107043094705134 | 0 | reference | post hoc |

## 4. 편향 시뮬레이션

### 4.1 매개변수 표

`results/stage9/sim_parameters.tsv`. free-living의 코사인 모양과 φ̄는 `assumption` 열에 있다.

| scenario | clock | phi_bar | A | cosine_denominator | beta_cos | beta_sin | sigma | sigma_u | total_sd | sd_diff | mean_delta | ci_high | morning_sd | assumption | unit | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| lab | TAms_chrono |  |  |  | -0.00570219212781147 | 0.0017218057888685 | 0.00806869915269403 | 0.0121288640391942 | 0.0145675409351731 | 0.0114108637724482 |  |  |  | f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24) from L1 Day 1 model M | score | post hoc |
| free-living | TAms_chrono | 12.616768998085 | 0.0169908846438399 | 0.708531970630542 |  |  |  |  | 0.0843795698011575 | 0.0877911523385128 | 0.0120385849794561 | 0.0351220917680009 | 0.0843795698011575 | f(t) = A*cos(2*pi*(t-phi_bar)/24). phi_bar is F1. A makes f(8)-f(20) equal the L3 mean delta, or its 95 percent upper limit in the upper scenario. The cosine shape and phi_bar are assumptions. | score | post hoc |
| free-living upper | TAms_chrono | 12.616768998085 | 0.0495702286189638 | 0.708531970630542 |  |  |  |  | 0.0843795698011575 | 0.0877911523385128 | 0.0120385849794561 | 0.0351220917680009 | 0.0843795698011575 | f(t) = A*cos(2*pi*(t-phi_bar)/24). phi_bar is F1. A makes f(8)-f(20) equal the L3 mean delta, or its 95 percent upper limit in the upper scenario. The cosine shape and phi_bar are assumptions. | score | post hoc |
| lab | TAms_mort |  |  |  | -0.0180759345568866 | 0.00244799198548219 | 0.0320527438270783 | 0.0529323718072866 | 0.061880646182705 | 0.0453294250315246 |  |  |  | f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24) from L1 Day 1 model M | score | post hoc |
| free-living | TAms_mort | 11.8599407935094 | 0.0778277650638591 | 1.06282342006586 |  |  |  |  | 0.387677419351225 | 0.438922406692367 | 0.0827171714412529 | 0.198125935133753 | 0.387677419351225 | f(t) = A*cos(2*pi*(t-phi_bar)/24). phi_bar is F1. A makes f(8)-f(20) equal the L3 mean delta, or its 95 percent upper limit in the upper scenario. The cosine shape and phi_bar are assumptions. | score | post hoc |
| free-living upper | TAms_mort | 11.8599407935094 | 0.186414724584707 | 1.06282342006586 |  |  |  |  | 0.387677419351225 | 0.438922406692367 | 0.0827171714412529 | 0.198125935133753 | 0.387677419351225 | f(t) = A*cos(2*pi*(t-phi_bar)/24). phi_bar is F1. A makes f(8)-f(20) equal the L3 mean delta, or its 95 percent upper limit in the upper scenario. The cosine shape and phi_bar are assumptions. | score | post hoc |
| lab | RA_default |  |  |  | -2.72907946012988 | 0.292727588701374 | 7.66228637689109 | 8.81161942253776 | 11.6771259036261 | 10.836109312986 |  |  |  | f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24) from L1 Day 1 model M | year | post hoc |
| free-living | RA_default | 11.2423863477089 | 1.727696501972 | 1.32168621726045 |  |  |  |  | 7.86497441621665 | 6.9894052136286 | 2.28347265426549 | 4.12124307188773 | 7.86497441621665 | f(t) = A*cos(2*pi*(t-phi_bar)/24). phi_bar is F1. A makes f(8)-f(20) equal the L3 mean delta, or its 95 percent upper limit in the upper scenario. The cosine shape and phi_bar are assumptions. | year | post hoc |
| free-living upper | RA_default | 11.2423863477089 | 3.11817057488131 | 1.32168621726045 |  |  |  |  | 7.86497441621665 | 6.9894052136286 | 2.28347265426549 | 4.12124307188773 | 7.86497441621665 | f(t) = A*cos(2*pi*(t-phi_bar)/24). phi_bar is F1. A makes f(8)-f(20) equal the L3 mean delta, or its 95 percent upper limit in the upper scenario. The cosine shape and phi_bar are assumptions. | year | post hoc |
| lab | RA_Peters |  |  |  | 1.89871848646548 | 2.50953110211271 | 8.63023837847397 | 14.2382197881013 | 16.649562072426 | 12.2050001613507 |  |  |  | f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24) from L1 Day 1 model M | year | post hoc |
| free-living | RA_Peters | 7.56604107260516 | 1.44045321375948 | 1.98710660382491 |  |  |  |  | 17.2404170486569 | 10.7770588259086 | 2.86233409356227 | 5.69601726028016 | 17.2404170486569 | f(t) = A*cos(2*pi*(t-phi_bar)/24). phi_bar is F1. A makes f(8)-f(20) equal the L3 mean delta, or its 95 percent upper limit in the upper scenario. The cosine shape and phi_bar are assumptions. | year | post hoc |
| free-living upper | RA_Peters | 7.56604107260516 | 2.86648801293101 | 1.98710660382491 |  |  |  |  | 17.2404170486569 | 10.7770588259086 | 2.86233409356227 | 5.69601726028016 | 17.2404170486569 | f(t) = A*cos(2*pi*(t-phi_bar)/24). phi_bar is F1. A makes f(8)-f(20) equal the L3 mean delta, or its 95 percent upper limit in the upper scenario. The cosine shape and phi_bar are assumptions. | year | post hoc |
| lab | RA_GTEx |  |  |  | -2.56734872788383 | -5.08540366219671 | 26.9797237185506 | 13.6720409220887 | 30.2461600026282 | 38.1550911918533 |  |  |  | f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24) from L1 Day 1 model M | year | post hoc |
| free-living | RA_GTEx | 11.518486825228 | 0.589920304949881 | 1.20982922265909 |  |  |  |  | 8.15872706186014 | 7.60125419068026 | 0.713702823968326 | 2.71235072911775 | 8.15872706186014 | f(t) = A*cos(2*pi*(t-phi_bar)/24). phi_bar is F1. A makes f(8)-f(20) equal the L3 mean delta, or its 95 percent upper limit in the upper scenario. The cosine shape and phi_bar are assumptions. | year | post hoc |
| free-living upper | RA_GTEx | 11.518486825228 | 2.24192859481131 | 1.20982922265909 |  |  |  |  | 8.15872706186014 | 7.60125419068026 | 0.713702823968326 | 2.71235072911775 | 8.15872706186014 | f(t) = A*cos(2*pi*(t-phi_bar)/24). phi_bar is F1. A makes f(8)-f(20) equal the L3 mean delta, or its 95 percent upper limit in the upper scenario. The cosine shape and phi_bar are assumptions. | year | post hoc |

### 4.2 설계 1 (b(s) 표, 위양성률 표, 최소 s 표, 몬테카를로 비교)

대응 t, 양측 α = 0.05. ncp = b / SD_diff × sqrt(n), df = n − 1.

b(s) (`results/stage9/sim_design1_bias.tsv`):

| scenario | clock | s | b | b_over_total_sd | total_sd | unit | title |
| --- | --- | --- | --- | --- | --- | --- | --- |
| lab | TAms_chrono | 0 | 0 | 0 | 0.0145675409351731 | score | post hoc |
| lab | TAms_chrono | 1 | 0.000907335652950247 | 0.0622847505277641 | 0.0145675409351731 | score | post hoc |
| lab | TAms_chrono | 2 | 0.00145692251692967 | 0.100011561554082 | 0.0145675409351731 | score | post hoc |
| lab | TAms_chrono | 3 | 0.00161130715539329 | 0.110609413254012 | 0.0145675409351731 | score | post hoc |
| lab | TAms_chrono | 4 | 0.00135996851036251 | 0.0933560795479825 | 0.0145675409351731 | score | post hoc |
| lab | TAms_chrono | 6 | -0.00026488327193883 | -0.0181831149895226 | 0.0145675409351731 | score | post hoc |
| lab | TAms_chrono | 8 | -0.00298225510708644 | -0.204719184957692 | 0.0145675409351731 | score | post hoc |
| lab | TAms_chrono | 12 | -0.00868444723489791 | -0.596150529011348 | 0.0145675409351731 | score | post hoc |
| lab | TAms_mort | 0 | 0 | 0 | 0.061880646182705 | score | post hoc |
| lab | TAms_mort | 1 | 0.00335461710855195 | 0.0542110872379599 | 0.061880646182705 | score | post hoc |
| lab | TAms_mort | 2 | 0.00572022399001831 | 0.0924396292360801 | 0.061880646182705 | score | post hoc |
| lab | TAms_mort | 3 | 0.00693560844477734 | 0.112080414032841 | 0.061880646182705 | score | post hoc |
| lab | TAms_mort | 4 | 0.006917944030755 | 0.111794954602276 | 0.061880646182705 | score | post hoc |
| lab | TAms_mort | 6 | 0.00327223200453611 | 0.0528797322974735 | 0.061880646182705 | score | post hoc |
| lab | TAms_mort | 8 | -0.00424004649537656 | -0.0685197514398549 | 0.061880646182705 | score | post hoc |
| lab | TAms_mort | 12 | -0.0223159810522631 | -0.360629412084262 | 0.061880646182705 | score | post hoc |
| lab | RA_default | 0 | 0 | 0 | 11.6771259036261 | year | post hoc |
| lab | RA_default | 1 | 0.518690997396996 | 0.0444194060831294 | 11.6771259036261 | year | post hoc |
| lab | RA_default | 2 | 0.891766677500593 | 0.0763686788050879 | 11.6771259036261 | year | post hoc |
| lab | RA_default | 3 | 1.09380254924844 | 0.0936705280285436 | 11.6771259036261 | year | post hoc |
| lab | RA_default | 4 | 1.11103020186099 | 0.0951458613215758 | 11.6771259036261 | year | post hoc |
| lab | RA_default | 6 | 0.599039088799219 | 0.0513002166580391 | 11.6771259036261 | year | post hoc |
| lab | RA_default | 8 | -0.507019056407903 | -0.0434198501063056 | 11.6771259036261 | year | post hoc |
| lab | RA_default | 12 | -3.23609851653779 | -0.277131422855763 | 11.6771259036261 | year | post hoc |
| lab | RA_Peters | 0 | 0 | 0 | 16.649562072426 | year | post hoc |
| lab | RA_Peters | 1 | -0.79204870022558 | -0.0475717437359704 | 16.649562072426 | year | post hoc |
| lab | RA_Peters | 2 | -1.61353133564192 | -0.0969113378852262 | 16.649562072426 | year | post hoc |
| lab | RA_Peters | 3 | -2.40846522220964 | -0.14465637064403 | 16.649562072426 | year | post hoc |
| lab | RA_Peters | 4 | -3.12267692924951 | -0.187553096932267 | 16.649562072426 | year | post hoc |
| lab | RA_Peters | 6 | -4.12306243775463 | -0.247637891003933 | 16.649562072426 | year | post hoc |
| lab | RA_Peters | 8 | -4.34663537203354 | -0.261066048051329 | 16.649562072426 | year | post hoc |
| lab | RA_Peters | 12 | -2.44791688556806 | -0.147025902238123 | 16.649562072426 | year | post hoc |
| lab | RA_GTEx | 0 | 0 | 0 | 30.2461600026282 | year | post hoc |
| lab | RA_GTEx | 1 | 1.33988067656597 | 0.0442991995165516 | 30.2461600026282 | year | post hoc |
| lab | RA_GTEx | 2 | 2.80110178364156 | 0.0926101621957352 | 30.2461600026282 | year | post hoc |
| lab | RA_GTEx | 3 | 4.28408351756162 | 0.141640575768605 | 30.2461600026282 | year | post hoc |
| lab | RA_GTEx | 4 | 5.68776312390269 | 0.18804909857676 | 30.2461600026282 | year | post hoc |
| lab | RA_GTEx | 6 | 7.88650544583827 | 0.260744023213293 | 30.2461600026282 | year | post hoc |
| lab | RA_GTEx | 8 | 8.80817751992155 | 0.291216389755135 | 30.2461600026282 | year | post hoc |
| lab | RA_GTEx | 12 | 6.24082879203771 | 0.206334582356749 | 30.2461600026282 | year | post hoc |
| free-living | TAms_chrono | 0 | 0 | 0 | 0.0843795698011575 | score | post hoc |
| free-living | TAms_chrono | 1 | 0.00390725675763172 | 0.0463057202927115 | 0.0843795698011575 | score | post hoc |
| free-living | TAms_chrono | 2 | 0.00713803558845414 | 0.084594358625851 | 0.0843795698011575 | score | post hoc |
| free-living | TAms_chrono | 3 | 0.00947216425426115 | 0.112256607571981 | 0.0843795698011575 | score | post hoc |
| free-living | TAms_chrono | 4 | 0.010750575743808 | 0.127407330579452 | 0.0843795698011575 | score | post hoc |
| free-living | TAms_chrono | 6 | 0.00986964324880958 | 0.11696721460026 | 0.0843795698011575 | score | post hoc |
| free-living | TAms_chrono | 8 | 0.00473128325407997 | 0.0560714313337856 | 0.0843795698011575 | score | post hoc |
| free-living | TAms_chrono | 12 | -0.0120385849794561 | -0.142671798491333 | 0.0843795698011575 | score | post hoc |
| free-living | TAms_mort | 0 | 0 | 0 | 0.387677419351225 | score | post hoc |
| free-living | TAms_mort | 1 | 0.015654412084611 | 0.0403799945604479 | 0.387677419351225 | score | post hoc |
| free-living | TAms_mort | 2 | 0.0274234825872286 | 0.07073788984956 | 0.387677419351225 | score | post hoc |
| free-living | TAms_mort | 3 | 0.0345051688024083 | 0.0890048454721773 | 0.387677419351225 | score | post hoc |
| free-living | TAms_mort | 4 | 0.0364168655176252 | 0.09393599859019 | 0.387677419351225 | score | post hoc |
| free-living | TAms_mort | 6 | 0.024570379097766 | 0.0633784117189089 | 0.387677419351225 | score | post hoc |
| free-living | TAms_mort | 8 | -0.00494172020300127 | -0.0127469900394797 | 0.387677419351225 | score | post hoc |
| free-living | TAms_mort | 12 | -0.0827171714412529 | -0.213365977259339 | 0.387677419351225 | score | post hoc |
| free-living | RA_default | 0 | 0 | 0 | 7.86497441621665 | year | post hoc |
| free-living | RA_default | 1 | 0.296701192571716 | 0.0377243684302333 | 7.86497441621665 | year | post hoc |
| free-living | RA_default | 2 | 0.495375245305969 | 0.0629849786014006 | 7.86497441621665 | year | post hoc |
| free-living | RA_default | 3 | 0.582482849833315 | 0.0740603616754689 | 7.86497441621665 | year | post hoc |
| free-living | RA_default | 4 | 0.552087766857339 | 0.0701957485988764 | 7.86497441621665 | year | post hoc |
| free-living | RA_default | 6 | 0.154941490303688 | 0.0197001899947973 | 7.86497441621665 | year | post hoc |
| free-living | RA_default | 8 | -0.589648560275406 | -0.0749714530615154 | 7.86497441621665 | year | post hoc |
| free-living | RA_default | 12 | -2.28347265426549 | -0.290334403320784 | 7.86497441621665 | year | post hoc |
| free-living | RA_Peters | 0 | 0 | 0 | 17.2404170486569 | year | post hoc |
| free-living | RA_Peters | 1 | -0.0910305529194535 | -0.00528006675607334 | 17.2404170486569 | year | post hoc |
| free-living | RA_Peters | 2 | -0.273389193215311 | -0.0158574582299103 | 17.2404170486569 | year | post hoc |
| free-living | RA_Peters | 3 | -0.534648480913311 | -0.0310113426725348 | 17.2404170486569 | year | post hoc |
| free-living | RA_Peters | 4 | -0.857004027308223 | -0.0497090079021601 | 17.2404170486569 | year | post hoc |
| free-living | RA_Peters | 6 | -1.59446537879272 | -0.092484153619528 | 17.2404170486569 | year | post hoc |
| free-living | RA_Peters | 8 | -2.28817107408936 | -0.132721329630922 | 17.2404170486569 | year | post hoc |
| free-living | RA_Peters | 12 | -2.86233409356227 | -0.166024643457523 | 17.2404170486569 | year | post hoc |
| free-living | RA_GTEx | 0 | 0 | 0 | 8.15872706186014 | year | post hoc |
| free-living | RA_GTEx | 1 | 0.109420271150594 | 0.0134114391033013 | 8.15872706186014 | year | post hoc |
| free-living | RA_GTEx | 2 | 0.187064897645943 | 0.0229281965467899 | 8.15872706186014 | year | post hoc |
| free-living | RA_GTEx | 3 | 0.227642526504202 | 0.0279017210378773 | 8.15872706186014 | year | post hoc |
| free-living | RA_GTEx | 4 | 0.228387859376383 | 0.0279930751016338 | 8.15872706186014 | year | post hoc |
| free-living | RA_GTEx | 6 | 0.112896430966786 | 0.0138375055459995 | 8.15872706186014 | year | post hoc |
| free-living | RA_GTEx | 8 | -0.12846355260778 | -0.0157455386892782 | 8.15872706186014 | year | post hoc |
| free-living | RA_GTEx | 12 | -0.713702823968326 | -0.087477227581824 | 8.15872706186014 | year | post hoc |
| free-living upper | TAms_chrono | 0 | 0 | 0 | 0.0843795698011575 | score | post hoc |
| free-living upper | TAms_chrono | 1 | 0.0113992658304002 | 0.135095093009633 | 0.0843795698011575 | score | post hoc |
| free-living upper | TAms_chrono | 2 | 0.020824934276642 | 0.246800668997442 | 0.0843795698011575 | score | post hoc |
| free-living upper | TAms_chrono | 3 | 0.0276346616107675 | 0.327504177562048 | 0.0843795698011575 | score | post hoc |
| free-living upper | TAms_chrono | 4 | 0.0313643761685627 | 0.37170580796363 | 0.0843795698011575 | score | post hoc |
| free-living upper | TAms_chrono | 6 | 0.0287942907321473 | 0.341247185782077 | 0.0843795698011575 | score | post hoc |
| free-living upper | TAms_chrono | 8 | 0.0138033302845623 | 0.163586165668895 | 0.0843795698011575 | score | post hoc |
| free-living upper | TAms_chrono | 12 | -0.0351220917680009 | -0.416239284589468 | 0.0843795698011575 | score | post hoc |
| free-living upper | TAms_mort | 0 | 0 | 0 | 0.387677419351225 | score | post hoc |
| free-living upper | TAms_mort | 1 | 0.037495782063044 | 0.0967190251260775 | 0.387677419351225 | score | post hoc |
| free-living upper | TAms_mort | 2 | 0.065685311012813 | 0.169432904095206 | 0.387677419351225 | score | post hoc |
| free-living upper | TAms_mort | 3 | 0.0826475170367793 | 0.213186306221005 | 0.387677419351225 | score | post hoc |
| free-living upper | TAms_mort | 4 | 0.0872264538257805 | 0.224997509454518 | 0.387677419351225 | score | post hoc |
| free-living upper | TAms_mort | 6 | 0.0588514966181245 | 0.151805324944156 | 0.387677419351225 | score | post hoc |
| free-living upper | TAms_mort | 8 | -0.011836513741096 | -0.0305318627040604 | 0.387677419351225 | score | post hoc |
| free-living upper | TAms_mort | 12 | -0.198125935133753 | -0.511058744317157 | 0.387677419351225 | score | post hoc |
| free-living upper | RA_default | 0 | 0 | 0 | 7.86497441621665 | year | post hoc |
| free-living upper | RA_default | 1 | 0.535490421583496 | 0.0680854626150313 | 7.86497441621665 | year | post hoc |
| free-living upper | RA_default | 2 | 0.894060103539363 | 0.11367616170447 | 7.86497441621665 | year | post hoc |
| free-living upper | RA_default | 3 | 1.05127311460673 | 0.133665166467563 | 7.86497441621665 | year | post hoc |
| free-living upper | RA_default | 4 | 0.99641564788813 | 0.126690259263099 | 7.86497441621665 | year | post hoc |
| free-living upper | RA_default | 6 | 0.2796405475972 | 0.0355551757448332 | 7.86497441621665 | year | post hoc |
| free-living upper | RA_default | 8 | -1.06420588805574 | -0.135309516819466 | 7.86497441621665 | year | post hoc |
| free-living upper | RA_default | 12 | -4.12124307188773 | -0.52399955216513 | 7.86497441621665 | year | post hoc |
| free-living upper | RA_Peters | 0 | 0 | 0 | 17.2404170486569 | year | post hoc |
| free-living upper | RA_Peters | 1 | -0.181149922997545 | -0.0105072819576402 | 17.2404170486569 | year | post hoc |
| free-living upper | RA_Peters | 2 | -0.544041859694461 | -0.0315561890503596 | 17.2404170486569 | year | post hoc |
| free-living upper | RA_Peters | 3 | -1.06394532431213 | -0.0617122730447532 | 17.2404170486569 | year | post hoc |
| free-living upper | RA_Peters | 4 | -1.70542975491797 | -0.098920446651888 | 17.2404170486569 | year | post hoc |
| free-living upper | RA_Peters | 6 | -3.1729707370461 | -0.184042574381533 | 17.2404170486569 | year | post hoc |
| free-living upper | RA_Peters | 8 | -4.55343838505805 | -0.264114166856119 | 17.2404170486569 | year | post hoc |
| free-living upper | RA_Peters | 12 | -5.69601726028016 | -0.330387440408461 | 17.2404170486569 | year | post hoc |
| free-living upper | RA_GTEx | 0 | 0 | 0 | 8.15872706186014 | year | post hoc |
| free-living upper | RA_GTEx | 1 | 0.415839957848824 | 0.0509687301334989 | 8.15872706186014 | year | post hoc |
| free-living upper | RA_GTEx | 2 | 0.710919999869342 | 0.0871361420083167 | 8.15872706186014 | year | post hoc |
| free-living upper | RA_GTEx | 3 | 0.865130908840683 | 0.106037486274659 | 8.15872706186014 | year | post hoc |
| free-living upper | RA_GTEx | 4 | 0.867963466162026 | 0.106384667556722 | 8.15872706186014 | year | post hoc |
| free-living upper | RA_GTEx | 6 | 0.429050728908343 | 0.0525879497690319 | 8.15872706186014 | year | post hoc |
| free-living upper | RA_GTEx | 8 | -0.488211898396849 | -0.059839224268098 | 8.15872706186014 | year | post hoc |
| free-living upper | RA_GTEx | 12 | -2.71235072911775 | -0.332447783649641 | 8.15872706186014 | year | post hoc |

위양성률 (`results/stage9/sim_design1_fpr.tsv`):

| scenario | clock | s | n | df | ncp | sd_diff | fpr | alpha | test | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| lab | TAms_chrono | 0 | 20 | 19 | 0 | 0.0114108637724482 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 0 | 50 | 49 | 0 | 0.0114108637724482 | 0.05 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 0 | 100 | 99 | 0 | 0.0114108637724482 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 0 | 200 | 199 | 0 | 0.0114108637724482 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 1 | 20 | 19 | 0.355602211868422 | 0.0114108637724482 | 0.0631864142601811 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 1 | 50 | 49 | 0.562256465248986 | 0.0114108637724482 | 0.0854837073950546 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 1 | 100 | 99 | 0.795150718687073 | 0.0114108637724482 | 0.123500360802521 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 1 | 200 | 199 | 1.12451293049797 | 0.0114108637724482 | 0.201247817096868 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 2 | 20 | 19 | 0.570995824815787 | 0.0114108637724482 | 0.0843605467596723 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 2 | 50 | 49 | 0.90282367043219 | 0.0114108637724482 | 0.143463818982711 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 2 | 100 | 99 | 1.27678547915666 | 0.0114108637724482 | 0.243987431390136 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 2 | 200 | 199 | 1.80564734086438 | 0.0114108637724482 | 0.435332563863717 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 3 | 20 | 19 | 0.631502120074508 | 0.0114108637724482 | 0.0921749689025298 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 3 | 50 | 49 | 0.998492523330293 | 0.0114108637724482 | 0.16495475888083 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 3 | 100 | 99 | 1.41208166842183 | 0.0114108637724482 | 0.287597392213601 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 3 | 200 | 199 | 1.99698504666059 | 0.0114108637724482 | 0.510960785116913 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 4 | 20 | 19 | 0.532997693614085 | 0.0114108637724482 | 0.0798768323859802 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 4 | 50 | 49 | 0.842743349718545 | 0.0114108637724482 | 0.131136922231633 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 4 | 100 | 99 | 1.1918190747717 | 0.0114108637724482 | 0.218631450792717 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 4 | 200 | 199 | 1.68548669943709 | 0.0114108637724482 | 0.388875172001528 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 6 | 20 | 19 | -0.103812825035717 | 0.0114108637724482 | 0.0511161347151495 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 6 | 50 | 49 | -0.164142488724709 | 0.0114108637724482 | 0.0529727825190948 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 6 | 100 | 99 | -0.232132533716156 | 0.0114108637724482 | 0.0560759033662722 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 6 | 200 | 199 | -0.328284977449418 | 0.0114108637724482 | 0.0623147376193568 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 8 | 20 | 19 | -1.16880286692976 | 0.0114108637724482 | 0.198835980252331 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 8 | 50 | 49 | -1.84803959761637 | 0.0114108637724482 | 0.441185123471937 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 8 | 100 | 99 | -2.61352266275159 | 0.0114108637724482 | 0.735043683625658 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 8 | 200 | 199 | -3.69607919523274 | 0.0114108637724482 | 0.95712175100079 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 12 | 20 | 19 | -3.40360112108769 | 0.0114108637724482 | 0.897461882744365 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 12 | 50 | 49 | -5.38156589466983 | 0.0114108637724482 | 0.999539518647072 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 12 | 100 | 99 | -7.61068347504657 | 0.0114108637724482 | 0.999999987656793 | 0.05 | paired t | post hoc |
| lab | TAms_chrono | 12 | 200 | 199 | -10.7631317893397 | 0.0114108637724482 | 1 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 0 | 20 | 19 | 0 | 0.0453294250315246 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 0 | 50 | 49 | 0 | 0.0453294250315246 | 0.05 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 0 | 100 | 99 | 0 | 0.0453294250315246 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 0 | 200 | 199 | 0 | 0.0453294250315246 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 1 | 20 | 19 | 0.330961704808267 | 0.0453294250315246 | 0.0614112196218065 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 1 | 50 | 49 | 0.523296402743213 | 0.0453294250315246 | 0.0806668559836599 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 1 | 100 | 99 | 0.740052869900505 | 0.0453294250315246 | 0.113440374094912 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 1 | 200 | 199 | 1.04659280548643 | 0.0453294250315246 | 0.180543538425135 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 2 | 20 | 19 | 0.56434908138855 | 0.0453294250315246 | 0.0835527580137878 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 2 | 50 | 49 | 0.892314246305779 | 0.0453294250315246 | 0.141242303583916 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 2 | 100 | 99 | 1.26192290902436 | 0.0453294250315246 | 0.239434275920206 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 2 | 200 | 199 | 1.78462849261156 | 0.0453294250315246 | 0.42712011294894 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 3 | 20 | 19 | 0.684257165717777 | 0.0453294250315246 | 0.0996714581243873 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 3 | 50 | 49 | 1.08190557447973 | 0.0453294250315246 | 0.185540121244118 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 3 | 100 | 99 | 1.53004553663629 | 0.0453294250315246 | 0.328504596537995 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 3 | 200 | 199 | 2.16381114895946 | 0.0453294250315246 | 0.576695372362695 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 4 | 20 | 19 | 0.68251442000639 | 0.0453294250315246 | 0.0994135818344268 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 4 | 50 | 49 | 1.07915005156449 | 0.0453294250315246 | 0.184832953886639 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 4 | 100 | 99 | 1.52614863875813 | 0.0453294250315246 | 0.327114973390419 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 4 | 200 | 199 | 2.15830010312899 | 0.0453294250315246 | 0.574547315782136 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 6 | 20 | 19 | 0.322833708797522 | 0.0453294250315246 | 0.060854278390722 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 6 | 50 | 49 | 0.510444912639854 | 0.0453294250315246 | 0.0791577993948348 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 6 | 100 | 99 | 0.721878118299631 | 0.0453294250315246 | 0.110291398402535 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 6 | 200 | 199 | 1.02088982527971 | 0.0453294250315246 | 0.174046398540755 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 8 | 20 | 19 | -0.418316896136587 | 0.0453294250315246 | 0.0682973906596859 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 8 | 50 | 49 | -0.661417087761852 | 0.0453294250315246 | 0.0994035681895954 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 8 | 100 | 99 | -0.935385015898127 | 0.0453294250315246 | 0.152603308397877 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 8 | 200 | 199 | -1.3228341755237 | 0.0453294250315246 | 0.260472282554318 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 12 | 20 | 19 | -2.20166263228596 | 0.0453294250315246 | 0.551643576693723 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 12 | 50 | 49 | -3.48113427865269 | 0.0453294250315246 | 0.926772454456908 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 12 | 100 | 99 | -4.92306730931252 | 0.0453294250315246 | 0.998220675887366 | 0.05 | paired t | post hoc |
| lab | TAms_mort | 12 | 200 | 199 | -6.96226855730539 | 0.0453294250315246 | 0.999999662490002 | 0.05 | paired t | post hoc |
| lab | RA_default | 0 | 20 | 19 | 0 | 10.836109312986 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| lab | RA_default | 0 | 50 | 49 | 0 | 10.836109312986 | 0.05 | 0.05 | paired t | post hoc |
| lab | RA_default | 0 | 100 | 99 | 0 | 10.836109312986 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| lab | RA_default | 0 | 200 | 199 | 0 | 10.836109312986 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| lab | RA_default | 1 | 20 | 19 | 0.214067299617754 | 10.836109312986 | 0.0547560323447824 | 0.05 | paired t | post hoc |
| lab | RA_default | 1 | 50 | 49 | 0.338470119676897 | 10.836109312986 | 0.0627113431457961 | 0.05 | paired t | post hoc |
| lab | RA_default | 1 | 100 | 99 | 0.478669033705112 | 10.836109312986 | 0.0761162164305704 | 0.05 | paired t | post hoc |
| lab | RA_default | 1 | 200 | 199 | 0.676940239353793 | 10.836109312986 | 0.103395402694777 | 0.05 | paired t | post hoc |
| lab | RA_default | 2 | 20 | 19 | 0.36803816819581 | 10.836109312986 | 0.0641320621512402 | 0.05 | paired t | post hoc |
| lab | RA_default | 2 | 50 | 49 | 0.58191943868745 | 10.836109312986 | 0.0880537013306602 | 0.05 | paired t | post hoc |
| lab | RA_default | 2 | 100 | 99 | 0.822958362400331 | 10.836109312986 | 0.128871735025222 | 0.05 | paired t | post hoc |
| lab | RA_default | 2 | 200 | 199 | 1.1638388773749 | 10.836109312986 | 0.212262360045288 | 0.05 | paired t | post hoc |
| lab | RA_default | 3 | 20 | 19 | 0.451419745489464 | 10.836109312986 | 0.0713412747987299 | 0.05 | paired t | post hoc |
| lab | RA_default | 3 | 50 | 49 | 0.713757288260114 | 10.836109312986 | 0.107723851162107 | 0.05 | paired t | post hoc |
| lab | RA_default | 3 | 100 | 99 | 1.0094052373001 | 10.836109312986 | 0.169984918325999 | 0.05 | paired t | post hoc |
| lab | RA_default | 3 | 200 | 199 | 1.42751457652023 | 10.836109312986 | 0.295189302972448 | 0.05 | paired t | post hoc |
| lab | RA_default | 4 | 20 | 19 | 0.458529714800727 | 10.836109312986 | 0.0720264888028669 | 0.05 | paired t | post hoc |
| lab | RA_default | 4 | 50 | 49 | 0.724999136818859 | 10.836109312986 | 0.109599407170171 | 0.05 | paired t | post hoc |
| lab | RA_default | 4 | 100 | 99 | 1.02530361199802 | 10.836109312986 | 0.173898595093265 | 0.05 | paired t | post hoc |
| lab | RA_default | 4 | 200 | 199 | 1.44999827363772 | 10.836109312986 | 0.302925118912336 | 0.05 | paired t | post hoc |
| lab | RA_default | 6 | 20 | 19 | 0.247227502980122 | 10.836109312986 | 0.0563494534228084 | 0.05 | paired t | post hoc |
| lab | RA_default | 6 | 50 | 49 | 0.390901004826626 | 10.836109312986 | 0.0669940054276074 | 0.05 | paired t | post hoc |
| lab | RA_default | 6 | 100 | 99 | 0.552817502571085 | 10.836109312986 | 0.0849831792441364 | 0.05 | paired t | post hoc |
| lab | RA_default | 6 | 200 | 199 | 0.781802009653252 | 10.836109312986 | 0.121707509909681 | 0.05 | paired t | post hoc |
| lab | RA_default | 8 | 20 | 19 | -0.209250210249764 | 10.836109312986 | 0.0545438342423184 | 0.05 | paired t | post hoc |
| lab | RA_default | 8 | 50 | 49 | -0.330853632629183 | 10.836109312986 | 0.0621418687982579 | 0.05 | paired t | post hoc |
| lab | RA_default | 8 | 100 | 99 | -0.467897694424596 | 10.836109312986 | 0.0749393544150181 | 0.05 | paired t | post hoc |
| lab | RA_default | 8 | 200 | 199 | -0.661707265258365 | 10.836109312986 | 0.100969135952303 | 0.05 | paired t | post hoc |
| lab | RA_default | 12 | 20 | 19 | -1.33555985010098 | 10.836109312986 | 0.245285486437345 | 0.05 | paired t | post hoc |
| lab | RA_default | 12 | 50 | 49 | -2.11170553889608 | 10.836109312986 | 0.543923556909383 | 0.05 | paired t | post hoc |
| lab | RA_default | 12 | 100 | 99 | -2.98640261284523 | 10.836109312986 | 0.840702615659858 | 0.05 | paired t | post hoc |
| lab | RA_default | 12 | 200 | 199 | -4.22341107779217 | 10.836109312986 | 0.98755098136202 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 0 | 20 | 19 | 0 | 12.2050001613507 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 0 | 50 | 49 | 0 | 12.2050001613507 | 0.05 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 0 | 100 | 99 | 0 | 12.2050001613507 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 0 | 200 | 199 | 0 | 12.2050001613507 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 1 | 20 | 19 | -0.290221173581493 | 12.2050001613507 | 0.0587618353472536 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 1 | 50 | 49 | -0.458879966862302 | 12.2050001613507 | 0.073498122997734 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 1 | 100 | 99 | -0.648954272637984 | 12.2050001613507 | 0.0984968984355024 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 1 | 200 | 199 | -0.917759933724603 | 12.2050001613507 | 0.149666666350289 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 2 | 20 | 19 | -0.591227480970732 | 12.2050001613507 | 0.0868808622353443 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 2 | 50 | 49 | -0.934812727575685 | 12.2050001613507 | 0.15039583284696 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 2 | 100 | 99 | -1.32202483761652 | 12.2050001613507 | 0.258143902617622 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 2 | 200 | 199 | -1.86962545515137 | 12.2050001613507 | 0.46049127895925 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 3 | 20 | 19 | -0.882505839755582 | 12.2050001613507 | 0.133621387414665 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 3 | 50 | 49 | -1.39536425101361 | 12.2050001613507 | 0.277399980950267 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 3 | 100 | 99 | -1.97334304823402 | 12.2050001613507 | 0.497739451786296 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 3 | 200 | 199 | -2.79072850202721 | 12.2050001613507 | 0.793114739572147 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 4 | 20 | 19 | -1.1442061111451 | 12.2050001613507 | 0.192495443613259 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 4 | 50 | 49 | -1.80914871195115 | 12.2050001613507 | 0.426203130253274 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 4 | 100 | 99 | -2.55852264479113 | 12.2050001613507 | 0.71690816142208 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 4 | 200 | 199 | -3.6182974239023 | 12.2050001613507 | 0.949582302269054 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 6 | 20 | 19 | -1.51076571313623 | 12.2050001613507 | 0.300076280112285 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 6 | 50 | 49 | -2.38873033219952 | 12.2050001613507 | 0.648685341517698 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 6 | 100 | 99 | -3.37817483264855 | 12.2050001613507 | 0.917013864698574 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 6 | 200 | 199 | -4.77746066439904 | 12.2050001613507 | 0.997399590130472 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 8 | 20 | 19 | -1.59268693761271 | 12.2050001613507 | 0.327525121743652 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 8 | 50 | 49 | -2.51825916122732 | 12.2050001613507 | 0.694519420544765 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 8 | 100 | 99 | -3.56135625937798 | 12.2050001613507 | 0.941400127226017 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 8 | 200 | 199 | -5.03651832245465 | 12.2050001613507 | 0.9988638278368 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 12 | 20 | 19 | -0.896961652935202 | 12.2050001613507 | 0.136456986874614 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 12 | 50 | 49 | -1.41822089855235 | 12.2050001613507 | 0.284917722449649 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 12 | 100 | 99 | -2.00566722917369 | 12.2050001613507 | 0.510502050576665 | 0.05 | paired t | post hoc |
| lab | RA_Peters | 12 | 200 | 199 | -2.83644179710469 | 12.2050001613507 | 0.805867198098819 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 0 | 20 | 19 | 0 | 38.1550911918533 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 0 | 50 | 49 | 0 | 38.1550911918533 | 0.05 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 0 | 100 | 99 | 0 | 38.1550911918533 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 0 | 200 | 199 | 0 | 38.1550911918533 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 1 | 20 | 19 | 0.157046631574013 | 38.1550911918533 | 0.0525564826986047 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 1 | 50 | 49 | 0.248312527315597 | 38.1550911918533 | 0.0568186751606831 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 1 | 100 | 99 | 0.351166943836857 | 38.1550911918533 | 0.0639670671009017 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 1 | 200 | 199 | 0.496625054631195 | 38.1550911918533 | 0.0784207072045597 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 2 | 20 | 19 | 0.328315504141989 | 38.1550911918533 | 0.0612283419227413 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 2 | 50 | 49 | 0.519112392117566 | 38.1550911918533 | 0.0801712234933819 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 2 | 100 | 99 | 0.734135785328601 | 38.1550911918533 | 0.112405964768813 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 2 | 200 | 199 | 1.03822478423513 | 38.1550911918533 | 0.178410016635914 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 3 | 20 | 19 | 0.502134927073615 | 38.1550911918533 | 0.0764735337029574 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 3 | 50 | 49 | 0.793945031137586 | 38.1550911918533 | 0.121790403128394 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 3 | 100 | 99 | 1.1228078308135 | 38.1550911918533 | 0.199274647655761 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 3 | 200 | 199 | 1.58789006227517 | 38.1550911918533 | 0.352270725135481 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 4 | 20 | 19 | 0.66665939394628 | 38.1550911918533 | 0.0970996701442222 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 4 | 50 | 49 | 1.05408105420886 | 38.1550911918533 | 0.17848406176494 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 4 | 100 | 99 | 1.49069572270269 | 38.1550911918533 | 0.314587634021115 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 4 | 200 | 199 | 2.10816210841771 | 38.1550911918533 | 0.55491142101082 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 6 | 20 | 19 | 0.924372697375808 | 38.1550911918533 | 0.141968543270409 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 6 | 50 | 49 | 1.46156156529055 | 38.1550911918533 | 0.299443392049913 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 6 | 100 | 99 | 2.06696018787715 | 38.1550911918533 | 0.534664018196221 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 6 | 200 | 199 | 2.9231231305811 | 38.1550911918533 | 0.828696373686176 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 8 | 20 | 19 | 1.03240134289787 | 38.1550911918533 | 0.16540194757304 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 8 | 50 | 49 | 1.63236985148689 | 38.1550911918533 | 0.359771075334203 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 8 | 100 | 99 | 2.30851958278172 | 38.1550911918533 | 0.62784325620127 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 8 | 200 | 199 | 3.26473970297378 | 38.1550911918533 | 0.901297560672325 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 12 | 20 | 19 | 0.731483897903189 | 38.1550911918533 | 0.106927251553383 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 12 | 50 | 49 | 1.15657759455607 | 38.1550911918533 | 0.205397560751847 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 12 | 100 | 99 | 1.63564772015805 | 38.1550911918533 | 0.367032666723466 | 0.05 | paired t | post hoc |
| lab | RA_GTEx | 12 | 200 | 199 | 2.31315518911214 | 38.1550911918533 | 0.633841009348404 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 0 | 20 | 19 | 0 | 0.0877911523385128 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 0 | 50 | 49 | 0 | 0.0877911523385128 | 0.05 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 0 | 100 | 99 | 0 | 0.0877911523385128 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 0 | 200 | 199 | 0 | 0.0877911523385128 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 1 | 20 | 19 | 0.199038091718434 | 0.0877911523385128 | 0.0541101106753013 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 1 | 50 | 49 | 0.314706855481874 | 0.0877911523385128 | 0.0609785315846476 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 1 | 100 | 99 | 0.445062703194256 | 0.0877911523385128 | 0.0725369652298078 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 1 | 200 | 199 | 0.629413710963748 | 0.0877911523385128 | 0.0960202570700029 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 2 | 20 | 19 | 0.363615976700063 | 0.0877911523385128 | 0.0637919627809757 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 2 | 50 | 49 | 0.574927339999457 | 0.0877911523385128 | 0.0871290994664452 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 2 | 100 | 99 | 0.81307004160632 | 0.0877911523385128 | 0.126939048781768 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 2 | 200 | 199 | 1.14985467999891 | 0.0877911523385128 | 0.20830284166928 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 3 | 20 | 19 | 0.482517944061203 | 0.0877911523385128 | 0.0744207319485956 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 3 | 50 | 49 | 0.762927857567559 | 0.0877911523385128 | 0.116159241613777 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 3 | 100 | 99 | 1.07894292328429 | 0.0877911523385128 | 0.187568364077656 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 3 | 200 | 199 | 1.52585571513512 | 0.0877911523385128 | 0.329688542370074 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 4 | 20 | 19 | 0.547641021220985 | 0.0877911523385128 | 0.0815662234917211 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 4 | 50 | 49 | 0.865896483599459 | 0.0877911523385128 | 0.135780212462896 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 4 | 100 | 99 | 1.22456255071753 | 0.0877911523385128 | 0.228207831988373 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 4 | 200 | 199 | 1.73179296719892 | 0.0877911523385128 | 0.40662614456986 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 6 | 20 | 19 | 0.502765771496284 | 0.0877911523385128 | 0.0765409687796574 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 6 | 50 | 49 | 0.79494248375001 | 0.0877911523385128 | 0.121975480981727 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 6 | 100 | 99 | 1.12421844182582 | 0.0877911523385128 | 0.199658881346407 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 6 | 200 | 199 | 1.58988496750002 | 0.0877911523385128 | 0.353006412928443 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 8 | 20 | 19 | 0.241014514449838 | 0.0877911523385128 | 0.0560332379986568 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 8 | 50 | 49 | 0.381077407410527 | 0.0877911523385128 | 0.0661432265297414 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 8 | 100 | 99 | 0.538924837873944 | 0.0877911523385128 | 0.0832195045568319 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 8 | 200 | 199 | 0.762154814821053 | 0.0877911523385128 | 0.118061581181479 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 12 | 20 | 19 | -0.613253013542293 | 0.0877911523385128 | 0.0897303089522834 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 12 | 50 | 49 | -0.969638152377865 | 0.0877911523385128 | 0.158232925354165 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 12 | 100 | 99 | -1.37127542568717 | 0.0877911523385128 | 0.274048265998102 | 0.05 | paired t | post hoc |
| free-living | TAms_chrono | 12 | 200 | 199 | -1.93927630475573 | 0.0877911523385128 | 0.488063635945646 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 0 | 20 | 19 | 0 | 0.438922406692367 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 0 | 50 | 49 | 0 | 0.438922406692367 | 0.05 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 0 | 100 | 99 | 0 | 0.438922406692367 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 0 | 200 | 199 | 0 | 0.438922406692367 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 1 | 20 | 19 | 0.159501219510619 | 0.438922406692367 | 0.0526371466172445 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 1 | 50 | 49 | 0.252193571614022 | 0.438922406692367 | 0.0570343658274856 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 1 | 100 | 99 | 0.356655569319861 | 0.438922406692367 | 0.0644105802961965 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 1 | 200 | 199 | 0.504387143228044 | 0.438922406692367 | 0.079329075152652 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 2 | 20 | 19 | 0.279415087085329 | 0.438922406692367 | 0.0581185650802226 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 2 | 50 | 49 | 0.44179404390197 | 0.438922406692367 | 0.0717614956373536 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 2 | 100 | 99 | 0.62479112866182 | 0.438922406692367 | 0.0948838465904079 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 2 | 200 | 199 | 0.88358808780394 | 0.438922406692367 | 0.142190482277128 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 3 | 20 | 19 | 0.35156967081595 | 0.438922406692367 | 0.0628869446221191 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 3 | 50 | 49 | 0.555880458007015 | 0.438922406692367 | 0.0846704093412028 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 3 | 100 | 99 | 0.786133682771688 | 0.438922406692367 | 0.121801022888275 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 3 | 200 | 199 | 1.11176091601403 | 0.438922406692367 | 0.197756842529212 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 4 | 20 | 19 | 0.371047755062326 | 0.438922406692367 | 0.0643659413442625 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 4 | 50 | 49 | 0.586678013344611 | 0.438922406692367 | 0.0886897294589715 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 4 | 100 | 99 | 0.829688003218053 | 0.438922406692367 | 0.13020134817678 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 4 | 200 | 199 | 1.17335602668922 | 0.438922406692367 | 0.214983842162857 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 6 | 20 | 19 | 0.250345104546244 | 0.438922406692367 | 0.0565112020373312 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 6 | 50 | 49 | 0.395830365719552 | 0.438922406692367 | 0.0674293555010592 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 6 | 100 | 99 | 0.559788671599693 | 0.438922406692367 | 0.0858860579331722 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 6 | 200 | 199 | 0.791660731439105 | 0.438922406692367 | 0.123574468933049 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 8 | 20 | 19 | -0.0503506867328359 | 0.438922406692367 | 0.0502624239706588 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 8 | 50 | 49 | -0.0796114259146917 | 0.438922406692367 | 0.0506983567759234 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 8 | 100 | 99 | -0.112587558248418 | 0.438922406692367 | 0.0514253438738309 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 8 | 200 | 199 | -0.159222851829383 | 0.438922406692367 | 0.0528811794094634 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 12 | 20 | 19 | -0.842796883590324 | 0.438922406692367 | 0.126084745060927 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 12 | 50 | 49 | -1.33257887851861 | 0.438922406692367 | 0.25728219063653 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 12 | 100 | 99 | -1.88455112293294 | 0.438922406692367 | 0.462736941880292 | 0.05 | paired t | post hoc |
| free-living | TAms_mort | 12 | 200 | 199 | -2.66515775703721 | 0.438922406692367 | 0.75563113531206 | 0.05 | paired t | post hoc |
| free-living | RA_default | 0 | 20 | 19 | 0 | 6.9894052136286 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living | RA_default | 0 | 50 | 49 | 0 | 6.9894052136286 | 0.05 | 0.05 | paired t | post hoc |
| free-living | RA_default | 0 | 100 | 99 | 0 | 6.9894052136286 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living | RA_default | 0 | 200 | 199 | 0 | 6.9894052136286 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| free-living | RA_default | 1 | 20 | 19 | 0.189842773545871 | 6.9894052136286 | 0.0537383077067782 | 0.05 | paired t | post hoc |
| free-living | RA_default | 1 | 50 | 49 | 0.300167780864256 | 6.9894052136286 | 0.059981975350593 | 0.05 | paired t | post hoc |
| free-living | RA_default | 1 | 100 | 99 | 0.424501346685666 | 6.9894052136286 | 0.0704809122936988 | 0.05 | paired t | post hoc |
| free-living | RA_default | 1 | 200 | 199 | 0.600335561728512 | 6.9894052136286 | 0.091789545788872 | 0.05 | paired t | post hoc |
| free-living | RA_default | 2 | 20 | 19 | 0.316963372137846 | 6.9894052136286 | 0.0604608457222924 | 0.05 | paired t | post hoc |
| free-living | RA_default | 2 | 50 | 49 | 0.501163095401573 | 6.9894052136286 | 0.0780924258311043 | 0.05 | paired t | post hoc |
| free-living | RA_default | 2 | 100 | 99 | 0.708751646477785 | 6.9894052136286 | 0.108069209549194 | 0.05 | paired t | post hoc |
| free-living | RA_default | 2 | 200 | 199 | 1.00232619080315 | 6.9894052136286 | 0.169457768349104 | 0.05 | paired t | post hoc |
| free-living | RA_default | 3 | 20 | 19 | 0.372698737058602 | 6.9894052136286 | 0.0644950743757666 | 0.05 | paired t | post hoc |
| free-living | RA_default | 3 | 50 | 49 | 0.589288445086693 | 6.9894052136286 | 0.0890409718198075 | 0.05 | paired t | post hoc |
| free-living | RA_default | 3 | 100 | 99 | 0.833379711191354 | 6.9894052136286 | 0.130935658522455 | 0.05 | paired t | post hoc |
| free-living | RA_default | 3 | 200 | 199 | 1.17857689017339 | 6.9894052136286 | 0.216485940366855 | 0.05 | paired t | post hoc |
| free-living | RA_default | 4 | 20 | 19 | 0.353250595295837 | 6.9894052136286 | 0.0630113491753043 | 0.05 | paired t | post hoc |
| free-living | RA_default | 4 | 50 | 49 | 0.558538232972603 | 6.9894052136286 | 0.085008232473079 | 0.05 | paired t | post hoc |
| free-living | RA_default | 4 | 100 | 99 | 0.789892344173759 | 6.9894052136286 | 0.122506852152566 | 0.05 | paired t | post hoc |
| free-living | RA_default | 4 | 200 | 199 | 1.11707646594521 | 6.9894052136286 | 0.199207185206759 | 0.05 | paired t | post hoc |
| free-living | RA_default | 6 | 20 | 19 | 0.0991385373332227 | 6.9894052136286 | 0.0510178273947433 | 0.05 | paired t | post hoc |
| free-living | RA_default | 6 | 50 | 49 | 0.15675179093531 | 6.9894052136286 | 0.0527106793940144 | 0.05 | paired t | post hoc |
| free-living | RA_default | 6 | 100 | 99 | 0.221680508666987 | 6.9894052136286 | 0.0555393348148979 | 0.05 | paired t | post hoc |
| free-living | RA_default | 6 | 200 | 199 | 0.313503581870619 | 6.9894052136286 | 0.0612238963382533 | 0.05 | paired t | post hoc |
| free-living | RA_default | 8 | 20 | 19 | -0.377283681031904 | 6.9894052136286 | 0.0648567868341098 | 0.05 | paired t | post hoc |
| free-living | RA_default | 8 | 50 | 49 | -0.596537878036642 | 6.9894052136286 | 0.0900250875671678 | 0.05 | paired t | post hoc |
| free-living | RA_default | 8 | 100 | 99 | -0.843631957588686 | 6.9894052136286 | 0.132993192016735 | 0.05 | paired t | post hoc |
| free-living | RA_default | 8 | 200 | 199 | -1.19307575607328 | 6.9894052136286 | 0.220691279496174 | 0.05 | paired t | post hoc |
| free-living | RA_default | 12 | 20 | 19 | -1.46106855265548 | 6.9894052136286 | 0.28396372233818 | 0.05 | paired t | post hoc |
| free-living | RA_default | 12 | 50 | 49 | -2.31015222201849 | 6.9894052136286 | 0.619730795644241 | 0.05 | paired t | post hoc |
| free-living | RA_default | 12 | 100 | 99 | -3.26704860352489 | 6.9894052136286 | 0.898886373312157 | 0.05 | paired t | post hoc |
| free-living | RA_default | 12 | 200 | 199 | -4.62030444403698 | 6.9894052136286 | 0.995829391854732 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 0 | 20 | 19 | 0 | 10.7770588259086 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 0 | 50 | 49 | 0 | 10.7770588259086 | 0.05 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 0 | 100 | 99 | 0 | 10.7770588259086 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 0 | 200 | 199 | 0 | 10.7770588259086 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 1 | 20 | 19 | -0.037774778377927 | 10.7770588259086 | 0.0501476952821078 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 1 | 50 | 49 | -0.0597271688911651 | 10.7770588259086 | 0.0503929974849888 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 1 | 100 | 99 | -0.0844669722880341 | 10.7770588259086 | 0.0508019522566729 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 1 | 200 | 199 | -0.11945433778233 | 10.7770588259086 | 0.0516204442046837 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 2 | 20 | 19 | -0.113447802451189 | 10.7770588259086 | 0.0513331006133491 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 2 | 50 | 49 | -0.179376725643295 | 10.7770588259086 | 0.0535514269042833 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 2 | 100 | 99 | -0.253676998178825 | 10.7770588259086 | 0.0572610581829044 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 2 | 200 | 199 | -0.35875345128659 | 10.7770588259086 | 0.0647262494219743 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 3 | 20 | 19 | -0.221862080684777 | 10.7770588259086 | 0.0551097449521849 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 3 | 50 | 49 | -0.350794750693973 | 10.7770588259086 | 0.0636610523486039 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 3 | 100 | 99 | -0.496098694040705 | 10.7770588259086 | 0.0780800617596416 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 3 | 200 | 199 | -0.701589501387946 | 10.7770588259086 | 0.107446699604341 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 4 | 20 | 19 | -0.35562935918014 | 10.7770588259086 | 0.0631884421976279 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 4 | 50 | 49 | -0.562299388917676 | 10.7770588259086 | 0.0854892158110812 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 4 | 100 | 99 | -0.795211421921481 | 10.7770588259086 | 0.123511871193457 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 4 | 200 | 199 | -1.12459877783535 | 10.7770588259086 | 0.201271453210905 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 6 | 20 | 19 | -0.661652317639619 | 10.7770588259086 | 0.0963809597307172 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 6 | 50 | 49 | -1.0461641714352 | 10.7770588259086 | 0.176510891384548 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 6 | 100 | 99 | -1.47949955971247 | 10.7770588259086 | 0.310675776670756 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 6 | 200 | 199 | -2.0923283428704 | 10.7770588259086 | 0.54868006544638 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 8 | 20 | 19 | -0.949518073245025 | 10.7770588259086 | 0.14717953074545 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 8 | 50 | 49 | -1.50131989547443 | 10.7770588259086 | 0.313065929531339 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 8 | 100 | 99 | -2.1231869576405 | 10.7770588259086 | 0.5567208081714 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 8 | 200 | 199 | -3.00263979094887 | 10.7770588259086 | 0.848061755568744 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 12 | 20 | 19 | -1.18777742812977 | 10.7770588259086 | 0.203819166970809 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 12 | 50 | 49 | -1.87804101311351 | 10.7770588259086 | 0.452802878574807 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 12 | 100 | 99 | -2.65595107143804 | 10.7770588259086 | 0.7486189501103 | 0.05 | paired t | post hoc |
| free-living | RA_Peters | 12 | 200 | 199 | -3.75608202622703 | 10.7770588259086 | 0.962292364948273 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 0 | 20 | 19 | 0 | 7.60125419068026 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 0 | 50 | 49 | 0 | 7.60125419068026 | 0.05 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 0 | 100 | 99 | 0 | 7.60125419068026 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 0 | 200 | 199 | 0 | 7.60125419068026 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 1 | 20 | 19 | 0.0643765247869682 | 7.60125419068026 | 0.0504290334802551 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 1 | 50 | 49 | 0.101788223086553 | 7.60125419068026 | 0.0511419271502734 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 1 | 100 | 99 | 0.143950285578861 | 7.60125419068026 | 0.0523313158352863 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 1 | 200 | 199 | 0.203576446173105 | 7.60125419068026 | 0.0547150660217363 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 2 | 20 | 19 | 0.110058108003604 | 7.60125419068026 | 0.0512545689168173 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 2 | 50 | 49 | 0.174017148130098 | 7.60125419068026 | 0.0533419559243584 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 2 | 100 | 99 | 0.246097410971073 | 7.60125419068026 | 0.0568319380426834 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 2 | 200 | 199 | 0.348034296260197 | 7.60125419068026 | 0.0638527620739189 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 3 | 20 | 19 | 0.13393162524608 | 7.60125419068026 | 0.0518585430744937 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 3 | 50 | 49 | 0.211764493252862 | 7.60125419068026 | 0.0549537799570122 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 3 | 100 | 99 | 0.299480218387263 | 7.60125419068026 | 0.0601364956539231 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 3 | 200 | 199 | 0.423528986505724 | 7.60125419068026 | 0.0705884708193121 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 4 | 20 | 19 | 0.134370135767174 | 7.60125419068026 | 0.0518707467857135 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 4 | 50 | 49 | 0.212457839265164 | 7.60125419068026 | 0.0549863672110379 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 4 | 100 | 99 | 0.300460757721278 | 7.60125419068026 | 0.0602033664225269 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 4 | 200 | 199 | 0.424915678530328 | 7.60125419068026 | 0.0707249843474191 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 6 | 20 | 19 | 0.0664216950850983 | 7.60125419068026 | 0.0504567339028281 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 6 | 50 | 49 | 0.105021921259061 | 7.60125419068026 | 0.0512156894064303 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 6 | 100 | 99 | 0.148523425391043 | 7.60125419068026 | 0.0524820199038106 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 6 | 200 | 199 | 0.210043842518122 | 7.60125419068026 | 0.0550203123667855 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 8 | 20 | 19 | -0.0755804842349074 | 7.60125419068026 | 0.0505914223516578 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 8 | 50 | 49 | -0.119503238420378 | 7.60125419068026 | 0.0515744039049044 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 8 | 100 | 99 | -0.169003100521604 | 7.60125419068026 | 0.0532151051073816 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 8 | 200 | 199 | -0.239006476840756 | 7.60125419068026 | 0.0565059043789996 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 12 | 20 | 19 | -0.419901240004164 | 7.60125419068026 | 0.0684375998513179 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 12 | 50 | 49 | -0.663922155371084 | 7.60125419068026 | 0.0997863791296272 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 12 | 100 | 99 | -0.938927716485764 | 7.60125419068026 | 0.153403550467705 | 0.05 | paired t | post hoc |
| free-living | RA_GTEx | 12 | 200 | 199 | -1.32784431074217 | 7.60125419068026 | 0.26208275521468 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 0 | 20 | 19 | 0 | 0.0877911523385128 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 0 | 50 | 49 | 0 | 0.0877911523385128 | 0.05 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 0 | 100 | 99 | 0 | 0.0877911523385128 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 0 | 200 | 199 | 0 | 0.0877911523385128 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 1 | 20 | 19 | 0.580685698077653 | 0.0877911523385128 | 0.085556070823447 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 1 | 50 | 49 | 0.918144705305122 | 0.0877911523385128 | 0.146751953084888 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 1 | 100 | 99 | 1.29845269446355 | 0.0877911523385128 | 0.250712329044246 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 1 | 200 | 199 | 1.83628941061024 | 0.0877911523385128 | 0.447355322121188 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 2 | 20 | 19 | 1.06083511673212 | 0.0877911523385128 | 0.172020820725417 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 2 | 50 | 49 | 1.67732759538205 | 0.0877911523385128 | 0.376348143419861 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 2 | 100 | 99 | 2.37209943393195 | 0.0877911523385128 | 0.651397825599452 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 2 | 200 | 199 | 3.3546551907641 | 0.0877911523385128 | 0.915969093258857 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 3 | 20 | 19 | 1.40772686656654 | 0.0877911523385128 | 0.267159905250382 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 3 | 50 | 49 | 2.2258116108811 | 0.0877911523385128 | 0.587893341892061 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 3 | 100 | 99 | 3.14777296739555 | 0.0877911523385128 | 0.876390625189227 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 3 | 200 | 199 | 4.45162322176219 | 0.0877911523385128 | 0.993245844484405 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 4 | 20 | 19 | 1.5977208480954 | 0.0877911523385128 | 0.329245514196255 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 4 | 50 | 49 | 2.52621847255868 | 0.0877911523385128 | 0.697247892406599 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 4 | 100 | 99 | 3.57261242540993 | 0.0877911523385128 | 0.942692064346363 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 4 | 200 | 199 | 5.05243694511737 | 0.0877911523385128 | 0.998922360158253 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 6 | 20 | 19 | 1.46679909594326 | 0.0877911523385128 | 0.285799688374672 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 6 | 50 | 49 | 2.31921300652827 | 0.0877911523385128 | 0.623107531494511 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 6 | 100 | 99 | 3.27986248786435 | 0.0877911523385128 | 0.901112965291502 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 6 | 200 | 199 | 4.63842601305653 | 0.0877911523385128 | 0.996045963830237 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 8 | 20 | 19 | 0.703150238036518 | 0.0877911523385128 | 0.1025121380667 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 8 | 50 | 49 | 1.11177814474248 | 0.0877911523385128 | 0.193324145938467 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 8 | 100 | 99 | 1.57229173064481 | 0.0877911523385128 | 0.343722539359023 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 8 | 200 | 199 | 2.22355628948496 | 0.0877911523385128 | 0.59982387198649 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 12 | 20 | 19 | -1.78914122011777 | 0.0877911523385128 | 0.397098534897742 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 12 | 50 | 49 | -2.82888065563241 | 0.0877911523385128 | 0.791914602128481 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 12 | 100 | 99 | -4.00064138953024 | 0.0877911523385128 | 0.977334953531161 | 0.05 | paired t | post hoc |
| free-living upper | TAms_chrono | 12 | 200 | 199 | -5.65776131126481 | 0.0877911523385128 | 0.999878874895024 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 0 | 20 | 19 | 0 | 0.438922406692367 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 0 | 50 | 49 | 0 | 0.438922406692367 | 0.05 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 0 | 100 | 99 | 0 | 0.438922406692367 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 0 | 200 | 199 | 0 | 0.438922406692367 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 1 | 20 | 19 | 0.382040726488807 | 0.438922406692367 | 0.0652368940292871 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 1 | 50 | 49 | 0.604059427325026 | 0.438922406692367 | 0.0910596510200633 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 1 | 100 | 99 | 0.854269034602377 | 0.438922406692367 | 0.135156371789253 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 1 | 200 | 199 | 1.20811885465005 | 0.438922406692367 | 0.225106740845356 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 2 | 20 | 19 | 0.669260982389576 | 0.438922406692367 | 0.0974753776611499 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 2 | 50 | 49 | 1.05819452671645 | 0.438922406692367 | 0.179515330883713 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 2 | 100 | 99 | 1.49651305131138 | 0.438922406692367 | 0.316628701375136 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 2 | 200 | 199 | 2.1163890534329 | 0.438922406692367 | 0.558143915683614 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 3 | 20 | 19 | 0.842087181916586 | 0.438922406692367 | 0.125953412806379 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 3 | 50 | 49 | 1.33145674164448 | 0.438922406692367 | 0.256929956150289 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 3 | 100 | 99 | 1.88296418174672 | 0.438922406692367 | 0.46211323527659 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 3 | 200 | 199 | 2.66291348328897 | 0.438922406692367 | 0.75492949113348 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 4 | 20 | 19 | 0.888741505180869 | 0.438922406692367 | 0.134838534609944 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 4 | 50 | 49 | 1.40522370374894 | 0.438922406692367 | 0.280630474441533 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 4 | 100 | 99 | 1.9872864200099 | 0.438922406692367 | 0.503245123163426 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 4 | 200 | 199 | 2.81044740749788 | 0.438922406692367 | 0.798675381800686 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 6 | 20 | 19 | 0.599631939537589 | 0.438922406692367 | 0.0879551274714954 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 6 | 50 | 49 | 0.948101343361578 | 0.438922406692367 | 0.153350627821641 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 6 | 100 | 99 | 1.34081777828609 | 0.438922406692367 | 0.264153289481583 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 6 | 200 | 199 | 1.89620268672316 | 0.438922406692367 | 0.470996552501596 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 8 | 20 | 19 | -0.120601039902032 | 0.438922406692367 | 0.0515066685129088 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 8 | 50 | 49 | -0.190686987137636 | 0.438922406692367 | 0.0540145052869898 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 8 | 100 | 99 | -0.269672123378109 | 0.438922406692367 | 0.0582100766069771 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 8 | 200 | 199 | -0.381373974275272 | 0.438922406692367 | 0.0666592981382367 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 12 | 20 | 19 | -2.0186850901658 | 0.438922406692367 | 0.482656204144616 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 12 | 50 | 49 | -3.19182138177315 | 0.438922406692367 | 0.878739347829163 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 12 | 100 | 99 | -4.51391708677602 | 0.438922406692367 | 0.993959411876006 | 0.05 | paired t | post hoc |
| free-living upper | TAms_mort | 12 | 200 | 199 | -6.3836427635463 | 0.438922406692367 | 0.999994400768756 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 0 | 20 | 19 | 0 | 6.9894052136286 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 0 | 50 | 49 | 0 | 6.9894052136286 | 0.05 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 0 | 100 | 99 | 0 | 6.9894052136286 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 0 | 200 | 199 | 0 | 6.9894052136286 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 1 | 20 | 19 | 0.342630866965884 | 6.9894052136286 | 0.0622356175566892 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 1 | 50 | 49 | 0.541746968145169 | 6.9894052136286 | 0.0829025209659391 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 1 | 100 | 99 | 0.766145909725403 | 6.9894052136286 | 0.118108114779567 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 1 | 200 | 199 | 1.08349393629034 | 6.9894052136286 | 0.190161009430081 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 2 | 20 | 19 | 0.572059883890071 | 6.9894052136286 | 0.0844907906351923 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 2 | 50 | 49 | 0.904506095552044 | 6.9894052136286 | 0.143822022932789 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 2 | 100 | 99 | 1.27916478757884 | 6.9894052136286 | 0.244720873466827 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 2 | 200 | 199 | 1.80901219110409 | 6.9894052136286 | 0.436650024535665 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 3 | 20 | 19 | 0.672651842418556 | 6.9894052136286 | 0.0979674076985867 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 3 | 50 | 49 | 1.06355594717565 | 6.9894052136286 | 0.180865661810513 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 3 | 100 | 99 | 1.50409524483837 | 6.9894052136286 | 0.319297635687451 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 3 | 200 | 199 | 2.1271118943513 | 6.9894052136286 | 0.562351289582971 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 4 | 20 | 19 | 0.637551566813693 | 6.9894052136286 | 0.093002158977407 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 4 | 50 | 49 | 1.00805753847014 | 6.9894052136286 | 0.167228512243033 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 4 | 100 | 99 | 1.42560864255692 | 6.9894052136286 | 0.292160280086502 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 4 | 200 | 199 | 2.01611507694029 | 6.9894052136286 | 0.518547285173767 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 6 | 20 | 19 | 0.178926605220526 | 6.9894052136286 | 0.0533199413334082 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 6 | 50 | 49 | 0.282907803249318 | 6.9894052136286 | 0.058861395434856 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 6 | 100 | 99 | 0.400092052256364 | 6.9894052136286 | 0.0681711836695155 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 6 | 200 | 199 | 0.565815606498635 | 6.9894052136286 | 0.0870426859867784 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 8 | 20 | 19 | -0.680926812801786 | 6.9894052136286 | 0.0991792718000354 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 8 | 50 | 49 | -1.07663982416637 | 6.9894052136286 | 0.184190339995486 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 8 | 100 | 99 | -1.52259864112707 | 6.9894052136286 | 0.325851194283827 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 8 | 200 | 199 | -2.15327964833275 | 6.9894052136286 | 0.572588557513331 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 12 | 20 | 19 | -2.63695675923096 | 6.9894052136286 | 0.706084228787634 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 12 | 50 | 49 | -4.16939472527303 | 6.9894052136286 | 0.983265750938127 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 12 | 100 | 99 | -5.89641456736797 | 6.9894052136286 | 0.999947442185748 | 0.05 | paired t | post hoc |
| free-living upper | RA_default | 12 | 200 | 199 | -8.33878945054607 | 6.9894052136286 | 0.999999999883728 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 0 | 20 | 19 | 0 | 10.7770588259086 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 0 | 50 | 49 | 0 | 10.7770588259086 | 0.05 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 0 | 100 | 99 | 0 | 10.7770588259086 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 0 | 200 | 199 | 0 | 10.7770588259086 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 1 | 20 | 19 | -0.0751714449154847 | 10.7770588259086 | 0.0505850359160679 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 1 | 50 | 49 | -0.118856490469408 | 10.7770588259086 | 0.0515573927842651 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 1 | 100 | 99 | -0.168088460797905 | 10.7770588259086 | 0.0531803331862056 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 1 | 200 | 199 | -0.237712980938815 | 10.7770588259086 | 0.0564354114342904 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 2 | 20 | 19 | -0.225760033518172 | 10.7770588259086 | 0.0552914264023348 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 2 | 50 | 49 | -0.35695795527669 | 10.7770588259086 | 0.0641490805203699 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 2 | 100 | 99 | -0.504814781549264 | 10.7770588259086 | 0.0790897928854998 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 2 | 200 | 199 | -0.71391591055338 | 10.7770588259086 | 0.10953084957189 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 3 | 20 | 19 | -0.441503402354198 | 10.7770588259086 | 0.0704041954325544 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 3 | 50 | 49 | -0.698078173076506 | 10.7770588259086 | 0.105160343219653 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 3 | 100 | 99 | -0.987231619961427 | 10.7770588259086 | 0.164632639187805 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 3 | 200 | 199 | -1.39615634615301 | 10.7770588259086 | 0.284559949055394 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 4 | 20 | 19 | -0.707698997370162 | 10.7770588259086 | 0.103208406194863 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 4 | 50 | 49 | -1.11897036475361 | 10.7770588259086 | 0.195230274161639 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 4 | 100 | 99 | -1.58246306572813 | 10.7770588259086 | 0.347426820892062 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 4 | 200 | 199 | -2.23794072950722 | 10.7770588259086 | 0.605343841322999 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 6 | 20 | 19 | -1.31668173539074 | 10.7770588259086 | 0.239733826252351 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 6 | 50 | 49 | -2.08185661868893 | 10.7770588259086 | 0.53231000749369 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 6 | 100 | 99 | -2.94418986506608 | 10.7770588259086 | 0.830350811773905 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 6 | 200 | 199 | -4.16371323737786 | 10.7770588259086 | 0.985503233573858 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 8 | 20 | 19 | -1.88953181472279 | 10.7770588259086 | 0.434154684554725 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 8 | 50 | 49 | -2.98761212293765 | 10.7770588259086 | 0.83363215145601 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 8 | 100 | 99 | -4.22512158336869 | 10.7770588259086 | 0.986920214135335 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 8 | 200 | 199 | -5.97522424587529 | 10.7770588259086 | 0.99996643323844 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 12 | 20 | 19 | -2.36366563470525 | 10.7770588259086 | 0.611502075744109 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 12 | 50 | 49 | -3.73728351636807 | 10.7770588259086 | 0.955726862890518 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 12 | 100 | 99 | -5.28531703528113 | 10.7770588259086 | 0.999468847058461 | 0.05 | paired t | post hoc |
| free-living upper | RA_Peters | 12 | 200 | 199 | -7.47456703273614 | 10.7770588259086 | 0.999999978515176 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 0 | 20 | 19 | 0 | 7.60125419068026 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 0 | 50 | 49 | 0 | 7.60125419068026 | 0.05 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 0 | 100 | 99 | 0 | 7.60125419068026 | 0.0500000000000001 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 0 | 200 | 199 | 0 | 7.60125419068026 | 0.0500000000000002 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 1 | 20 | 19 | 0.244656050221471 | 7.60125419068026 | 0.0562175873572773 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 1 | 50 | 49 | 0.386835181020195 | 7.60125419068026 | 0.066639167217347 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 1 | 100 | 99 | 0.547067559401811 | 7.60125419068026 | 0.0842474705265262 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 1 | 200 | 199 | 0.77367036204039 | 7.60125419068026 | 0.120186461704449 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 2 | 20 | 19 | 0.418263988124761 | 7.60125419068026 | 0.0682927179499681 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 2 | 50 | 49 | 0.661333432849932 | 7.60125419068026 | 0.0993908111500333 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 2 | 100 | 99 | 0.935266709987131 | 7.60125419068026 | 0.152576640013233 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 2 | 200 | 199 | 1.32266686569986 | 7.60125419068026 | 0.260418593984969 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 3 | 20 | 19 | 0.508992719642445 | 7.60125419068026 | 0.0772113674879379 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 3 | 50 | 49 | 0.804788153256826 | 7.60125419068026 | 0.123815704183232 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 3 | 100 | 99 | 1.138142321173 | 7.60125419068026 | 0.203477382579869 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 3 | 200 | 199 | 1.60957630651365 | 7.60125419068026 | 0.360298076273529 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 4 | 20 | 19 | 0.510659231657911 | 7.60125419068026 | 0.0773922527220469 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 4 | 50 | 49 | 0.80742314011528 | 7.60125419068026 | 0.124312316083823 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 4 | 100 | 99 | 1.1418687553249 | 7.60125419068026 | 0.204507249653092 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 4 | 200 | 199 | 1.61484628023056 | 7.60125419068026 | 0.362258531847358 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 6 | 20 | 19 | 0.252428499710264 | 7.60125419068026 | 0.0566204399309389 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 6 | 50 | 49 | 0.399124502711794 | 7.60125419068026 | 0.0677234343484366 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 6 | 100 | 99 | 0.564447284810437 | 7.60125419068026 | 0.0864960992060771 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 6 | 200 | 199 | 0.798249005423589 | 7.60125419068026 | 0.12483607017901 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 8 | 20 | 19 | -0.287235491631909 | 7.60125419068026 | 0.0585816165624799 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 8 | 50 | 49 | -0.454159189197533 | 7.60125419068026 | 0.0730114270622248 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 8 | 100 | 99 | -0.64227808483952 | 7.60125419068026 | 0.097484008806531 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 8 | 200 | 199 | -0.908318378395067 | 7.60125419068026 | 0.147570935515682 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 12 | 20 | 19 | -1.59578944657964 | 7.60125419068026 | 0.328584996900607 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 12 | 50 | 49 | -2.52316465862563 | 7.60125419068026 | 0.696202327770711 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 12 | 100 | 99 | -3.56829368032884 | 7.60125419068026 | 0.942199051849431 | 0.05 | paired t | post hoc |
| free-living upper | RA_GTEx | 12 | 200 | 199 | -5.04632931725126 | 7.60125419068026 | 0.998900237467291 | 0.05 | paired t | post hoc |

n = 100에서 위양성률이 0.10을 넘는 가장 작은 s (`results/stage9/sim_design1_smallest_s.tsv`):

| scenario | clock | n | threshold | smallest_s | status | title |
| --- | --- | --- | --- | --- | --- | --- |
| lab | TAms_chrono | 100 | 0.1 | 1 | exceeds | post hoc |
| lab | TAms_mort | 100 | 0.1 | 1 | exceeds | post hoc |
| lab | RA_default | 100 | 0.1 | 2 | exceeds | post hoc |
| lab | RA_Peters | 100 | 0.1 | 2 | exceeds | post hoc |
| lab | RA_GTEx | 100 | 0.1 | 2 | exceeds | post hoc |
| free-living | TAms_chrono | 100 | 0.1 | 2 | exceeds | post hoc |
| free-living | TAms_mort | 100 | 0.1 | 3 | exceeds | post hoc |
| free-living | RA_default | 100 | 0.1 | 2 | exceeds | post hoc |
| free-living | RA_Peters | 100 | 0.1 | 4 | exceeds | post hoc |
| free-living | RA_GTEx | 100 | 0.1 | 12 | exceeds | post hoc |
| free-living upper | TAms_chrono | 100 | 0.1 | 1 | exceeds | post hoc |
| free-living upper | TAms_mort | 100 | 0.1 | 1 | exceeds | post hoc |
| free-living upper | RA_default | 100 | 0.1 | 1 | exceeds | post hoc |
| free-living upper | RA_Peters | 100 | 0.1 | 3 | exceeds | post hoc |
| free-living upper | RA_GTEx | 100 | 0.1 | 2 | exceeds | post hoc |

몬테카를로, lab, n = 50, s = 3, 2000회, 시드 20261019 (`results/stage9/sim_design1_montecarlo.tsv`):

| scenario | clock | n | s | n_rep | n_reject | n_failed | mc_fpr | analytical_fpr | seed | generator | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| lab | TAms_chrono | 50 | 3 | 2000 | 320 | 0 | 0.16 | 0.16495475888083 | 20261019 | model M: f(t) + random intercept + residual; paired difference | post hoc |
| lab | TAms_mort | 50 | 3 | 2000 | 400 | 0 | 0.2 | 0.185540121244118 | 20261019 | model M: f(t) + random intercept + residual; paired difference | post hoc |
| lab | RA_default | 50 | 3 | 2000 | 214 | 0 | 0.107 | 0.107723851162107 | 20261019 | model M: f(t) + random intercept + residual; paired difference | post hoc |
| lab | RA_Peters | 50 | 3 | 2000 | 553 | 0 | 0.2765 | 0.277399980950267 | 20261019 | model M: f(t) + random intercept + residual; paired difference | post hoc |
| lab | RA_GTEx | 50 | 3 | 2000 | 216 | 0 | 0.108 | 0.121790403128394 | 20261019 | model M: f(t) + random intercept + residual; paired difference | post hoc |

### 4.3 설계 2

두 표본 t, 같은 분산, 양측 α = 0.05. SD는 총 SD. ncp = b / total_sd × sqrt(n/2), df = 2n − 2. `results/stage9/sim_design2_fpr.tsv`.

| scenario | clock | s | n_per_group | df | ncp | total_sd | b | fpr | alpha | test | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| lab | TAms_chrono | 0 | 25 | 48 | 0 | 0.0145675409351731 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 0 | 50 | 98 | 0 | 0.0145675409351731 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 0 | 100 | 198 | 0 | 0.0145675409351731 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 0 | 200 | 398 | 0 | 0.0145675409351731 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 1 | 25 | 48 | 0.220209847313472 | 0.0145675409351731 | 0.000907335652950247 | 0.0553536469715115 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 1 | 50 | 98 | 0.311423752638821 | 0.0145675409351731 | 0.000907335652950247 | 0.0609640585198261 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 1 | 100 | 198 | 0.440419694626944 | 0.0145675409351731 | 0.000907335652950247 | 0.0722818312153772 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 1 | 200 | 398 | 0.622847505277641 | 0.0145675409351731 | 0.000907335652950247 | 0.0952688059333402 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 2 | 25 | 48 | 0.353594266859734 | 0.0145675409351731 | 0.00145692251692967 | 0.0638701759974917 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 2 | 50 | 98 | 0.500057807770408 | 0.0145675409351731 | 0.00145692251692967 | 0.0785306771376699 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 2 | 100 | 198 | 0.707188533719469 | 0.0145675409351731 | 0.00145692251692967 | 0.108385649626457 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 2 | 200 | 398 | 1.00011561554082 | 0.0145675409351731 | 0.00145692251692967 | 0.169509303863389 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 3 | 25 | 48 | 0.391063330874884 | 0.0145675409351731 | 0.00161130715539329 | 0.0669941747883293 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 3 | 50 | 98 | 0.553047066270059 | 0.0145675409351731 | 0.00161130715539329 | 0.0850056537731027 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 3 | 100 | 198 | 0.782126661749769 | 0.0145675409351731 | 0.00161130715539329 | 0.121764995965562 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 3 | 200 | 398 | 1.10609413254012 | 0.0145675409351731 | 0.00161130715539329 | 0.196944459480502 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 4 | 25 | 48 | 0.330063584566846 | 0.0145675409351731 | 0.00135996851036251 | 0.0620735857418233 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 4 | 50 | 98 | 0.466780397739913 | 0.0145675409351731 | 0.00135996851036251 | 0.0748138957388647 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 4 | 100 | 198 | 0.660127169133692 | 0.0145675409351731 | 0.00135996851036251 | 0.100718309824148 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 4 | 200 | 398 | 0.933560795479825 | 0.0145675409351731 | 0.00135996851036251 | 0.153739439557103 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 6 | 25 | 48 | -0.0642870195609309 | 0.0145675409351731 | -0.00026488327193883 | 0.0504549388055142 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 6 | 50 | 98 | -0.0909155749476129 | 0.0145675409351731 | -0.00026488327193883 | 0.0509289634307085 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 6 | 100 | 198 | -0.128574039121862 | 0.0145675409351731 | -0.00026488327193883 | 0.0518775132877407 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 6 | 200 | 398 | -0.181831149895226 | 0.0145675409351731 | -0.00026488327193883 | 0.0537776849577626 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 8 | 25 | 48 | -0.723791619612834 | 0.0145675409351731 | -0.00298225510708644 | 0.109346463720358 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 8 | 50 | 98 | -1.02359592478846 | 0.0145675409351731 | -0.00298225510708644 | 0.173449941700561 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 8 | 100 | 198 | -1.44758323922567 | 0.0145675409351731 | -0.00298225510708644 | 0.302077464992125 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 8 | 200 | 398 | -2.04719184957692 | 0.0145675409351731 | -0.00298225510708644 | 0.532821969078047 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 12 | 25 | 48 | -2.10771040835936 | 0.0145675409351731 | -0.00868444723489791 | 0.542028257099216 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 12 | 50 | 98 | -2.98075264505674 | 0.0145675409351731 | -0.00868444723489791 | 0.83926898837844 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 12 | 100 | 198 | -4.21542081671872 | 0.0145675409351731 | -0.00868444723489791 | 0.987288902110044 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_chrono | 12 | 200 | 398 | -5.96150529011348 | 0.0145675409351731 | -0.00868444723489791 | 0.999966557931584 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 0 | 25 | 48 | 0 | 0.061880646182705 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 0 | 50 | 98 | 0 | 0.061880646182705 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 0 | 100 | 198 | 0 | 0.061880646182705 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 0 | 200 | 398 | 0 | 0.061880646182705 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 1 | 25 | 48 | 0.191665137007285 | 0.061880646182705 | 0.00335461710855195 | 0.054052564757886 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 1 | 50 | 98 | 0.2710554361898 | 0.061880646182705 | 0.00335461710855195 | 0.0582932719380374 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 1 | 100 | 198 | 0.38333027401457 | 0.061880646182705 | 0.00335461710855195 | 0.0668313820335917 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 1 | 200 | 398 | 0.542110872379599 | 0.061880646182705 | 0.00335461710855195 | 0.0841223863904338 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 2 | 25 | 48 | 0.326823443416012 | 0.061880646182705 | 0.00572022399001831 | 0.0618361401994675 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 2 | 50 | 98 | 0.4621981461804 | 0.061880646182705 | 0.00572022399001831 | 0.0743230750514935 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 2 | 100 | 198 | 0.653646886832025 | 0.061880646182705 | 0.00572022399001831 | 0.0997066790874509 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 2 | 200 | 398 | 0.924396292360801 | 0.061880646182705 | 0.00572022399001831 | 0.151657280085991 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 3 | 25 | 48 | 0.39626410400409 | 0.061880646182705 | 0.00693560844477734 | 0.0674534714870174 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 3 | 50 | 98 | 0.560402070164206 | 0.061880646182705 | 0.00693560844477734 | 0.0859588138901046 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 3 | 100 | 198 | 0.792528208008179 | 0.061880646182705 | 0.00693560844477734 | 0.123736249096962 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 3 | 200 | 398 | 1.12080414032841 | 0.061880646182705 | 0.00693560844477734 | 0.200974036484743 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 4 | 25 | 48 | 0.395254852508557 | 0.061880646182705 | 0.006917944030755 | 0.0673638502392566 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 4 | 50 | 98 | 0.558974773011379 | 0.061880646182705 | 0.006917944030755 | 0.0857728036351863 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 4 | 100 | 198 | 0.790509705017114 | 0.061880646182705 | 0.006917944030755 | 0.12335153033894 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 4 | 200 | 398 | 1.11794954602276 | 0.061880646182705 | 0.006917944030755 | 0.200187904884313 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 6 | 25 | 48 | 0.186958086474364 | 0.061880646182705 | 0.00327223200453611 | 0.0538555075966922 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 6 | 50 | 98 | 0.264398661487367 | 0.061880646182705 | 0.00327223200453611 | 0.057889098144843 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 6 | 100 | 198 | 0.373916172948728 | 0.061880646182705 | 0.00327223200453611 | 0.0660077179438767 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 6 | 200 | 398 | 0.528797322974735 | 0.061880646182705 | 0.00327223200453611 | 0.0824412115253661 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 8 | 25 | 48 | -0.242253904441691 | 0.061880646182705 | -0.00424004649537656 | 0.0564834155322713 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 8 | 50 | 98 | -0.342598757199275 | 0.061880646182705 | -0.00424004649537656 | 0.0632861974404418 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 8 | 100 | 198 | -0.484507808883381 | 0.061880646182705 | -0.00424004649537656 | 0.0770308967000197 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 8 | 200 | 398 | -0.685197514398549 | 0.061880646182705 | -0.00424004649537656 | 0.105006365410268 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 12 | 25 | 48 | -1.2750175139005 | 0.061880646182705 | -0.0223159810522631 | 0.239398926096095 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 12 | 50 | 98 | -1.80314706042131 | 0.061880646182705 | -0.0223159810522631 | 0.430829220892027 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 12 | 100 | 198 | -2.55003502780099 | 0.061880646182705 | -0.0223159810522631 | 0.718265113096993 | 0.05 | two-sample t, equal variance | post hoc |
| lab | TAms_mort | 12 | 200 | 398 | -3.60629412084262 | 0.061880646182705 | -0.0223159810522631 | 0.949248954787359 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 0 | 25 | 48 | 0 | 11.6771259036261 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 0 | 50 | 98 | 0 | 11.6771259036261 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 0 | 100 | 198 | 0 | 11.6771259036261 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 0 | 200 | 398 | 0 | 11.6771259036261 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 1 | 25 | 48 | 0.157046316288299 | 11.6771259036261 | 0.518690997396996 | 0.0527186591607556 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 1 | 50 | 98 | 0.222097030415647 | 11.6771259036261 | 0.518690997396996 | 0.0555591319761699 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 1 | 100 | 198 | 0.314092632576598 | 11.6771259036261 | 0.518690997396996 | 0.0612658287853699 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 1 | 200 | 398 | 0.444194060831294 | 11.6771259036261 | 0.518690997396996 | 0.0727820947688259 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 2 | 25 | 48 | 0.270004053266675 | 11.6771259036261 | 0.891766677500593 | 0.0580611704195056 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 2 | 50 | 98 | 0.381843394025439 | 11.6771259036261 | 0.891766677500593 | 0.066533612306331 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 2 | 100 | 198 | 0.54000810653335 | 11.6771259036261 | 0.891766677500593 | 0.0836864556466729 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 2 | 200 | 398 | 0.763686788050879 | 11.6771259036261 | 0.891766677500593 | 0.11868168922199 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 3 | 25 | 48 | 0.331175327831539 | 11.6771259036261 | 1.09380254924844 | 0.0621556107070903 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 3 | 50 | 98 | 0.468352640142718 | 11.6771259036261 | 1.09380254924844 | 0.0749834710422724 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 3 | 100 | 198 | 0.662350655663077 | 11.6771259036261 | 1.09380254924844 | 0.101067872924509 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 3 | 200 | 398 | 0.936705280285436 | 11.6771259036261 | 1.09380254924844 | 0.15445885866032 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 4 | 25 | 48 | 0.336391418711605 | 11.6771259036261 | 1.11103020186099 | 0.0625442362321558 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 4 | 50 | 98 | 0.475729306607879 | 11.6771259036261 | 1.11103020186099 | 0.0757870565394282 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 4 | 100 | 198 | 0.672782837423211 | 11.6771259036261 | 1.11103020186099 | 0.102724735193889 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 4 | 200 | 398 | 0.951458613215758 | 11.6771259036261 | 1.11103020186099 | 0.157868293473596 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 6 | 25 | 48 | 0.181373655376193 | 11.6771259036261 | 0.599039088799219 | 0.0536281306102746 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 6 | 50 | 98 | 0.256501083290195 | 11.6771259036261 | 0.599039088799219 | 0.0574228481334654 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 6 | 100 | 198 | 0.362747310752385 | 11.6771259036261 | 0.599039088799219 | 0.0650579440628883 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 6 | 200 | 398 | 0.513002166580391 | 11.6771259036261 | 0.599039088799219 | 0.0805038154743072 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 8 | 25 | 48 | -0.153512352241361 | 11.6771259036261 | -0.507019056407903 | 0.052597493327987 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 8 | 50 | 98 | -0.217099250531528 | 11.6771259036261 | -0.507019056407903 | 0.0553109833590944 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 8 | 100 | 198 | -0.307024704482721 | 11.6771259036261 | -0.507019056407903 | 0.0607614618079939 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 8 | 200 | 398 | -0.434198501063056 | 11.6771259036261 | -0.507019056407903 | 0.0717568345734548 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 12 | 25 | 48 | -0.979807541905932 | 11.6771259036261 | -3.23609851653779 | 0.160484690962897 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 12 | 50 | 98 | -1.38565711427881 | 11.6771259036261 | -3.23609851653779 | 0.278740233160922 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 12 | 100 | 198 | -1.95961508381187 | 11.6771259036261 | -3.23609851653779 | 0.496113513415064 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_default | 12 | 200 | 398 | -2.77131422855763 | 11.6771259036261 | -3.23609851653779 | 0.789492394882309 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 0 | 25 | 48 | 0 | 16.649562072426 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 0 | 50 | 98 | 0 | 16.649562072426 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 0 | 100 | 198 | 0 | 16.649562072426 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 0 | 200 | 398 | 0 | 16.649562072426 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 1 | 25 | 48 | -0.168191512942867 | 16.649562072426 | -0.79204870022558 | 0.053118969739538 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 1 | 50 | 98 | -0.237858718679852 | 16.649562072426 | -0.79204870022558 | 0.0563792201522185 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 1 | 100 | 198 | -0.336383025885733 | 16.649562072426 | -0.79204870022558 | 0.0629335920890122 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 1 | 200 | 398 | -0.475717437359704 | 16.649562072426 | -0.79204870022558 | 0.0761752920479935 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 2 | 25 | 48 | -0.342633320962521 | 16.649562072426 | -1.61353133564192 | 0.0630174874035266 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 2 | 50 | 98 | -0.484556689426131 | 16.649562072426 | -1.61353133564192 | 0.0767659710099606 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 2 | 100 | 198 | -0.685266641925043 | 16.649562072426 | -1.61353133564192 | 0.104743844337468 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 2 | 200 | 398 | -0.969113378852262 | 16.649562072426 | -1.61353133564192 | 0.162021906337018 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 3 | 25 | 48 | -0.511437503121142 | 16.649562072426 | -2.40846522220964 | 0.0792485711162898 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 3 | 50 | 98 | -0.723281853220151 | 16.649562072426 | -2.40846522220964 | 0.110519312303115 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 3 | 100 | 198 | -1.02287500624228 | 16.649562072426 | -2.40846522220964 | 0.174536017948851 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 3 | 200 | 398 | -1.4465637064403 | 16.649562072426 | -2.40846522220964 | 0.302949633730386 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 4 | 25 | 48 | -0.663100333366721 | 16.649562072426 | -3.12267692924951 | 0.099618950806567 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 4 | 50 | 98 | -0.937765484661337 | 16.649562072426 | -3.12267692924951 | 0.153119583321264 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 4 | 100 | 198 | -1.32620066673344 | 16.649562072426 | -3.12267692924951 | 0.261543396224426 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 4 | 200 | 398 | -1.87553096932267 | 16.649562072426 | -3.12267692924951 | 0.464620940319167 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 6 | 25 | 48 | -0.87553216003808 | 16.649562072426 | -4.12306243775463 | 0.13767797614513 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 6 | 50 | 98 | -1.23818945501966 | 16.649562072426 | -4.12306243775463 | 0.232229025030839 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 6 | 100 | 198 | -1.75106432007616 | 16.649562072426 | -4.12306243775463 | 0.414057504400519 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 6 | 200 | 398 | -2.47637891003933 | 16.649562072426 | -4.12306243775463 | 0.695131556716286 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 8 | 25 | 48 | -0.923007864573339 | 16.649562072426 | -4.34663537203354 | 0.14772511824988 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 8 | 50 | 98 | -1.30533024025664 | 16.649562072426 | -4.34663537203354 | 0.252827393142204 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 8 | 100 | 198 | -1.84601572914668 | 16.649562072426 | -4.34663537203354 | 0.45116468407676 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 8 | 200 | 398 | -2.61066048051329 | 16.649562072426 | -4.34663537203354 | 0.740341964441599 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 12 | 25 | 48 | -0.519815062413235 | 16.649562072426 | -2.44791688556806 | 0.0802289655420244 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 12 | 50 | 98 | -0.735129511190614 | 16.649562072426 | -2.44791688556806 | 0.112566331045418 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 12 | 100 | 198 | -1.03963012482647 | 16.649562072426 | -2.44791688556806 | 0.178760632973903 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_Peters | 12 | 200 | 398 | -1.47025902238123 | 16.649562072426 | -2.44791688556806 | 0.311222762020996 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 0 | 25 | 48 | 0 | 30.2461600026282 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 0 | 50 | 98 | 0 | 30.2461600026282 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 0 | 100 | 198 | 0 | 30.2461600026282 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 0 | 200 | 398 | 0 | 30.2461600026282 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 1 | 25 | 48 | 0.156621321896447 | 30.2461600026282 | 1.33988067656597 | 0.0527039409615525 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 1 | 50 | 98 | 0.221495997582758 | 30.2461600026282 | 1.33988067656597 | 0.0555289871078398 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 1 | 100 | 198 | 0.313242643792895 | 30.2461600026282 | 1.33988067656597 | 0.0612045517920359 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 1 | 200 | 398 | 0.442991995165516 | 30.2461600026282 | 1.33988067656597 | 0.0726575097930531 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 2 | 25 | 48 | 0.327426368476952 | 30.2461600026282 | 2.80110178364156 | 0.0618801421219042 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 2 | 50 | 98 | 0.463050810978676 | 30.2461600026282 | 2.80110178364156 | 0.0744140235255015 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 2 | 100 | 198 | 0.654852736953904 | 30.2461600026282 | 2.80110178364156 | 0.0998941161592696 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 2 | 200 | 398 | 0.926101621957352 | 30.2461600026282 | 2.80110178364156 | 0.152043084979774 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 3 | 25 | 48 | 0.500775058085738 | 30.2461600026282 | 4.28408351756162 | 0.0780249916072955 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 3 | 50 | 98 | 0.708202878843026 | 30.2461600026282 | 4.28408351756162 | 0.107965474699451 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 3 | 100 | 198 | 1.00155011617148 | 30.2461600026282 | 4.28408351756162 | 0.169261851324697 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 3 | 200 | 398 | 1.41640575768605 | 30.2461600026282 | 4.28408351756162 | 0.29256970421323 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 4 | 25 | 48 | 0.664853963998222 | 30.2461600026282 | 5.68776312390269 | 0.099887266826021 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 4 | 50 | 98 | 0.940245492883799 | 30.2461600026282 | 5.68776312390269 | 0.15368083128265 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 4 | 100 | 198 | 1.32970792799644 | 30.2461600026282 | 5.68776312390269 | 0.262672673073346 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 4 | 200 | 398 | 1.8804909857676 | 30.2461600026282 | 5.68776312390269 | 0.466586275799588 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 6 | 25 | 48 | 0.92186933483991 | 30.2461600026282 | 7.88650544583827 | 0.147477579267601 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 6 | 50 | 98 | 1.30372011606647 | 30.2461600026282 | 7.88650544583827 | 0.252321836625227 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 6 | 100 | 198 | 1.84373866967982 | 30.2461600026282 | 7.88650544583827 | 0.450268288744589 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 6 | 200 | 398 | 2.60744023213293 | 30.2461600026282 | 7.88650544583827 | 0.739299619635145 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 8 | 25 | 48 | 1.0296054199426 | 30.2461600026282 | 8.80817751992155 | 0.172329749535056 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 8 | 50 | 98 | 1.45608194877567 | 30.2461600026282 | 8.80817751992155 | 0.302515201339528 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 8 | 100 | 198 | 2.0592108398852 | 30.2461600026282 | 8.80817751992155 | 0.535590070569922 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 8 | 200 | 398 | 2.91216389755135 | 30.2461600026282 | 8.80817751992155 | 0.827713480624163 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 12 | 25 | 48 | 0.729502911888758 | 30.2461600026282 | 6.24082879203771 | 0.110308806790912 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 12 | 50 | 98 | 1.03167291178375 | 30.2461600026282 | 6.24082879203771 | 0.175458646310955 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 12 | 100 | 198 | 1.45900582377752 | 30.2461600026282 | 6.24082879203771 | 0.306037962922577 | 0.05 | two-sample t, equal variance | post hoc |
| lab | RA_GTEx | 12 | 200 | 398 | 2.06334582356749 | 30.2461600026282 | 6.24082879203771 | 0.539222569252544 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 0 | 25 | 48 | 0 | 0.0843795698011575 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 0 | 50 | 98 | 0 | 0.0843795698011575 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 0 | 100 | 198 | 0 | 0.0843795698011575 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 0 | 200 | 398 | 0 | 0.0843795698011575 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 1 | 25 | 48 | 0.163715444133519 | 0.0843795698011575 | 0.00390725675763172 | 0.0529548797984718 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 1 | 50 | 98 | 0.231528601463558 | 0.0843795698011575 | 0.00390725675763172 | 0.056043014585381 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 1 | 100 | 198 | 0.327430888267038 | 0.0843795698011575 | 0.00390725675763172 | 0.0622497022989846 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 1 | 200 | 398 | 0.463057202927115 | 0.0843795698011575 | 0.00390725675763172 | 0.0747833127467541 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 2 | 25 | 48 | 0.29908622317233 | 0.0843795698011575 | 0.00713803558845414 | 0.0599015981375636 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 2 | 50 | 98 | 0.422971793129255 | 0.0843795698011575 | 0.00713803558845414 | 0.070327927486344 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 2 | 100 | 198 | 0.598172446344659 | 0.0843795698011575 | 0.00713803558845414 | 0.0914812630325457 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 2 | 200 | 398 | 0.84594358625851 | 0.0843795698011575 | 0.00713803558845414 | 0.134722836207532 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 3 | 25 | 48 | 0.396887042235726 | 0.0843795698011575 | 0.00947216425426115 | 0.0675089063047283 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 3 | 50 | 98 | 0.561283037859908 | 0.0843795698011575 | 0.00947216425426115 | 0.0860738751965078 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 3 | 100 | 198 | 0.793774084471452 | 0.0843795698011575 | 0.00947216425426115 | 0.12397423268144 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 3 | 200 | 398 | 1.12256607571981 | 0.0843795698011575 | 0.00947216425426115 | 0.201460256327024 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 4 | 25 | 48 | 0.450452937128034 | 0.0843795698011575 | 0.010750575743808 | 0.0726142176529522 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 4 | 50 | 98 | 0.637036652897261 | 0.0843795698011575 | 0.010750575743808 | 0.0966871595208289 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 4 | 100 | 198 | 0.900905874256068 | 0.0843795698011575 | 0.010750575743808 | 0.14593686612444 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 4 | 200 | 398 | 1.27407330579452 | 0.0843795698011575 | 0.010750575743808 | 0.246039834557784 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 6 | 25 | 48 | 0.413541553101729 | 0.0843795698011575 | 0.00986964324880958 | 0.0690244425257341 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 6 | 50 | 98 | 0.584836073001299 | 0.0843795698011575 | 0.00986964324880958 | 0.0892211329159162 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 6 | 100 | 198 | 0.827083106203458 | 0.0843795698011575 | 0.00986964324880958 | 0.13048532140117 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 6 | 200 | 398 | 1.1696721460026 | 0.0843795698011575 | 0.00986964324880958 | 0.214739527026052 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 8 | 25 | 48 | 0.198242446634778 | 0.0843795698011575 | 0.00473128325407997 | 0.054336204913862 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 8 | 50 | 98 | 0.280357156668928 | 0.0843795698011575 | 0.00473128325407997 | 0.0588751875540401 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 8 | 100 | 198 | 0.396484893269557 | 0.0843795698011575 | 0.00473128325407997 | 0.0680178145061677 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 8 | 200 | 398 | 0.560714313337856 | 0.0843795698011575 | 0.00473128325407997 | 0.0865456065106364 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 12 | 25 | 48 | -0.504420980986511 | 0.0843795698011575 | -0.0120385849794561 | 0.0784403364085932 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 12 | 50 | 98 | -0.713358992456665 | 0.0843795698011575 | -0.0120385849794561 | 0.108832254054911 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 12 | 100 | 198 | -1.00884196197302 | 0.0843795698011575 | -0.0120385849794561 | 0.171052333040761 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_chrono | 12 | 200 | 398 | -1.42671798491333 | 0.0843795698011575 | -0.0120385849794561 | 0.296099780490004 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 0 | 25 | 48 | 0 | 0.387677419351225 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 0 | 50 | 98 | 0 | 0.387677419351225 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 0 | 100 | 198 | 0 | 0.387677419351225 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 0 | 200 | 398 | 0 | 0.387677419351225 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 1 | 25 | 48 | 0.142764839889843 | 0.387677419351225 | 0.015654412084611 | 0.0522460471645631 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 1 | 50 | 98 | 0.201899972802239 | 0.387677419351225 | 0.015654412084611 | 0.0545914186017243 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 1 | 100 | 198 | 0.285529679779686 | 0.387677419351225 | 0.015654412084611 | 0.0592996765818354 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 1 | 200 | 398 | 0.403799945604479 | 0.387677419351225 | 0.015654412084611 | 0.0687878786863757 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 2 | 25 | 48 | 0.250096207997255 | 0.387677419351225 | 0.0274234825872286 | 0.0569116903092861 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 2 | 50 | 98 | 0.3536894492478 | 0.387677419351225 | 0.0274234825872286 | 0.0641671872795934 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 2 | 100 | 198 | 0.500192415994509 | 0.387677419351225 | 0.0274234825872286 | 0.0788349155333796 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 2 | 200 | 398 | 0.7073788984956 | 0.387677419351225 | 0.0274234825872286 | 0.108710303725491 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 3 | 25 | 48 | 0.314679648959187 | 0.387677419351225 | 0.0345051688024083 | 0.0609675676809457 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 3 | 50 | 98 | 0.445024227360887 | 0.387677419351225 | 0.0345051688024083 | 0.0725284986100671 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 3 | 100 | 198 | 0.629359297918374 | 0.387677419351225 | 0.0345051688024083 | 0.0960098494589039 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 3 | 200 | 398 | 0.890048454721774 | 0.387677419351225 | 0.0345051688024083 | 0.144046790327845 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 4 | 25 | 48 | 0.332113908003267 | 0.387677419351225 | 0.0364168655176252 | 0.0622250799209513 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 4 | 50 | 98 | 0.46967999295095 | 0.387677419351225 | 0.0364168655176252 | 0.075127098064471 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 4 | 100 | 198 | 0.664227816006533 | 0.387677419351225 | 0.0364168655176252 | 0.101363966459141 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 4 | 200 | 398 | 0.9393599859019 | 0.387677419351225 | 0.0364168655176252 | 0.155068208943327 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 6 | 25 | 48 | 0.224076523536367 | 0.387677419351225 | 0.024570379097766 | 0.0555439247184076 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 6 | 50 | 98 | 0.316892058594544 | 0.387677419351225 | 0.024570379097766 | 0.0613549645655218 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 6 | 100 | 198 | 0.448153047072734 | 0.387677419351225 | 0.024570379097766 | 0.0730806492811701 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 6 | 200 | 398 | 0.633784117189089 | 0.387677419351225 | 0.024570379097766 | 0.0969052346003626 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 8 | 25 | 48 | -0.0450674154831674 | 0.387677419351225 | -0.00494172020300127 | 0.0502235478574016 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 8 | 50 | 98 | -0.0637349501973985 | 0.387677419351225 | -0.00494172020300127 | 0.0504564071783159 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 8 | 100 | 198 | -0.0901348309663348 | 0.387677419351225 | -0.00494172020300127 | 0.0509221690977382 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 8 | 200 | 398 | -0.127469900394797 | 0.387677419351225 | -0.00494172020300127 | 0.0518544037585341 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 12 | 25 | 48 | -0.754362646972868 | 0.387677419351225 | -0.0827171714412529 | 0.114592081455277 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 12 | 50 | 98 | -1.0668298862967 | 0.387677419351225 | -0.0827171714412529 | 0.184391524159916 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 12 | 100 | 198 | -1.50872529394574 | 0.387677419351225 | -0.0827171714412529 | 0.323545835512233 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | TAms_mort | 12 | 200 | 398 | -2.13365977259339 | 0.387677419351225 | -0.0827171714412529 | 0.566944614922377 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 0 | 25 | 48 | 0 | 7.86497441621665 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 0 | 50 | 98 | 0 | 7.86497441621665 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 0 | 100 | 198 | 0 | 7.86497441621665 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 0 | 200 | 398 | 0 | 7.86497441621665 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 1 | 25 | 48 | 0.133375783664988 | 7.86497441621665 | 0.296701192571716 | 0.051959998338423 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 1 | 50 | 98 | 0.188621842151166 | 7.86497441621665 | 0.296701192571716 | 0.0540059721497429 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 1 | 100 | 198 | 0.266751567329977 | 7.86497441621665 | 0.296701192571716 | 0.0581111817353313 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 1 | 200 | 398 | 0.377243684302333 | 7.86497441621665 | 0.296701192571716 | 0.066376850470101 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 2 | 25 | 48 | 0.2226855274097 | 7.86497441621665 | 0.495375245305969 | 0.0554750880813367 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 2 | 50 | 98 | 0.314924893007003 | 7.86497441621665 | 0.495375245305969 | 0.0612135376526829 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 2 | 100 | 198 | 0.4453710548194 | 7.86497441621665 | 0.495375245305969 | 0.0727916124014737 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 2 | 200 | 398 | 0.629849786014006 | 7.86497441621665 | 0.495375245305969 | 0.0963130497711159 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 3 | 25 | 48 | 0.261842919789262 | 7.86497441621665 | 0.582482849833315 | 0.057579123117579 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 3 | 50 | 98 | 0.370301808377345 | 7.86497441621665 | 0.582482849833315 | 0.065540898509767 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 3 | 100 | 198 | 0.523685839578524 | 7.86497441621665 | 0.582482849833315 | 0.0816502458328177 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 3 | 200 | 398 | 0.740603616754689 | 7.86497441621665 | 0.582482849833315 | 0.114495031513086 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 4 | 25 | 48 | 0.248179449223658 | 7.86497441621665 | 0.552087766857339 | 0.0568057354091894 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 4 | 50 | 98 | 0.350978742994382 | 7.86497441621665 | 0.552087766857339 | 0.0639491952825565 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 4 | 100 | 198 | 0.496358898447316 | 7.86497441621665 | 0.552087766857339 | 0.0783884169876189 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 4 | 200 | 398 | 0.701957485988764 | 7.86497441621665 | 0.552087766857339 | 0.107793349171666 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 6 | 25 | 48 | 0.0696506896799226 | 7.86497441621665 | 0.154941490303688 | 0.0505340453488954 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 6 | 50 | 98 | 0.0985009499739863 | 7.86497441621665 | 0.154941490303688 | 0.0510905497203941 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 6 | 100 | 198 | 0.139301379359845 | 7.86497441621665 | 0.154941490303688 | 0.0522043105007556 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 6 | 200 | 398 | 0.197001899947973 | 7.86497441621665 | 0.154941490303688 | 0.0544360812972682 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 8 | 25 | 48 | -0.265064114276033 | 7.86497441621665 | -0.589648560275406 | 0.0577675887928374 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 8 | 50 | 98 | -0.374857265307577 | 7.86497441621665 | -0.589648560275406 | 0.0659289634038538 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 8 | 100 | 198 | -0.530128228552065 | 7.86497441621665 | -0.589648560275406 | 0.0824460606499854 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 8 | 200 | 398 | -0.749714530615154 | 7.86497441621665 | -0.589648560275406 | 0.11613105709449 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 12 | 25 | 48 | -1.02648712699938 | 7.86497441621665 | -2.28347265426549 | 0.171570080420318 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 12 | 50 | 98 | -1.45167201660392 | 7.86497441621665 | -2.28347265426549 | 0.300999135070376 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 12 | 100 | 198 | -2.05297425399876 | 7.86497441621665 | -2.28347265426549 | 0.533124133697149 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_default | 12 | 200 | 398 | -2.90334403320784 | 7.86497441621665 | -2.28347265426549 | 0.825458568391616 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 0 | 25 | 48 | 0 | 17.2404170486569 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 0 | 50 | 98 | 0 | 17.2404170486569 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 0 | 100 | 198 | 0 | 17.2404170486569 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 0 | 200 | 398 | 0 | 17.2404170486569 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 1 | 25 | 48 | -0.0186678550416856 | 17.2404170486569 | -0.0910305529194535 | 0.050038351732797 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 1 | 50 | 98 | -0.0264003337803667 | 17.2404170486569 | -0.0910305529194535 | 0.050078291704646 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 1 | 100 | 198 | -0.0373357100833712 | 17.2404170486569 | -0.0910305529194535 | 0.0501581509732633 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 1 | 200 | 398 | -0.0528006675607334 | 17.2404170486569 | -0.0910305529194535 | 0.0503178797134135 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 2 | 25 | 48 | -0.05606458123376 | 17.2404170486569 | -0.273389193215311 | 0.0503459830031648 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 2 | 50 | 98 | -0.0792872911495515 | 17.2404170486569 | -0.273389193215311 | 0.0507064321453169 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 2 | 100 | 198 | -0.11212916246752 | 17.2404170486569 | -0.273389193215311 | 0.0514275638269902 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 2 | 200 | 398 | -0.158574582299103 | 17.2404170486569 | -0.273389193215311 | 0.0528715846163812 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 3 | 25 | 48 | -0.109641653487246 | 17.2404170486569 | -0.534648480913311 | 0.0513239951554073 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 3 | 50 | 98 | -0.155056713362674 | 17.2404170486569 | -0.534648480913311 | 0.0527050070865339 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 3 | 100 | 198 | -0.219283306974491 | 17.2404170486569 | -0.534648480913311 | 0.0554728840814137 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 3 | 200 | 398 | -0.310113426725348 | 17.2404170486569 | -0.534648480913311 | 0.0610344653512057 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 4 | 25 | 48 | -0.175747882868365 | 17.2404170486569 | -0.857004027308223 | 0.0534061002705143 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 4 | 50 | 98 | -0.2485450395108 | 17.2404170486569 | -0.857004027308223 | 0.0569676775456934 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 4 | 100 | 198 | -0.351495765736731 | 17.2404170486569 | -0.857004027308223 | 0.0641311541987358 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 4 | 200 | 398 | -0.497090079021601 | 17.2404170486569 | -0.857004027308223 | 0.0786145698079597 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 6 | 25 | 48 | -0.326980860883333 | 17.2404170486569 | -1.59446537879272 | 0.0618476206223033 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 6 | 50 | 98 | -0.46242076809764 | 17.2404170486569 | -1.59446537879272 | 0.0743468038492738 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 6 | 100 | 198 | -0.653961721766666 | 17.2404170486569 | -1.59446537879272 | 0.0997555813755188 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 6 | 200 | 398 | -0.92484153619528 | 17.2404170486569 | -1.59446537879272 | 0.151757937240817 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 8 | 25 | 48 | -0.469240760950599 | 17.2404170486569 | -2.28817107408936 | 0.0745642936749357 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 8 | 50 | 98 | -0.663606648154608 | 17.2404170486569 | -2.28817107408936 | 0.100748828377309 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 8 | 100 | 198 | -0.938481521901198 | 17.2404170486569 | -2.28817107408936 | 0.154341388602704 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 8 | 200 | 398 | -1.32721329630922 | 17.2404170486569 | -2.28817107408936 | 0.262915522993251 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 12 | 25 | 48 | -0.586985756164467 | 17.2404170486569 | -2.86233409356227 | 0.0886986743969436 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 12 | 50 | 98 | -0.830123217287616 | 17.2404170486569 | -2.86233409356227 | 0.130271349927137 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 12 | 100 | 198 | -1.17397151232893 | 17.2404170486569 | -2.86233409356227 | 0.215152341994119 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_Peters | 12 | 200 | 398 | -1.66024643457523 | 17.2404170486569 | -2.86233409356227 | 0.380818715586047 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 0 | 25 | 48 | 0 | 8.15872706186014 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 0 | 50 | 98 | 0 | 8.15872706186014 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 0 | 100 | 198 | 0 | 8.15872706186014 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 0 | 200 | 398 | 0 | 8.15872706186014 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 1 | 25 | 48 | 0.0474165976770739 | 8.15872706186014 | 0.109420271150594 | 0.0502474641619519 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 1 | 50 | 98 | 0.0670571955165066 | 8.15872706186014 | 0.109420271150594 | 0.0505052436500125 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 1 | 100 | 198 | 0.0948331953541479 | 8.15872706186014 | 0.109420271150594 | 0.0510208738364377 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 1 | 200 | 398 | 0.134114391033013 | 8.15872706186014 | 0.109420271150594 | 0.0520530133249603 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 2 | 25 | 48 | 0.0810634162930657 | 8.15872706186014 | 0.187064897645943 | 0.0507234809240671 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 2 | 50 | 98 | 0.11464098273395 | 8.15872706186014 | 0.187064897645943 | 0.051477563760721 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 2 | 100 | 198 | 0.162126832586131 | 8.15872706186014 | 0.187064897645943 | 0.0529872782255009 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 2 | 200 | 398 | 0.229281965467899 | 8.15872706186014 | 0.187064897645943 | 0.0560145244022503 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 3 | 25 | 48 | 0.098647480763292 | 8.15872706186014 | 0.227642526504202 | 0.051071620331205 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 3 | 50 | 98 | 0.139508605189387 | 8.15872706186014 | 0.227642526504202 | 0.0521890453815843 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 3 | 100 | 198 | 0.197294961526584 | 8.15872706186014 | 0.227642526504202 | 0.0544276143973131 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 3 | 200 | 398 | 0.279017210378773 | 8.15872706186014 | 0.227642526504202 | 0.0589218021397372 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 4 | 25 | 48 | 0.0989704661531477 | 8.15872706186014 | 0.228387859376383 | 0.0510786536762272 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 4 | 50 | 98 | 0.139965375508169 | 8.15872706186014 | 0.228387859376383 | 0.0522034223679787 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 4 | 100 | 198 | 0.197940932306295 | 8.15872706186014 | 0.228387859376383 | 0.0544567318102879 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 4 | 200 | 398 | 0.279930751016338 | 8.15872706186014 | 0.228387859376383 | 0.0589806217044097 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 6 | 25 | 48 | 0.0489229700314134 | 8.15872706186014 | 0.112896430966786 | 0.0502634398027136 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 6 | 50 | 98 | 0.0691875277299973 | 8.15872706186014 | 0.112896430966786 | 0.0505378663169399 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 6 | 100 | 198 | 0.0978459400628267 | 8.15872706186014 | 0.112896430966786 | 0.0510868117091086 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 6 | 200 | 398 | 0.138375055459995 | 8.15872706186014 | 0.112896430966786 | 0.0521857038390387 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 8 | 25 | 48 | -0.0556688859031189 | 8.15872706186014 | -0.12846355260778 | 0.0503411154352735 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 8 | 50 | 98 | -0.0787276934463911 | 8.15872706186014 | -0.12846355260778 | 0.0506964913400638 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 8 | 100 | 198 | -0.111337771806238 | 8.15872706186014 | -0.12846355260778 | 0.0514074667578257 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 8 | 200 | 398 | -0.157455386892782 | 8.15872706186014 | -0.12846355260778 | 0.05283112466775 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 12 | 25 | 48 | -0.309278704112533 | 8.15872706186014 | -0.713702823968326 | 0.060592086531255 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 12 | 50 | 98 | -0.43738613790912 | 8.15872706186014 | -0.713702823968326 | 0.0717531246304581 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 12 | 100 | 198 | -0.618557408225066 | 8.15872706186014 | -0.713702823968326 | 0.0944135925968895 | 0.05 | two-sample t, equal variance | post hoc |
| free-living | RA_GTEx | 12 | 200 | 398 | -0.87477227581824 | 8.15872706186014 | -0.713702823968326 | 0.140760221108298 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 0 | 25 | 48 | 0 | 0.0843795698011575 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 0 | 50 | 98 | 0 | 0.0843795698011575 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 0 | 100 | 198 | 0 | 0.0843795698011575 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 0 | 200 | 398 | 0 | 0.0843795698011575 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 1 | 25 | 48 | 0.477633281860695 | 0.0843795698011575 | 0.0113992658304002 | 0.0754623127529894 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 1 | 50 | 98 | 0.675475465048166 | 0.0843795698011575 | 0.0113992658304002 | 0.102620483357812 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 1 | 100 | 198 | 0.95526656372139 | 0.0843795698011575 | 0.0113992658304002 | 0.158212945289015 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 1 | 200 | 398 | 1.35095093009633 | 0.0843795698011575 | 0.0113992658304002 | 0.270649308997418 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 2 | 25 | 48 | 0.872572133247337 | 0.0843795698011575 | 0.020824934276642 | 0.137070233502759 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 2 | 50 | 98 | 1.23400334498721 | 0.0843795698011575 | 0.020824934276642 | 0.230978127687953 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 2 | 100 | 198 | 1.74514426649467 | 0.0843795698011575 | 0.020824934276642 | 0.411766108245685 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 2 | 200 | 398 | 2.46800668997442 | 0.0843795698011575 | 0.020824934276642 | 0.692200557700429 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 3 | 25 | 48 | 1.15790212410524 | 0.0843795698011575 | 0.0276346616107675 | 0.20563010881495 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 3 | 50 | 98 | 1.63752088781024 | 0.0843795698011575 | 0.0276346616107675 | 0.367669140909428 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 3 | 100 | 198 | 2.31580424821048 | 0.0843795698011575 | 0.0276346616107675 | 0.634811132674576 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 3 | 200 | 398 | 3.27504177562048 | 0.0843795698011575 | 0.0276346616107675 | 0.904421736915134 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 4 | 25 | 48 | 1.31417848708754 | 0.0843795698011575 | 0.0313643761685627 | 0.251370921396815 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 4 | 50 | 98 | 1.85852903981815 | 0.0843795698011575 | 0.0313643761685627 | 0.452449784440617 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 4 | 100 | 198 | 2.62835697417507 | 0.0843795698011575 | 0.0313643761685627 | 0.743969698644312 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 4 | 200 | 398 | 3.7170580796363 | 0.0843795698011575 | 0.0313643761685627 | 0.959777520201215 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 6 | 25 | 48 | 1.20649099563666 | 0.0843795698011575 | 0.0287942907321473 | 0.219261087468037 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 6 | 50 | 98 | 1.70623592891038 | 0.0843795698011575 | 0.0287942907321473 | 0.393539023642497 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 6 | 100 | 198 | 2.41298199127332 | 0.0843795698011575 | 0.0287942907321473 | 0.670508740582778 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 6 | 200 | 398 | 3.41247185782077 | 0.0843795698011575 | 0.0287942907321473 | 0.925667403047117 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 8 | 25 | 48 | 0.57836443526391 | 0.0843795698011575 | 0.0138033302845623 | 0.0875507220514583 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 8 | 50 | 98 | 0.817930828344477 | 0.0843795698011575 | 0.0138033302845623 | 0.127870083481858 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 8 | 100 | 198 | 1.15672887052782 | 0.0843795698011575 | 0.0138033302845623 | 0.210235341196208 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 8 | 200 | 398 | 1.63586165668895 | 0.0843795698011575 | 0.0138033302845623 | 0.371601283980646 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 12 | 25 | 48 | -1.47162810364725 | 0.0843795698011575 | -0.0351220917680009 | 0.302658114344387 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 12 | 50 | 98 | -2.08119642294734 | 0.0843795698011575 | -0.0351220917680009 | 0.540179382843203 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 12 | 100 | 198 | -2.9432562072945 | 0.0843795698011575 | -0.0351220917680009 | 0.833724741017844 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_chrono | 12 | 200 | 398 | -4.16239284589468 | 0.0843795698011575 | -0.0351220917680009 | 0.985823411804794 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 0 | 25 | 48 | 0 | 0.387677419351225 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 0 | 50 | 98 | 0 | 0.387677419351225 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 0 | 100 | 198 | 0 | 0.387677419351225 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 0 | 200 | 398 | 0 | 0.387677419351225 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 1 | 25 | 48 | 0.341953392682007 | 0.387677419351225 | 0.037495782063044 | 0.0629655025265563 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 1 | 50 | 98 | 0.483595125630388 | 0.387677419351225 | 0.037495782063044 | 0.0766584226756659 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 1 | 100 | 198 | 0.683906785364015 | 0.387677419351225 | 0.037495782063044 | 0.104521976388782 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 1 | 200 | 398 | 0.967190251260775 | 0.387677419351225 | 0.037495782063044 | 0.161565565407718 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 2 | 25 | 48 | 0.599035777209249 | 0.387677419351225 | 0.065685311012813 | 0.0903333730080287 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 2 | 50 | 98 | 0.847164520476028 | 0.387677419351225 | 0.065685311012813 | 0.133691282961963 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 2 | 100 | 198 | 1.1980715544185 | 0.387677419351225 | 0.065685311012813 | 0.222143187749677 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 2 | 200 | 398 | 1.69432904095206 | 0.387677419351225 | 0.065685311012813 | 0.393817565550933 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 3 | 25 | 48 | 0.753727413924923 | 0.387677419351225 | 0.0826475170367793 | 0.114480717660514 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 3 | 50 | 98 | 1.06593153110503 | 0.387677419351225 | 0.0826475170367793 | 0.184159443357939 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 3 | 100 | 198 | 1.50745482784985 | 0.387677419351225 | 0.0826475170367793 | 0.323093188979162 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 3 | 200 | 398 | 2.13186306221005 | 0.387677419351225 | 0.0826475170367793 | 0.566239711257506 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 4 | 25 | 48 | 0.795486323426871 | 0.387677419351225 | 0.0872264538257805 | 0.122015671335646 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 4 | 50 | 98 | 1.12498754727259 | 0.387677419351225 | 0.0872264538257805 | 0.199838025383745 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 4 | 100 | 198 | 1.59097264685374 | 0.387677419351225 | 0.0872264538257805 | 0.353393392127528 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 4 | 200 | 398 | 2.24997509454518 | 0.387677419351225 | 0.0872264538257805 | 0.612029918163356 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 6 | 25 | 48 | 0.536712873441199 | 0.387677419351225 | 0.0588514966181245 | 0.0822575167961802 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 6 | 50 | 98 | 0.759026624720779 | 0.387677419351225 | 0.0588514966181245 | 0.116803757605808 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 6 | 100 | 198 | 1.0734257468824 | 0.387677419351225 | 0.0588514966181245 | 0.187496373525081 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 6 | 200 | 398 | 1.51805324944156 | 0.387677419351225 | 0.0588514966181245 | 0.32820880449204 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 8 | 25 | 48 | -0.107946435801489 | 0.387677419351225 | -0.011836513741096 | 0.0512833383912046 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 8 | 50 | 98 | -0.152659313520302 | 0.387677419351225 | -0.011836513741096 | 0.0526218765382301 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 8 | 100 | 198 | -0.215892871602977 | 0.387677419351225 | -0.011836513741096 | 0.0553044307077604 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 8 | 200 | 398 | -0.305318627040604 | 0.387677419351225 | -0.011836513741096 | 0.060693837148047 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 12 | 25 | 48 | -1.80686551845672 | 0.387677419351225 | -0.198125935133753 | 0.425037727825456 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 12 | 50 | 98 | -2.55529372158578 | 0.387677419351225 | -0.198125935133753 | 0.715739009582134 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 12 | 100 | 198 | -3.61373103691344 | 0.387677419351225 | -0.198125935133753 | 0.949099487813791 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | TAms_mort | 12 | 200 | 398 | -5.11058744317157 | 0.387677419351225 | -0.198125935133753 | 0.999150221648051 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 0 | 25 | 48 | 0 | 7.86497441621665 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 0 | 50 | 98 | 0 | 7.86497441621665 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 0 | 100 | 198 | 0 | 7.86497441621665 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 0 | 200 | 398 | 0 | 7.86497441621665 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 1 | 25 | 48 | 0.240718461576559 | 7.86497441621665 | 0.535490421583496 | 0.056401184589331 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 1 | 50 | 98 | 0.340427313075156 | 7.86497441621665 | 0.535490421583496 | 0.0631170869674865 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 1 | 100 | 198 | 0.481436923153118 | 7.86497441621665 | 0.535490421583496 | 0.0766847457787555 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 1 | 200 | 398 | 0.680854626150313 | 7.86497441621665 | 0.535490421583496 | 0.104295945551753 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 2 | 25 | 48 | 0.401905924002446 | 7.86497441621665 | 0.894060103539363 | 0.0679588229630331 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 2 | 50 | 98 | 0.56838080852235 | 7.86497441621665 | 0.894060103539363 | 0.0870078829229117 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 2 | 100 | 198 | 0.803811848004893 | 7.86497441621665 | 0.894060103539363 | 0.125906224104752 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 2 | 200 | 398 | 1.1367616170447 | 7.86497441621665 | 0.894060103539363 | 0.205405325235684 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 3 | 25 | 48 | 0.472577728088214 | 7.86497441621665 | 1.05127311460673 | 0.0749193625455018 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 3 | 50 | 98 | 0.668325832337817 | 7.86497441621665 | 1.05127311460673 | 0.101488775235673 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 3 | 100 | 198 | 0.945155456176429 | 7.86497441621665 | 1.05127311460673 | 0.155872109425963 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 3 | 200 | 398 | 1.33665166467563 | 7.86497441621665 | 1.05127311460673 | 0.265976409673589 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 4 | 25 | 48 | 0.447917707176096 | 7.86497441621665 | 0.99641564788813 | 0.0723574393098166 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 4 | 50 | 98 | 0.633451296315496 | 7.86497441621665 | 0.99641564788813 | 0.0961526258009999 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 4 | 100 | 198 | 0.895835414352193 | 7.86497441621665 | 0.99641564788813 | 0.144830603170745 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 4 | 200 | 398 | 1.26690259263099 | 7.86497441621665 | 0.99641564788813 | 0.243810146585272 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 6 | 25 | 48 | 0.125706529377255 | 7.86497441621665 | 0.2796405475972 | 0.0517408451067976 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 6 | 50 | 98 | 0.177775878724166 | 7.86497441621665 | 0.2796405475972 | 0.0535575741917017 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 6 | 100 | 198 | 0.25141305875451 | 7.86497441621665 | 0.2796405475972 | 0.0572014184269416 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 6 | 200 | 398 | 0.355551757448332 | 7.86497441621665 | 0.2796405475972 | 0.0645330782037177 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 8 | 25 | 48 | -0.478391384510598 | 7.86497441621665 | -1.06420588805574 | 0.0755442519177847 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 8 | 50 | 98 | -0.676547584097329 | 7.86497441621665 | -1.06420588805574 | 0.102791297945255 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 8 | 100 | 198 | -0.956782769021195 | 7.86497441621665 | -1.06420588805574 | 0.158566223607093 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 8 | 200 | 398 | -1.35309516819466 | 7.86497441621665 | -1.06420588805574 | 0.271353697023813 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 12 | 25 | 48 | -1.85261818337339 | 7.86497441621665 | -4.12124307188773 | 0.442656297988794 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 12 | 50 | 98 | -2.61999776082565 | 7.86497441621665 | -4.12124307188773 | 0.737054171753636 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 12 | 100 | 198 | -3.70523636674678 | 7.86497441621665 | -4.12124307188773 | 0.957937862950945 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_default | 12 | 200 | 398 | -5.2399955216513 | 7.86497441621665 | -4.12124307188773 | 0.999457204668619 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 0 | 25 | 48 | 0 | 17.2404170486569 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 0 | 50 | 98 | 0 | 17.2404170486569 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 0 | 100 | 198 | 0 | 17.2404170486569 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 0 | 200 | 398 | 0 | 17.2404170486569 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 1 | 25 | 48 | -0.0371488516204322 | 17.2404170486569 | -0.181149922997545 | 0.0501518857456207 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 1 | 50 | 98 | -0.0525364097882009 | 17.2404170486569 | -0.181149922997545 | 0.0503100838776626 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 1 | 100 | 198 | -0.0742977032408643 | 17.2404170486569 | -0.181149922997545 | 0.0506264665946228 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 1 | 200 | 398 | -0.105072819576402 | 17.2404170486569 | -0.181149922997545 | 0.0512595446110096 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 2 | 25 | 48 | -0.11156797632957 | 17.2404170486569 | -0.544041859694461 | 0.0513709659824255 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 2 | 50 | 98 | -0.157780945251798 | 17.2404170486569 | -0.544041859694461 | 0.0528010530917548 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 2 | 100 | 198 | -0.22313595265914 | 17.2404170486569 | -0.544041859694461 | 0.0556675289045848 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 2 | 200 | 398 | -0.315561890503596 | 17.2404170486569 | -0.544041859694461 | 0.0614281287007591 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 3 | 25 | 48 | -0.218185833761904 | 17.2404170486569 | -1.06394532431213 | 0.055255382853614 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 3 | 50 | 98 | -0.308561365223766 | 17.2404170486569 | -1.06394532431213 | 0.0607622160383058 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 3 | 100 | 198 | -0.436371667523808 | 17.2404170486569 | -1.06394532431213 | 0.07186946903625 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 3 | 200 | 398 | -0.617122730447532 | 17.2404170486569 | -1.06394532431213 | 0.0944243138207476 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 4 | 25 | 48 | -0.349736593127761 | 17.2404170486569 | -1.70542975491797 | 0.0635669258356088 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 4 | 50 | 98 | -0.49460223325944 | 17.2404170486569 | -1.70542975491797 | 0.0779029419210641 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 4 | 100 | 198 | -0.699473186255521 | 17.2404170486569 | -1.70542975491797 | 0.107089928918182 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 4 | 200 | 398 | -0.98920446651888 | 17.2404170486569 | -1.70542975491797 | 0.166846075415291 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 6 | 25 | 48 | -0.65068876186106 | 17.2404170486569 | -3.1729707370461 | 0.0977415154199243 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 6 | 50 | 98 | -0.920212871907667 | 17.2404170486569 | -3.1729707370461 | 0.149192027236919 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 6 | 100 | 198 | -1.30137752372212 | 17.2404170486569 | -3.1729707370461 | 0.253625843606511 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 6 | 200 | 398 | -1.84042574381533 | 17.2404170486569 | -3.1729707370461 | 0.45073921546197 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 8 | 25 | 48 | -0.933784591956984 | 17.2404170486569 | -4.55343838505805 | 0.150084224028432 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 8 | 50 | 98 | -1.32057083428059 | 17.2404170486569 | -4.55343838505805 | 0.257640263984831 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 8 | 100 | 198 | -1.86756918391397 | 17.2404170486569 | -4.55343838505805 | 0.459661511532151 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 8 | 200 | 398 | -2.64114166856119 | 17.2404170486569 | -4.55343838505805 | 0.750100655592225 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 12 | 25 | 48 | -1.16809599765845 | 17.2404170486569 | -5.69601726028016 | 0.208444393244674 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 12 | 50 | 98 | -1.65193720204231 | 17.2404170486569 | -5.69601726028016 | 0.373050439047628 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 12 | 100 | 198 | -2.33619199531689 | 17.2404170486569 | -5.69601726028016 | 0.642410459387534 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_Peters | 12 | 200 | 398 | -3.30387440408461 | 17.2404170486569 | -5.69601726028016 | 0.909213624270223 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 0 | 25 | 48 | 0 | 8.15872706186014 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 0 | 50 | 98 | 0 | 8.15872706186014 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 0 | 100 | 198 | 0 | 8.15872706186014 | 0 | 0.0500000000000002 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 0 | 200 | 398 | 0 | 8.15872706186014 | 0 | 0.05 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 1 | 25 | 48 | 0.180201673529321 | 8.15872706186014 | 0.415839957848824 | 0.0535812947611972 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 1 | 50 | 98 | 0.254843650667494 | 8.15872706186014 | 0.415839957848824 | 0.0573268233209833 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 1 | 100 | 198 | 0.360403347058642 | 8.15872706186014 | 0.415839957848824 | 0.0648623900252189 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 1 | 200 | 398 | 0.509687301334989 | 8.15872706186014 | 0.415839957848824 | 0.080105076040522 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 2 | 25 | 48 | 0.308072784502574 | 8.15872706186014 | 0.710919999869342 | 0.0605091575766872 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 2 | 50 | 98 | 0.435680710041583 | 8.15872706186014 | 0.710919999869342 | 0.0715819095536957 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 2 | 100 | 198 | 0.616145569005147 | 8.15872706186014 | 0.710919999869342 | 0.0940611999813695 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 2 | 200 | 398 | 0.871361420083167 | 8.15872706186014 | 0.710919999869342 | 0.140034662761428 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 3 | 25 | 48 | 0.374899128023934 | 8.15872706186014 | 0.865130908840683 | 0.0656066901456799 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 3 | 50 | 98 | 0.530187431373294 | 8.15872706186014 | 0.865130908840683 | 0.0821280768945353 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 3 | 100 | 198 | 0.749798256047867 | 8.15872706186014 | 0.865130908840683 | 0.115816174775354 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 3 | 200 | 398 | 1.06037486274659 | 8.15872706186014 | 0.865130908840683 | 0.184762679376251 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 4 | 25 | 48 | 0.376126599218174 | 8.15872706186014 | 0.867963466162026 | 0.0657099306185587 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 4 | 50 | 98 | 0.531923337783612 | 8.15872706186014 | 0.867963466162026 | 0.0823420963769134 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 4 | 100 | 198 | 0.752253198436349 | 8.15872706186014 | 0.867963466162026 | 0.116258479110633 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 4 | 200 | 398 | 1.06384667556722 | 8.15872706186014 | 0.867963466162026 | 0.185669429687523 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 6 | 25 | 48 | 0.1859264794519 | 8.15872706186014 | 0.429050728908343 | 0.0538129806634243 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 6 | 50 | 98 | 0.26293974884516 | 8.15872706186014 | 0.429050728908343 | 0.0578018850601331 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 6 | 100 | 198 | 0.3718529589038 | 8.15872706186014 | 0.429050728908343 | 0.0658300282940916 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 6 | 200 | 398 | 0.525879497690319 | 8.15872706186014 | 0.429050728908343 | 0.0820786537922872 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 8 | 25 | 48 | -0.211563606304574 | 8.15872706186014 | -0.488211898396849 | 0.0549402951812097 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 8 | 50 | 98 | -0.29919612134049 | 8.15872706186014 | -0.488211898396849 | 0.0601151438929141 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 8 | 100 | 198 | -0.423127212609147 | 8.15872706186014 | -0.488211898396849 | 0.0705479908751206 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 8 | 200 | 398 | -0.59839224268098 | 8.15872706186014 | -0.488211898396849 | 0.0917192436421061 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 12 | 25 | 48 | -1.1753804110455 | 8.15872706186014 | -2.71235072911775 | 0.210470317817763 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 12 | 50 | 98 | -1.6622389182482 | 8.15872706186014 | -2.71235072911775 | 0.376911498859321 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 12 | 100 | 198 | -2.35076082209099 | 8.15872706186014 | -2.71235072911775 | 0.64780670309624 | 0.05 | two-sample t, equal variance | post hoc |
| free-living upper | RA_GTEx | 12 | 200 | 398 | -3.32447783649641 | 8.15872706186014 | -2.71235072911775 | 0.912527015390056 | 0.05 | two-sample t, equal variance | post hoc |

### 4.4 설계 3

lab. 사전·사후 시각은 각각 독립 Uniform(8, 8+W). δ = 0.3 × sqrt(2) × σ를 사후 표본에 더한다. n = 50. 2000회, 시드 20261019. `power` = n_reject / 2000. `bias` = mean(estimates) − δ. `results/stage9/sim_design3_power.tsv`.

| scenario | clock | W | n | delta | n_rep | n_reject | n_failed | power | mean_estimate | bias | seed | time_draw | effect | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| lab | TAms_chrono | 0 | 50 | 0.00342325913173446 | 2000 | 1120 | 0 | 0.56 | 0.00346809829972338 | 4.48391679889163e-05 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | TAms_mort | 0 | 50 | 0.0135988275094574 | 2000 | 1092 | 0 | 0.546 | 0.013540425547066 | -5.84019623913225e-05 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_default | 0 | 50 | 3.2508327938958 | 2000 | 1126 | 0 | 0.563 | 3.28951210614491 | 0.0386793122491103 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_Peters | 0 | 50 | 3.6615000484052 | 2000 | 1097 | 0 | 0.5485 | 3.67906047567822 | 0.0175604272730201 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_GTEx | 0 | 50 | 11.446527357556 | 2000 | 1138 | 0 | 0.569 | 11.699695272672 | 0.253167915116043 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | TAms_chrono | 2 | 50 | 0.00342325913173446 | 2000 | 1072 | 0 | 0.536 | 0.00343495172760903 | 1.16925958745752e-05 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | TAms_mort | 2 | 50 | 0.0135988275094574 | 2000 | 1076 | 0 | 0.538 | 0.0135255901963005 | -7.32373131568809e-05 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_default | 2 | 50 | 3.2508327938958 | 2000 | 1091 | 0 | 0.5455 | 3.23135312822552 | -0.0194796656702789 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_Peters | 2 | 50 | 3.6615000484052 | 2000 | 1092 | 0 | 0.546 | 3.61458208162861 | -0.0469179667765909 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_GTEx | 2 | 50 | 11.446527357556 | 2000 | 1104 | 0 | 0.552 | 11.5982867974976 | 0.151759439941641 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | TAms_chrono | 4 | 50 | 0.00342325913173446 | 2000 | 1041 | 0 | 0.5205 | 0.00337586315070141 | -4.73959810330495e-05 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | TAms_mort | 4 | 50 | 0.0135988275094574 | 2000 | 1065 | 0 | 0.5325 | 0.0133797359600323 | -0.000219091549425118 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_default | 4 | 50 | 3.2508327938958 | 2000 | 1072 | 0 | 0.536 | 3.21634535436132 | -0.0344874395344781 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_Peters | 4 | 50 | 3.6615000484052 | 2000 | 1069 | 0 | 0.5345 | 3.62836081353738 | -0.0331392348678246 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_GTEx | 4 | 50 | 11.446527357556 | 2000 | 1089 | 0 | 0.5445 | 11.5700792214264 | 0.123551863870455 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | TAms_chrono | 8 | 50 | 0.00342325913173446 | 2000 | 1048 | 0 | 0.524 | 0.00340600792475218 | -1.72512069822784e-05 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | TAms_mort | 8 | 50 | 0.0135988275094574 | 2000 | 1089 | 0 | 0.5445 | 0.0136335692839545 | 3.47417744971188e-05 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_default | 8 | 50 | 3.2508327938958 | 2000 | 1106 | 0 | 0.553 | 3.28196342198294 | 0.0311306280871482 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_Peters | 8 | 50 | 3.6615000484052 | 2000 | 1034 | 0 | 0.517 | 3.64157541834391 | -0.0199246300612885 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |
| lab | RA_GTEx | 8 | 50 | 11.446527357556 | 2000 | 1078 | 0 | 0.539 | 11.2763638696472 | -0.170163487908811 | 20261019 | independent Uniform(8, 8+W) | delta added to the post sample | post hoc |

## 5. 완화

### 5.1 코사이너 자료 구성 보정 표

MCPcounter 8개 점수를 참가자 안에서 중심화해 모형 M의 고정 효과로 넣었다. 보정 모형의 LRT는 cos와 sin을 빼고 구성 공변량은 남긴다. L1은 조성이 비유한인 샘플을 두 모형에서 함께 뺐다. `results/stage9/composition_cosinor.tsv`.

| n_subjects | n_samples | two_a_before | r2_within_before | lrt_p_before | peak_before | two_a_after | r2_within_after | lrt_p_after | peak_after | two_a_ratio | reml_converged_after | ml_converged_after | null_converged_after | dataset | clock | n_dropped | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 24 | 221 | 0.0345210353298852 | 0.362193176953973 | 5.49745798912395e-20 | 14.2335397784943 | 0.0132221784395503 | 0.0422755862812658 | 0.0138250307973468 | 14.1083433112827 | 0.383018015340455 | TRUE | TRUE | TRUE | GSE39445 | TAms_chrono | 0 | post hoc |
| 24 | 221 | 0.148755555495313 | 0.28699925739554 | 3.45851149483147e-15 | 13.4288594848278 | 0.0488314470574515 | 0.0214541628839459 | 0.1182067715791 | 13.2871444504737 | 0.328266375631196 | TRUE | TRUE | TRUE | GSE39445 | TAms_mort | 0 | post hoc |
| 24 | 221 | 1.98498899424706 | 0.00680181474094321 | 0.571779550797631 | 12.0064473932843 | 6.16223397454987 | 0.0339481999890606 | 0.0241636234971642 | 0.761657320658042 | 3.10441719949551 | TRUE | TRUE | TRUE | GSE39445 | RA_default | 0 | post hoc |
| 24 | 221 | 3.86882060477724 | 0.0237318052396572 | 0.0822398090507771 | 6.84503499428336 | 3.97630516530413 | 0.0292955331944914 | 0.0466761034079478 | 8.98987208955333 | 1.02778225498338 | TRUE | TRUE | TRUE | GSE39445 | RA_Peters | 0 | post hoc |
| 24 | 221 | 4.56221470455666 | 0.00480935568643748 | 0.643663244864561 | 13.7137336093387 | 9.54901600524505 | 0.00970435211021647 | 0.346688262360935 | 23.9149808479063 | 2.09306589532221 | TRUE | TRUE | TRUE | GSE39445 | RA_GTEx | 0 | post hoc |
| 14 | 130 | 0.00943066206760042 | 0.0286364463384637 | 0.181165188487945 | 11.4579688514277 | 0.00634984246734852 | 0.0118840365758605 | 0.492172715280185 | 5.14452008511751 | 0.673318842498213 | TRUE | TRUE | TRUE | GSE56931 | TAms_chrono | 0 | post hoc |
| 14 | 130 | 0.0629468359700225 | 0.0821416660164724 | 0.00673707816031792 | 12.0688456403715 | 0.0192213217427161 | 0.00856317090970282 | 0.597649222517841 | 8.27948252129182 | 0.305358028668351 | TRUE | TRUE | TRUE | GSE56931 | TAms_mort | 0 | post hoc |
| 14 | 130 | 4.61694071235948 | 0.0658778147622046 | 0.0187399469939367 | 9.88701727265618 | 2.99011335055806 | 0.0183751807714839 | 0.333505411967642 | 4.20023629641222 | 0.647639538136926 | TRUE | TRUE | TRUE | GSE56931 | RA_default | 0 | post hoc |
| 14 | 130 | 6.37997357324778 | 0.145491804603583 | 0.00010511653009091 | 4.54192615314065 | 3.96015059605792 | 0.0464277146029072 | 0.0614925932925041 | 4.98299949901876 | 0.620715830652254 | TRUE | TRUE | TRUE | GSE56931 | RA_Peters | 0 | post hoc |
| 14 | 130 | 8.4129561513483 | 0.0366339113390549 | 0.0977660921595171 | 13.2671292169983 | 8.23140216627795 | 0.0208364265492177 | 0.266428527393462 | 14.8114452451722 | 0.978419715756957 | TRUE | TRUE | TRUE | GSE56931 | RA_GTEx | 0 | post hoc |
| 11 | 153 | 0.0418702301338288 | 0.265399481367344 | 2.96936394371969e-10 | 12.188050742382 | 0.0150799665852315 | 0.0537996633334796 | 0.019570886047184 | 8.72787340418235 | 0.360159629814113 | TRUE | TRUE | TRUE | GSE113883 | TAms_chrono | 0 | post hoc |
| 11 | 153 | 0.102562103238955 | 0.158916853730887 | 4.44466652817908e-06 | 10.0600990696032 | 0.0430937000835042 | 0.0318094069732833 | 0.100054897837397 | 8.04872492523851 | 0.420171766398961 | TRUE | TRUE | TRUE | GSE113883 | TAms_mort | 0 | post hoc |
| 11 | 153 | 6.770165944749 | 0.134425171735901 | 3.24675832907807e-05 | 11.812631972311 | 4.46359984000769 | 0.0654036674780254 | 0.0078792275682327 | 6.92721808426326 | 0.659304347402251 | TRUE | TRUE | TRUE | GSE113883 | RA_default | 0 | post hoc |
| 11 | 153 | 3.1551856815567 | 0.0352688623430064 | 0.0773443137299791 | 11.8362208570625 | 1.72992981067056 | 0.00632996044783596 | 0.63271309500557 | 15.1906614130393 | 0.548281459561217 | TRUE | TRUE | TRUE | GSE113883 | RA_Peters | 0 | post hoc |
| 11 | 153 | 4.34369415553422 | 0.0315491728014815 | 0.104384139695734 | 6.16828794535857 | 5.67674067646209 | 0.0447939538479821 | 0.0392477906194442 | 6.78164274664445 | 1.30689235318961 | TRUE | TRUE | TRUE | GSE113883 | RA_GTEx | 0 | post hoc |
| 22 | 146 | 0.0120084580513382 | 0.242268321685822 | 2.41897966602237e-08 | 10.975433301276 | 0.00595043705460161 | 0.0833385451329599 | 0.00382091898160207 | 10.6645085901725 | 0.49552049306934 | TRUE | TRUE | TRUE | L1 | TAms_chrono | 1 | post hoc |
| 22 | 146 | 0.0361347645463202 | 0.154487919960743 | 3.24731259384485e-05 | 11.4080740851573 | 0.0121361812197449 | 0.0324896199329729 | 0.124248747645643 | 4.30504825720064 | 0.335858870871784 | TRUE | TRUE | TRUE | L1 | TAms_mort | 1 | post hoc |
| 22 | 146 | 5.53144269979541 | 0.0675877700825838 | 0.0123776996539389 | 11.6484216527209 | 1.02719888607263 | 0.00181318387756768 | 0.880526515294065 | 3.88877116410112 | 0.185701803638067 | TRUE | TRUE | TRUE | L1 | RA_default | 1 | post hoc |
| 22 | 146 | 6.14429124215941 | 0.0633669376082697 | 0.0177414850485553 | 3.52037349612089 | 4.6370661312965 | 0.0319030439610186 | 0.136552737089083 | 3.55144198542743 | 0.754695041061693 | TRUE | TRUE | TRUE | L1 | RA_Peters | 1 | post hoc |
| 22 | 146 | 11.9030037723114 | 0.0274636847418713 | 0.193593128358903 | 16.2006789645578 | 11.9402218377118 | 0.0242788505428779 | 0.185300814470306 | 21.7030089077681 | 1.00312677926618 | TRUE | TRUE | TRUE | L1 | RA_GTEx | 1 | post hoc |

### 5.2 L3 구성 보정 표

Δclock을 Δcomposition 8개와 절편에 회귀했다. 절편이 구성 보정 mean Δ다. 시각 열은 H5의 mean Δ − c다. `results/stage9/composition_l3.tsv`.

| clock | n | mean_delta_before | ci_low_before | ci_high_before | p_before | mean_delta_composition | ci_low_composition | ci_high_composition | p_composition | df_residual | mean_delta_time | ci_low_time | ci_high_time | time_source | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TAms_chrono | 58 | 0.0120385849794561 | -0.0110449218090888 | 0.0351220917680009 | 0.300740570646651 | 0.0166461557079569 | -0.0107039008919362 | 0.04399621230785 | 0.227143866381081 | 49 | -0.00708671521880421 | -0.0301702220073491 | 0.0159967915697406 | H5 mean_delta_minus_c | post hoc |
| TAms_mort | 58 | 0.0827171714412529 | -0.0326915922512471 | 0.198125935133753 | 0.156685969793716 | 0.0828733287438787 | -0.0489965626069634 | 0.214743220094721 | 0.212597308464689 | 49 | -0.00528640391088912 | -0.120695167603389 | 0.110122359781611 | H5 mean_delta_minus_c | post hoc |
| RA_default | 58 | 2.28347265426549 | 0.445702236643247 | 4.12124307188773 | 0.0157896851797883 | 1.81652100332312 | -0.385596454242882 | 4.01863846088912 | 0.103768923991725 | 49 | -1.38502667791645 | -3.2227970955387 | 0.452743739705789 | H5 mean_delta_minus_c | post hoc |
| RA_Peters | 58 | 2.86233409356227 | 0.0286509268443877 | 5.69601726028016 | 0.047803206296958 | 2.02172266582033 | -1.46068313169132 | 5.50412846333198 | 0.248992062313986 | 49 |  |  |  | H5 c is not defined for this clock | post hoc |
| RA_GTEx | 58 | 0.713702823968326 | -1.2849450811811 | 2.71235072911775 | 0.477486960820646 | 0.159401851182051 | -2.201328211144 | 2.5201319135081 | 0.89262141551924 | 49 |  |  |  | H5 c is not defined for this clock | post hoc |

3.2 회귀 계수 (`results/stage9/composition_l3_betas.tsv`):

| clock | term | beta | title |
| --- | --- | --- | --- |
| TAms_chrono | (Intercept) | 0.0166461557079569 | post hoc |
| TAms_chrono | T_cells | -0.0374183026012343 | post hoc |
| TAms_chrono | CD8_T_cells | -0.0469392673757116 | post hoc |
| TAms_chrono | Cytotoxic_lymphocytes | 0.00882508030585047 | post hoc |
| TAms_chrono | NK_cells | -0.0454040803928078 | post hoc |
| TAms_chrono | B_lineage | 0.0413805996968943 | post hoc |
| TAms_chrono | Monocytic_lineage | 0.0576095055931148 | post hoc |
| TAms_chrono | Myeloid_dendritic_cells | 0.135332514659117 | post hoc |
| TAms_chrono | Neutrophils | -0.0747162683642148 | post hoc |
| TAms_mort | (Intercept) | 0.0828733287438787 | post hoc |
| TAms_mort | T_cells | -0.17931309104282 | post hoc |
| TAms_mort | CD8_T_cells | -0.482022433555404 | post hoc |
| TAms_mort | Cytotoxic_lymphocytes | 0.473649587738716 | post hoc |
| TAms_mort | NK_cells | -0.514736150034478 | post hoc |
| TAms_mort | B_lineage | 0.216123068469115 | post hoc |
| TAms_mort | Monocytic_lineage | -0.112613490139108 | post hoc |
| TAms_mort | Myeloid_dendritic_cells | 0.558939320380493 | post hoc |
| TAms_mort | Neutrophils | -0.397244991291884 | post hoc |
| RA_default | (Intercept) | 1.81652100332312 | post hoc |
| RA_default | T_cells | -0.0202761719337729 | post hoc |
| RA_default | CD8_T_cells | 1.25913415601617 | post hoc |
| RA_default | Cytotoxic_lymphocytes | -9.42340086079396 | post hoc |
| RA_default | NK_cells | 8.9218013508008 | post hoc |
| RA_default | B_lineage | -1.03602913428704 | post hoc |
| RA_default | Monocytic_lineage | -5.73134156082354 | post hoc |
| RA_default | Myeloid_dendritic_cells | -1.38637398888501 | post hoc |
| RA_default | Neutrophils | -5.44541420196543 | post hoc |
| RA_Peters | (Intercept) | 2.02172266582033 | post hoc |
| RA_Peters | T_cells | 0.0540740482654746 | post hoc |
| RA_Peters | CD8_T_cells | -1.76783983645659 | post hoc |
| RA_Peters | Cytotoxic_lymphocytes | -15.7994115877294 | post hoc |
| RA_Peters | NK_cells | 11.0289989651614 | post hoc |
| RA_Peters | B_lineage | -0.811500183308397 | post hoc |
| RA_Peters | Monocytic_lineage | -2.1037673554493 | post hoc |
| RA_Peters | Myeloid_dendritic_cells | -9.22906524309916 | post hoc |
| RA_Peters | Neutrophils | -8.37920631736404 | post hoc |
| RA_GTEx | (Intercept) | 0.159401851182051 | post hoc |
| RA_GTEx | T_cells | -0.801161720563749 | post hoc |
| RA_GTEx | CD8_T_cells | -3.14307051844255 | post hoc |
| RA_GTEx | Cytotoxic_lymphocytes | -5.99788027543367 | post hoc |
| RA_GTEx | NK_cells | 2.28379421460835 | post hoc |
| RA_GTEx | B_lineage | -0.644616817478975 | post hoc |
| RA_GTEx | Monocytic_lineage | -1.47199497026536 | post hoc |
| RA_GTEx | Myeloid_dendritic_cells | -8.63485863380512 | post hoc |
| RA_GTEx | Neutrophils | 1.69949376947723 | post hoc |

### 5.3 나이 연관 표

아침 보정 = 점수 − β̂ᵀ(Δcomposition / 2), 저녁 보정 = 점수 + β̂ᵀ(Δcomposition / 2). `results/stage9/age_correlation_l3.tsv`.

| clock | morning_r_before | morning_r_after | evening_r_before | evening_r_after | person_mean_r | n | max_abs_person_mean_change | title |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TAms_chrono | 0.271730437294783 | 0.25049954468322 | 0.344473699826594 | 0.379420830334003 | 0.357659887896734 | 58 | 2.77555756156289e-17 | post hoc |
| TAms_mort | 0.459242940208403 | 0.454368246593447 | 0.462942410027015 | 0.504173016355484 | 0.542417587627288 | 58 | 1.11022302462516e-16 | post hoc |
| RA_default | 0.446171503104569 | 0.484858641020211 | 0.517349079713589 | 0.499271914644175 | 0.536867692830739 | 58 | 7.105427357601e-15 | post hoc |
| RA_Peters | 0.716354841847154 | 0.731704823753873 | 0.709903447516876 | 0.70297198458081 | 0.740682571652942 | 58 | 1.4210854715202e-14 | post hoc |
| RA_GTEx | 0.295924783031412 | 0.306451158604095 | 0.482587756628502 | 0.495235502983008 | 0.438835111423925 | 58 | 7.105427357601e-15 | post hoc |

### 5.4 요약 표

코사이너 행의 효과는 2A다. L3 행의 효과는 mean Δ다. `results/stage9/mitigation_summary.tsv`.

| dataset | clock | effect | effect_before | effect_composition | effect_time | two_a_ratio | title |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE39445 | TAms_chrono | two_a | 0.0345210353298852 | 0.0132221784395503 |  | 0.383018015340455 | post hoc |
| GSE39445 | TAms_mort | two_a | 0.148755555495313 | 0.0488314470574515 |  | 0.328266375631196 | post hoc |
| GSE39445 | RA_default | two_a | 1.98498899424706 | 6.16223397454987 |  | 3.10441719949551 | post hoc |
| GSE39445 | RA_Peters | two_a | 3.86882060477724 | 3.97630516530413 |  | 1.02778225498338 | post hoc |
| GSE39445 | RA_GTEx | two_a | 4.56221470455666 | 9.54901600524505 |  | 2.09306589532221 | post hoc |
| GSE56931 | TAms_chrono | two_a | 0.00943066206760042 | 0.00634984246734852 |  | 0.673318842498213 | post hoc |
| GSE56931 | TAms_mort | two_a | 0.0629468359700225 | 0.0192213217427161 |  | 0.305358028668351 | post hoc |
| GSE56931 | RA_default | two_a | 4.61694071235948 | 2.99011335055806 |  | 0.647639538136926 | post hoc |
| GSE56931 | RA_Peters | two_a | 6.37997357324778 | 3.96015059605792 |  | 0.620715830652254 | post hoc |
| GSE56931 | RA_GTEx | two_a | 8.4129561513483 | 8.23140216627795 |  | 0.978419715756957 | post hoc |
| GSE113883 | TAms_chrono | two_a | 0.0418702301338288 | 0.0150799665852315 |  | 0.360159629814113 | post hoc |
| GSE113883 | TAms_mort | two_a | 0.102562103238955 | 0.0430937000835042 |  | 0.420171766398961 | post hoc |
| GSE113883 | RA_default | two_a | 6.770165944749 | 4.46359984000769 |  | 0.659304347402251 | post hoc |
| GSE113883 | RA_Peters | two_a | 3.1551856815567 | 1.72992981067056 |  | 0.548281459561217 | post hoc |
| GSE113883 | RA_GTEx | two_a | 4.34369415553422 | 5.67674067646209 |  | 1.30689235318961 | post hoc |
| L1 | TAms_chrono | two_a | 0.0120084580513382 | 0.00595043705460161 |  | 0.49552049306934 | post hoc |
| L1 | TAms_mort | two_a | 0.0361347645463202 | 0.0121361812197449 |  | 0.335858870871784 | post hoc |
| L1 | RA_default | two_a | 5.53144269979541 | 1.02719888607263 |  | 0.185701803638067 | post hoc |
| L1 | RA_Peters | two_a | 6.14429124215941 | 4.6370661312965 |  | 0.754695041061693 | post hoc |
| L1 | RA_GTEx | two_a | 11.9030037723114 | 11.9402218377118 |  | 1.00312677926618 | post hoc |
| L3 | TAms_chrono | mean_delta | 0.0120385849794561 | 0.0166461557079569 | -0.00708671521880421 |  | post hoc |
| L3 | TAms_mort | mean_delta | 0.0827171714412529 | 0.0828733287438787 | -0.00528640391088912 |  | post hoc |
| L3 | RA_default | mean_delta | 2.28347265426549 | 1.81652100332312 | -1.38502667791645 |  | post hoc |
| L3 | RA_Peters | mean_delta | 2.86233409356227 | 2.02172266582033 |  |  | post hoc |
| L3 | RA_GTEx | mean_delta | 0.713702823968326 | 0.159401851182051 |  |  | post hoc |

## 6. 그림 (파일 목록, 원자료와 n), docx 확인 결과

- `figures/Fig6_simulation.pdf`, `figures/Fig6_simulation.png`. PNG 3600×4200, info dpi `(299.9994, 299.9994)`.
- `figures/Fig7_mitigation.pdf`, `figures/Fig7_mitigation.png`. PNG 3600×2400, info dpi `(299.9994, 299.9994)`.

`figures/figure_sources.tsv`의 Fig 6·Fig 7 행:

| figure | panel | n_subjects | n_samples | role | source |
| --- | --- | --- | --- | --- | --- |
| Fig6 | A b(s); B false positive rate at s=2 and 4; C design 3 power | design 1 n 20-200; design 2 n per group 25-200; design 3 n 50 | design 3 Monte Carlo 2000 | post hoc | results/stage9/sim_design1_bias.tsv; results/stage9/sim_design1_fpr.tsv; results/stage9/sim_design2_fpr.tsv; results/stage9/sim_design3_power.tsv |
| Fig7 | A cosinor 2A before and after; B L3 mean delta | GSE39445 24/221; GSE56931 14/130; GSE113883 11/153; L1 22/146; L3 58 | L3 116; L1 drops non-finite composition | post hoc | results/stage9/composition_cosinor.tsv; results/stage9/composition_l3.tsv; results/confirm/h5.tsv; results/confirm/family_P.tsv; results/confirm/family_S1.tsv |

지시서 8 docx 확인. `/usr/bin/soffice` (LibreOffice 7.3) `--headless --norestore --convert-to pdf`. 변환 종료 상태 0. 경고: `failed to launch javaldx - java may not function correctly`. `results/stage9/docx_check.tsv`.

| file | converter | exit_status | pages | note | title |
| --- | --- | --- | --- | --- | --- |
| results/stage8/tables/Table1.docx | LibreOffice 7.3 soffice --headless --convert-to pdf | 0 | 2 | PDF written. pdftotext shows the title and dataset ids. Cell text wraps across lines. Warning: failed to launch javaldx. | post hoc |
| results/stage8/tables/Table2.docx | LibreOffice 7.3 soffice --headless --convert-to pdf | 0 | 32 | PDF written. Warning: failed to launch javaldx. | post hoc |
| results/stage8/tables/Table3.docx | LibreOffice 7.3 soffice --headless --convert-to pdf | 0 | 34 | PDF written. Warning: failed to launch javaldx. | post hoc |
| results/stage8/tables/TableS1.docx | LibreOffice 7.3 soffice --headless --convert-to pdf | 0 | 7 | PDF written. Warning: failed to launch javaldx. | post hoc |
| results/stage8/tables/TableS2.docx | LibreOffice 7.3 soffice --headless --convert-to pdf | 0 | 36 | PDF written. Warning: failed to launch javaldx. | post hoc |
| results/stage8/tables/TableS3.docx | LibreOffice 7.3 soffice --headless --convert-to pdf | 0 | 2 | PDF written. Warning: failed to launch javaldx. | post hoc |

## 7. 판단이 필요한 사항

- 설계 2의 ncp 식은 지시서에 없다. 지시한 같은 분산 두 표본 t에 맞춰 ncp = b / total_sd × sqrt(n/2), df = 2n − 2로 계산했다.
- `prereg_gate_pass`는 `in_prereg_window`와 `lrt_p_below_0.01`을 함께 요구한다. 창 판정 열은 `in_prereg_window`다.

## 8. 지시서와 다르게 한 것

없음.

## 9. 산출물 목록, 커밋 해시

- Part A–C `50e5f00580c580d77fa67215b9dd097528ca8268` (2026-10-19 16:33:32 +0900). `scripts/56_part9_analysis.R`. `results/stage9/`의 L2 표 4개, 시뮬레이션 표 7개, 완화 표 5개.
- 그림 `c31ddf74158a59a282e8f3546780b50081d74d55` (2026-10-19 16:33:39 +0900). `scripts/57_part9_figures.py`, `figures/Fig6_simulation.pdf`, `figures/Fig6_simulation.png`, `figures/Fig7_mitigation.pdf`, `figures/Fig7_mitigation.png`, `figures/figure_sources.tsv`, `results/stage9/docx_check.tsv`.
- 보고서: 이 보고서를 넣는 커밋. `scripts/58_part9_report.py`, `reports/TOD_보고서_9.md`.
