# Preregistration: Time-of-day dependence of transcriptomic aging clocks in human blood (project TOD)

Version 1.0. This document is frozen before any expression value of the three confirmation datasets is read. Development analyses (TOD stages 1–5) are complete. Their outputs are referenced by path and SHA-256 in Section 10.

---

## 1. Questions

- **Q1.** Do transcriptomic aging clock scores vary with blood-sampling time within the same person?
- **Q2.** How much of that variation is explained by daily changes in blood cell composition?
- **Q3.** Does the morning–evening difference depend on age?
- **Q4.** Does misalignment of sleep with the circadian system, or a shift of the sleep schedule, change the clock rhythm?
- **Q5.** Does a time-of-day correction estimated in development data transfer to an independent cohort?

---

## 2. Data

**Development (fully analysed):** GSE39445 (Agilent microarray, leukocytes, 24 participants, extended wakefulness), GSE56931 (Affymetrix microarray, whole blood, 14 participants, baseline including normal sleep), GSE113883 (RNA-seq TPM, whole blood, 11 participants, sampling every 2 h for 29 h).

**Confirmation (locked).** Only sample metadata, platform annotation, and supplementary file names and sizes have been read.

| Code | Dataset | Technology | Tissue | Design |
|---|---|---|---|---|
| L1 | GSE48113 | Microarray, GPL15331 | Leukocytes | 22 participants, forced desynchrony. Day 1: sleep in phase with melatonin. Day 4: sleep out of phase |
| L2 | GSE122725 | Microarray, GPL24324 (Clariom S HT) | PBMC | 11 participants. Period 1: day-oriented schedule. Period 2: night-oriented schedule after a 10-h sleep delay. Bright-light (3) and dim-light (8) groups |
| L3 | GSE295758 | RNA-seq, raw transcript counts | Whole blood | 58 participants aged 20–75 with chronic low back pain. One morning (07:00–09:00) and one evening (19:00–21:00) sample each |

**Exclusions fixed before unlocking:**
- L1: GSM1168837 (missing sampling time) and GSM1168735 (Day 1 sample labelled out of phase).
- L3: subject8 and subject16 were sequenced in two batches. Only their batch-1 pair enters the primary analyses. Their batch-2 pair is used only for a technical-replicate check.

---

## 3. Time axes

| Dataset | Time axis |
|---|---|
| L1 | Clock time of sampling (`time sample taken`), hours in [0, 24) |
| L2 | Hours since habitual wake (`relclocktime`), modulo 24 |
| L3 | Morning = 8.0 h, evening = 20.0 h (window midpoints) |

---

## 4. Clocks

**Primary clocks.** Selected in development because their peak times lay within 2 h of their circular mean in all three development datasets (Table F1).
- **TAms_chrono:** multispecies chronological composite of Tyshkovskiy et al. (2026), file `EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl` (Zenodo 10.5281/zenodo.22166800, version 2; md5 `7dc8a470f5db20a75d7052829efbd37a`).
- **TAms_mort:** multispecies mortality composite, file `EN_Mortality_Multispecies_Multitissue_scaleddiff.pkl` (md5 `87307d7f5572d9932382912d3fd9571c`).
- **RA_default:** RNAAgeCalc 1.22.0, tissue `blood`, default signature (DESeq2), `stype = "all"`. Linear score: intercept + coefficients × log2 expression, genes absent from the data dropped.

**Secondary clocks:** RA_Peters and RA_GTEx (RNAAgeCalc 1.22.0, linear score as above).

**Module-level outcomes (secondary):** the reproducible outcomes in Table F3. They are one multispecies module clock and module components of the TAms composites.

**Scoring of TAms and module outcomes** follows the code of the authors' earlier revision (repository `clock-specificity`, commit in Table F5):
- Preprocessing: `array_to_diff` for microarrays, `counts_to_diff` for RNA-seq counts.
- Reference samples: all analysed samples of the dataset. Scores are computed once per dataset.
- Score: sum of effective weights × scaled-diff values over the selected non-zero genes, without intercept.
- Components: per-module partial sums of the same score, using the gene sets of the 23 rodent modules. Features outside every module form the component "outside".

---

## 5. Preprocessing

| Dataset | Steps |
|---|---|
| L1 | Series-matrix values. Log judgement by the GEO2R rule. Probes → genes by `GENE_SYMBOL`, single-symbol probes only, highest-mean probe per gene |
| L2 | Series-matrix values. Log judgement by the GEO2R rule. Probesets → SYMBOL with `clariomshumanhttranscriptcluster.db` 8.8.0, highest-mean probeset per symbol |
| L3 | Raw transcript counts (`GSE295758_rawTranscriptCounts.csv.gz`) summed to genes with the gene identifier column of the file. If the file has no gene column, transcripts are mapped to genes with GENCODE release 38 (GRCh38.p13). RNAAgeCalc converts counts to FPKM with its internal gene lengths; RA linear scores use log2(FPKM + 1). TAms and module outcomes use `counts_to_diff`. Gate genes use log2(CPM + 1) |

**Composition:** MCPcounter 1.2.0 scores of T cells, CD8 T cells, cytotoxic lymphocytes, NK cells, B lineage, monocytic lineage, myeloid dendritic cells, and neutrophils. Neutrophils are omitted in L2 (PBMC).

**No sample is removed for quality.** Samples whose median lies more than 3 MAD from the median of sample medians are flagged and used only in a sensitivity analysis.

---

## 6. Quality gates

Gates are applied before any hypothesis test. A dataset that fails its gate is not tested and is reported as failing.

| Dataset | Gate |
|---|---|
| L1 (Day 1) | PER1 and PER3 cosinor LRT p < 0.01, and both peaks within 04:00–12:00 |
| L2 (period 1) | PER1 and PER3 cosinor LRT p < 0.01, and both peaks within 18–26 h after habitual wake (modulo 24) |
| L3 | PER1 morning − evening > 0, paired t-test p < 0.01 |

Rationale: in development, these criteria passed in every dataset whose time variable was correct. Criteria based on ARNTL were unstable and are not used.

---

## 7. Hypotheses, tests and decision rules

**Model M (L1, L2):** `score ~ cos(2πt/24) + sin(2πt/24) + (1 | participant)` in lme4. Estimates by REML. Tests by ML likelihood ratio (df 2) against `score ~ 1 + (1 | participant)`.

**Effect sizes:** 2A (peak-to-trough), R²_within, 2A/σ, peak time. Intervals from a participant-level bootstrap, 2,000 resamples, seed 20261015.

### H1 — primary family P (Holm, family-wise α = 0.05)

| Test | Dataset | Clocks | Test |
|---|---|---|---|
| P1–P3 | L1, Day 1 | TAms_chrono, TAms_mort, RA_default | Model M |
| P4–P6 | L2, period 1 | TAms_chrono, TAms_mort, RA_default | Model M |
| P7–P9 | L3 | TAms_chrono, TAms_mort, RA_default | Paired two-sided t-test of Δ = morning − evening |

- A clock is **replicated** if its test is significant after Holm correction in at least two of the three datasets.
- Secondary checks (reported without multiplicity correction):
  - L1: the peak lies within ±3 h of the development circular mean φ̄ (Table F1).
  - L3: the sign of mean Δ equals the sign predicted from development (Table F2).

### H1, secondary clocks — family S1 (Holm across 6 tests)
RA_Peters and RA_GTEx in L1, L2 and L3, tested as in family P.

### H2 — composition (estimation; no p-value family)
Computed for each primary clock × dataset whose family-P test is significant.
- **L1, L2:** projection mediation `p_med = 1 − (r_a · r_u) / |r_u|²`. Here r_u and r_a are the (β_cos, β_sin) of model M without and with the composition scores as fixed covariates. 95% bootstrap interval from 2,000 participant resamples, refitting both models in each resample.
- **L3:** across participants, `Δclock = α + Σ_k γ_k Δcomp_k + e`. `p_med = 1 − α / mean(Δclock)`. 95% bootstrap interval from 2,000 participant resamples.
- H2 is supported for a clock × dataset if the lower bound of the interval is above 0.

### H3 — age — family S2 (Holm across the 3 primary clocks)
- L3: `Δclock ~ age (per 10 years) + sex + opioid group (factor 0/1/2)`. Two-sided test of the age coefficient.
- Expected direction: a smaller morning–evening difference with age, that is, a negative coefficient when mean Δ > 0. The observed direction is reported.

### H4 — misalignment — family S3 (Holm across 6 tests)
- **L1:** Day 1 and Day 4 together, `score ~ (cos + sin) × day + (1 | participant)`, t = clock time. ML likelihood ratio test of the two interaction terms (df 2).
- **L2:** period 1 and period 2 together, the same model with t = hours since habitual wake.
- Reported: amplitude and peak by condition, and the mean difference between conditions. The light groups of L2 are described only (3 vs 8 participants).

### H5 — transfer of a correction (estimation)
- For each primary clock, `c = f(8) − f(20)`, where `f(t) = β_cos cos(2πt/24) + β_sin sin(2πt/24)` from the GSE113883 fit of model M (Table F4).
- In L3, report `mean(Δ) − c` with a 95% interval and the ratio `mean(Δ) / c`.
- The correction is called **calibrated** for a clock if the interval of `mean(Δ) − c` includes 0 and the family-P test of that clock in L3 is significant.

### Module-level outcomes — family S4 (descriptive)
- **L1:** model M for each outcome in Table F3. An outcome is **confirmed** if LRT p < 0.05 and its peak lies within ±3 h of its development φ̄. Report the number confirmed out of the list.
- **L2:** model M on the wake-referenced axis; LRT p only.
- **L3:** Δ for each outcome. Sign agreement with the development prediction for outcomes with |d_exp| ≥ 0.2 (Table F3).
- All other module clocks and components: model M or Δ, with BH-adjusted q within dataset × clock, reported as exploratory.

---

## 8. Sensitivity analyses
Pre-specified and reported, with no inferential weight.
- Excluding flagged samples (3 MAD).
- Random slopes for cos and sin.
- L3: batch-adjusted counts (`GSE295758_combatVstTranscriptCounts.csv.gz`) for the primary clocks.
- L3: batch-2 replicate pairs of subject8 and subject16 compared with their batch-1 pairs.
- L3: regression of Δ on sex and opioid group, testing the intercept.

---

## 9. Deviations
- Any change after unlocking is reported in a deviation table with its reason and the affected outputs. The preregistered analysis is always reported first.
- Format adaptations needed to read the locked files (column names, identifier formats) are made without inspecting values and are listed.

---

## 10. Frozen development parameters
Filled from the development outputs before freezing. Each table lists its source file and SHA-256.

- **Table F1.** Primary and secondary clocks: φ̄ (h), maximum circular distance, and per development dataset 2A, σ, σ_u, peak.
- **Table F2.** L3 predictions: d_exp per development dataset, median d_exp, predicted sign.
- **Table F3.** Reproducible module-level outcomes: φ̄, per-dataset q and peak, median d_exp.
- **Table F4.** H5 constants from GSE113883: β_cos, β_sin, c.
- **Table F5.** Software versions, code commits, clock files and md5, development output files and SHA-256.
- **Table F6.** Power computed from development parameters.


## Appendix. Tables F1-F6

### Table F1

| clock | role | phi_bar | max_circ_dist | dataset | two_a | sigma | sigma_u | peak_h | lrt_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TAms_chrono | primary | 12.616768998085 | 1.61677078040933 | GSE39445 | 0.0345210353298852 | 0.0171269287358502 | 0.0276429862743754 | 14.2335397784943 | 5.49745798912395e-20 |
| TAms_chrono | primary | 12.616768998085 | 1.61677078040933 | GSE56931 | 0.00943066206760042 | 0.0197534498990921 | 0.0331124542924689 | 11.4579688514277 | 0.181165188487945 |
| TAms_chrono | primary | 12.616768998085 | 1.61677078040933 | GSE113883 | 0.0418702301338288 | 0.0253117943835499 | 0.0416315613478055 | 12.188050742382 | 2.96936394371969e-10 |
| TAms_mort | primary | 11.8599407935094 | 1.79984172390618 | GSE39445 | 0.148755555495313 | 0.0860831456543238 | 0.161562476919705 | 13.4288594848278 | 3.45851149483147e-15 |
| TAms_mort | primary | 11.8599407935094 | 1.79984172390618 | GSE56931 | 0.0629468359700225 | 0.0740143176007128 | 0.177501208368076 | 12.0688456403715 | 0.00673707816031792 |
| TAms_mort | primary | 11.8599407935094 | 1.79984172390618 | GSE113883 | 0.102562103238955 | 0.0897840656082817 | 0.184575482929547 | 10.0600990696032 | 4.44466652817908e-06 |
| RA_default | primary | 11.2423863477089 | 1.35536907505411 | GSE39445 | 1.98498899424841 | 9.34019861854301 | 8.81869510215541 | 12.0064473932836 | 0.571779550797501 |
| RA_default | primary | 11.2423863477089 | 1.35536907505411 | GSE56931 | 4.61694071236019 | 6.68803391836383 | 6.90600678789726 | 9.88701727265479 | 0.0187399469939356 |
| RA_default | primary | 11.2423863477089 | 1.35536907505411 | GSE113883 | 6.770165944749 | 6.28931435261543 | 6.86875908783177 | 11.812631972311 | 3.24675832907807e-05 |
| RA_Peters | secondary | 7.56604107260516 | 4.27017978445731 | GSE39445 | 3.86882060477993 | 9.12159876450608 | 14.613080804435 | 6.84503499428439 | 0.0822398090507864 |
| RA_Peters | secondary | 7.56604107260516 | 4.27017978445731 | GSE56931 | 6.37997357324785 | 5.58777765048951 | 13.9973116863226 | 4.54192615314062 | 0.000105116530090887 |
| RA_Peters | secondary | 7.56604107260516 | 4.27017978445731 | GSE113883 | 3.1551856815567 | 6.04743769494992 | 13.6399297750127 | 11.8362208570625 | 0.0773443137299791 |
| RA_GTEx | secondary | 11.518486825228 | 5.35019887986941 | GSE39445 | 4.5622147045564 | 25.1603887545092 | 16.2536695399875 | 13.7137336093391 | 0.643663244864414 |
| RA_GTEx | secondary | 11.518486825228 | 5.35019887986941 | GSE56931 | 8.41295616809185 | 14.1878999046968 | 11.0606886507247 | 13.2671292133393 | 0.0977660921594727 |
| RA_GTEx | secondary | 11.518486825228 | 5.35019887986941 | GSE113883 | 4.34369415553422 | 8.81987355742966 | 10.4317015658613 | 6.16828794535857 | 0.104384139695734 |

Source files:

- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/clock_summary_long.tsv` SHA-256 `bb2be685a97ab6a2141ff49c7412902b3d936bfbd803a3c3ee48be00301b9dad`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/clock_summary_phase.tsv` SHA-256 `38077898da876ea43de2166396e270ded605aa11c393727872b3f3dadd3d7d7b`

### Table F2

| clock | role | d_GSE39445 | d_GSE56931 | d_GSE113883 | d_exp_median | predicted_sign |
| --- | --- | --- | --- | --- | --- | --- |
| TAms_chrono | primary | -0.087085971883512 | 0.208443813199012 | 0.534281736707176 | 0.208443813199012 | + |
| TAms_mort | primary | 0.182025471981377 | 0.291250817026961 | 0.693084284817919 | 0.291250817026961 | + |
| RA_default | primary | 0.0749177677231186 | 0.429771444888433 | 0.412448894924069 | 0.412448894924069 | + |
| RA_Peters | secondary | 0.286305297017999 | 0.498487573681539 | 0.197988187071947 | 0.286305297017999 | + |
| RA_GTEx | secondary | 0.00960009963837758 | 0.0799546478017167 | 0.308963171982023 | 0.0799546478017167 | + |

Source files:

- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/gse295758_expected_d.tsv` SHA-256 `b6c1ef3c4db6202a467ae3636d0de3aa7e55b47bf766bcd67ece61ca538c7e98`

### Table F3

| clock | family | module | parent_clock | phi_bar | q_GSE39445 | q_GSE56931 | q_GSE113883 | peak_GSE39445 | peak_GSE56931 | peak_GSE113883 | d_exp_median | predicted_sign |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TM_turquoise_mortality | TM | turquoise |  | 11.5724011994242 | 1.7030021051141e-07 | 0.00657277863416067 | 3.693684980855e-06 | 12.6215119266462 | 11.1568084467969 | 10.9484328837053 | 0.467000081226947 | + |
| CMP_TAms_chrono_turquoise | CMP | turquoise | TAms_chrono | 13.2837328415609 | 2.42711785115231e-08 | 0.00626983142167074 | 0.0107579084724153 | 13.1432263174321 | 12.1422987334514 | 14.5730749724508 | 0.20086732450169 | + |
| CMP_TAms_chrono_yellow | CMP | yellow | TAms_chrono | 12.6128269162576 | 0.00155470474860942 | 4.69738653218747e-08 | 9.00796184650676e-09 | 13.045845902978 | 11.2503142784199 | 13.523129022035 | 0.1496592336641 | + |
| CMP_TAms_chrono_orange | CMP | orange | TAms_chrono | 1.07332762742524 | 4.17596855270051e-09 | 0.00616665634539307 | 9.86133709312738e-06 | 2.44005820482784 | 0.0303979642477727 | 0.765259718398529 | -0.280140424913859 | - |
| CMP_TAms_chrono_brown | CMP | brown | TAms_chrono | 1.32505139804729 | 2.43965028951066e-08 | 0.0359875447558034 | 2.67476690758748e-05 | 0.91282893412168 | 0.84038477850989 | 2.22823394895771 | -0.168830899708354 | - |
| CMP_TAms_chrono_red | CMP | red | TAms_chrono | 0.781691026210673 | 6.00047155091236e-06 | 0.00206822511320336 | 0.00113591207810166 | 2.08110231751181 | 23.2474694413108 | 1.00061811727099 | -0.169417336437209 | - |
| CMP_TAms_mort_turquoise | CMP | turquoise | TAms_mort | 13.003811149926 | 6.77796115376162e-09 | 0.000125036076009668 | 3.88973233356303e-07 | 13.2360241903645 | 11.9113647277402 | 13.8564175468774 | 0.1888870282562 | + |
| CMP_TAms_mort_darkred | CMP | darkred | TAms_mort | 14.9666596497573 | 0.00299942175721486 | 0.0275912785277717 | 9.48016174502148e-08 | 17.0651005937284 | 13.327427235357 | 14.5608059570561 | -0.167328334678284 | - |

Source files:

- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/reproducible_modules.tsv` SHA-256 `6b6c9a1b27cc34b060358fc656d51e3ba329c3a3acdac3c3ee536c6c2ab4dadd`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/gse295758_expected_d.tsv` SHA-256 `b6c1ef3c4db6202a467ae3636d0de3aa7e55b47bf766bcd67ece61ca538c7e98`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/h1_point_estimates.tsv` SHA-256 `60f24af5ff25abce59570d1f79d391adbcd9093adad2f8faf3f90c34e2c5c2f1`

### Table F4

| clock | dataset | beta_cos | beta_sin | formula | c |
| --- | --- | --- | --- | --- | --- |
| TAms_chrono | GSE113883 | -0.0209097494754014 | -0.00103025227051267 | c = f(8) - f(20); f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24); c = -beta_cos + beta_sin*sqrt(3) | 0.019125300198260260499475555777287663859727086340733 |
| TAms_mort | GSE113883 | -0.0448086053054778 | 0.0249386275840794 | c = f(8) - f(20); f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24); c = -beta_cos + beta_sin*sqrt(3) | 0.088003575352142004086395806850464311947738779118951 |
| RA_default | GSE113883 | -3.38101122840144 | 0.165981334106488 | c = f(8) - f(20); f(t) = beta_cos*cos(2*pi*t/24) + beta_sin*sin(2*pi*t/24); c = -beta_cos + beta_sin*sqrt(3) | 3.6684993321819421764583897146097149703131816543407 |

Source files:

- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/h1_point_estimates.tsv` SHA-256 `60f24af5ff25abce59570d1f79d391adbcd9093adad2f8faf3f90c34e2c5c2f1`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage4/rnaseq/h1_point_estimates.tsv` SHA-256 `f0ad35b89eb5b0a97979b809448bde321cd97c4ec658af6cae87cb835e6ff85e`

### Table F5

| item | value |
| --- | --- |
| R | 4.5.3 |
| lme4 | 2.0.6 |
| RNAAgeCalc | 1.22.0 |
| MCPcounter | 1.2.0 |
| impute | 1.84.0 |
| limma | 3.66.0 |
| clariomshumanhttranscriptcluster.db | 8.8.0 |
| org.Hs.eg.db | 3.22.0 |
| python | 3.12.7 |
| numpy | 1.26.4 |
| pandas | 2.2.3 |
| scikit-learn | 1.3.2 |
| scipy | 1.13.1 |
| joblib | 1.6.0 |
| clock_python | /HDD8T1/WORK/GeroSci/tod4_clock_env |
| clock_specificity_commit | 9d109499f0bf47a7869cd7b4bcffe311b8be9b75 |
| tAge_tag | 1.5.0 |
| md5_TAms_chrono | 7dc8a470f5db20a75d7052829efbd37a |
| md5_TAms_mort | 87307d7f5572d9932382912d3fd9571c |
| sha256_clock_summary_long | bb2be685a97ab6a2141ff49c7412902b3d936bfbd803a3c3ee48be00301b9dad |
| path_clock_summary_long | /HDD8T1/WORK/GeroSci/TOD/results/stage5/clock_summary_long.tsv |
| sha256_clock_summary_phase | 38077898da876ea43de2166396e270ded605aa11c393727872b3f3dadd3d7d7b |
| path_clock_summary_phase | /HDD8T1/WORK/GeroSci/TOD/results/stage5/clock_summary_phase.tsv |
| sha256_expected_d | b6c1ef3c4db6202a467ae3636d0de3aa7e55b47bf766bcd67ece61ca538c7e98 |
| path_expected_d | /HDD8T1/WORK/GeroSci/TOD/results/stage5/gse295758_expected_d.tsv |
| sha256_reproducible | 6b6c9a1b27cc34b060358fc656d51e3ba329c3a3acdac3c3ee536c6c2ab4dadd |
| path_reproducible | /HDD8T1/WORK/GeroSci/TOD/results/stage5/reproducible_modules.tsv |
| sha256_h1_stage5 | 60f24af5ff25abce59570d1f79d391adbcd9093adad2f8faf3f90c34e2c5c2f1 |
| path_h1_stage5 | /HDD8T1/WORK/GeroSci/TOD/results/stage5/h1_point_estimates.tsv |
| sha256_h1_stage4_rnaseq | f0ad35b89eb5b0a97979b809448bde321cd97c4ec658af6cae87cb835e6ff85e |
| path_h1_stage4_rnaseq | /HDD8T1/WORK/GeroSci/TOD/results/stage4/rnaseq/h1_point_estimates.tsv |
| sha256_g2_clocks | 9984f8ad52ecab1499a5d12e36e2923c03d0acf675a193f5416dd9b3cf4181ff |
| path_g2_clocks | /HDD8T1/WORK/GeroSci/TOD/results/stage4/rnaseq/g2_clocks.tsv |
| sha256_g2_tams | 2fd4c6b63eed1cd60c3adf4b265803e7ef660e670769823d4df39e8aa64e257a |
| path_g2_tams | /HDD8T1/WORK/GeroSci/TOD/results/stage5/g2_tams.tsv |
| sha256_power | e508d01b1b3f1ab9893c0e0332be9df90140ac49defb3faeee071c3525bdfb37 |
| path_power | /HDD8T1/WORK/GeroSci/TOD/results/stage5/locked_design_power.tsv |

Source files:

- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/clock_summary_long.tsv` SHA-256 `bb2be685a97ab6a2141ff49c7412902b3d936bfbd803a3c3ee48be00301b9dad`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/clock_summary_phase.tsv` SHA-256 `38077898da876ea43de2166396e270ded605aa11c393727872b3f3dadd3d7d7b`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/gse295758_expected_d.tsv` SHA-256 `b6c1ef3c4db6202a467ae3636d0de3aa7e55b47bf766bcd67ece61ca538c7e98`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/reproducible_modules.tsv` SHA-256 `6b6c9a1b27cc34b060358fc656d51e3ba329c3a3acdac3c3ee536c6c2ab4dadd`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/h1_point_estimates.tsv` SHA-256 `60f24af5ff25abce59570d1f79d391adbcd9093adad2f8faf3f90c34e2c5c2f1`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage4/rnaseq/h1_point_estimates.tsv` SHA-256 `f0ad35b89eb5b0a97979b809448bde321cd97c4ec658af6cae87cb835e6ff85e`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage4/rnaseq/g2_clocks.tsv` SHA-256 `9984f8ad52ecab1499a5d12e36e2923c03d0acf675a193f5416dd9b3cf4181ff`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/g2_tams.tsv` SHA-256 `2fd4c6b63eed1cd60c3adf4b265803e7ef660e670769823d4df39e8aa64e257a`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/locked_design_power.tsv` SHA-256 `e508d01b1b3f1ab9893c0e0332be9df90140ac49defb3faeee071c3525bdfb37`

### Table F6

| block | clock | design | param_dataset | d_exp_median | power_alpha_0.05 | power_alpha_0.01 | frac_p_lt_0.05 | frac_p_lt_0.01 | n_rep | n_converged |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| L3 | TAms_chrono | GSE295758 | median_of_three | 0.208443813199012 | 0.345077590733327 | 0.150545763753431 |  |  |  |  |
| L3 | TAms_mort | GSE295758 | median_of_three | 0.291250817026961 | 0.587339593367421 | 0.336502748601378 |  |  |  |  |
| L3 | RA_default | GSE295758 | median_of_three | 0.412448894924069 | 0.870300205086949 | 0.682114842365414 |  |  |  |  |
| L3 | RA_Peters | GSE295758 | median_of_three | 0.286305297017999 | 0.572873321105643 | 0.323267563613161 |  |  |  |  |
| L3 | RA_GTEx | GSE295758 | median_of_three | 0.0799546478017167 | 0.0919710023400413 | 0.0243768273749048 |  |  |  |  |
| L1L2 | RA_GTEx | GSE48113 | GSE39445 |  |  |  | 0.0935 | 0.027 | 2000 | 2000 |
| L1L2 | RA_GTEx | GSE48113 | GSE56931 |  |  |  | 0.579 | 0.3245 | 2000 | 2000 |
| L1L2 | RA_GTEx | GSE48113 | GSE113883 |  |  |  | 0.468 | 0.2485 | 2000 | 2000 |
| L1L2 | RA_Peters | GSE48113 | GSE39445 |  |  |  | 0.361 | 0.1545 | 2000 | 2000 |
| L1L2 | RA_Peters | GSE48113 | GSE56931 |  |  |  | 0.998 | 0.968 | 2000 | 2000 |
| L1L2 | RA_Peters | GSE48113 | GSE113883 |  |  |  | 0.477 | 0.244 | 2000 | 2000 |
| L1L2 | RA_default | GSE48113 | GSE39445 |  |  |  | 0.112 | 0.0375 | 2000 | 2000 |
| L1L2 | RA_default | GSE48113 | GSE56931 |  |  |  | 0.773 | 0.5335 | 2000 | 2000 |
| L1L2 | RA_default | GSE48113 | GSE113883 |  |  |  | 0.9795 | 0.931 | 2000 | 2000 |
| L1L2 | TAms_chrono | GSE48113 | GSE39445 |  |  |  | 1 | 1 | 2000 | 2000 |
| L1L2 | TAms_chrono | GSE48113 | GSE56931 |  |  |  | 0.431 | 0.201 | 2000 | 2000 |
| L1L2 | TAms_chrono | GSE48113 | GSE113883 |  |  |  | 1 | 1 | 2000 | 2000 |
| L1L2 | TAms_mort | GSE48113 | GSE39445 |  |  |  | 1 | 1 | 2000 | 2000 |
| L1L2 | TAms_mort | GSE48113 | GSE56931 |  |  |  | 0.8865 | 0.7175 | 2000 | 2000 |
| L1L2 | TAms_mort | GSE48113 | GSE113883 |  |  |  | 0.9945 | 0.9775 | 2000 | 2000 |
| L1L2 | RA_GTEx | GSE122725 | GSE39445 |  |  |  | 0.0755 | 0.0225 | 2000 | 2000 |
| L1L2 | RA_GTEx | GSE122725 | GSE56931 |  |  |  | 0.3395 | 0.1505 | 2000 | 2000 |
| L1L2 | RA_GTEx | GSE122725 | GSE113883 |  |  |  | 0.1995 | 0.067 | 2000 | 2000 |
| L1L2 | RA_Peters | GSE122725 | GSE39445 |  |  |  | 0.1635 | 0.059 | 2000 | 2000 |
| L1L2 | RA_Peters | GSE122725 | GSE56931 |  |  |  | 0.83 | 0.6185 | 2000 | 2000 |
| L1L2 | RA_Peters | GSE122725 | GSE113883 |  |  |  | 0.2995 | 0.1155 | 2000 | 2000 |
| L1L2 | RA_default | GSE122725 | GSE39445 |  |  |  | 0.0975 | 0.03 | 2000 | 2000 |
| L1L2 | RA_default | GSE122725 | GSE56931 |  |  |  | 0.456 | 0.22 | 2000 | 2000 |
| L1L2 | RA_default | GSE122725 | GSE113883 |  |  |  | 0.842 | 0.644 | 2000 | 2000 |
| L1L2 | TAms_chrono | GSE122725 | GSE39445 |  |  |  | 1 | 0.9985 | 2000 | 2000 |
| L1L2 | TAms_chrono | GSE122725 | GSE56931 |  |  |  | 0.237 | 0.092 | 2000 | 2000 |
| L1L2 | TAms_chrono | GSE122725 | GSE113883 |  |  |  | 0.997 | 0.976 | 2000 | 2000 |
| L1L2 | TAms_mort | GSE122725 | GSE39445 |  |  |  | 0.9975 | 0.988 | 2000 | 2000 |
| L1L2 | TAms_mort | GSE122725 | GSE56931 |  |  |  | 0.6175 | 0.3815 | 2000 | 2000 |
| L1L2 | TAms_mort | GSE122725 | GSE113883 |  |  |  | 0.8635 | 0.6735 | 2000 | 2000 |

Source files:

- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/gse295758_expected_d.tsv` SHA-256 `b6c1ef3c4db6202a467ae3636d0de3aa7e55b47bf766bcd67ece61ca538c7e98`
- `/HDD8T1/WORK/GeroSci/TOD/results/stage5/locked_design_power.tsv` SHA-256 `e508d01b1b3f1ab9893c0e0332be9df90140ac49defb3faeee071c3525bdfb37`

