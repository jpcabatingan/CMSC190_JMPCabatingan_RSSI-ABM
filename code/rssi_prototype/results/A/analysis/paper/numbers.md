# Window A seeded runs: numbers

All values are computed from the raw CSVs by analysis.py.

## 1. Main comparison (recomputed)

| config | n | mean | sd | lo | hi | ever | end | chem | bio | cane | maxD |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chemical (n_half none) | 2 | 30,013 | 1,024 | 29,289 | 30,738 | 0.192 | 0.078 | 34,882 | 0 | 61.63 | 0.103 |
| hybrid (n_half 10) | 2 | 42,379 | 4,789 | 38,992 | 45,765 | 0.249 | 0.110 | 0 | 52,891 | 61.62 | 0.111 |
| hybrid (n_half none) | 2 | 42,777 | 145 | 42,675 | 42,880 | 0.258 | 0.111 | 0 | 52,902 | 61.62 | 0.116 |
| hybrid (n_half 20) | 2 | 42,789 | 1,021 | 42,067 | 43,511 | 0.256 | 0.111 | 0 | 53,128 | 61.62 | 0.115 |
| biological | 2 | 44,123 | 4,959 | 40,616 | 47,629 | 0.263 | 0.115 | 0 | 54,588 | 61.62 | 0.116 |
| chemical (n_half 20) | 2 | 754,011 | 17,092 | 741,924 | 766,097 | 0.340 | 0.275 | 85,602 | 0 | 61.37 | 0.956 |
| chemical (n_half 10) | 2 | 1,919,710 | 110,037 | 1,841,902 | 1,997,518 | 0.491 | 0.438 | 125,474 | 0 | 60.92 | 0.983 |
| none | 2 | 10,763,962 | 222,266 | 10,606,796 | 10,921,127 | 0.965 | 0.963 | 0 | 0 | 57.24 | 1.039 |

Check against summary.md section 1: 8 of 8 configurations matched by name; largest difference in mean load = 0 K-days (rounding only if 0 or 1).

## 2. Replicate variability

Main comparison, n = 2 per configuration. CV = SD / mean; relative range = (max - min) / mean.

| config | cv_load | rr_load | rr_ever | rr_apps |
|---|---|---|---|---|
| biological | 11.2% | 15.9% | 17.2% | 15.7% |
| chemical (n_half none) | 3.4% | 4.8% | 8.1% | 5.9% |
| chemical (n_half 10) | 5.7% | 8.1% | 8.3% | 8.0% |
| chemical (n_half 20) | 2.3% | 3.2% | 5.3% | 3.6% |
| hybrid (n_half none) | 0.3% | 0.5% | 0.5% | 0.4% |
| hybrid (n_half 10) | 11.3% | 16.0% | 13.6% | 15.8% |
| hybrid (n_half 20) | 2.4% | 3.4% | 3.2% | 1.9% |
| none | 2.1% | 2.9% | 0.8% | nan% |

Pooled within-pair SD of ln(load) by strategy, from all 20 replicate pairs (main 8, wind 8, saturating 4). For small values this is about the CV.

| strategy | pairs | sd_ln |
|---|---|---|
| none | 4 | 0.020 |
| chemical | 6 | 0.045 |
| biological | 4 | 0.089 |
| hybrid | 6 | 0.070 |

Baseline pool for pest load (main at n_half 20 + saturating-damage runs + OAT case 0; damage does not affect the pest, so these are replicates of one process):

| strategy | size | mean | std | min | max | cv |
|---|---|---|---|---|---|---|
| none | 5 | 10,879,737 | 189,567 | 10,606,796 | 11,137,303 | 1.7% |
| chemical | 5 | 775,970 | 40,820 | 741,924 | 846,833 | 5.3% |
| biological | 5 | 42,653 | 3,184 | 39,987 | 47,629 | 7.5% |
| hybrid | 5 | 44,407 | 3,098 | 41,118 | 47,688 | 7.0% |

Adjacent configurations in the main ranking: ratio of mean loads and whether it exceeds the noise band exp(2 x sqrt(s_a^2/2 + s_b^2/2)) from the pooled SDs (n = 2 each).

| lower load | higher load | ratio | noise band | verdict |
|---|---|---|---|---|
| chemical (n_half none) | hybrid (n_half 10) | 1.41 | 1.12 | distinct |
| hybrid (n_half 10) | hybrid (n_half none) | 1.01 | 1.15 | tie |
| hybrid (n_half none) | hybrid (n_half 20) | 1.00 | 1.15 | tie |
| hybrid (n_half 20) | biological | 1.03 | 1.17 | tie |
| biological | chemical (n_half 20) | 17.09 | 1.15 | distinct |
| chemical (n_half 20) | chemical (n_half 10) | 2.55 | 1.09 | distinct |
| chemical (n_half 10) | none | 5.61 | 1.07 | distinct |

## 3. Efficacy decline: where chemical would match biological/hybrid

Chemical mean load: no decline 30,013; n_half 20 754,011 (25.1x); n_half 10 1,919,710 (64.0x).
Biological and hybrid baseline runs pooled (n = 14): mean 43,258, range 38,992-47,688.
Ratios to that pool: chemical no decline 0.69; n_half 20 17.4x; n_half 10 44.4x.

Break-even between 'no decline' (decline rate 1/n_half = 0) and n_half 20 (rate 0.05), two interpolations over only two points:
- log-linear in load: rate 0.0057, n_half about 176
- linear in load: rate 0.00091, n_half about 1,093
The two differ by a factor of several: the data only bound the break-even at n_half > 20.

## 4. Hybrid vs biological

Hybrid runs in all files: 31; with any chemical application: 3 (oat case 13: 435 chemical, 145,106 biological, oat case 14: 209,538 chemical, 0 biological, oat case 16: 40,896 chemical, 44,732 biological).
Max cell density in the 9 baseline hybrid runs: 0.107-0.116 K; farms are 2x2 cells and Level 3 needs a FARM MEAN of 0.11 K, so one cell at that density gives a farm mean near 0.029 K if its neighbours are clean.

Consistency check: with bio:chem ratio 1.0 (case 10) biological is the same process as chemical without decline: biological case 10 = 28,925 vs chemical no decline = 30,013 (ratio 0.96).

Baseline, no chemicals used: hybrid n = 9 mean 43,594 (CV 6.7%); biological n = 5 mean 42,653 (CV 7.5%); ratio hybrid/biological 1.022; difference of mean ln(load) = +0.022, SE 0.040.

OAT cases where hybrid could differ (single runs):

| case | assumption | hybrid chem apps | hybrid bio apps | hybrid/biological | hybrid/chemical | hybrid max_D |
|---|---|---|---|---|---|---|
| 11 | level thresholds x0.5 | 0 | 27,856 | 0.96 | 0.02 | 0.07 |
| 12 | level thresholds x2 | 0 | 79,977 | 0.94 | 0.11 | 0.23 |
| 13 | trigger T = 2 | 435 | 145,106 | 0.96 | 0.32 | 0.43 |
| 14 | trigger T = 3 | 209,538 | 0 | 2.03 | 1.00 | 0.92 |
| 15 | re-treatment 7 d | 0 | 63,887 | 0.91 | 0.02 | 0.11 |
| 16 | re-treatment 28 d | 40,896 | 44,732 | 1.28 | 0.62 | 0.92 |

## 5. One-at-a-time screen (1 run per cell)

Noise check, OAT case 0 (single run) over the pooled baseline mean: none 1.02; chemical 0.99; biological 0.94; hybrid 1.07.

Tie rule for two single runs: |ln ratio| < 2 x sqrt(s_a^2 + s_b^2) with the pooled SDs of section 2 (biological vs hybrid: within a factor 1.25). Baseline ranking (pooled means): biological ~ hybrid < chemical < none.

Load ratio to the pooled baseline; (ns) = within 2 SD of replicate noise.

| case | assumption | none | chemical | biological | hybrid | ranking | changes? |
|---|---|---|---|---|---|---|---|
| 1 | settlement survival 0.01 | 0.04 | 0.01 | 0.27 | 0.23 | chemical < biological ~ hybrid < none | YES |
| 2 | settlement survival 0.10 | 1.32 | 3.34 | 2.04 | 2.00 | biological ~ hybrid < chemical < none | no |
| 3 | neighbour transmission x0.25 | 0.45 | 0.28 | 0.24 | 0.24 | biological ~ hybrid < chemical < none | no |
| 4 | neighbour transmission x4 | 1.70 | 7.28 | 5.01 | 4.80 | biological ~ hybrid < chemical < none | no |
| 5 | long jump 0 | 0.75 | 0.95 (ns) | 0.81 | 0.71 | biological ~ hybrid < chemical < none | no |
| 6 | long jump 0.2 | 1.29 | 1.23 | 1.49 | 1.46 | biological ~ hybrid < chemical < none | no |
| 7 | parasitism 0 | 1.31 | 2.27 | 1.50 | 1.42 | biological ~ hybrid < chemical < none | no |
| 8 | parasitism 0.668 | 0.48 | 0.07 | 0.54 | 0.49 | biological ~ hybrid < chemical < none | no |
| 9 | bio:chem ratio 0.8 | 0.98 (ns) | 1.05 (ns) | 1.33 | 1.25 | biological ~ hybrid < chemical < none | no |
| 10 | bio:chem ratio 1.0 | 1.03 (ns) | 0.95 (ns) | 0.68 | 0.70 | biological ~ hybrid < chemical < none | no |
| 11 | level thresholds x0.5 | 1.02 (ns) | 0.82 | 0.28 | 0.26 | biological ~ hybrid < chemical < none | no |
| 12 | level thresholds x2 | 0.98 (ns) | 1.53 | 3.15 | 2.86 | biological ~ hybrid < chemical < none | no |
| 13 | trigger T = 2 | 0.96 (ns) | 2.76 | 16.55 | 15.22 | biological ~ hybrid < chemical < none | no |
| 14 | trigger T = 3 | 1.04 (ns) | 3.69 | 33.07 | 64.50 | biological < chemical ~ hybrid < none | YES |
| 15 | re-treatment 7 d | 1.04 (ns) | 2.15 | 0.89 (ns) | 0.78 | biological ~ hybrid < chemical < none | no |
| 16 | re-treatment 28 d | 0.99 (ns) | 2.82 | 24.84 | 30.57 | biological < hybrid < chemical < none | YES |
| 17 | wind bias 0.5 | 1.02 (ns) | 0.91 (ns) | 1.02 (ns) | 0.96 (ns) | biological ~ hybrid < chemical < none | no |
| 18 | wind bias 1.0 | 0.91 | 0.89 | 0.93 (ns) | 0.88 (ns) | biological ~ hybrid < chemical < none | no |

Ranking unchanged in 15 of 18 cases. 'none' is last in 18 of 18. Chemical (n_half 20) is worse than biological in 17 of 18.

Influence = largest |log10 load ratio| over the low and high case (1.0 = tenfold).

| assumption | none | chemical | biological | hybrid | max_control |
|---|---|---|---|---|---|
| settlement survival | 1.42 | 1.96 | 0.57 | 0.64 | 1.96 |
| trigger level | 0.02 | 0.57 | 1.52 | 1.81 | 1.81 |
| re-treatment interval | 0.02 | 0.45 | 1.40 | 1.49 | 1.49 |
| parasitism | 0.32 | 1.16 | 0.27 | 0.31 | 1.16 |
| neighbour transmission | 0.35 | 0.86 | 0.70 | 0.68 | 0.86 |
| level thresholds | 0.01 | 0.18 | 0.55 | 0.59 | 0.59 |
| long jump | 0.13 | 0.09 | 0.17 | 0.16 | 0.17 |
| bio:chem ratio | 0.01 | 0.02 | 0.17 | 0.15 | 0.17 |
| wind bias | 0.04 | 0.05 | 0.03 | 0.06 | 0.06 |

## 6. Secondary outputs (main comparison means)

Farms = 10,611 farmer agents; applications are farm-level events; per farm per year averages over ALL farms, treated or not.

| config | ever | end | apps_per_farm_yr | averted_per_app | reduction_vs_none |
|---|---|---|---|---|---|
| chemical (n_half none) | 0.192 | 0.078 | 1.10 | 308 | 99.72% |
| hybrid (n_half 10) | 0.249 | 0.110 | 1.66 | 203 | 99.61% |
| hybrid (n_half none) | 0.258 | 0.111 | 1.66 | 203 | 99.60% |
| hybrid (n_half 20) | 0.256 | 0.111 | 1.67 | 202 | 99.60% |
| biological | 0.263 | 0.115 | 1.71 | 196 | 99.59% |
| chemical (n_half 20) | 0.340 | 0.275 | 2.69 | 117 | 93.00% |
| chemical (n_half 10) | 0.491 | 0.438 | 3.94 | 70 | 82.17% |
| none | 0.965 | 0.963 | 0.00 | nan | 0.00% |

Cane per harvested ha: none 57.18-57.31; all controlled configurations 60.90-61.63 (Y0 = 61.64; set by the assumed damage slope).

## 7. Checks

Wind direction at bias 1.0 (n = 2 each): ratio 315/135 vs noise band exp(2 x s).

| strategy | load 135 | load 315 | ratio | band | verdict |
|---|---|---|---|---|---|
| none | 9,995,631 | 9,683,729 | 0.97 | 1.04 | no difference |
| chemical | 640,478 | 666,400 | 1.04 | 1.09 | no difference |
| biological | 39,537 | 38,132 | 0.96 | 1.20 | no difference |
| hybrid | 36,674 | 40,837 | 1.11 | 1.15 | no difference |

Saturating vs linear damage (n = 2 each). Pest load cannot differ by construction (damage does not act on the pest); only cane can.

| strategy | load linear | load saturating | cane linear | cane saturating |
|---|---|---|---|---|
| none | 10,763,962 | 10,866,730 | 57.24 | 57.60 |
| chemical | 754,011 | 802,936 | 61.37 | 61.37 |
| biological | 44,123 | 42,516 | 61.62 | 61.62 |
| hybrid | 42,789 | 44,385 | 61.62 | 61.62 |

## 8. Earlier island runs (random start over 1% of cells, superseded design; context only)

| strategy | n_half | size | mean |
|---|---|---|---|
| chemical | 0 | 2 | 821,714 |
| none | 20 | 2 | 151,844,929 |

none / chemical (no decline): island random start 185x; window A seeded 359x. Different start and extent: not comparable beyond the sign.

## 9. Chemical at slower decline rates (window A, 2 runs each)

Biological/hybrid pool: mean 43,258, range 38,992-47,688.

| n_half | size | mean | min | max | ratio_to_bio_hybrid |
|---|---|---|---|---|---|
| 160 | 2 | 34,073 | 33,085 | 35,061 | 0.79 |
| 80 | 2 | 38,761 | 35,381 | 42,141 | 0.90 |
| 40 | 2 | 168,728 | 160,878 | 176,579 | 3.90 |

## 10. Hybrid test: chemical WITHOUT decline, conditions where farms reach Level 3 (window A, 2 runs each)

Case 14 = treat only at Level 3; case 16 = re-treat every 28 days.

| oat_case | strategy | n | mean | lo | hi | chem | bio |
|---|---|---|---|---|---|---|---|
| 14 | biological | 2 | 1,378,064 | 1,350,026 | 1,406,101 | 0 | 161,452 |
| 14 | chemical | 2 | 1,322,895 | 1,309,574 | 1,336,216 | 142,009 | 0 |
| 14 | hybrid | 2 | 1,220,375 | 1,173,803 | 1,266,947 | 131,110 | 0 |
| 16 | biological | 2 | 1,119,909 | 1,114,429 | 1,125,389 | 0 | 90,630 |
| 16 | chemical | 2 | 262,998 | 258,022 | 267,974 | 64,860 | 0 |
| 16 | hybrid | 2 | 693,075 | 672,924 | 713,226 | 38,023 | 51,265 |

Case 14: hybrid/biological 0.89; hybrid/chemical 0.92; biological/chemical 1.04. Noise band for 2-run means: about 1.17.

Case 16: hybrid/biological 0.62; hybrid/chemical 2.64; biological/chemical 4.26. Noise band for 2-run means: about 1.17.

## 11. Island hindcast, model day 731 = 2 Oct 2026 (seeded 22 May 2025)

SRA FOI Table I, 2 Oct 2026: validated 72,103 ha, reported 122,010 ha. 1 cell = 1 ha; 356,869 mapped cane cells.

| run | n | ha ever infested, 2 Oct 2026 | ratio to SRA validated | LGUs infested | Spearman | Spearman non-seed |
|---|---|---|---|---|---|---|
| hybrid, baseline | 2 | 2,908 | 0.04 | 35 | 0.09 | 0.10 |
| none, baseline | 1 | 21,275 | 0.30 | 42 | -0.20 | -0.29 |
| none, long jump 0.2 | 1 | 44,072 | 0.61 | 50 | 0.18 | 0.22 |
| none, transmission x4 | 1 | 105,535 | 1.46 | 49 | 0.00 | -0.04 |
| none, both | 1 | 253,532 | 3.52 | 56 | 0.02 | -0.02 |
| none, transmission x2.5 | 1 | 58,136 | 0.81 | 48 | -0.04 | -0.09 |
| none, transmission x3 | 1 | 70,430 | 0.98 | 51 | -0.12 | -0.18 |

Tuning (no control): transmission x3 gives 0.98 of the validated area; log-log interpolation between x3 and x4 puts the match at about x3.05 [inference; single runs]. One free parameter fitted to one number, with control ignored: tuned, not calibrated.

## 12. Trajectory: model area ever infested (ha, island) vs dated official figures

Model value = weekly series at the week containing the date. Observed metrics and scopes differ between rows (see column); early figures reflect detection catching up and SRA attributes part of the mid-2026 rise to validation effort. Sources are in OBS in analysis.py; the manuscript cites only the 22 May 2025, 11 Sep 2025 and 2 Oct 2026 figures.

| date | model day | observed ha | metric; scope | hybrid, baseline | none, baseline | none, long jump 0.2 | none, transmission x4 | none, both | none, transmission x2.5 | none, transmission x3 |
|---|---|---|---|---|---|---|---|---|---|---|
| 2025-05-22 | 233 | 87 | affected, at declaration (= the model seed); Negros | 87 (seed) | 87 (seed) | 87 (seed) | 87 (seed) | 87 (seed) | 87 (seed) | 87 (seed) |
| 2025-06-19 | 261 | 1,505 | infested (SRA statement; date approximate); Negros Island Region | 144 | 154 | 167 | 436 | 403 | 261 | 312 |
| 2025-07-24 | 296 | 2,876 | infested (Res. 771); Negros Occidental only | 293 | 492 | 546 | 1,460 | 1,423 | 949 | 1,103 |
| 2025-09-11 | 345 | 4,923 | affected (SRA status report); Negros Island | 472 | 1,162 | 1,176 | 3,832 | 3,737 | 2,457 | 2,861 |
| 2026-06-26 | 633 | 16,019 | validated; Negros Island | 1,994 | 11,746 | 19,174 | 53,240 | 122,774 | 30,195 | 35,210 |
| 2026-08-14 | 682 | 58,374 | validated; Negros Island | 2,400 | 15,831 | 29,494 | 76,481 | 191,543 | 42,222 | 50,342 |
| 2026-10-02 | 731 | 72,103 | validated (tuning / hindcast point); Negros Island | 2,876 | 20,933 | 43,062 | 103,782 | 250,551 | 57,165 | 69,153 |

Observed multiplication between dates vs model (no control):

| period | days | observed | hybrid, baseline | none, baseline | none, long jump 0.2 | none, transmission x4 | none, both | none, transmission x2.5 | none, transmission x3 |
|---|---|---|---|---|---|---|---|---|---|
| 2025-05-22 to 2025-09-11 | 112 | 56.6x | 5.4x | 13.4x | 13.5x | 44.0x | 43.0x | 28.2x | 32.9x |
| 2025-09-11 to 2026-06-26 | 288 | 3.3x | 4.2x | 10.1x | 16.3x | 13.9x | 32.9x | 12.3x | 12.3x |
| 2026-06-26 to 2026-10-02 | 98 | 4.5x | 1.4x | 1.8x | 2.2x | 1.9x | 2.0x | 1.9x | 2.0x |

## 13. Symmetry test: the biological agent also halves in efficacy every n_half_bio applications per farm

References: chemical n_half 20 = 754,011; chemical n_half 80 = 38,761; biological without decline (n = 5) = 42,653. Hybrid keeps separate application counts for its two agents.

| n_half_bio | strategy | n | mean | lo | hi | chem | bio | vs chemical, same decline | vs biological, no decline |
|---|---|---|---|---|---|---|---|---|---|
| 20 | biological | 2 | 1,205,918 | 1,185,914 | 1,225,922 | 0 | 116,590 | 1.60 | 28.3 |
| 20 | hybrid | 2 | 384,126 | 365,398 | 402,853 | 24,526 | 86,507 | 0.51 | 9.0 |
| 80 | biological | 2 | 75,619 | 72,436 | 78,802 | 0 | 69,679 | 1.95 | 1.8 |
| 80 | hybrid | 2 | 67,082 | 63,623 | 70,541 | 250 | 67,086 | 1.73 | 1.6 |

n_half_bio 20: hybrid/biological 0.32; hybrid/chemical 0.51.

n_half_bio 80: hybrid/biological 0.89; hybrid/chemical 1.73.

