# Seeded runs: analysis report (15 Nov 2026, updated with the extra runs)

Scope: 116 + 18 runs on window A, plus 16 island runs for the hindcast and island comparison (Cadiz, Manapla, Victorias), seeded start, 1,095 days after seeding. Every number is [computed] by `analysis.py` from the raw CSVs and listed in `numbers.md`, unless tagged otherwise. Loads are cumulative pest density (K-days) summed over cells.

## Findings that hold

1. **Any control beats none by a wide margin.** No control: 10.8 M K-days, 96.5% of cane ever infested. Every controlled configuration cuts load by 82% to 99.7%. "None" ranks last in 18 of 18 sensitivity cases.
2. **Chemical control depends on the efficacy-decline assumption more than on anything else in the main comparison.** Without decline it has the lowest load (30,013). With efficacy halving every 20 applications per farm the load is 25x higher (754,011), and every 10 applications 64x higher (1,919,710). Applications rise as efficacy falls: 1.10, 2.69 and 3.94 per farm per year (averaged over all 10,611 farms).
3. **Biological and hybrid are indistinguishable at baseline.** Pooling runs that are replicates of the same process (see Variability): hybrid 43,594 (n = 9), biological 42,653 (n = 5), ratio 1.02, difference in mean ln(load) +0.022 with SE 0.040.
4. **Three tiers are separated by more than replicate noise:** chemical without decline (0.69x the biological/hybrid pool) < biological ≈ hybrid < chemical with decline (17x and 44x the pool) < none. Gaps between tiers are 1.4x to 17x; the noise band for a 2-run comparison is 1.07x to 1.17x.
5. **For biological and hybrid, when and how often farmers act matters more than any biological parameter.** Waiting for Level 2 raises biological load 16.6x, waiting for Level 3 33x, and a 28-day re-treatment interval 25x. Wind bias, wind direction, long-jump rate and the bio:chem efficacy ratio (0.8 to 1.0) move loads by less than 1.5x.

## What the data do not support (push-back on the earlier summary)

- **"Biological beats chemical" is conditional, not a result.** In the model, biological kill is 0.9x the chemical first-application kill and never declines; only chemicals decline. Chemical with decline losing to biological follows from that asymmetry. Defensible claim: *if* fungal efficacy holds at 0.8 to 1.0 of chemical efficacy and does not decline, a biological programme matches chemical control and is not exposed to efficacy loss. Whether Metarhizium efficacy holds in the field is not in the model's evidence (lab and Egyptian field analogues only).
- **"Hybrid behaves as biological" is a consequence of the rules at baseline.** With the trigger at Level 1 and 14-day re-treatment, biological treatment keeps farm mean density far below Level 3 (0.11 K): the highest single-cell density in baseline hybrid runs is 0.107 to 0.116 K, about 0.03 K as a 2x2 farm mean [inference from max_D; farm means are not saved]. Hybrid used chemicals in 3 of 31 baseline-design runs, all where treatment is delayed.
- **Hybrid was then tested where it can act** (chemical without decline, 2 runs each; noise band 1.17x). Treating only at Level 3: chemical, biological and hybrid tie (ratios 0.89 to 1.04). Re-treating every 28 days: chemical 262,998 < hybrid 693,075 < biological 1,119,909; hybrid is 0.62x biological but 2.6x chemical. With biological efficacy held constant hybrid is never the best strategy; it sits between its two components when both are used. It becomes the best when both agents decline quickly (see the symmetry test below).
- **The saturating-damage check says nothing about pest load.** Damage does not act on the pest in the model, so those 8 runs are extra baseline replicates. It only shows the cane ranking is unchanged, and cane differs trivially (61.4 to 61.6 t/ha controlled, 57.2 none) because of the assumed slope.
- **The break-even decline rate is now located on window A.** Chemical halving every 160 or 80 applications stays within or below the biological/hybrid band (0.79x and 0.90x the pool); every 40 applications it is 3.9x worse. Chemical matches biological only if efficacy halves no faster than about every 80 applications per farm. This holds for window A at baseline settings only.

## Robustness (one-at-a-time, 1 run per cell)

Baseline ranking: biological ≈ hybrid < chemical (n_half 20) < none. Tie rule: two single runs tie when their ratio is inside the replicate noise band (1.25x for biological vs hybrid). The ranking holds in 15 of 18 cases.

| Case | Assumption | Ranking | Change |
|---|---|---|---|
| 1 | settlement survival 0.01 | chemical < biological ≈ hybrid < none | yes |
| 2-13, 15, 17, 18 | all others | biological ≈ hybrid < chemical < none | no |
| 14 | trigger at Level 3 | biological < chemical ≈ hybrid < none | yes |
| 16 | re-treatment every 28 d | biological < hybrid < chemical < none | yes |

Chemical (n_half 20) is worse than biological in 17 of 18 cases. The exception, case 1, has outbreaks about 26x smaller than baseline (no control); noise there is probably larger than the baseline band, so treat it as "ranking unclear when the pest barely establishes" [inference].

Influence (largest log10 change in load): settlement survival 1.96 (chemical), trigger level 1.81 (hybrid), re-treatment interval 1.49, parasitism 1.16 (chemical), neighbour transmission 0.86, level thresholds 0.59, long jump 0.17, bio:chem ratio 0.17, wind bias 0.06. Settlement survival and transmission have no Philippine data; they set the scale of every result.

## Variability

- With n = 2, per-configuration SDs are unreliable (CV 0.3% to 11.3% for configurations that are the same process). Pooled over all 20 replicate pairs, the SD of ln(load) is 0.020 (none), 0.045 (chemical), 0.089 (biological), 0.070 (hybrid).
- Because damage does not affect the pest and hybrid used no chemicals, the baseline has more replicates than designed: n = 5 per strategy (CV 1.7%, 5.3%, 7.5%, 7.0%), and 9 for hybrid across n_half settings.
- Two replicates are enough for the tier ranking. They are not enough to rank biological against hybrid, which is why the pooled comparison is used.
- The single case-0 run sits 0.94x to 1.07x from the pooled baseline: the size of error to expect in every OAT cell.
- Wind direction (135 vs 315 at bias 1.0): ratios 0.96 to 1.11, all inside the noise bands. No detectable effect.

## Hindcast against the SRA record (island, seeded 22 May 2025, read on 2 Oct 2026)

SRA FOI Table I: 72,103 ha validated, 122,010 ha reported.

| Run | Model ha ever infested | Ratio to validated | Spearman vs SRA by LGU (non-seed LGUs) |
|---|---|---|---|
| hybrid, baseline (n = 2) | 2,908 | 0.04 | 0.09 (0.10) |
| none, baseline | 21,275 | 0.30 | -0.20 (-0.29) |
| none, long jump 0.2 | 44,072 | 0.61 | 0.18 (0.22) |
| none, transmission x2.5 | 58,136 | 0.81 | -0.04 (-0.09) |
| none, transmission x3 | 70,430 | 0.98 | -0.12 (-0.18) |
| none, transmission x4 | 105,535 | 1.46 | 0.00 (-0.04) |
| none, both | 253,532 | 3.52 | 0.02 (-0.02) |

- At baseline the model under-predicts the infested area even with no control (0.30x validated).
- **Tuned to one point:** with no control, neighbour transmission x3 (0.30 per K per day) gives 0.98 of the validated area. This is one free parameter fitted to one number with control ignored; because the real outbreak was treated, the rate needed with control is higher. Tuned, not calibrated.
- **Weak check on dates not used for tuning** (x3 run vs official figures, mixed metrics): 19 Jun 2025 312 vs 1,505; 24 Jul 2025 1,103 vs 2,876 (Negros Occidental only); 11 Sep 2025 2,861 vs 4,923; 26 Jun 2026 35,210 vs 16,019; 14 Aug 2026 50,342 vs 58,374. The tuned curve is 0.2x to 0.6x the 2025 figures and 2.2x the June 2026 figure. The official series is not a smooth spread curve (x57 in 112 days, x3.3 in 288, x4.5 in 98); no single model curve follows it.
- No run reproduces which municipalities are worst hit (rank correlations -0.29 to 0.22). Magnitude is reachable; trajectory shape and spatial pattern are not.
- A calibration target for SP2, not a validation.

## If the biological agent also loses efficacy (window A, 2 runs each)

| Fungus halves every | Strategy | Load | Chemical apps | Biological apps |
|---|---|---|---|---|
| 20 applications | biological | 1,205,918 | 0 | 116,590 |
| 20 applications | hybrid (chemical arm n_half 20) | 384,126 | 24,526 | 86,507 |
| 80 applications | biological | 75,619 | 0 | 69,679 |
| 80 applications | hybrid (chemical arm n_half 20) | 67,082 | 250 | 67,086 |

- With the same decline rate, biological is worse than chemical: 1.60x at n_half 20 (chemical 754,011) and 1.95x at n_half 80 (chemical 38,761). Biological's lead in the main comparison comes from the no-decline assumption.
- When both agents decline every 20 applications, hybrid is the best strategy: 0.32x biological and 0.51x chemical. Hybrid counts applications of its two agents separately, so switching spreads the decline over two agents [inference from the model code]. The model has no resistance biology; do not call this resistance management.
- At n_half 80 hybrid and biological tie (ratio 0.89; noise band 1.17).

## What not to claim

1. Absolute losses, tonnes, or income. Yield differences come from an assumed damage slope.
2. That biological control is superior in the field. The model assumes non-declining fungal efficacy.
3. That hybrid is redundant or superior in general. Its rank depends on whether biological efficacy holds.
4. Island-scale conclusions. Window A is 3 of 57 LGUs; the earlier island runs used a different start and are not comparable (none/chemical 185x there vs 359x here: same sign only).
5. Calibrated spread, a resistance mechanism (decline is a swept per-farm half-life, not genetics), or effects of delivery by field access.

## Suggested Results paragraph (about 190 words)

On the three-municipality window, all control strategies reduced cumulative pest load by at least 82% relative to no control (10.8 million K-days; 96% of cane cells ever infested), and the island-wide runs gave the same ranking. Chemical control without efficacy decline gave the lowest load (3.0 x 10^4), followed by biological and hybrid control, which were indistinguishable (4.3 x 10^4; ratio 1.02). Chemical control stayed level with biological control only when its efficacy halved no faster than about every 80 applications per farm; at 20 and 10 applications its load was 17 and 44 times higher. This advantage of biological control depends on the assumption that its efficacy does not decline: under an equal decline it performed 1.6 to 1.95 times worse than chemical control, and the hybrid strategy then gave the lowest load. The baseline ranking held in 15 of 18 one-at-a-time cases; outcomes were most sensitive to crawler settlement survival, the infestation level that triggers treatment, and the re-treatment interval. Seeded from the 22 May 2025 declaration, the baseline model under-predicted the SRA validated area of 2 October 2026; a threefold neighbour-transmission rate matched the area but not its trajectory or spatial pattern.

## Defence-slide headlines

1. Chemical control is best only if it never loses efficacy.
2. Acting at Level 1 matters more than which agent is used.
3. Hybrid adds nothing while biological efficacy holds; it is the best strategy when both agents lose efficacy quickly.
4. Baseline spread is too slow against the 2 Oct 2026 SRA record; transmission x3 matches the area but not the trajectory or the map.
