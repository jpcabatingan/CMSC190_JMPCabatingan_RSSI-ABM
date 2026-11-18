# Batch results

Folders: `main/` (batch_main_chemical, batch_main_nochemical), `oat/` (batch_oat), `wind/`
(batch_wind_direction), `sanity/` (sanity_cohort, sanity_no_pest). Purpose of each:
`docs/experiments.md`.

Batch CSVs: one row per run, no header, appended. Column order is in `COLUMNS.csv`; in
`wind/batch_wind_direction.csv` the second column is wind_dir, not oat_case. Move old CSVs out
before re-running a batch. Sanity CSVs are overwritten each run and have their own header.

OAT cases (oat_case): 0 baseline; 1-2 settlement survival 0.01 / 0.10; 3-4 neighbour
transmission x0.25 / x4; 5-6 long-jump 0 / 0.2; 7-8 parasitism 0 / 0.668; 9-10 bio:chemical
ratio 0.8 / 1.0; 11-12 level thresholds x0.5 / x2; 13-14 trigger T 2 / 3; 15-16 re-treatment
7 / 28 d; 17-18 wind bias 0.5 / 1.0.

## 8 Nov 2026: island-wide runs
The model now runs on all of Negros (includes/landscape_negros.asc). Replicates were cut to fit the
20 Nov deadline: main 10, series 10, OAT 3, wind 3, saturating 3 (repeat: facet of each batch
experiment). Batch experiments run with gui = false (no recolouring, no cars). Stylised-grid results
are in _archive_stylised/. Window A inputs (includes/*_A.*) are kept as a fallback landscape.
New file: series/lgu_<tag>.csv, one row per LGU per run: tag, seed, lgu index, LGU name,
cane cells, cells ever infested, SRA validated ha (FOI Table I, 2 Oct 2026; -1 = no row).
