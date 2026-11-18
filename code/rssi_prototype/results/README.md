# Batch results

Output CSVs of the runs reported in the manuscript and IEEE paper. All runs use the seeded start
(87 cane cells on day 233, 22 May 2025) and run 1,095 days after seeding.

| Folder | Contents |
|---|---|
| `A/main/` | Window A: main comparison (`batch_main`), chemical decline 40/80/160 (`batch_decline`), hybrid test (`batch_hybridtest`), biological decline (`batch_biodecline`), saturating damage (`batch_saturating`) |
| `A/oat/` | Window A one-at-a-time sensitivity (`batch_oat`), one run per case |
| `A/wind/` | Window A wind direction, from the north-east (135) and south-west (315) |
| `A/analysis/fable/` | Report, numbers and figures used in the paper |
| `seeded/main/` | Island main comparison |
| `seeded/series/`, `seeded/hindcast/` | Island weekly series and per-municipality counts on day 731 (2 Oct 2026) |
| `seeded/analysis/` | Island summary |
| `sanity/` | Life-table and no-pest checks, 11 Nov 2026, run before the seeded start was added |

Batch CSVs: one row per run, no header. Column order is in `COLUMNS.csv`; in `A/wind/` the second
column is wind_dir, not oat_case; in `A/main/batch_biodecline.csv` the first column (biodec20, biodec80)
gives the biological decline setting. Sanity CSVs have their own header. The seed column reads 1.0 in
every row because the seed was not logged, so single runs cannot be reproduced exactly.

OAT cases (oat_case): 0 baseline; 1-2 settlement survival 0.01 / 0.10; 3-4 neighbour transmission
x0.25 / x4; 5-6 long-jump 0 / 0.2; 7-8 parasitism 0 / 0.668; 9-10 bio:chemical ratio 0.8 / 1.0;
11-12 level thresholds x0.5 / x2; 13-14 trigger level 2 / 3; 15-16 re-treatment 7 / 28 d;
17-18 wind bias 0.5 / 1.0.

Series files `series/lgu_<tag>.csv`: one row per municipality per run: tag, seed, index, name, cane
cells, cells ever infested, SRA validated ha (FOI Table I, 2 Oct 2026; -1 = no row).
