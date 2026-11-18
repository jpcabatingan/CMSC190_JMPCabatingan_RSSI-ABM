/**
 * Name:        rssi_prototype
 * Description: Initial prototype of the RSSI (Pulvinaria tenuivalvata) management model on a
 *              Negros Island cane landscape. Implements docs/model-odd-d.md with the values in
 *              docs/parameter-table.md (decisions D1-D16, 3 Nov 2026).
 *
 *              8 Nov 2026: moved from the stylised 50 x 50 grid onto a real landscape, window A
 *              (Cadiz City, Manapla, City of Victorias, Negros Occidental): 1-ha cells in UTM 51N,
 *              cane from Di Tommaso et al. 2024 (GEDI-S2 sugarcane map, 2019-2022), municipal
 *              boundaries from OCHA COD-AB (valid 13 Feb 2025), roads from OpenStreetMap (ODbL,
 *              downloaded 8 Nov 2026). Inputs built by SP/GIS/scripts/03-04; see includes/.
 *              8 Nov 2026 (later): extended to ALL of Negros Island (57 LGUs, both provinces; J.M.):
 *              1,327 x 2,224 cells of 1 ha, 356,869 cane cells; roads OSM trunk/primary/secondary island-wide.
 *              Inputs: SP/GIS/scripts/06-07. Island-wide RSSI validated area by LGU (SRA FOI Table I)
 *              is carried on the LGU layer for display and for the per-LGU validation output.
 *              Non-cane cells hold no crop and no pest. Farms are 2 x 2-cell blocks of cane cells
 *              (1-4 ha) [assumption: no farm-boundary data were found].
 *              13 Nov 2026: seeded start (seed_mode "may2025"): 87 cane cells seeded on day 233 (22 May 2025)
 *              across the 11 LGUs first reported; replaces the random start used before 13 Nov.
 *
 *              Shows mechanics only: daily clock, crop calendar, stage-structured pest per cell
 *              (63 daily age classes), parasitism, neighbour (optionally wind-biased) and long-jump spread, cumulative
 *              damage, weekly farmer monitoring, four strategies with a common trigger.
 *              Outputs are relative comparisons between strategies under stated assumptions, not forecasts.
 *
 *              Status: written 3 Nov 2026; actions given "()" for the 2026 GAML grammar (4 Nov). 4 Nov: farm-level
 *              outlines, treatment markers, road strip, compare_strategies experiment added. 4 Nov 01:00:
 *              sanity_cohort, sanity_no_pest, batch_main_*, batch_oat added. 4 Nov 11:00: sanity checks
 *              passed; main batch dry run (320 runs). Code audited against the ODD+D (notes/code-audit-
 *              2026-11-04.md): saturating damage switch, max density output, weekly series added.
 *
 * Author: Joanne Maryz Cabatingan
 */

model rssi_prototype

global {

    // ---------------- Simulation settings ----------------
    int horizon <- 1095;              // 3 crop years, 1 tick = 1 day (D5)
    int farm_side <- 2;               // farms = cane cells in each 2 x 2-cell block (1-4 ha) [assumption]

    // ---------------- Landscape: all of Negros Island (8 Nov 2026) ----------------
    // landscape_negros.asc: value = lgu index x 10 + cane fraction of the 1-ha cell; -9999 outside Negros
    // (LGU index 1-57: SP/GIS/derived/negros/lgu_index.csv, ordered by PSGC code)
    string landscape <- "negros";     // "negros" (all of Negros) or "A" (Cadiz, Manapla, Victorias window; sensitivity runs, 13 Nov)
    grid_file landscape_file <- grid_file("../includes/landscape_" + landscape + ".asc");
    file lgu_shp <- shape_file("../includes/lgus_" + landscape + ".shp");
    file road_shp <- shape_file("../includes/roads_" + landscape + ".shp");
    geometry shape <- envelope(landscape_file);
    float cane_min <- 0.5;            // a cell is cane if >= this share of its 10 m pixels is mapped cane [assumption]
    int gap_cells <- 0;               // 0 = neighbour spread only between adjacent cane cells; 1 = may jump one non-cane cell [assumption]
    list<cell> cane_cells;
    // Event-driven bookkeeping (8 Nov 2026, island scale): same results as the full daily scans,
    // only cells where something happens are visited. Verified identical on window A (seed 42).
    map<int, list<cell>> cal;         // cane cells by plant_day mod crop_length: who plants/harvests today
    list<cell> roadside_cells;
    list<cell> marked_cells;
    list<cell> act;                   // cells with D > 0 at the end of today
    list<cell> act_prev;
    int n_ever <- 0;
    list<cell> timed_cells;           // cells with a treatment still running
    list<farmer> to_treat;            // farmers at or above the trigger level since the last weekly check
    // shared colour objects (one per colour, not one per cell: memory at island scale)
    rgb c_out <- rgb(232, 238, 243);
    rgb c_land <- rgb(236, 232, 222);
    rgb c_tan <- rgb(222, 205, 160);
    rgb c_green <- rgb(120, 180, 90);
    bool gui <- true;                 // false in batch experiments: no recolouring, no cars
    int n_cane <- 1;
    bool show_lgu_names <- true;      // display only

    // ---------------- Pest life history (A1 Guerrero et al. 2026) ----------------
    // Daily age classes 0..62. Integer stage spans inside A1 Table 2 ranges, sum 63 d:
    // egg 0-1 | I1 2-7 | I2 8-13 | I3 14-18 | young adult 19-31 | mature 32-48 | senescing 49-62
    int n_ages <- 63;
    int nymph_first <- 2;
    int adult_first <- 19;
    int young_adult_last <- 31;
    int lay_first <- 38;              // adult age 19 (A1 Table 3; age counted from the moult to adult, A1 Methods p.176 'Oviposition in RSSI')
    int lay_last <- 53;               // adult age 34
    float eggs_per_day <- 197.0 / 16.0;   // 197 eggs over 16 laying days [computed]

    // ---------------- Pest assumptions ----------------
    float s_settle <- 0.05;           // egg -> settled crawler, lumped [assumption], swept 0.01-0.10
    float K <- 1.0;                   // capacity; densities are in units of K

    // ---------------- Spread [assumption; Burgonio: no rate data] ----------------
    float beta_nb <- 0.1;             // neighbour transmission, PLACEHOLDER: in a Python mirror (not GAMA)
                                      // no-control spread saturated the landscape in about 1 year
    float p_jump <- 0.05;             // long jump via planting material, per planting x infested share
    float inoculum <- 0.01;           // density added per spread event (K)
    float infested_th <- 0.01;        // cell counts as infested at D >= this (K)
    // Wind bias on neighbour spread (Doc Arnejo's suggestion, 7 Oct). No RSSI wind data:
    // a swept sensitivity switch, not a calibrated wind model. 0 = no wind effect.
    float wind_bias <- 0.0;           // 0-1; weight = max(0, 1 + wind_bias x cos(angle to wind))
    float wind_dir <- 135.0;          // degrees the wind blows TOWARD: 0 = right (+x), 90 = down (+y); 135 = toward SW, i.e. wind FROM the NE (PAGASA Dumaguete normals, A22). Only acts when wind_bias > 0
    // Roadside exposure (Burgonio: roadside plantations worst hit; likely traffic and cane
    // transport, not wind). Real roads (OSM trunk, primary, secondary); cane cells within
    // roadside_m of one get their long-jump probability multiplied. 1 = no effect [assumption], swept {1, 3}.
    bool show_road <- true;
    float roadside_m <- 100.0;        // buffer around trunk/primary/secondary roads [assumption]
    float road_mult <- 1.0;
    int n_cars <- 150;                 // cosmetic only: vehicles moving on the real road network, no effect on the model
    graph road_net;
    bool show_legend <- true;         // display only: legend box, bottom-right of each map
    float road_daily <- 0.0;          // daily chance a roadside cell is inoculated x infested share
                                      // (traffic-borne crawlers) [assumption]; 0 = no effect

    // ---------------- Crop (SRA FOI; A9) ----------------
    float Y0 <- 61.64;                // TC/ha, pre-outbreak Negros mean [computed, FOI Table II]
    float recovery <- 1.61;           // LKg/TC, FOI Table III NEGROS CY2025-26
    int planting_window <- 243;       // 1 Oct - 31 May (A9), uniform [assumption]
    int crop_length <- 365;           // harvest at 12 months, then ratoon [assumption]
    float h_harv <- 0.9;              // share of pests removed at harvest [assumption]

    // ---------------- Damage (Alves et al. 2026 structure) ----------------
    float L_ref <- 0.3;               // loss after a full season at capacity [assumption], 0.1-0.5
    float L_max <- 1.0;               // maximum loss, swept 0.5-1.0
    bool dmg_saturating <- false;     // DM1 robustness check: loss = L_max x (1 - exp(-beta x pest_days))

    // ---------------- Control ----------------
    float k_chem_nymph <- 0.743;      // median % reduction, 1st application, 6 compounds, A4 p.3 [computed]; A4 gives no spray-to-count interval: applying it over 14 d is the model's assumption (D9)
    float k_chem_adult <- 0.549;
    float r_bc <- 0.9;                // Metarhizium = r_bc x chemical, 0.80-1.00 (A6 field; A5 lab)
    int treat_days <- 14;             // both agents act over 14 days (D9)
    int n_half <- 20;                 // chemical kill halves after n_half uses; 0 = no decline (D16)
    int n_half_bio <- 0;              // same decline for the biological agent; 0 = none (baseline). Symmetry test only (15 Nov 2026) [assumption]
    int r_interval <- 14;             // re-treatment interval (A4 p.3 example rotation)

    // ---------------- Natural enemies (A2, A3) ----------------
    float p_para <- 0.3485;           // {0, 0.3485, 0.668} (D15)
    int young_adult_days <- 13;
    float h_para <- 0.0;              // daily hazard, set in init

    // ---------------- Farmer (D8, D9) ----------------
    string strategy <- "hybrid";      // none / chemical / biological / hybrid
    int trigger_T <- 1;               // same trigger for every strategy
    // Level thresholds (13 Nov 2026, J.M. OK): PGNO Provincial Biosecurity Task Force infestation levels,
    // Level 1 Low 1-5, Level 2 Moderate 6-10, Level 3 High >10 insects per leaf (rssi.mis-negros.com,
    // data as of 8 Nov 2026). The model's density is relative to K, so K must be given in insects per
    // leaf: k_leaf = 100 [assumption], chosen so that Level 1 (1 per leaf) equals the infested
    // threshold infested_th = 0.01 K. Swept by OAT cases 11 (k_leaf 200 = thresholds x0.5) and
    // 12 (k_leaf 50 = thresholds x2). Applied to the farm mean density; PGNO counts one leaf per lot.
    float k_leaf <- 100.0;            // insects per leaf at carrying capacity K [assumption]
    float th1 <- 1.0 / k_leaf;        // Level 1: >= 1 insect per leaf
    float th2 <- 6.0 / k_leaf;        // Level 2: >= 6
    float th3 <- 11.0 / k_leaf;       // Level 3: > 10
    int monitor_every <- 7;

    // ---------------- Initial condition [assumption] ----------------
    float init_share <- 0.01;         // share of cells seeded
    // Seeded start (13 Nov 2026, after the mock defence): "random" = init_share of cane cells anywhere;
    // "may2025" = the SRA's first declared infestation, 22 May 2025: 87.04 ha in 11 Negros Occidental
    // LGUs (Inquirer 2025; Watchmen Daily Journal 27 May 2025; PNA 26 May 2025, as SRA figures).
    // No per-LGU split was published for that date, so the 87 cells are split equally across the 11
    // LGUs and placed at random within each LGU's cane [assumption].
    string seed_mode <- "random";
    int seed_ha <- 87;
    // Model day 0 = 1 Oct (start of the planting window, A9), so 22 May = day 233 [computed]. The
    // landscape runs pest-free until then, so all cane planted 1 Oct-21 May is standing when the
    // seeds go in. Hindcast end, 2 Oct 2026 = day 233 + 498 = 731.
    int seed_day <- 233;
    int n_seeded <- 0;               // cells actually seeded (check: = seed_ha)
    list<string> seed_lgu_names <- ["Silay City", "City of Talisay", "City of Victorias", "Cadiz City",
        "Enrique B. Magalona (Saravia)", "Manapla", "Toboso", "Bago City", "La Carlota City", "La Castellana", "Murcia"];
    float init_density <- 0.05;       // K, as third instars

    // ---------------- Sanity check and one-at-a-time sensitivity ----------------
    bool cohort_test <- false;        // true only in experiment sanity_cohort
    float eggs_laid_total <- 0.0;     // eggs laid by the single test cohort (per female)
    int cohort_death_day <- -1;       // cycle at which the test cohort is gone
    float cohort_after_young_adult <- -1.0;  // survivors at the end of the young-adult stage
    int oat_case <- 0;                // 0 = baseline; 1-18 = one parameter low/high (see apply_oat)
    int harvest_count <- 0;
    string sanity_tag <- "";          // set by sanity experiments; non-empty = write a summary row

    // ---------------- Outputs ----------------
    float share_infested <- 0.0;
    float share_ever <- 0.0;
    float mean_density <- 0.0;
    float cum_kdays <- 0.0;
    int apps_chem <- 0;
    int apps_bio <- 0;
    float cane_total <- 0.0;          // t cane harvested, whole landscape
    float sugar_total <- 0.0;         // LKg
    float inf_share_prev <- 0.0;      // infested share at the end of the previous day
    float run_seed <- 0.0;            // this simulation's random seed, logged in batch rows
    float max_D <- 0.0;               // highest cell density reached in the run (ODD pattern 1: should stay near K)
    string batch_file <- "";          // batch experiments: end-of-run row appended here (relative to models/)
    string batch_label <- "";         // first column of that row (main, oat, wind, sat, long)
    bool save_series <- false;        // write a weekly time series row to ../results/series/
    string series_tag <- "";
    string series_dir <- "../results/series/";   // seeded runs: ../results/seeded/series/
    int lgu_day <- -1;                // extra per-LGU snapshot on this day (seeded: 731 = 2 Oct 2026)

    init {
        run_seed <- seed;
        do apply_oat;
        h_para <- 1.0 - (1.0 - p_para) ^ (1.0 / young_adult_days);

        // landscape: decode grid values; only cane cells take part in the model
        ask cell {
            if (grid_value > 0.0) {
                lgu_id <- int(grid_value / 10.0);
                cane_frac <- grid_value - lgu_id * 10.0;
                is_cane <- cane_frac >= cane_min;
            }
            color <- lgu_id = 0 ? c_out : c_land;   // outside Negros LGUs / non-cane land
        }
        cane_cells <- cell where (each.is_cane);
        n_cane <- max([1, length(cane_cells)]);
        ask cane_cells {
            plant_day <- rnd(planting_window - 1);
            nbs <- (gap_cells = 0 ? neighbors : (self neighbors_at (1 + gap_cells))) where (each.is_cane);
            color <- c_tan;
        }
        create lgu from: lgu_shp with: [name::string(read("adm3_name")), idx::int(read("idx")), valid_ha::float(read("valid_ha"))];
        // lgus_A has no idx field; its records are in grid-code order (1 Cadiz, 2 Manapla, 3 Victorias:
        // checked 13 Nov by matching polygon and cell-code centroids)
        if (lgu all_match (each.idx = 0)) { loop i from: 0 to: length(lgu) - 1 { lgu[i].idx <- i + 1; } }
        create road from: road_shp with: [hwy::string(read("hwy"))];
        if (cohort_test) {
            // one cell, one cohort of density 1, no recruitment, spread, control or harvest
            ask cane_cells { plant_day <- 999999; }
            ask first(cane_cells) { plant_day <- 0; }
        } else {
            if (seed_mode = "random") {
                ask (int(init_share * n_cane)) among cane_cells {
                    seeded <- true;
                }
            }
        }
        ask road where (each.hwy in ["trunk", "primary", "secondary"]) {
            ask (cell overlapping (self + roadside_m)) where (each.is_cane) { roadside <- true; }
        }
        roadside_cells <- cane_cells where (each.roadside);
        loop i from: 0 to: length(cane_cells) - 1 { cane_cells[i].ci <- i; }
        cal <- cane_cells group_by (each.plant_day mod crop_length);
        // cosmetic traffic on the real road network (no model state)
        road_net <- as_edge_graph(road);
        // Vehicles removed from the map (J.M., 8 Nov 2026); uncomment to show them again.
        // create car number: (gui ? n_cars : 0) {
        //     location <- any_location_in(one_of(road));
        //     col <- one_of([rgb(250, 250, 250), rgb(255, 210, 0), rgb(40, 40, 40), rgb(220, 30, 30)]);
        // }

        map<point, list<cell>> blocks <- cane_cells group_by ({each.grid_x div farm_side, each.grid_y div farm_side});
        loop b over: blocks.values {
            create farmer {
                my_cells <- b;
                ask my_cells { owner <- myself; }
                shape <- union(my_cells collect each.shape);
                location <- shape.location;
                lgu_id <- first(my_cells).lgu_id;
            }
        }
    }

    // One-at-a-time sensitivity: each case moves one assumption to its low or high value
    // (docs/parameter-table.md); everything else stays at baseline.
    action apply_oat() {
        switch oat_case {
            match 1  { s_settle <- 0.01; }
            match 2  { s_settle <- 0.10; }
            match 3  { beta_nb <- beta_nb * 0.25; }
            match 4  { beta_nb <- beta_nb * 4.0; }
            match 5  { p_jump <- 0.0; }
            match 6  { p_jump <- 0.2; }
            match 7  { p_para <- 0.0; }
            match 8  { p_para <- 0.668; }
            match 9  { r_bc <- 0.8; }
            match 10 { r_bc <- 1.0; }
            match 11 { th1 <- th1 * 0.5; th2 <- th2 * 0.5; th3 <- th3 * 0.5; }   // = k_leaf 200
            match 12 { th1 <- th1 * 2.0; th2 <- th2 * 2.0; th3 <- th3 * 2.0; }   // = k_leaf 50
            match 13 { trigger_T <- 2; }
            match 14 { trigger_T <- 3; }
            match 15 { r_interval <- 7; }
            match 16 { r_interval <- 28; }
            match 17 { wind_bias <- 0.5; }
            match 18 { wind_bias <- 1.0; }
        }
    }

    reflex track_cohort when: cohort_test {
        float alive <- empty(first(cane_cells).pop) ? 0.0 : sum(first(cane_cells).pop);
        if (cycle = young_adult_last) { cohort_after_young_adult <- alive; }
        if (cohort_death_day < 0 and cycle > 0 and alive <= 0.0) { cohort_death_day <- cycle; }
    }

    // One global reflex fixes the ODD schedule order (steps 1-9).
    // Seeded start: seed_ha planted cells, split equally across the 11 LGUs of 22 May 2025
    action seed_may2025() {
        list<int> ids <- (lgu where (each.name in seed_lgu_names)) collect each.idx;
        list<cell> pool <- cane_cells where (each.planted and (ids contains each.lgu_id));
        list<cell> chosen <- [];
        int per <- int(seed_ha / max([1, length(ids)]));
        loop i over: ids { chosen <- chosen + (per among (pool where (each.lgu_id = i))); }
        chosen <- chosen + ((seed_ha - length(chosen)) among (pool - chosen));
        ask chosen { do add_pest(a: 14, x: init_density); }
        n_seeded <- length(chosen);
    }

    reflex daily_step {
        inf_share_prev <- share_infested;

        // 1 Calendar and crop: only the cells whose planting or harvest day is today
        list<cell> today <- cal[cycle mod crop_length];
        if (today != nil) { ask today { do calendar_crop; } }
        if (seed_mode = "may2025" and cycle = seed_day) { do seed_may2025; }

        // 2-4 Pest development, reproduction, mortality (occupied cells only)
        ask cane_cells where (each.occ) {
            do age_pests;
            do reproduce;
            do mortality;
            do update_density;
        }
        // treatment timers run down every day in every treated cell, occupied or not (fix 4 Nov)
        ask timed_cells {
            if (chem_left > 0) { chem_left <- chem_left - 1; }
            if (bio_left > 0) { bio_left <- bio_left - 1; }
        }
        timed_cells <- timed_cells where (each.chem_left > 0 or each.bio_left > 0);

        // 5 Spread (synchronous: reads the snapshot D_prev). Only planted cane cells next to a
        // source can receive; they are visited in cane_cells order, as in the full scan.
        list<cell> dp <- cane_cells where (each.D > 0.0 or each.D_prev > 0.0);
        ask dp { D_prev <- D; }
        ask dp where (each.D_prev > 0.0) {
            ask nbs where (each.planted and not each.marked) { marked <- true; marked_cells <+ self; }
        }
        ask (marked_cells sort_by (each.ci)) { marked <- false; do receive_spread; }
        marked_cells <- [];
        if (road_daily > 0.0) {
            ask roadside_cells where (each.planted) {
                if (flip(min([1.0, road_daily * inf_share_prev]))) {
                    do add_pest(a: nymph_first, x: inoculum);
                    D <- D + inoculum;
                }
            }
        }

        // 6 Damage (cells with D = 0 would add nothing)
        act <- cane_cells where (each.D > 0.0);
        ask act where (each.planted) { do damage; }

        // 7 Farmer monitoring, every 7 days
        // (farmers with no pest on any cell and level 0 stay at level 0: skipping them changes nothing)
        if (cycle mod monitor_every = 0) {
            ask remove_duplicates((act collect each.owner) + (farmer where (each.level > 0))) { do monitor; }
            to_treat <- (strategy = "none") ? [] : (farmer where (each.level >= trigger_T));
        }

        // 8 Control (order of treatments does not matter: no random draws, counts are additive)
        ask to_treat { do control; }

        // 9 Record
        share_infested <- (act count (each.D >= infested_th)) / n_cane;
        share_ever <- n_ever / n_cane;
        mean_density <- (act sum_of (each.D)) / n_cane;
        max_D <- max([max_D, empty(act) ? 0.0 : act max_of (each.D)]);
        if (gui) {
            ask remove_duplicates(act + act_prev + (today = nil ? [] : today)) { do recolour; }
        }
        act_prev <- act;
    }

    // One row per run at the last day, written by the model itself so that it also works in headless
    // batch runs (experiment-level reflexes did not run headless in GAMA 2026.04.9, 8 Nov 2026).
    // Columns: results/COLUMNS.csv (for label "wind" the second column is wind_dir, not oat_case).
    reflex save_batch_row when: batch_file != "" and cycle = horizon - 1 {
        save (batch_label + "," + (batch_label = "wind" ? wind_dir : oat_case) + "," + strategy + "," + n_half + "," + run_seed + ","
            + share_ever + "," + share_infested + "," + cum_kdays + "," + apps_chem + "," + apps_bio + "," + cane_total + ","
            + harvest_count + "," + max_D + "\n") to: batch_file format: "text" rewrite: false;
    }

    // Weekly time series (for the time-scale check, ODD pattern 3). One row per run per week.
    reflex save_weekly when: save_series and cycle mod 7 = 0 {
        save (series_tag + "," + strategy + "," + n_half + "," + oat_case + "," + run_seed + "," + cycle + ","
            + share_infested + "," + share_ever + "," + mean_density + "," + apps_chem + "," + apps_bio + ","
            + (farmer count (each.level = 1)) + "," + (farmer count (each.level = 2)) + "," + (farmer count (each.level = 3)) + "," + length(farmer) + "\n")
            to: series_dir + "series_" + series_tag + ".csv" format: "text" rewrite: false;
    }

    // Per-LGU end-of-run output (time-scale run only): cane cells and cells ever infested per LGU, with
    // SRA's validated RSSI area for the same LGU (FOI Table I, 2 Oct 2026) for a spatial comparison.
    reflex save_lgu when: save_series and (cycle = horizon - 1 or cycle = lgu_day) {
        string rows <- "";
        loop l over: lgu {
            list<cell> cc <- cane_cells where (each.lgu_id = l.idx);
            rows <- rows + series_tag + "," + run_seed + "," + l.idx + "," + l.name + "," + length(cc) + ","
                + (cc count (each.ever_infested)) + "," + l.valid_ha + "\n";
        }
        save rows to: series_dir + "lgu_" + series_tag + (cycle = lgu_day ? "_d" + cycle : "") + ".csv" format: "text" rewrite: false;
    }

    // Sanity experiments write one summary row at the end to ../results/sanity/<tag>.csv
    reflex save_sanity when: sanity_tag != "" and cycle = horizon - 1 {
        if (sanity_tag = "cohort") {
            save ("p_para,eggs_per_female,cohort_gone_after_day,left_after_young_adult\n" + p_para + "," + eggs_laid_total + "," + cohort_death_day + "," + cohort_after_young_adult)
                to: "../results/sanity/sanity_cohort_p" + p_para + ".csv" format: "text" rewrite: true;
        } else {
            save ("share_ever,apps,cane_per_harvested_ha,harvests\n" + share_ever + "," + (apps_chem + apps_bio) + "," + (harvest_count = 0 ? 0.0 : cane_total / harvest_count) + "," + harvest_count)
                to: "../results/sanity/sanity_no_pest.csv" format: "text" rewrite: true;
        }
    }

    reflex stop when: cycle >= horizon {
        do pause;
    }
}

grid cell file: landscape_file neighbors: 8 use_individual_shapes: false use_regular_agents: false schedules: [] {
    int lgu_id <- 0;                  // 1 Cadiz City, 2 Manapla, 3 City of Victorias, 0 outside
    float cane_frac <- 0.0;
    bool is_cane <- false;
    list<cell> nbs;                   // cane neighbours used for spread (gap_cells)
    int ci <- -1;                     // position in cane_cells (keeps visiting order)
    bool occ <- false;                // = sum(pop) > 0
    bool marked <- false;
    int plant_day <- 0;
    bool planted <- false;
    bool seeded <- false;
    bool roadside <- false;
    int crop_age <- 0;
    list<float> pop;                  // density by daily age class (K)
    float D <- 0.0;                   // density of ages 2-62 (feeding, infesting stages)
    float D_prev <- 0.0;
    bool ever_infested <- false;
    float pest_days <- 0.0;           // this crop
    float cum_kdays_cell <- 0.0;      // whole run
    float loss <- 0.0;
    farmer owner;

    // active treatments: daily hazards and days left
    int chem_left <- 0;
    float hc_nymph <- 0.0;
    float hc_adult <- 0.0;
    int bio_left <- 0;
    float hb_nymph <- 0.0;
    float hb_adult <- 0.0;

    action calendar_crop() {
        // called only on days with cycle mod crop_length = plant_day mod crop_length
        if (not planted) {
            if (cycle = plant_day) {
                planted <- true;
                crop_age <- 0;
                if (seeded) { do add_pest(a: 14, x: init_density); }
                if (cohort_test) { do add_pest(a: 0, x: 1.0); }
                if (flip(min([1.0, p_jump * inf_share_prev * (roadside ? road_mult : 1.0)]))) { do add_pest(a: nymph_first, x: inoculum); }
                occ <- (not empty(pop)) and sum(pop) > 0.0;
            }
        } else {
            if (cycle > plant_day) {
                float cane <- Y0 * (1.0 - loss);
                cane_total <- cane_total + cane;
                harvest_count <- harvest_count + 1;
                sugar_total <- sugar_total + cane * recovery;
                if (not empty(pop)) { loop a from: 0 to: n_ages - 1 { pop[a] <- pop[a] * (1.0 - h_harv); } }
                crop_age <- 0;
                loss <- 0.0;
                pest_days <- 0.0;
                // ratoon start: long jump via planting material
                if (flip(min([1.0, p_jump * inf_share_prev * (roadside ? road_mult : 1.0)]))) { do add_pest(a: nymph_first, x: inoculum); }
                occ <- (not empty(pop)) and sum(pop) > 0.0;
            }
        }
    }

    // add density x to age class a; the 63-class age list exists only in cells that have had pests
    action add_pest(int a, float x) {
        if (empty(pop)) { pop <- list_with(n_ages, 0.0); }
        pop[a] <- pop[a] + x;
        occ <- true;
    }

    action age_pests() {
        // shift every class one day older; age 62 dies at the end of day 63
        pop <- [0.0] + copy_between(pop, 0, n_ages - 1);
    }

    action reproduce() {
        float layers <- sum(copy_between(pop, lay_first, lay_last + 1));
        float dens <- sum(copy_between(pop, nymph_first, n_ages));
        if (cohort_test) {
            eggs_laid_total <- eggs_laid_total + layers * eggs_per_day;
        } else {
            float settled <- layers * eggs_per_day * s_settle * max([0.0, 1.0 - dens / K]);
            pop[0] <- pop[0] + settled;
        }
    }

    action mortality() {
        // parasitism on young adults (D3)
        if (h_para > 0.0) {
            loop a from: adult_first to: young_adult_last { pop[a] <- pop[a] * (1.0 - h_para); }
        }
        // treatments: same stages and time course for both agents (D9)
        if (chem_left > 0 or bio_left > 0) {
            float sn <- (chem_left > 0 ? (1.0 - hc_nymph) : 1.0) * (bio_left > 0 ? (1.0 - hb_nymph) : 1.0);
            float sa <- (chem_left > 0 ? (1.0 - hc_adult) : 1.0) * (bio_left > 0 ? (1.0 - hb_adult) : 1.0);
            loop a from: nymph_first to: adult_first - 1 { pop[a] <- pop[a] * sn; }
            loop a from: adult_first to: n_ages - 1 { pop[a] <- pop[a] * sa; }
        }
    }

    action update_density() {
        D <- sum(copy_between(pop, nymph_first, n_ages));
        if (D < 1.0e-9 and sum(pop) < 1.0e-9) {
            pop <- [];                // release the age list (memory); recreated by add_pest when needed
            D <- 0.0;
        }
        occ <- (not empty(pop)) and sum(pop) > 0.0;
    }

    action receive_spread() {
        // weighted by wind: neighbours upwind of this cell push more crawlers here
        float pressure <- 0.0;
        loop nb over: nbs {
            if (nb.D_prev > 0.0) {
                float ang <- atan2(float(grid_y - nb.grid_y), float(grid_x - nb.grid_x));  // degrees
                pressure <- pressure + nb.D_prev * max([0.0, 1.0 + wind_bias * cos(ang - wind_dir)]);
            }
        }
        if (pressure > 0.0 and flip(1.0 - exp(-beta_nb * pressure))) {
            do add_pest(a: nymph_first, x: inoculum);
            D <- D + inoculum;
        }
    }

    action damage() {
        pest_days <- pest_days + D;
        cum_kdays_cell <- cum_kdays_cell + D;
        cum_kdays <- cum_kdays + D;
        loss <- dmg_saturating ? L_max * (1.0 - exp(-(L_ref / crop_length) * pest_days))
                               : min([L_max, (L_ref / crop_length) * pest_days]);
        if (D >= infested_th and not ever_infested) { ever_infested <- true; n_ever <- n_ever + 1; }
    }

    action recolour() {
        if (not planted) {
            color <- c_tan;
        } else if (D < infested_th) {
            color <- c_green;
        } else {
            int r <- int(150 + 105 * min([1.0, D]));
            int g <- int(150 * (1.0 - min([1.0, D])));
            color <- rgb(r, g, 40);
        }
    }

    aspect default {
        draw shape color: color;
    }

    // active treatment markers drawn over the map: magenta = chemical, cyan = Metarhizium
    aspect treat {
        if (chem_left > 0) {
            draw circle(shape.width * 0.32) color: rgb(255, 0, 200) border: #white;
        } else if (bio_left > 0) {
            draw circle(shape.width * 0.32) color: rgb(0, 200, 255) border: #white;
        }
    }
}

species farmer schedules: [] {
    list<cell> my_cells;
    int lgu_id <- 0;
    int level <- 0;
    int last_treat <- -10000;         // cycle of the last treatment (same rule as the old days_since_treat counter)
    int n_chem <- 0;
    int n_bio <- 0;

    // Farm infestation level (none drawn at level 0): the farm is filled in its level colour and a
    // marker of fixed screen size is drawn on it, so levels stay visible when the whole island is shown.
    // Level 1 yellow, 2 orange, 3 black. A smaller dot inside shows an active treatment:
    // magenta = chemical, cyan = biological.
    aspect level_outline {
        if (level > 0) {
            rgb c <- level = 1 ? rgb(255, 215, 0) : (level = 2 ? rgb(255, 130, 0) : rgb(20, 20, 20));
            draw shape color: c;
            draw circle(3.5 #px) at: location color: c border: (level = 3 ? #white : rgb(60, 60, 60));
            if (my_cells one_matches (each.chem_left > 0)) {
                draw circle(1.8 #px) at: location color: rgb(255, 0, 200);
            } else if (my_cells one_matches (each.bio_left > 0)) {
                draw circle(1.8 #px) at: location color: rgb(0, 200, 255);
            }
        }
    }

    action monitor() {
        list<cell> live <- my_cells where (each.planted);
        float dmean <- empty(live) ? 0.0 : live mean_of (each.D);
        level <- dmean >= th3 ? 3 : (dmean >= th2 ? 2 : (dmean >= th1 ? 1 : 0));
    }

    action control() {
        if (strategy != "none" and level >= trigger_T and (cycle - last_treat) >= r_interval) {
            do treat;
        }
    }

    action treat() {
        string agent_used <- strategy;
        if (strategy = "hybrid") {
            agent_used <- level <= 2 ? "biological" : "chemical";
        }

        if (agent_used = "chemical") {
            n_chem <- n_chem + 1;
            apps_chem <- apps_chem + 1;
            float factor <- n_half > 0 ? 0.5 ^ ((n_chem - 1) / float(n_half)) : 1.0;
            float hn <- 1.0 - (1.0 - k_chem_nymph * factor) ^ (1.0 / treat_days);
            float ha <- 1.0 - (1.0 - k_chem_adult * factor) ^ (1.0 / treat_days);
            ask my_cells where (each.planted) {
                if (chem_left = 0 and bio_left = 0) { timed_cells <+ self; }
                chem_left <- treat_days;
                hc_nymph <- hn;
                hc_adult <- ha;
            }
        } else {
            n_bio <- n_bio + 1;
            apps_bio <- apps_bio + 1;
            float fb <- n_half_bio > 0 ? 0.5 ^ ((n_bio - 1) / float(n_half_bio)) : 1.0;
            float hn <- 1.0 - (1.0 - r_bc * k_chem_nymph * fb) ^ (1.0 / treat_days);
            float ha <- 1.0 - (1.0 - r_bc * k_chem_adult * fb) ^ (1.0 / treat_days);
            ask my_cells where (each.planted) {
                if (chem_left = 0 and bio_left = 0) { timed_cells <+ self; }
                bio_left <- treat_days;
                hb_nymph <- hn;
                hb_adult <- ha;
            }
        }
        last_treat <- cycle;
    }
}

// Cosmetic only: cars moving along the road to illustrate traffic-borne spread.
// They carry no model state and do not infest cells; road_daily does that.
species car skills: [moving] {
    rgb col <- #white;

    reflex drive when: show_road and gui {
        do wander on: road_net speed: 1500.0;  // metres per tick; cosmetic only
    }

    aspect default {
        if (show_road) {
            draw circle(450) color: col border: #black;
        }
    }
}

species lgu schedules: [] {
    int idx;                          // = lgu_id of its cells
    float valid_ha;                   // SRA FOI Table I, validated RSSI area, 2 Oct 2026 (-1 = no row)
    aspect default {
        draw shape.contour + 60 color: rgb(60, 60, 60);
        if (show_lgu_names) { draw name at: shape.location color: rgb(40, 40, 40) font: font("Arial", 9, #plain); }
    }
}

species road schedules: [] {
    string hwy;
    aspect default {
        if (show_road) {
            draw shape + (hwy in ["trunk", "primary"] ? 140 : 80) color: rgb(95, 95, 95);
        }
    }
}

experiment prototype type: gui {
    parameter "Strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"] category: "Farmer";
    parameter "Trigger level T" var: trigger_T min: 1 max: 3 category: "Farmer";
    parameter "Level 1 threshold (K)" var: th1 min: 0.001 max: 1.0 category: "Farmer";
    parameter "Level 2 threshold (K)" var: th2 min: 0.001 max: 1.0 category: "Farmer";
    parameter "Level 3 threshold (K)" var: th3 min: 0.001 max: 1.0 category: "Farmer";

    parameter "Settlement survival" var: s_settle min: 0.0 max: 1.0 category: "Pest";
    parameter "Parasitism p (young adults)" var: p_para among: [0.0, 0.3485, 0.668] category: "Pest";
    parameter "Neighbour transmission" var: beta_nb min: 0.0 max: 1.0 category: "Spread";
    parameter "Long-jump probability" var: p_jump min: 0.0 max: 1.0 category: "Spread";
    parameter "Wind bias (0 = none)" var: wind_bias min: 0.0 max: 1.0 category: "Spread";
    parameter "Wind direction (deg, toward)" var: wind_dir min: 0.0 max: 359.0 category: "Spread";
    parameter "Show road" var: show_road category: "Spread";
    parameter "Roadside long-jump multiplier (1 = none)" var: road_mult min: 1.0 max: 5.0 category: "Spread";
    parameter "Roadside daily inoculation (0 = none)" var: road_daily min: 0.0 max: 0.5 category: "Spread";

    parameter "Efficacy half-life (applications; 0 = none)" var: n_half among: [0, 20, 10] category: "Control";
    parameter "Bio : chemical efficacy ratio" var: r_bc min: 0.8 max: 1.0 category: "Control";
    parameter "Re-treatment interval (d)" var: r_interval min: 7 max: 28 category: "Control";

    parameter "Loss after a season at capacity" var: L_ref min: 0.1 max: 0.5 category: "Damage";
    parameter "Maximum loss" var: L_max among: [0.5, 1.0] category: "Damage";

    output {
        monitor "Day" value: cycle;
        monitor "Crop year" value: 1 + int(cycle / 365);
        monitor "Share infested now" value: share_infested with_precision 3;
        monitor "Share ever infested" value: share_ever with_precision 3;
        monitor "Mean density (K)" value: mean_density with_precision 4;
        monitor "Cumulative K-days" value: cum_kdays with_precision 1;
        monitor "Chemical applications" value: apps_chem;
        monitor "Biological applications" value: apps_bio;
        monitor "Cane harvested (t)" value: cane_total with_precision 0;

        display "Landscape" type: 2d {
            grid cell;
            species farmer aspect: level_outline;
            species road;
            species lgu;
            // species car;   // vehicles removed (8 Nov 2026)
            overlay position: {0, 0} size: {215 #px, 222 #px} background: rgb(255, 255, 255) transparency: 0.1 border: rgb(90, 90, 90) visible: show_legend {
                float x1 <- 12 #px; float x2 <- 30 #px; float dy <- 17 #px; float y0 <- 18 #px;
                font f <- font("Arial", 10, #plain);
                draw "Legend" at: {x1, y0} color: #black font: font("Arial", 11, #bold);
                draw square(10 #px) at: {x1 + 5 #px, y0 + dy - 4 #px} color: rgb(222, 205, 160);
                draw "cane, not yet planted" at: {x2, y0 + dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 2 * dy - 4 #px} color: rgb(120, 180, 90);
                draw "cane, no RSSI" at: {x2, y0 + 2 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 3 * dy - 4 #px} color: rgb(155, 142, 40);
                draw square(10 #px) at: {x1 + 16 #px, y0 + 3 * dy - 4 #px} color: rgb(255, 0, 40);
                draw "    RSSI density low to high" at: {x2, y0 + 3 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 4 * dy - 4 #px} color: rgb(236, 232, 222) border: rgb(180, 180, 180);
                draw "other land (not cane)" at: {x2, y0 + 4 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 5 * dy - 4 #px} color: rgb(232, 238, 243) border: rgb(180, 180, 180);
                draw "outside model area" at: {x2, y0 + 5 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 6 * dy - 4 #px} color: rgb(255, 215, 0) border: rgb(60, 60, 60);
                draw circle(5 #px) at: {x1 + 19 #px, y0 + 6 * dy - 4 #px} color: rgb(255, 130, 0) border: rgb(60, 60, 60);
                draw circle(5 #px) at: {x1 + 33 #px, y0 + 6 * dy - 4 #px} color: rgb(20, 20, 20) border: #white;
                draw "farm level 1 / 2 / 3" at: {x2 + 22 #px, y0 + 6 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 7 * dy - 4 #px} color: rgb(255, 0, 200);
                draw "dot in marker: chemical" at: {x2, y0 + 7 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 8 * dy - 4 #px} color: rgb(0, 200, 255);
                draw "dot in marker: biological" at: {x2, y0 + 8 * dy} color: #black font: f;
                draw line([{x1, y0 + 9 * dy - 4 #px}, {x1 + 12 #px, y0 + 9 * dy - 4 #px}]) + (1.5 #px) color: rgb(95, 95, 95);
                draw "road (OpenStreetMap)" at: {x2, y0 + 9 * dy} color: #black font: f;
                draw line([{x1, y0 + 10 * dy - 4 #px}, {x1 + 12 #px, y0 + 10 * dy - 4 #px}]) + (1.2 #px) color: rgb(60, 60, 60);
                draw "municipal boundary" at: {x2, y0 + 10 * dy} color: #black font: f;
                draw "1 cell = 1 ha; 1 tick = 1 day" at: {x1, y0 + 11 * dy + 6 #px} color: rgb(90, 90, 90) font: font("Arial", 9, #italic);
            }
        }

        display "Spread" type: 2d refresh: every(1 #cycle) {
            chart "Area infested" type: series {
                data "infested now" value: share_infested color: #red;
                data "ever infested" value: share_ever color: #darkred;
            }	
        }

        display "Control" type: 2d refresh: every(7 #cycle) {
            chart "Applications" type: series {
                data "chemical" value: apps_chem color: #orange;
                data "biological" value: apps_bio color: #blue;
            }
        }
    }
}

// Four strategies side by side: same seed, same landscape and start, one map each.
// The two sliders are live (reflex sync_sliders): move them while the run is going and the change
// acts from the next day. demo_road and demo_wind have the values built in. The map label shows
// the road and wind values actually in use.
experiment compare_strategies type: gui {
    float seed_value <- 42.0;
    float road_demo <- 0.0;           // roadside daily inoculation for all four maps (0 = no road effect)
    float wind_demo <- 0.0;           // wind bias for all four maps (0 = no wind)
    parameter "Roadside daily inoculation (all maps)" var: road_demo min: 0.0 max: 0.5;
    parameter "Wind bias (all maps)" var: wind_demo min: 0.0 max: 1.0;

    // Live sliders: every step, copy the slider values into all four simulations, so a change
    // acts from the next day without reloading.
    reflex sync_sliders {
        ask simulations {
            road_daily <- myself.road_demo;
            wind_bias <- myself.wind_demo;
        }
    }

    action _init_() {
        create simulation with: [seed::seed_value, strategy::"none", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, strategy::"chemical", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, strategy::"biological", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, strategy::"hybrid", road_daily::road_demo, wind_bias::wind_demo];
    }

    output {
        layout #split;
        display "Map" type: 2d {
            grid cell;
            species farmer aspect: level_outline;
            graphics "label" {
                draw strategy + "  |  day " + cycle + "  |  ever infested " + (share_ever with_precision 2)
                    + "  |  road " + road_daily + "  wind " + wind_bias
                    at: {world.shape.width * 0.02, world.shape.height * 0.04} color: #black font: font("Arial", 14, #bold);
            }
            species road;
            species lgu;
            // species car;   // vehicles removed (8 Nov 2026)
            overlay position: {0, 0} size: {215 #px, 222 #px} background: rgb(255, 255, 255) transparency: 0.1 border: rgb(90, 90, 90) visible: show_legend {
                float x1 <- 12 #px; float x2 <- 30 #px; float dy <- 17 #px; float y0 <- 18 #px;
                font f <- font("Arial", 10, #plain);
                draw "Legend" at: {x1, y0} color: #black font: font("Arial", 11, #bold);
                draw square(10 #px) at: {x1 + 5 #px, y0 + dy - 4 #px} color: rgb(222, 205, 160);
                draw "cane, not yet planted" at: {x2, y0 + dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 2 * dy - 4 #px} color: rgb(120, 180, 90);
                draw "cane, no RSSI" at: {x2, y0 + 2 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 3 * dy - 4 #px} color: rgb(155, 142, 40);
                draw square(10 #px) at: {x1 + 16 #px, y0 + 3 * dy - 4 #px} color: rgb(255, 0, 40);
                draw "    RSSI density low to high" at: {x2, y0 + 3 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 4 * dy - 4 #px} color: rgb(236, 232, 222) border: rgb(180, 180, 180);
                draw "other land (not cane)" at: {x2, y0 + 4 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 5 * dy - 4 #px} color: rgb(232, 238, 243) border: rgb(180, 180, 180);
                draw "outside model area" at: {x2, y0 + 5 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 6 * dy - 4 #px} color: rgb(255, 215, 0) border: rgb(60, 60, 60);
                draw circle(5 #px) at: {x1 + 19 #px, y0 + 6 * dy - 4 #px} color: rgb(255, 130, 0) border: rgb(60, 60, 60);
                draw circle(5 #px) at: {x1 + 33 #px, y0 + 6 * dy - 4 #px} color: rgb(20, 20, 20) border: #white;
                draw "farm level 1 / 2 / 3" at: {x2 + 22 #px, y0 + 6 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 7 * dy - 4 #px} color: rgb(255, 0, 200);
                draw "dot in marker: chemical" at: {x2, y0 + 7 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 8 * dy - 4 #px} color: rgb(0, 200, 255);
                draw "dot in marker: biological" at: {x2, y0 + 8 * dy} color: #black font: f;
                draw line([{x1, y0 + 9 * dy - 4 #px}, {x1 + 12 #px, y0 + 9 * dy - 4 #px}]) + (1.5 #px) color: rgb(95, 95, 95);
                draw "road (OpenStreetMap)" at: {x2, y0 + 9 * dy} color: #black font: f;
                draw line([{x1, y0 + 10 * dy - 4 #px}, {x1 + 12 #px, y0 + 10 * dy - 4 #px}]) + (1.2 #px) color: rgb(60, 60, 60);
                draw "municipal boundary" at: {x2, y0 + 10 * dy} color: #black font: f;
                draw "1 cell = 1 ha; 1 tick = 1 day" at: {x1, y0 + 11 * dy + 6 #px} color: rgb(90, 90, 90) font: font("Arial", 9, #italic);
            }
        }
    }
}

experiment video_A type: gui {
    // For the presentation video (Doc Arnejo, 13 Nov): window A (Cadiz, Manapla, Victorias), seeded
    // start on day 233 (22 May) in those 3 of the 11 first-declared LGUs (87 cells, 29 each).
    // Days 0-232 have no pests and run fast. Record with GAMA's own movie export or any screen recorder.
    parameter "Landscape" var: landscape init: "A" among: ["A"] category: "Map";
    string seed_demo <- "may2025";
    parameter "Seeding" var: seed_demo among: ["may2025", "random"] category: "Map";
    float seed_value <- 42.0;
    float road_demo <- 0.0;           // roadside daily inoculation for all four maps (0 = no road effect)
    float wind_demo <- 0.0;           // wind bias for all four maps (0 = no wind)
    parameter "Roadside daily inoculation (all maps)" var: road_demo min: 0.0 max: 0.5;
    parameter "Wind bias (all maps)" var: wind_demo min: 0.0 max: 1.0;

    // Live sliders: every step, copy the slider values into all four simulations, so a change
    // acts from the next day without reloading.
    reflex sync_sliders {
        ask simulations {
            road_daily <- myself.road_demo;
            wind_bias <- myself.wind_demo;
        }
    }

    action _init_() {
        create simulation with: [seed::seed_value, landscape::"A", seed_mode::seed_demo, strategy::"none", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, landscape::"A", seed_mode::seed_demo, strategy::"chemical", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, landscape::"A", seed_mode::seed_demo, strategy::"biological", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, landscape::"A", seed_mode::seed_demo, strategy::"hybrid", road_daily::road_demo, wind_bias::wind_demo];
    }

    output {
        layout #split;
        display "Map" type: 2d {
            grid cell;
            species farmer aspect: level_outline;
            graphics "label" {
                draw strategy + "  |  day " + cycle + "  |  ever infested " + (share_ever with_precision 2)
                    + "  |  road " + road_daily + "  wind " + wind_bias
                    at: {world.shape.width * 0.02, world.shape.height * 0.04} color: #black font: font("Arial", 14, #bold);
            }
            species road;
            species lgu;
            // species car;   // vehicles removed (8 Nov 2026)
            overlay position: {0, 0} size: {215 #px, 222 #px} background: rgb(255, 255, 255) transparency: 0.1 border: rgb(90, 90, 90) visible: show_legend {
                float x1 <- 12 #px; float x2 <- 30 #px; float dy <- 17 #px; float y0 <- 18 #px;
                font f <- font("Arial", 10, #plain);
                draw "Legend" at: {x1, y0} color: #black font: font("Arial", 11, #bold);
                draw square(10 #px) at: {x1 + 5 #px, y0 + dy - 4 #px} color: rgb(222, 205, 160);
                draw "cane, not yet planted" at: {x2, y0 + dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 2 * dy - 4 #px} color: rgb(120, 180, 90);
                draw "cane, no RSSI" at: {x2, y0 + 2 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 3 * dy - 4 #px} color: rgb(155, 142, 40);
                draw square(10 #px) at: {x1 + 16 #px, y0 + 3 * dy - 4 #px} color: rgb(255, 0, 40);
                draw "    RSSI density low to high" at: {x2, y0 + 3 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 4 * dy - 4 #px} color: rgb(236, 232, 222) border: rgb(180, 180, 180);
                draw "other land (not cane)" at: {x2, y0 + 4 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 5 * dy - 4 #px} color: rgb(232, 238, 243) border: rgb(180, 180, 180);
                draw "outside model area" at: {x2, y0 + 5 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 6 * dy - 4 #px} color: rgb(255, 215, 0) border: rgb(60, 60, 60);
                draw circle(5 #px) at: {x1 + 19 #px, y0 + 6 * dy - 4 #px} color: rgb(255, 130, 0) border: rgb(60, 60, 60);
                draw circle(5 #px) at: {x1 + 33 #px, y0 + 6 * dy - 4 #px} color: rgb(20, 20, 20) border: #white;
                draw "farm level 1 / 2 / 3" at: {x2 + 22 #px, y0 + 6 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 7 * dy - 4 #px} color: rgb(255, 0, 200);
                draw "dot in marker: chemical" at: {x2, y0 + 7 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 8 * dy - 4 #px} color: rgb(0, 200, 255);
                draw "dot in marker: biological" at: {x2, y0 + 8 * dy} color: #black font: f;
                draw line([{x1, y0 + 9 * dy - 4 #px}, {x1 + 12 #px, y0 + 9 * dy - 4 #px}]) + (1.5 #px) color: rgb(95, 95, 95);
                draw "road (OpenStreetMap)" at: {x2, y0 + 9 * dy} color: #black font: f;
                draw line([{x1, y0 + 10 * dy - 4 #px}, {x1 + 12 #px, y0 + 10 * dy - 4 #px}]) + (1.2 #px) color: rgb(60, 60, 60);
                draw "municipal boundary" at: {x2, y0 + 10 * dy} color: #black font: f;
                draw "1 cell = 1 ha; 1 tick = 1 day" at: {x1, y0 + 11 * dy + 6 #px} color: rgb(90, 90, 90) font: font("Arial", 9, #italic);
            }
        }
    }
}


// Demo: roadside exposure ON (road_daily 0.5, the slider maximum), wind off. Watch the "none" map:
// cells along the road become infested ahead of the rest. Illustration of a swept assumption
// (Burgonio: roadside plantations worst hit); not a calibrated rate.
experiment demo_road type: gui {
    float seed_value <- 42.0;
    float road_demo <- 0.5;           // roadside daily inoculation for all four maps (0 = no road effect)
    float wind_demo <- 0.0;           // wind bias for all four maps (0 = no wind)

    action _init_() {
        create simulation with: [seed::seed_value, strategy::"none", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, strategy::"chemical", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, strategy::"biological", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, strategy::"hybrid", road_daily::road_demo, wind_bias::wind_demo];
    }

    output {
        layout #split;
        display "Map" type: 2d {
            grid cell;
            species farmer aspect: level_outline;
            graphics "label" {
                draw strategy + "  |  day " + cycle + "  |  ever infested " + (share_ever with_precision 2)
                    + "  |  road " + road_daily + "  wind " + wind_bias
                    at: {world.shape.width * 0.02, world.shape.height * 0.04} color: #black font: font("Arial", 14, #bold);
            }
            species road;
            species lgu;
            // species car;   // vehicles removed (8 Nov 2026)
            overlay position: {0, 0} size: {215 #px, 222 #px} background: rgb(255, 255, 255) transparency: 0.1 border: rgb(90, 90, 90) visible: show_legend {
                float x1 <- 12 #px; float x2 <- 30 #px; float dy <- 17 #px; float y0 <- 18 #px;
                font f <- font("Arial", 10, #plain);
                draw "Legend" at: {x1, y0} color: #black font: font("Arial", 11, #bold);
                draw square(10 #px) at: {x1 + 5 #px, y0 + dy - 4 #px} color: rgb(222, 205, 160);
                draw "cane, not yet planted" at: {x2, y0 + dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 2 * dy - 4 #px} color: rgb(120, 180, 90);
                draw "cane, no RSSI" at: {x2, y0 + 2 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 3 * dy - 4 #px} color: rgb(155, 142, 40);
                draw square(10 #px) at: {x1 + 16 #px, y0 + 3 * dy - 4 #px} color: rgb(255, 0, 40);
                draw "    RSSI density low to high" at: {x2, y0 + 3 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 4 * dy - 4 #px} color: rgb(236, 232, 222) border: rgb(180, 180, 180);
                draw "other land (not cane)" at: {x2, y0 + 4 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 5 * dy - 4 #px} color: rgb(232, 238, 243) border: rgb(180, 180, 180);
                draw "outside model area" at: {x2, y0 + 5 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 6 * dy - 4 #px} color: rgb(255, 215, 0) border: rgb(60, 60, 60);
                draw circle(5 #px) at: {x1 + 19 #px, y0 + 6 * dy - 4 #px} color: rgb(255, 130, 0) border: rgb(60, 60, 60);
                draw circle(5 #px) at: {x1 + 33 #px, y0 + 6 * dy - 4 #px} color: rgb(20, 20, 20) border: #white;
                draw "farm level 1 / 2 / 3" at: {x2 + 22 #px, y0 + 6 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 7 * dy - 4 #px} color: rgb(255, 0, 200);
                draw "dot in marker: chemical" at: {x2, y0 + 7 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 8 * dy - 4 #px} color: rgb(0, 200, 255);
                draw "dot in marker: biological" at: {x2, y0 + 8 * dy} color: #black font: f;
                draw line([{x1, y0 + 9 * dy - 4 #px}, {x1 + 12 #px, y0 + 9 * dy - 4 #px}]) + (1.5 #px) color: rgb(95, 95, 95);
                draw "road (OpenStreetMap)" at: {x2, y0 + 9 * dy} color: #black font: f;
                draw line([{x1, y0 + 10 * dy - 4 #px}, {x1 + 12 #px, y0 + 10 * dy - 4 #px}]) + (1.2 #px) color: rgb(60, 60, 60);
                draw "municipal boundary" at: {x2, y0 + 10 * dy} color: #black font: f;
                draw "1 cell = 1 ha; 1 tick = 1 day" at: {x1, y0 + 11 * dy + 6 #px} color: rgb(90, 90, 90) font: font("Arial", 9, #italic);
            }
        }
    }
}

// Demo: wind bias ON (wind_bias 1, wind from the NE, wind_dir 135), road effect off. Watch the
// "none" map: patches grow toward the lower left (SW) and barely toward the upper right (NE).
// Illustration of a swept assumption (Doc Arnejo's suggestion); not a calibrated wind model.
experiment demo_wind type: gui {
    float seed_value <- 42.0;
    float road_demo <- 0.0;           // roadside daily inoculation for all four maps (0 = no road effect)
    float wind_demo <- 1.0;           // wind bias for all four maps (0 = no wind)

    action _init_() {
        create simulation with: [seed::seed_value, strategy::"none", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, strategy::"chemical", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, strategy::"biological", road_daily::road_demo, wind_bias::wind_demo];
        create simulation with: [seed::seed_value, strategy::"hybrid", road_daily::road_demo, wind_bias::wind_demo];
    }

    output {
        layout #split;
        display "Map" type: 2d {
            grid cell;
            species farmer aspect: level_outline;
            graphics "label" {
                draw strategy + "  |  day " + cycle + "  |  ever infested " + (share_ever with_precision 2)
                    + "  |  road " + road_daily + "  wind " + wind_bias
                    at: {world.shape.width * 0.02, world.shape.height * 0.04} color: #black font: font("Arial", 14, #bold);
            }
            species road;
            species lgu;
            // species car;   // vehicles removed (8 Nov 2026)
            overlay position: {0, 0} size: {215 #px, 222 #px} background: rgb(255, 255, 255) transparency: 0.1 border: rgb(90, 90, 90) visible: show_legend {
                float x1 <- 12 #px; float x2 <- 30 #px; float dy <- 17 #px; float y0 <- 18 #px;
                font f <- font("Arial", 10, #plain);
                draw "Legend" at: {x1, y0} color: #black font: font("Arial", 11, #bold);
                draw square(10 #px) at: {x1 + 5 #px, y0 + dy - 4 #px} color: rgb(222, 205, 160);
                draw "cane, not yet planted" at: {x2, y0 + dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 2 * dy - 4 #px} color: rgb(120, 180, 90);
                draw "cane, no RSSI" at: {x2, y0 + 2 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 3 * dy - 4 #px} color: rgb(155, 142, 40);
                draw square(10 #px) at: {x1 + 16 #px, y0 + 3 * dy - 4 #px} color: rgb(255, 0, 40);
                draw "    RSSI density low to high" at: {x2, y0 + 3 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 4 * dy - 4 #px} color: rgb(236, 232, 222) border: rgb(180, 180, 180);
                draw "other land (not cane)" at: {x2, y0 + 4 * dy} color: #black font: f;
                draw square(10 #px) at: {x1 + 5 #px, y0 + 5 * dy - 4 #px} color: rgb(232, 238, 243) border: rgb(180, 180, 180);
                draw "outside model area" at: {x2, y0 + 5 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 6 * dy - 4 #px} color: rgb(255, 215, 0) border: rgb(60, 60, 60);
                draw circle(5 #px) at: {x1 + 19 #px, y0 + 6 * dy - 4 #px} color: rgb(255, 130, 0) border: rgb(60, 60, 60);
                draw circle(5 #px) at: {x1 + 33 #px, y0 + 6 * dy - 4 #px} color: rgb(20, 20, 20) border: #white;
                draw "farm level 1 / 2 / 3" at: {x2 + 22 #px, y0 + 6 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 7 * dy - 4 #px} color: rgb(255, 0, 200);
                draw "dot in marker: chemical" at: {x2, y0 + 7 * dy} color: #black font: f;
                draw circle(5 #px) at: {x1 + 5 #px, y0 + 8 * dy - 4 #px} color: rgb(0, 200, 255);
                draw "dot in marker: biological" at: {x2, y0 + 8 * dy} color: #black font: f;
                draw line([{x1, y0 + 9 * dy - 4 #px}, {x1 + 12 #px, y0 + 9 * dy - 4 #px}]) + (1.5 #px) color: rgb(95, 95, 95);
                draw "road (OpenStreetMap)" at: {x2, y0 + 9 * dy} color: #black font: f;
                draw line([{x1, y0 + 10 * dy - 4 #px}, {x1 + 12 #px, y0 + 10 * dy - 4 #px}]) + (1.2 #px) color: rgb(60, 60, 60);
                draw "municipal boundary" at: {x2, y0 + 10 * dy} color: #black font: f;
                draw "1 cell = 1 ha; 1 tick = 1 day" at: {x1, y0 + 11 * dy + 6 #px} color: rgb(90, 90, 90) font: font("Arial", 9, #italic);
            }
        }
    }
}

// ---------------- Sanity checks ----------------
// 1. Life table: one cohort in one cell, nothing else on. Expect: the cohort is gone after 63
//    days (monitor reads 63), eggs per female = 197 x (1 - p) (197 at p = 0), and the share left
//    after the young-adult stage = 1 - p (0.6515 at p = 0.3485).
experiment sanity_cohort type: gui {
    float p_test <- 0.0;
    parameter "Parasitism p (0, 0.3485, 0.668)" var: p_test among: [0.0, 0.3485, 0.668];
    action _init_() {
        create simulation with: [cohort_test::true, strategy::"none", beta_nb::0.0, p_jump::0.0,
            road_daily::0.0, crop_length::999999, p_para::p_test, horizon::80, sanity_tag::"cohort", gui::false];
    }
    output {
        monitor "Day" value: cycle;
        monitor "Cohort alive (start = 1)" value: (empty(first(cane_cells).pop) ? 0.0 : sum(first(cane_cells).pop)) with_precision 4;
        monitor "Eggs laid per female (expect 197 x (1 - p))" value: eggs_laid_total with_precision 2;
        monitor "Cohort gone after day (expect 63)" value: cohort_death_day;
        monitor "Left after young adults (expect 1 - p)" value: cohort_after_young_adult with_precision 4;
    }
}

// 1b. Same life-table check as a one-click batch: three deterministic runs (p = 0, 0.3485, 0.668),
//     one summary file each in ../results/sanity/. Use this one: the GUI version above does not pick
//     up a changed p unless the experiment is reloaded.
experiment sanity_cohort_all type: batch repeat: 1 keep_seed: true parallel: 1 until: cycle >= horizon {
    parameter "gui" var: gui among: [false];
    parameter "Parasitism p" var: p_para among: [0.0, 0.3485, 0.668];
    parameter "cohort" var: cohort_test among: [true];
    parameter "strategy" var: strategy among: ["none"];
    parameter "beta_nb" var: beta_nb among: [0.0];
    parameter "p_jump" var: p_jump among: [0.0];
    parameter "road_daily" var: road_daily among: [0.0];
    parameter "crop_length" var: crop_length among: [999999];
    parameter "horizon" var: horizon among: [80];
    parameter "tag" var: sanity_tag among: ["cohort"];
    method exhaustive;
}

// 2. No pest: nothing seeded, nothing spreads. Expect: share infested stays 0, no applications,
//    loss 0, cane per harvested hectare = Y0 (61.64).
experiment sanity_no_pest type: gui {
    action _init_() {
        create simulation with: [init_share::0.0, p_jump::0.0, road_daily::0.0, strategy::"hybrid", sanity_tag::"no_pest", gui::false];
    }
    output {
        monitor "Day" value: cycle;
        monitor "Share ever infested (expect 0)" value: share_ever;
        monitor "Applications (expect 0)" value: apps_chem + apps_bio;
        monitor "Cane per harvested ha (expect 61.64)" value: harvest_count = 0 ? 0.0 : (cane_total / harvest_count) with_precision 2;
    }
}

// ---------------- Batch experiments ----------------
// parallel: N = simulations held in memory at once. Each holds the whole island (over 11 GB on the laptop, 8 Nov; measure M on the campus machine);
// set N = MEM / M rounded down (README_campus.md). If GAMA warns that memory is low, lower N.
// Rows are appended to ../results/<experiment folder>/*.csv, one per run (main/, oat/, wind/).
// Purpose and goal of each experiment: docs/experiments.md. Columns (no header row):
// experiment, oat_case, strategy, n_half, seed, share_ever, share_infested_end, cum_kdays,
// apps_chem, apps_bio, cane_total, harvest_count, max_D (max_D added 4 Nov; the 320 dry-run rows lack it).
// apps_* count farm-level treatment events (one per 4-ha farm), not hectares treated.

// Main experiment, part 1: strategies that use chemicals x efficacy decline (6 configs x 10 = 60 runs; island scale, 8 Nov 2026)
experiment batch_main_chemical type: batch repeat: 10 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "gui" var: gui among: [false];
    parameter "Strategy" var: strategy among: ["chemical", "hybrid"];
    parameter "Efficacy half-life" var: n_half among: [0, 20, 10];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

// Main experiment, part 2: strategies without chemicals (2 configs x 10 = 20 runs)
experiment batch_main_nochemical type: batch repeat: 10 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "gui" var: gui among: [false];
    parameter "Strategy" var: strategy among: ["none", "biological"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

// One-at-a-time sensitivity: baseline (case 0) + 18 cases x 4 strategies x 3 replicates = 228 runs (island scale, 8 Nov 2026),
// n_half fixed at 20.
experiment batch_oat type: batch repeat: 3 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "gui" var: gui among: [false];
    parameter "OAT case" var: oat_case among: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18];
    parameter "Strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

// Wind-direction check (docs/notes: weather-data note): does direction matter? On the real map it may.
// wind_bias 1, direction 135 (from NE, PAGASA Dumaguete annual) vs 315 (from SW), 4 strategies,
// 3 replicates = 24 runs. If outputs match within replicate spread, direction does not matter for SP1.
experiment batch_wind_direction type: batch repeat: 3 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "gui" var: gui among: [false];
    parameter "Wind bias" var: wind_bias among: [1.0];
    parameter "Wind direction" var: wind_dir among: [135.0, 315.0];
    parameter "Strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/wind/batch_wind_direction.csv"];
    parameter "label" var: batch_label among: ["wind"];
}

// Time-scale check (ODD pattern 3; D7, D13): weekly share ever infested under the hybrid strategy,
// 10 replicates, baseline parameters; also writes per-LGU ever-infested cells (results/series/lgu_*.csv). Compared later with the verified point (30.1% of planted
// area validated by 2 Oct 2026, A21) once the model's outbreak start date is defined.
experiment batch_series_hybrid type: batch repeat: 10 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "gui" var: gui among: [false];
    parameter "Strategy" var: strategy among: ["hybrid"];
    parameter "save" var: save_series among: [true];
    parameter "tag" var: series_tag among: ["hybrid_baseline"];
    method exhaustive;
}

// Damage-form check (DM1): main configurations at n_half = 20 with saturating damage, 3 replicates.
// Under linear damage the ranking by loss equals the ranking by pest load; this checks it under
// saturation. Rows go to ../results/main/batch_saturating.csv.
experiment batch_saturating type: batch repeat: 3 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "gui" var: gui among: [false];
    parameter "Strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    parameter "Saturating" var: dmg_saturating among: [true];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/main/batch_saturating.csv"];
    parameter "label" var: batch_label among: ["sat"];
}

// No-pest sanity check as a one-click batch (runs headless; same check as sanity_no_pest above).
// Expect results/sanity/sanity_no_pest.csv: share_ever 0, apps 0, cane per harvested ha 61.64.
experiment sanity_no_pest_all type: batch repeat: 1 keep_seed: true parallel: 1 until: cycle >= horizon {
    parameter "gui" var: gui among: [false];
    parameter "init_share" var: init_share among: [0.0];
    parameter "p_jump" var: p_jump among: [0.0];
    parameter "road_daily" var: road_daily among: [0.0];
    parameter "Strategy" var: strategy among: ["hybrid"];
    parameter "tag" var: sanity_tag among: ["no_pest"];
    method exhaustive;
}

// Longer-horizon hybrid check (Doc Arnejo, 7 Nov: the hybrid may separate from biological control
// only in a longer run). 2,190 d = 6 crop years, baseline parameters, 10 replicates. Weekly series to
// results/series/series_hybrid_long.csv (with per-LGU end rows), one end row per run to
// results/main/batch_long.csv (max_D shows whether any cell nears Level 3 densities).
experiment batch_long_hybrid type: batch repeat: 10 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "gui" var: gui among: [false];
    parameter "Strategy" var: strategy among: ["hybrid"];
    parameter "horizon" var: horizon among: [2190];
    parameter "save" var: save_series among: [true];
    parameter "tag" var: series_tag among: ["hybrid_long"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/main/batch_long.csv"];
    parameter "label" var: batch_label among: ["long"];
}

// ===== Window A (Cadiz, Manapla, Victorias): sensitivity, wind and damage-form checks at a scale a
// laptop can run. Seeded start (87 cells on day 233 in Cadiz, Manapla, Victorias), horizon 1,328 d.
// Run counts cut 13 Nov (Doc Arnejo: as few runs as acceptable): 2 replicates for main, wind and
// damage form; 1 per cell for the OAT screen, read against the replicate spread of A_main.
// Same design as the island batches, split so each JVM holds at most
// 12 runs. Results go to results/A/; analyse with: python analysis/analyse_results.py results/A

experiment A_main_chemical type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "strategy" var: strategy among: ["chemical", "hybrid"];
    parameter "n_half" var: n_half among: [0, 20, 10];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment A_main_nochemical type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "strategy" var: strategy among: ["none", "biological"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment A_oat_0 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [0];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_1 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [1];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_2 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [2];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_3 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [3];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_4 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [4];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_5 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [5];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_6 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [6];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_7 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [7];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_8 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [8];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_9 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [9];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_10 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [10];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_11 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [11];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_12 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [12];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_13 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [13];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_14 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [14];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_15 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [15];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_16 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [16];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_17 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [17];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_oat_18 type: batch repeat: 1 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [18];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/oat/batch_oat.csv"];
    parameter "label" var: batch_label among: ["oat"];
}

experiment A_wind_135 type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "wind_bias" var: wind_bias among: [1.0];
    parameter "wind_dir" var: wind_dir among: [135.0];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/wind/batch_wind_direction.csv"];
    parameter "label" var: batch_label among: ["wind"];
}

experiment A_wind_315 type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "wind_bias" var: wind_bias among: [1.0];
    parameter "wind_dir" var: wind_dir among: [315.0];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/wind/batch_wind_direction.csv"];
    parameter "label" var: batch_label among: ["wind"];
}

experiment A_saturating type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "dmg_saturating" var: dmg_saturating among: [true];
    parameter "strategy" var: strategy among: ["none", "chemical", "biological", "hybrid"];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/main/batch_saturating.csv"];
    parameter "label" var: batch_label among: ["sat"];
}

// ===== Seeded island runs (13 Nov 2026, after the mock defence; J.M. chose "everything seeded").
// Start: pest-free crop calendar from day 0 (1 Oct); 87 cells seeded on day 233 (22 May 2025) in the
// 11 LGUs of SRA's first declared infestation; then 1,095 days, so horizon = 233 + 1,095 = 1,328.
// Required set (13 Nov, run counts cut): hybrid, none, chemical n0, biological; the decline variants
// (chemical n20, n10) are optional on the island and covered on window A.
// One configuration x 2 replicates per experiment (one parallel pair per JVM: the 11-12 Nov campus
// run stopped every multi-pair batch after its first pair). Results: results/seeded/.
// isl_s_hybrid also writes the weekly series and a per-LGU snapshot on day 731 (2 Oct 2026) for the
// comparison with SRA FOI Table I (the hindcast).

experiment isl_s_hybrid type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["hybrid"];
    parameter "n_half" var: n_half among: [20];
    parameter "save" var: save_series among: [true];
    parameter "tag" var: series_tag among: ["hybrid_baseline"];
    parameter "dir" var: series_dir among: ["../results/seeded/series/"];
    parameter "lgu_day" var: lgu_day among: [731];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment isl_s_none type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["none"];
    parameter "n_half" var: n_half among: [20];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment isl_s_chemical_n0 type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["chemical"];
    parameter "n_half" var: n_half among: [0];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment isl_s_biological type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["biological"];
    parameter "n_half" var: n_half among: [20];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment isl_s_chemical_n20 type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["chemical"];
    parameter "n_half" var: n_half among: [20];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment isl_s_chemical_n10 type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["chemical"];
    parameter "n_half" var: n_half among: [10];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

// ===== Laptop versions (15 Nov 2026): same island seeded runs, ONE simulation per JVM (16 GB RAM).
// Run each twice for 2 replicates; rows append to the same files as isl_s_*.
experiment lap_s_hybrid type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["hybrid"];
    parameter "n_half" var: n_half among: [20];
    parameter "save" var: save_series among: [true];
    parameter "tag" var: series_tag among: ["hybrid_baseline"];
    parameter "dir" var: series_dir among: ["../results/seeded/series/"];
    parameter "lgu_day" var: lgu_day among: [731];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment lap_s_none type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["none"];
    parameter "n_half" var: n_half among: [20];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment lap_s_chemical_n0 type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["chemical"];
    parameter "n_half" var: n_half among: [0];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment lap_s_biological type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["biological"];
    parameter "n_half" var: n_half among: [20];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

// ===== Laptop follow-up (15 Nov 2026): no-control run with the weekly series (upper bound for the
// hindcast; tag none_baseline) and chemical with efficacy decline on the island.
experiment lap_s_none_series type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["none"];
    parameter "n_half" var: n_half among: [20];
    parameter "save" var: save_series among: [true];
    parameter "tag" var: series_tag among: ["none_baseline"];
    parameter "dir" var: series_dir among: ["../results/seeded/series/"];
    parameter "lgu_day" var: lgu_day among: [731];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment lap_s_chemical_n20 type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["chemical"];
    parameter "n_half" var: n_half among: [20];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

experiment lap_s_chemical_n10 type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "strategy" var: strategy among: ["chemical"];
    parameter "n_half" var: n_half among: [10];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/seeded/main/batch_main.csv"];
    parameter "label" var: batch_label among: ["main"];
}

// ===== Extras for the defence (15 Nov 2026) =====
// (1) Window A: where does declining chemical match biological? n_half 40, 80, 160 (2 replicates).
experiment A_decline type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "strategy" var: strategy among: ["chemical"];
    parameter "n_half" var: n_half among: [40, 80, 160];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/main/batch_decline.csv"];
    parameter "label" var: batch_label among: ["decline"];
}

// (2) Window A: is the hybrid rule ever better than biological? Conditions where farms reach Level 3
// (OAT 14: trigger at Level 3; OAT 16: re-treatment every 28 d) with chemical efficacy NOT declining.
experiment A_hybrid_test type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "oat_case" var: oat_case among: [14, 16];
    parameter "strategy" var: strategy among: ["chemical", "biological", "hybrid"];
    parameter "n_half" var: n_half among: [0];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/main/batch_hybridtest.csv"];
    parameter "label" var: batch_label among: ["hybtest"];
}

// (3) Island hindcast probes, no control, to day 735 only (2 Oct 2026 = day 731): can faster spread
// assumptions reach the SRA validated area? b4 = neighbour transmission x4 (OAT 4); j02 = long jump 0.2
// (OAT 6); b4j02 = both. One run each; weekly series and per-LGU snapshot in results/seeded/hindcast/.
experiment lap_h_b4 type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [735];
    parameter "strategy" var: strategy among: ["none"];
    parameter "oat_case" var: oat_case among: [4];
    parameter "save" var: save_series among: [true];
    parameter "tag" var: series_tag among: ["none_b4"];
    parameter "dir" var: series_dir among: ["../results/seeded/hindcast/"];
    parameter "lgu_day" var: lgu_day among: [731];
    method exhaustive;
}

experiment lap_h_j02 type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [735];
    parameter "strategy" var: strategy among: ["none"];
    parameter "oat_case" var: oat_case among: [6];
    parameter "save" var: save_series among: [true];
    parameter "tag" var: series_tag among: ["none_j02"];
    parameter "dir" var: series_dir among: ["../results/seeded/hindcast/"];
    parameter "lgu_day" var: lgu_day among: [731];
    method exhaustive;
}

experiment lap_h_b4j02 type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [735];
    parameter "strategy" var: strategy among: ["none"];
    parameter "oat_case" var: oat_case among: [4];
    parameter "p_jump" var: p_jump among: [0.2];
    parameter "save" var: save_series among: [true];
    parameter "tag" var: series_tag among: ["none_b4j02"];
    parameter "dir" var: series_dir among: ["../results/seeded/hindcast/"];
    parameter "lgu_day" var: lgu_day among: [731];
    method exhaustive;
}

// ===== Follow-up tests (15 Nov 2026, night) =====
// (4) Window A symmetry test: what if the biological agent loses efficacy at the same rate as the
// chemical? n_half_bio 20 and 80, biological and hybrid (hybrid's chemical arm at n_half 20), 2 replicates.
experiment A_biodecline_20 type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "strategy" var: strategy among: ["biological", "hybrid"];
    parameter "n_half" var: n_half among: [20];
    parameter "n_half_bio" var: n_half_bio among: [20];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/main/batch_biodecline.csv"];
    parameter "label" var: batch_label among: ["biodec20"];
}

experiment A_biodecline_80 type: batch repeat: 2 keep_seed: false parallel: 2 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["A"];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [1328];
    parameter "gui" var: gui among: [false];
    parameter "strategy" var: strategy among: ["biological", "hybrid"];
    parameter "n_half" var: n_half among: [20];
    parameter "n_half_bio" var: n_half_bio among: [80];
    method exhaustive;
    parameter "out" var: batch_file among: ["../results/A/main/batch_biodecline.csv"];
    parameter "label" var: batch_label among: ["biodec80"];
}

// (5) Island tuning probes, no control, to day 735: neighbour transmission x2.5 and x3 (beta_nb 0.25, 0.30),
// to find the rate that gives the SRA validated area on 2 Oct 2026 (day 731). Tuned to one point, not calibrated.
experiment lap_h_b25 type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [735];
    parameter "strategy" var: strategy among: ["none"];
    parameter "beta_nb" var: beta_nb among: [0.25];
    parameter "save" var: save_series among: [true];
    parameter "tag" var: series_tag among: ["none_b25"];
    parameter "dir" var: series_dir among: ["../results/seeded/hindcast/"];
    parameter "lgu_day" var: lgu_day among: [731];
    method exhaustive;
}

experiment lap_h_b30 type: batch repeat: 1 keep_seed: false parallel: 1 until: cycle >= horizon {
    parameter "landscape" var: landscape among: ["negros"];
    parameter "gui" var: gui among: [false];
    parameter "seed_mode" var: seed_mode among: ["may2025"];
    parameter "horizon" var: horizon among: [735];
    parameter "strategy" var: strategy among: ["none"];
    parameter "beta_nb" var: beta_nb among: [0.30];
    parameter "save" var: save_series among: [true];
    parameter "tag" var: series_tag among: ["none_b30"];
    parameter "dir" var: series_dir among: ["../results/seeded/hindcast/"];
    parameter "lgu_day" var: lgu_day among: [731];
    method exhaustive;
}

