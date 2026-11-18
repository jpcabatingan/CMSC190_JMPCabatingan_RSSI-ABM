# Development Set-up

The model is `rssi_prototype/models/rssi_prototype.gaml`. Experiments, inputs, run scripts, analysis scripts and the results reported in the paper are in `rssi_prototype/`.

## Prerequisites
* GAMA platform 2026.11.0 (https://gama-platform.org), which bundles its own Java
* RAM: about 11 GB free for one island simulation; window A runs need far less
* Python 3 with numpy, pandas, matplotlib and scipy (analysis only)
* Optional, to rebuild the landscape inputs: geopandas, rasterio, shapely, and the raw data listed in `rssi_prototype/gis-scripts/README.md` (not included; about 900 MB)

## Folder contents
* `models/rssi_prototype.gaml`: the model and all experiments
* `includes/`: landscape inputs for all of Negros (`*_negros.*`) and window A (`*_A.*`): 1-ha grid of municipality and cane fraction, municipal boundaries, roads
* `gis-scripts/`: scripts that built `includes/` from the sugarcane map (Di Tommaso et al. 2024), OCHA COD-AB boundaries and OpenStreetMap roads
* `run_batches.bat`, `run_batches.sh`: run batch experiments headless, one after another
* `analysis/`: summary tables and figures from the batch CSVs
* `results/`: output CSVs of the reported runs (`A/` window A, `seeded/` island), `run_log.txt`, and `A/analysis/paper/` (numbers and figures for the window A results). Batch CSVs have no header; column order is in `results/COLUMNS.csv`

## Instructions
1. Install GAMA 2026.11.0.
2. In GAMA, import `rssi_prototype` as a project (File > Import > Existing project) and open `models/rssi_prototype.gaml`.
3. To watch the model, run the GUI experiment `video_A` (window A, seeded start, four strategies side by side) or `prototype` (one strategy, with parameter sliders).
4. To rerun batches headless, edit `GAMA_DIR`, `MEM` and `OPTS` at the top of `run_batches.bat` (Windows) or `run_batches.sh`, then run for example:
```
run_batches.bat windowA     # window A main comparison, sensitivity, wind, damage form
run_batches.bat laptop      # island seeded runs, one simulation per process
run_batches.bat extras2     # biological decline test and island hindcast probes
```
Keep `OPTS=-hpc 1` on a 16 GB machine so that only one simulation is held in memory. Rows are appended to the CSVs in `results/`; move old CSVs out before rerunning a batch.
5. To rebuild the analysis:
```
python analysis/analyse_results.py results/seeded
cd results/A/analysis/paper && python analysis.py
```
