#!/bin/bash
# Runs the RSSI batch experiments headless (no GAMA window), one after another (Linux/macOS).
# Usage:  ./run_batches.sh                       (all experiments, in the order below)
#         ./run_batches.sh batch_main_chemical   (only the experiments named)
# Results go to results/<folder>/ next to this file; times go to results/run_log.txt.

# ===== edit these two lines for the machine =====
GAMA_DIR="$HOME/GAMA"        # folder that contains headless/gama-headless.sh
MEM="24g"
OPTS=""                       # e.g. "-hpc 1" = one simulation at a time (16 GB machines)
# ================================================

HERE="$(cd "$(dirname "$0")" && pwd)"
MODEL="$HERE/models/rssi_prototype.gaml"
LOG="$HERE/results/run_log.txt"
EXPS=(sanity_cohort_all sanity_no_pest_all batch_series_hybrid batch_main_chemical batch_main_nochemical
      batch_long_hybrid batch_oat batch_wind_direction batch_saturating)
[ $# -gt 0 ] && EXPS=("$@")
# 13 Nov (seeded start): "./run_batches.sh island" or "./run_batches.sh windowA"
[ "$1" = island ] && EXPS=(isl_s_hybrid isl_s_none isl_s_chemical_n0 isl_s_biological)
[ "$1" = windowA ] && EXPS=(A_main_chemical A_main_nochemical A_oat_0 A_oat_1 A_oat_2 A_oat_3 A_oat_4 A_oat_5 A_oat_6 A_oat_7 A_oat_8 A_oat_9 A_oat_10 A_oat_11 A_oat_12 A_oat_13 A_oat_14 A_oat_15 A_oat_16 A_oat_17 A_oat_18 A_wind_135 A_wind_315 A_saturating)

if [ ! -f "$GAMA_DIR/headless/gama-headless.sh" ]; then
  echo "Cannot find $GAMA_DIR/headless/gama-headless.sh. Edit GAMA_DIR at the top of this file."; exit 1
fi
echo "$(date) machine=$(hostname) MEM=$MEM experiments: ${EXPS[*]}" >> "$LOG"
WORK="$(mktemp -d)"; cd "$WORK"     # GAMA writes its temporary workspace here
for e in "${EXPS[@]}"; do
  echo "$(date) START $e" >> "$LOG"
  bash "$GAMA_DIR/headless/gama-headless.sh" -m "$MEM" $OPTS -batch "$e" "$MODEL"
  rc=$?
  echo "$(date) END $e exit code $rc" >> "$LOG"
done
echo "Done. See $LOG"
