@echo off
REM Runs the RSSI batch experiments headless (no GAMA window), one after another.
REM Usage:  run_batches.bat                      (all experiments, in the order below)
REM         run_batches.bat batch_main_chemical  (only the experiments named)
REM Results go to results\<folder>\ next to this file; times go to results\run_log.txt.
setlocal EnableDelayedExpansion

REM ===== edit these lines for the machine =====
set "GAMA_DIR=C:\RSSI\GAMA"
set "MEM=11g"
REM OPTS: extra headless options. "-hpc 1" = one simulation at a time (16 GB laptops); leave empty on 32 GB+
set "OPTS=-hpc 1"
REM ================================================

set "MODEL=%~dp0models\rssi_prototype.gaml"
set "LOG=%~dp0results\run_log.txt"
set "EXPS=sanity_cohort_all sanity_no_pest_all batch_series_hybrid batch_main_chemical batch_main_nochemical batch_long_hybrid batch_oat batch_wind_direction batch_saturating"
REM 13 Nov (seeded start): "run_batches.bat island" runs the missing island pairs (campus PC);
REM "run_batches.bat windowA" runs the window A sensitivity set (laptop is enough).
set "ISLAND=isl_s_hybrid isl_s_none isl_s_chemical_n0 isl_s_biological"
set "WINDOWA=A_main_chemical A_main_nochemical A_oat_0 A_oat_1 A_oat_2 A_oat_3 A_oat_4 A_oat_5 A_oat_6 A_oat_7 A_oat_8 A_oat_9 A_oat_10 A_oat_11 A_oat_12 A_oat_13 A_oat_14 A_oat_15 A_oat_16 A_oat_17 A_oat_18 A_wind_135 A_wind_315 A_saturating"
if not "%~1"=="" set "EXPS=%*"
if /I "%~1"=="island" set "EXPS=%ISLAND%"
if /I "%~1"=="windowA" set "EXPS=%WINDOWA%"
if /I "%~1"=="laptop2" set "EXPS=lap_s_none_series lap_s_chemical_n20 lap_s_chemical_n10 lap_s_chemical_n20 lap_s_chemical_n10"
if /I "%~1"=="extras" set "EXPS=lap_h_b4 lap_h_j02 lap_h_b4j02 A_decline A_hybrid_test"
if /I "%~1"=="extras2" set "EXPS=A_biodecline_20 A_biodecline_80 lap_h_b25 lap_h_b30"
REM "run_batches.bat laptop": island seeded runs, one simulation per JVM, hybrid first, each twice
if /I "%~1"=="laptop" set "EXPS=lap_s_hybrid lap_s_hybrid lap_s_none lap_s_chemical_n0 lap_s_biological lap_s_none lap_s_chemical_n0 lap_s_biological"

if not exist "%GAMA_DIR%\headless\gama-headless.bat" (
  echo Cannot find "%GAMA_DIR%\headless\gama-headless.bat". Edit GAMA_DIR at the top of this file.
  exit /b 1
)

echo %date% %time% machine=%COMPUTERNAME% MEM=%MEM% OPTS=%OPTS% experiments: %EXPS%>>"%LOG%"
pushd "%GAMA_DIR%\headless"
for %%E in (%EXPS%) do (
  echo !date! !time! START %%E>>"%LOG%"
  call gama-headless.bat -m %MEM% %OPTS% -batch %%E "%MODEL%"
  echo !date! !time! END %%E exit code !errorlevel!>>"%LOG%"
)
popd
echo Done. See "%LOG%".
endlocal
