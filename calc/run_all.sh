#!/usr/bin/env bash
# Rerun every vacuum-pops calculation, refresh its .txt, and name any script that fails.
# Each script carries its own assertions; the validations include planted failures.
cd "$(dirname "$0")"
fail=0; n=0
run() { n=$((n+1)); if "$@" > "$OUT" 2>&1; then echo "ok    $OUT"; else echo "FAIL  $OUT"; fail=$((fail+1)); fi; }
OUT=p1_kinematics.txt run Rscript p1_kinematics.R
S=.
OUT=p1_reference.txt run python3 p1_reference.py $S/p1_cases.csv
OUT=p2_validate.txt run Rscript p2_validate.R
OUT=p3_spectrum.txt run Rscript p3_spectrum.R 1500
OUT=p3b_highratio.txt run Rscript p3b_highratio.R
OUT=p3c_examples.txt run Rscript p3c_examples.R
OUT=p4_check3d.txt run Rscript p4_check3d.R
OUT=p4b_brute3d.txt run Rscript p4b_brute3d.R
OUT=p4c_converge.txt run Rscript p4c_converge.R
OUT=p5_observables.txt run Rscript p5_observables.R
OUT=p5b_multiplets.txt run Rscript p5b_multiplets.R
OUT=p6_hotgas.txt run Rscript p6_hotgas.R 1000
OUT=p7_mechanism.txt run Rscript p7_mechanism.R
OUT=p8_figure.txt run Rscript p8_figure.R
OUT=p9_loaded_validate.txt run Rscript p9_loaded_validate.R
# p10 solves 50 loaded collapses on 14 cores and takes about an hour; LOADED=0 skips it, and p11 then
# reads the p10_loaded.rds already on disk.
[ "${LOADED:-1}" = 0 ] || { OUT=p10_loaded.txt run Rscript p10_loaded.R 14; }
OUT=p11_loaded_summary.txt run Rscript p11_loaded_summary.R
OUT=p12_halo_directions.txt run Rscript p12_halo_directions.R
OUT=p15_schematic.txt run Rscript p15_schematic.R
OUT=p14_photons_per_proton.txt run Rscript p14_photons_per_proton.R
OUT=p13_famous_events.txt run Rscript p13_famous_events.R
echo "$n scripts run, $fail failed"
exit $fail
