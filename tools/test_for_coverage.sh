#!/usr/bin/env bash

set -e

make covered_tests
IMGTHLA_DIR=../foreign/IMGTHLA
LOG_FILE=covered_tests.log

echo testing parsing >> ${LOG_FILE}
_build/default/src/tests/alignment_parsing.exe >> ${LOG_FILE}

echo mas alignment >> ${LOG_FILE}
cp src/tests/mas_align_script.sh .
ln -s _build/default/src/tests/mas_align.exe . >> ${LOG_FILE}
./mas_align_script.sh >> ${LOG_FILE}
rm mas_align_script.sh
unlink mas_align.exe

echo testing round trip graph construction >> ${LOG_FILE}
cp src/tests/round_trips.sh .
ln -s _build/default/src/tests/round_trip.exe .
./round_trips.sh >> ${LOG_FILE}
rm ./round_trips.sh
unlink round_trip.exe

echo testing merging of A >> ${LOG_FILE}
time _build/default/src/tests/merged_sensible.exe A >> ${LOG_FILE}

echo testing merging of B >> ${LOG_FILE}
time _build/default/src/tests/merged_sensible.exe B >> ${LOG_FILE}

echo testing merging of C >> ${LOG_FILE}
time _build/default/src/tests/merged_sensible.exe C >> ${LOG_FILE}

echo testing allele differences between A_gen >> ${LOG_FILE}
time _build/default/src/tests/allele_distances.exe A_gen >> ${LOG_FILE}

echo imputing A >> ${LOG_FILE}
time _build/default/src/app/mhc2gpdf.exe $IMGTHLA_DIR/alignments/A_gen.txt --no-pdf >> ${LOG_FILE}

echo imputing B >> ${LOG_FILE}
time _build/default/src/app/mhc2gpdf.exe $IMGTHLA_DIR/alignments/B_gen.txt --no-pdf >> ${LOG_FILE}

echo imputing C >> ${LOG_FILE}
time _build/default/src/app/mhc2gpdf.exe $IMGTHLA_DIR/alignments/C_gen.txt --no-pdf >> ${LOG_FILE}

echo Testing relative PHMM calculations >> ${LOG_FILE}
time _build/default/src/tests/relative_phmm.exe >> ${LOG_FILE}
