
set -e

# TODO: pin this to a specific version!
git clone -b 3320 https://github.com/ANHIG/IMGTHLA.git IMGTHLA

# for test coverage
${HOME}/opam install ocveralls parmap oml core_bench qcheck

make
make apps

export IMGTHLA_DIR="IMGTHLA"

echo current dir:
ls

case "$TEST" in
  scraps)
    make scraps
    ;;
  parsing)
    make covered_tests
    echo testing parsing
    _build/default/src/tests/alignment_parsing.exe
    ;;
  mas)
    echo mas alignment
    make covered_tests
    cp src/tests/mas_align_script.sh .
    cp _build/default/src/tests/mas_align.exe .
    ./mas_align_script.sh
    ;;
  round)
    echo testing round trip graph construction
    make covered_tests
    cp src/tests/round_trips.sh .
    cp _build/default/src/tests/round_trip.exe .
    ./round_trips.sh
    ;;
  mergeA)
    echo testing merging of A
    make tests
    time _build/default/src/tests/merged_sensible.exe A > merged_A.log
    ;;
  mergeB)
    echo testing merging of B
    make tests
    time _build/default/src/tests/merged_sensible.exe B > merged_B.log
    ;;
  mergeC)
    echo testing merging of C
    make tests
    time _build/default/src/tests/merged_sensible.exe C > merged_C.log
    ;;
 alleleDiffA)
    echo testing allele differences between A_gen
    make covered_tests
    time _build/default/src/tests/allele_distances.exe A_gen > A_sim.csv
    ;;
  impute)
    echo imputing A
    time _build/default/src/app/mhc2gpdf.exe $IMGTHLA_DIR/alignments/A_gen.txt --no-pdf
    echo imputing B
    time _build/default/src/app/mhc2gpdf.exe $IMGTHLA_DIR/alignments/B_gen.txt --no-pdf
    echo imputing C
    time _build/default/src/app/mhc2gpdf.exe $IMGTHLA_DIR/alignments/C_gen.txt --no-pdf
    ;;
  relPhmm)
    echo Testing relative PHMM calculations
    time _build/default/src/tests/relative_phmm.exe > relative_phmm.log
    ;;
  *)
    ;;
esac
