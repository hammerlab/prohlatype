
set -e

# TODO: pin this to a specific version!
git clone https://github.com/ANHIG/IMGTHLA.git


eval `opam config env`
export OPAMYES="true"

make setup

# for test coverage
opam install bisect_ppx ocveralls

make build
make tools
make covered_tests

export IMGTHLA_DIR="IMGTHLA"

echo current dir:
ls

case "$TEST" in
  parsing)
    echo testing parsing
    ./test_parsing.native
    export BISECT="true"
    ;;
  mas)
    echo mas alignment
    cp src/scripts/mas_align_script.sh .
    ./mas_align_script.sh
    export BISECT="true"
    ;;
  round)
    echo testing round trip graph construction
    cp src/scripts/round_trip_tests.sh .
    ./round_trip_tests.sh
    export BISECT="true"
    ;;
  mergeA)
    echo testing merging of A
    time ./merged_sensible_test.native A
    export BISECT="true"
    ;;
  mergeB)
    echo testing merging of B
    time ./merged_sensible_test.native B
    export BISECT="true"
    ;;
  mergeC)
    echo testing merging of C
    time ./merged_sensible_test.native C
    export BISECT="true"
    ;;
 alleleDiffA)
    echo testing allele differences between A_gen
    time ./test_allele_distances.native A_gen > A_sim.csv
    export BISECT="true"
    ;;
  biologicalKmers)
    echo testing that biological Kmers are found
    cp src/scripts/biological_kmer_tests.sh .
    time ./biological_kmer_tests.sh
    export BISECT="true"
    ;;
  impute)
    echo imputing A
    time ./mhc2gpdf.native --alignment $IMGTHLA_DIR/alignments/A_gen.txt --no-pdf
    echo imputing B
    time ./mhc2gpdf.native --alignment $IMGTHLA_DIR/alignments/B_gen.txt --no-pdf
    echo imputing C
    time ./mhc2gpdf.native --alignment $IMGTHLA_DIR/alignments/C_gen.txt --no-pdf
    ;;
  *)
    ;;
esac

if [ -n "$TEST" -a -n "$BISECT" ] ; then
  ocveralls --repo_token $COVERALLSTOKEN --git --send bisect*.out
fi
