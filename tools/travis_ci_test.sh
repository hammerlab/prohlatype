
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
    ;;
  round)
    echo testing round trip graph construction
    cp src/scripts/round_trip_tests.sh .
    ./round_trip_tests.sh
    ;;
  mergeA)
    echo testing merging of A
    ./merged_sensible_test.native A
    ;;
  mergeB)
    echo testing merging of B
    ./merged_sensible_test.native B
    ;;
  mergeC)
    echo testing merging of C
    ./merged_sensible_test.native C
    ;;
  adjA)
    echo testing adjcent finding for A_gen only
    ./adjacents.native A_gen
    ;;
  *)
    ;;
esac

#echo testing same alignment
#./same_alignments_test.native

# Full adjacents calculation takes too long
#cp src/scripts/adjacent_tests.sh .
#./adjacent_tests.sh

ocveralls --repo_token $COVERALLSTOKEN --git --send bisect*.out
