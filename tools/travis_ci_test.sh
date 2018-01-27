
set -e

# TODO: pin this to a specific version!
git clone https://github.com/ANHIG/IMGTHLA.git

#eval `opam config env`
#export OPAMYES="true"

#opam install ppx_deriving.4.2.1 ppx_deriving_yojson.3.1 nonstd.0.0.3 sosa.0.3.0 ocamlgraph.1.8.8 cmdliner.1.0.2 biocaml.0.8.0 parany.3.0.0

# for test coverage
${HOME}/opam install bisect_ppx ocveralls

make
make apps

export IMGTHLA_DIR="IMGTHLA"

echo current dir:
ls

case "$TEST" in
  scraps)
    # TODO: Figure out a way to keep this in sync with jbuilder
    ${HOME}/opam install parmap oml core_bench
    make scraps
    ;;
  parsing)
    make covered_tests
    echo testing parsing
    _build/default/src/tests/alignment_parsing.exe
    export BISECT="true"
    ;;
  mas)
    echo mas alignment
    make covered_tests
    cp src/tests/mas_align_script.sh .
    cp _build/default/src/tests/mas_align.exe .
    ./mas_align_script.sh
    export BISECT="true"
    ;;
  round)
    echo testing round trip graph construction
    make covered_tests
    cp src/tests/round_trips.sh .
    cp _build/default/src/tests/round_trip.exe .
    ./round_trips.sh
    export BISECT="true"
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
    export BISECT="true"
    ;;
  impute)
    echo imputing A
    time _build/default/src/app/mhc2gpdf.exe $IMGTHLA_DIR/alignments/A_gen.txt --no-pdf
    echo imputing B
    time _build/default/src/app/mhc2gpdf.exe $IMGTHLA_DIR/alignments/B_gen.txt --no-pdf
    echo imputing C
    time _build/default/src/app/mhc2gpdf.exe $IMGTHLA_DIR/alignments/C_gen.txt --no-pdf
    ;;
  *)
    ;;
esac

if [ -n "$TEST" -a -n "$BISECT" ] ; then
  ocveralls --repo_token $COVERALLSTOKEN --git --send bisect*.out
fi

