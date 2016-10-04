
set -e

# TODO: pin this to a specific version!
git clone https://github.com/ANHIG/IMGTHLA.git

eval `opam config env`
export OPAMYES="true"

make setup
make build
make tools
make tests

export IMGTHLA_DIR="IMGTHLA"

echo current dir:
ls

echo testing parsing
./test_parsing.native

echo testing round trip graph construction
cp src/scripts/round_trip_test.sh .
./round_trip_test.sh

echo testing same alignment
./same_alignments_test.native

echo testing adjcent finding
cp src/scripts/adjacent_tests.sh .
./adjacent_tests.sh
