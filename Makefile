PACKAGES=unix ppx_deriving.std ppx_deriving_yojson nonstd sosa ocamlgraph cmdliner biocaml.unix parany
SETUP_PACKAGE_NAMES=ppx_deriving.4.2.1 ppx_deriving_yojson.3.1 nonstd.0.0.3 sosa.0.3.0 ocamlgraph.1.8.8 cmdliner.1.0.2 biocaml.0.8.0 parany.3.0.0
TOOLS=mhc2gpdf par_type multi_par align2fasta mpjson2tsv
TESTS=test_parsing round_trip benchmark_k merged_sensible_test mas_align_test test_allele_distances biological_kmers


.PHONY: default setup clean build tools tests covered_tests

default: build

all: build tools tests

build:
	jbuilder build

setup:
	opam install --deps-only

clean:
	jbuilder clean

apps:
	jbuilder build @apps  

scraps:
	jbuilder build @scraps  

tests:
	jbuilder build @tests

covered_tests:
	BISECT_ENABLE=Yes jbuilder build @tests

## Coverage

report_dir:
	mkdir report_dir

report: report_dir
	bisect-ppx-report -html report_dir $(shell ls -t bisect*.out | head -1)

clean_reports:
	rm -rf report_dir bisect*.out
