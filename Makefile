PACKAGES=unix ppx_deriving.std nonstd sosa ocamlgraph cmdliner biocaml.unix parany
SETUP_PACKAGE_NAMES=ocamlfind ocamlbuild ppx_deriving.4.1 nonstd.0.0.2 sosa.0.2.0 ocamlgraph.1.8.6 cmdliner.1.0.0 biocaml.0.6.0
TOOLS=mhc2gpdf par_type multi_par align2fasta allele_distances
TESTS=test_parsing round_trip benchmark_k merged_sensible_test mas_align_test test_allele_distances biological_kmers


.PHONY: default setup clean build tools tests covered_tests

default: build

all: build tools tests

build:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib prohlatype.cma

setup:
	opam install ocamlfind ocamlbuild $(SETUP_PACKAGE_NAMES)

#cli:
#	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib -I src/app cli.native

clean:
	ocamlbuild -clean

## Tests --- this might not scale

test_parsing:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts test_parsing.native

round_trip:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts round_trip.native

benchmark_k:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts benchmark_k.native

merged_sensible_test:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts merged_sensible_test.native

mas_align_test:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts mas_align_test.native

exon_edit_distances:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -package edist -I src/lib/ -I src/scripts exon_edit_distances.native

test_allele_distances:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts test_allele_distances.native

biological_kmers:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts biological_kmers.native

tests:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts $(foreach t, $(TESTS),$(t).native)

covered_tests:
	ocamlbuild -use-ocamlfind -package bisect_ppx $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts $(foreach t, $(TESTS),$(t).native)

## Tools:

mhc2gpdf:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app mhc2gpdf.native

par_type:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app par_type.native

multi_par:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app multi_par.native

explore:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -package oml -I src/lib/ -I src/app explore_alignment.native

align2fasta:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app align2fasta.native

allele_distances:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app allele_distances.native

time_phmm:
	corebuild -package core_bench -I src/app time_phmm.native

tools:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app $(foreach t, $(TOOLS),$(t).native)

## Coverage

report_dir:
	mkdir report_dir

report: report_dir
	bisect-ppx-report -html report_dir $(shell ls -t bisect*.out | head -1)

clean_reports:
	rm -rf report_dir bisect*.out
