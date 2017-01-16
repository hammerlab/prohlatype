PACKAGES=unix ppx_deriving.std nonstd sosa ocamlgraph cmdliner biocaml.unix
SETUP_PACKAGE_NAMES=ocamlfind ocamlbuild ppx_deriving nonstd sosa ocamlgraph cmdliner biocaml
TOOLS=mhc2gpdf type align2fasta allele_distances
TESTS=test_parsing round_trip same_alignments_test check_multiple adjacents benchmark_k merged_sensible_test mas_align_test test_allele_distances


.PHONY: default setup clean build tools tests covered_tests

default: build

build:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib prohlatype.cma

setup:
	opam install ocamlfind ocamlbuild $(SETUP_PACKAGE_NAMES)

#cli:
#	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib -I src/app cli.native

clean:
	ocamlbuild -clean

## Tests --- this might not scale

parsing:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts test_parsing.native

round_trip:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts round_trip.native

same_alignment:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts same_alignments_test.native

check_multiple:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts check_multiple.native

adjacents:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts adjacents.native

benchmark_k:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts benchmark_k.native

merged_tests:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts merged_sensible_test.native

mas_align_test:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts mas_align_test.native

exon_edit_distances:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -package edist -I src/lib/ -I src/scripts exon_edit_distances.native

test_allele_distances:
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts test_allele_distances.native

tests:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts $(foreach t, $(TESTS),$(t).native)

covered_tests:
	ocamlbuild -use-ocamlfind -package bisect_ppx $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts $(foreach t, $(TESTS),$(t).native)

## Tools:

mhc2gpdf:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app mhc2gpdf.native

type:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app type.native

explore:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -package oml -I src/lib/ -I src/app explore_alignment.native

align2fasta:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app align2fasta.native

allele_distances:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app allele_distances.native

tools:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app $(foreach t, $(TOOLS),$(t).native)

## Coverage

report_dir:
	mkdir report_dir

report: report_dir
	bisect-ppx-report -html report_dir $(shell ls -t bisect*.out | head -1)

clean_reports:
	rm -rf report_dir bisect*.out

## Throw Away Scripts
# ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -package future.unix -package biocaml -package biocaml.unix -I src/lib/ -I src/app type.native
