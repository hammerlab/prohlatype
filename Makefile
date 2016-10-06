PACKAGES=unix ppx_deriving.std nonstd sosa ocamlgraph cmdliner extlib biocaml.unix
SETUP_PACKAGE_NAMES=ocamlfind ocamlbuild ppx_deriving nonstd sosa ocamlgraph cmdliner extlib biocaml
TOOLS=mhc2gpdf type
TESTS=test_parsing round_trip same_alignments_test check_multiple adjacents benchmark_k

.PHONY: default setup clean build tools tests

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

tests:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts $(foreach t, $(TESTS),$(t).native)

## Tools:

mhc2gpdf:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app mhc2gpdf.native

type:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app type.native

explore:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -package oml -I src/lib/ -I src/app explore_alignment.native

tools:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app $(foreach t, $(TOOLS),$(t).native)

## Throw Away Scripts


# ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -package future.unix -package biocaml -package biocaml.unix -I src/lib/ -I src/app type.native

