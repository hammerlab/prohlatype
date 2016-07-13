PACKAGES=ppx_deriving.std nonstd sosa ocamlgraph cmdliner extlib

.PHONY: default setup clean build

default: build

build:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib prohlatype.cma

setup:
	opam install ocamlfind ocamlbuild $(PACKAGES)

cli:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib -I src/app cli.native

clean:
	ocamlbuild -clean

tests:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts test_parsing.native
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts round_trip.native
	ocamlbuild -use-ocamlfind -package unix $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts same_alignments_test.native


## Tools:

mhc2gpdf:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app mhc2gpdf.native

type:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/app type.native

tools: mhc2gpdf type

## Throw Away Scripts

benchmark_k:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts benchmark_k.native


# ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -package future.unix -package biocaml -package biocaml.unix -I src/lib/ -I src/app type.native

