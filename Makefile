PACKAGES=nonstd sosa ocamlgraph

.PHONY: default setup clean build

default: build

build:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib mas_parser.cmo to_graph.cmo

setup:
	opam install $(PACKAGES)

cli:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib -I src/app cli.native

clean:
	ocamlbuild -clean

test:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts/ test_parsing.native
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts/ test_graphing.native

