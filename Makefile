PACKAGES=sosa ocamlgraph

.PHONY: clean build

default: 
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib mas_parser.cmo to_graph.cmo

cli:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib -I src/app cli.native

clean:
	ocamlbuild -clean

test:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) -I src/lib/ -I src/scripts/ test_parsing.native

