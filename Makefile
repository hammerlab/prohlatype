PACKAGES=sosa ocamlgraph

.PHONY: clean build

default: 
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) to_graph.cmo

cli:
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) cli.native

clean:
	ocamlbuild -clean


