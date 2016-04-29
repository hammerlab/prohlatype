PACKAGES=sosa ocamlgraph

.PHONY: clean build

default: 
	ocamlbuild -use-ocamlfind $(foreach package, $(PACKAGES),-package $(package)) parser.cmo

clean:
	ocamlbuild -clean


