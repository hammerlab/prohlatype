.PHONY: default setup clean build tests covered_tests

default: build

build:
	dune build

setup:
	opam install --deps-only ./prohlatype.opam

clean:
	dune clean

apps:
	dune build @apps

scraps:
	dune build @scraps

tests:
	dune build @tests

all:
	dune build @apps @scraps @tests

release:
	patch -p1 < tools/static_patch &&\
	dune build @apps &&\
	zip -l prohlatype.zip _build/default/src/app/*.exe

## Coverage

covered_tests:
	make clean &&\
	BISECT_ENABLE=Yes dune build @tests @apps

report_dir:
	mkdir report_dir

report: report_dir
	cd _build/default && bisect-ppx-report -html ../../report_dir ../../bisect*.out && cd -

clean_reports:
	rm -rf report_dir bisect*.out
