.PHONY: default setup clean build tests covered_tests

default: build

build:
	jbuilder build

setup:
	opam install --deps-only ./prohlatype.opam

clean:
	jbuilder clean

apps:
	jbuilder build @apps

scraps:
	jbuilder build @scraps

tests:
	jbuilder build @tests

all:
	jbuilder build @apps @scraps @tests

release:
	patch -p1 < tools/static_patch &&\
	jbuilder build @apps &&\
	zip -l prohlatype.zip _build/default/src/app/*.exe

## Coverage

covered_tests:
	BISECT_ENABLE=Yes jbuilder build @tests

report_dir:
	mkdir report_dir

report: report_dir
	bisect-ppx-report -html report_dir $(shell ls -t bisect*.out | head -1)

clean_reports:
	rm -rf report_dir bisect*.out
