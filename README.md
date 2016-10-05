[![Build Status](https://travis-ci.org/hammerlab/prohlatype.svg?branch=master)](https://travis-ci.org/hammerlab/prohlatype)

HLA Typing 
----------

#### Status 

We are using the multiple allele alignment files from IMGT:

```
HLA-A Genomic Sequence Alignments
IPD-IMGT/HLA Release: 3.24.0
Sequences Aligned: 2016 April 15
Steven GE Marsh, Anthony Nolan Research Institute.
Please see http://hla.alleles.org/terms.html for terms of use.
                   
 gDNA              -300
                   |
 A*01:01:01:01     CAGGAGCAGA GGGGTCAGGG CGAAGTCCCA GGGCCCCAGG CGTGGCTCTC AGGGTCTCAG GCCCCGAAGG CGGTGTATGG ATTGGGGAGT CCCAGCCTTG 
 A*01:01:01:02N    ********** ********** ********** ********** ********** ********** ********** *********- ---------- ---------- 
 A*01:01:01:03     ********** ********** ********** ********** ********** ********** ********** ********** ********** ********** 
 ...

```

as a starting point. This repository cotains code to parses those files into graphs.

#### Utilities

The `mhc2gpdf` utility can be used to create arbitrary graphs.

```
$ ./mhc2gpdf.native -o ref_vs0204_vs0301 -f ../foreign/IMGTHLA/alignments/A_gen.txt -a A*02:04 -a A*03:01:01:01
```

creates [Reference vs A0204 vs 0301](regression_tests/ref_vs0204_vs0301.pdf).


#### Building and Getting started

Assuming that you have [opam](http://opam.ocaml.org/) installed on relatively
recent version of `OCaml`.

```
$ make setup
$ make mhc2gpdf
$ ./mhc2gpdf.native --help
```

