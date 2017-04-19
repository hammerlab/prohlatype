This branch contains the all of the previous graph alignment code before we switch
to using Parametric PHMMs.

[![Build Status](https://travis-ci.org/hammerlab/prohlatype.svg?branch=master)](https://travis-ci.org/hammerlab/prohlatype)
[![Coverage Status](https://coveralls.io/repos/hammerlab/prohlatype/badge.svg?branch=HEAD&service=github)](https://coveralls.io/github/hammerlab/prohlatype?branch=HEAD)

Probabilistic HLA Typing
------------------------

### In Brief

The goal of this project is to achieve HLA typing by inferring the posterior distribution of types, for a given gene, based upon reads using Bayes's Theorem.

#### Status

  1. Alignment file parsing -> **works!**

     We'll move away from this data source at some point as it is meant for human consumption.

  2. Graph representation of HLA alignment -> works pretty _well_ but open for
     improvement. In particular the following issues remain:

     a. Indexing the graph is presently inefficient and awkward. We're using a
        k-mer index to represent positions. In particular, they encode the full
        combinatorial variation in the alleles, but not the biologically
        relevant ones. Compressed text indices would probably work better and
        give us nicer properties such as being able to uniquely determine
        positions.

     b. At present there are accessory structures that don't fit neatly into
        the data structure (such as the [boundary array](src/lib/ref_graph.ml#L108))
        and superfluous computations that could be performed at compile time
        (such as the [adjacency calculations](src/lib/ref_graph.ml#L1445)).

  3. Alignment to the graph -> **works** pretty well.
     The algorithm probably doesn't have the cleanest
     [parameterization](src/lib/alignment.ml#L105) as this is a bit of a
     work-in-progress to determine what is actually necessary for downstream
     features.

  4. Inference (aka. typing), needs refinements to reduce variance.


### In Depth

#### Alignment file parsing and graph construction

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

as a starting point. These files are parsed by the [Mas_parser](src/lib/mas_parser.mli), and graph construction is performed in [Ref_graph](src/lib/ref_graph.ml)


##### mhc2gpdf

The `mhc2gpdf` utility can be used to create arbitrary graphs:

```
$ ./mhc2gpdf.native -f /path/to/IMGT/alignments/dir/A_nuc.txt -o demo --allele "A\*02:01:01:01" --allele-regex "A\\*11:" --no-reference
```

creates [demo.pdf](demo/demo.pdf).
This tool is described in [this](http://www.hammerlab.org/2016/11/01/graph-visualizations-of-mhc-alleles/)
blog post.


#### Building and Development started

Assuming that you have [opam](http://opam.ocaml.org/) installed on relatively
recent version of `OCaml`.

```
$ make setup                    # opam install dependencies
$ make                          # build library
$ make tools                    # build tools such as mhc2gpdf and type
$ make tests                    # build regression and unit tests
```
