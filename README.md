[![Build Status](https://travis-ci.org/hammerlab/prohlatype.svg?branch=master)](https://travis-ci.org/hammerlab/prohlatype)
[![Coverage Status](https://coveralls.io/repos/hammerlab/prohlatype/badge.svg?branch=HEAD&service=github)](https://coveralls.io/github/hammerlab/prohlatype?branch=HEAD)

Probabilistic HLA Typing
------------------------

Paper: [Prohlatype: A Probabilistic Framework for HLA Typing](https://doi.org/10.1101/244962)

This project provides a set of tools to calculate the full posterior
distribution of HLA types given read data.

Instead of:

```
	A1  	A2  	B1  	B2  	C1	    C2  	Reads	Objective
0	A*31:01	A*02:01	B*45:01	B*15:03	C*16:01	C*02:10	538.0	513.79

```

one can calculate:

| Allele 1      | Allele 2          | Log P       |     P |
|---------------|-------------------|-------------|-------|
|  A*02:05:01:01|	        A*30:114|	-23046.81 | 0.5000|
|  A*02:05:01:01|	      A*30:01:01|	-23046.81 | 0.5000|
|  A*02:05:01:01|	        A*30:106|   -23103.15 | 0.0000|
|  A*02:05:01:02|           A*30:114|   -23146.35 | 0.0000|
| ... | | |
|        B*07:36|	   B*57:03:01:02|	-13717.33 | 0.5000|
|        B*07:36|      B*57:03:01:01|	-13717.33 | 0.5000|
|        B*07:36|	      B*57:03:03|	-13804.74 | 0.0000|
|       B*27:157|	   B*57:03:01:02|	-13816.17 | 0.0000|
| ... | | |
|       C*06:103|	         C*18:10|	-11936.35 | 0.3338|
|       C*06:103|	         C*18:02|	-11936.36 | 0.3331|
|       C*06:103|	         C*18:01|	-11936.36 | 0.3331|
|       C*15:102|	         C*18:02|	-11951.72 | 0.0000|


## How:

### There are three options to obtain the software:

  1. If you are running on Linux, standalone binaries are available with each [release](https://github.com/hammerlab/prohlatype/releases).
  2. Use the linked Dockerfile (TODO).
  3. Build the software from source:

      a. Install [opam](https://opam.ocaml.org/).

      b. Make sure that the opam packages are up to date:

        $ opam update

      c. Make sure that you're on the relevant compiler:

        $ opam switch 4.05.0
        $ eval `opam config env`

      d. Get source:

        $ git clone https://github.com/hammerlab/prohlatype.git prohlatype
        $ cd prohlatype

      e. Install the dependent packages:

        $ make setup

      f. Build the programs (afterwards they'll be in `_build/default/src/apps`):

        $ make

### Make sure that you have [IMGT/HLA](https://github.com/ANHIG/IMGTHLA) available:

`$ git clone https://github.com/ANHIG/IMGTHLA.git imgthla`

### "Prohla"-typing:

  1. Create an imputed HLA reference sequence via `align2fasta`.
     This step makes sure that all alleles have sequence information that spans
     the entire locus. This way, reads that originate from a region for which
     we normally do **not** have sequence information will still align (in the
     next filtering step), albeit poorly:

        $ align2fasta path-to-imgthla/alignments -o imputed_hla_class_I.fasta

     This step needs to be performed only once, per each IMGT version.
     Run `$align2fasta --help` for further information.

  2. Filter your data against the reference, by first aligning. Ex:

        $ bwa mem imputed_hla_class_I.fasta ${SAMPLE}.fastq | \
            samtools view -F 4 -bT imputed_hla_class_I.fasta -o ${SAMPLE}.bam

     While fundamentally, the algorithms here are *alignment* based. They're
     too slow to run for all sequences. Sequences that do not originate from
     the HLA-region would just act as background noice.

  3. and then convert aligned reads back to FASTQ:

        $ samtools fastq ${SAMPLE}.bam > ${SAMPLE}_filtered.fastq

  4. Infer types:

        $ multi_par path-to-imgthla/aignments ${SAMPLE}_filtered.fast -o ${SAMPLE}_output.tsv

    See `$ multi_par --help` for further detail and optimizations.
