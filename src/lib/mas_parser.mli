(* Parse (or Lex, I don't want to get into the details of the exact difference)
   HLA alignment as provided by IMGT. The file looks a this:

---------------------------START FILE ------------------------------------------
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
 A*80:01:01:02     ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------


 gDNA              -200
                   |
 A*01:01:01:01     GGGATTCCCC AACTCCGCAG TTTCTTTTCT CCCTCTCCCA ACCTACGTAG GGTCCTTCAT CCTGGATACT CACGACGCGG ACCCAGTTCT CACTCCCATT
 A
   ...

 A*80:01:01:02     ********** ********** ********** ********** ********** ********** ********** ********** ********** *


Please see http://hla.alleles.org/terms.html for terms of use.
---------------------------END FILE --------------------------------------------

*)

(** Elements that describe alignment sequence.

    All positions are given reletive to reference. *)
type 'sr seq_elems =
  | Start of int * string
  (** Start of sequence, contains position and name *)

  | End of int
  (** End of a sequence.*)

  | Boundary of (int * int)
  (** Boundary of count and position.

      Exon boundaries are represented in the alignment files with '|', we
      preserve them for downstream tools to strip out. The boundary count
      starts at 0. *)

  | Nuc of (int * 'sr)
  (** Nucleotide sequence of position and sequence. *)

  | Gap of (int * int)
  (** Gap of position and length. *)

type result =
  { reference : string
  (** The name of the reference allele *)

  ; ref_elems : string seq_elems list
  (** The sequence elements of the reference. *)

  ; alt_elems : (string * string seq_elems list) list
  (** Seqeucen elements of alternative alleles in an associated list.*)
  }

(** Parse an input channel. *)
val from_in_channel : in_channel -> result

(** Parse an alignment file. *)
val from_file : string -> result
