(* Multiple Sequence Alignment

   .... files (as provided by IMGT) looks a this:

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

We provide facilities for parsing them and generally working with
representations of this data.  *)

(** We keep track of all {alignment_element}'s with regard to their position,
    character by character, in the alignment. *)
type position = int

(** Elements that describe alignment sequence.

    All positions are given relative the global alignment. *)
type 'sr sequence = { start : position; s : 'sr }
and gap = { gstart : position; length : int }
and boundary = { idx : int; pos : position }
and 'sr alignment_element =
  | Start of position
  (** Start of sequence, contains position and name *)

  | End of position
  (** End of a sequence.*)

  | Boundary of boundary
  (** Boundary of count and position.

      Exon boundaries are represented in the alignment files with '|', we
      preserve them for downstream tools to strip out. The boundary count
      starts at 0. *)

  | Sequence of 'sr sequence
  (** Nucleotide sequence of position and sequence. *)

  | Gap of gap
  (** Gap of position and length. *)

val start_position : 'a alignment_element -> position

val end_position : ('a -> int) -> 'a alignment_element -> position

(** [al_el_to_string] converts an alignment element to string. *)
val al_el_to_string : string alignment_element -> string

val is_sequence : 'a alignment_element -> bool
val is_gap      : 'a alignment_element -> bool
val is_end      : 'a alignment_element -> bool
val is_start    : 'a alignment_element -> bool

module Parser : sig

  (** [find_align_date] parses an in channel of an alignment file up to the
      alignment date, or returns None if it is not found. *)
  val find_align_date : in_channel -> string option

  type result =
    { align_date  : string
    (** When the sequences were aligned by IMGT. *)

    ; reference   : string
    (** The name of the reference allele *)

    ; ref_elems   : string alignment_element list
    (** The sequence elements of the reference. *)

    ; alt_elems   : (string * string alignment_element list) list
    (** Sequence elements of alternative alleles in an associated list.*)
    }

  (* Report invariant parsing violations to stdout. *)
  val report : bool ref

  (** Parse an input channel. *)
  val from_in_channel : in_channel -> result

  (** Parse an alignment file. *)
  val from_file : string -> result

  val sequence_length : 'a alignment_element list -> int

end (* Parser *)

module Boundaries : sig

  type marker =
    { index       : int                                    (** Which segment? *)
    ; position    : position             (** Position of the boundary marker. *)
    ; length      : int            (** Length to next boundary, or 1 + the total
                                             length of the following segment. *)
    ; seq_length  : int
    }

  val marker_to_string : marker -> string

  val before_start : marker -> bool

  val to_boundary : offset:int -> marker -> 'a alignment_element

  val bounded : string alignment_element list -> (marker * string) list

end (* Boundaries *)

val split_sequence : string sequence -> pos:position -> string sequence * string sequence

val split_gap : gap -> pos:position -> gap * gap

val allele_sequences : reference:string alignment_element list ->
  allele:string alignment_element list -> (Boundaries.marker * string) list

(** [allele_sequence boundary_character reference allele ()] will convert the
    reference and allele alignment elements into a string representing the
    allele sequence.

    @param boundary_char if specified will be inserted into the resulting
    sequence.*)
val allele_sequence : ?boundary_char:char -> reference:string alignment_element list ->
    allele:string alignment_element list -> unit -> string

val reference_sequence_from_ref_alignment_elements : ?boundary_char:char ->
  string alignment_element list -> string

(** Will return the sequence of the reference. *)
val reference_sequence : ?boundary_char:char -> Parser.result -> string

val split_by_boundaries_rev : ' a alignment_element list ->
  'a alignment_element list list

(* Segments are parts of a DNA sequence that have relevant biological
   interpretation, specifically: UTR, intron and exons. We want to measure
   things (such as string distance) on a per-segment basis. *)
module Segments : sig

  type relationship =
    | Missing
    | Partial of int    (* sequence length *)
    | Full of int       (* sequence length, might be > reference length *)

  type 'a t =
    { seq_length    : int
    ; mismatches    : int
    ; relationship  : 'a
    }

  val distances : reference:string alignment_element list
                -> allele:string alignment_element list
                -> relationship t list

  val distances_between : reference:string alignment_element list
                        -> allele1:string alignment_element list
                        -> allele2:string alignment_element list
                        -> (relationship * relationship) t list

end (* Segments *)
