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

val equal_position : position -> position -> bool
val compare_position : position -> position -> int
val pp_position : Format.formatter -> position -> unit

type boundary_label =
  | UTR5
  | UTR3
  | Exon of int
  | Intron of int

val equal_boundary_label : boundary_label -> boundary_label -> bool
val compare_boundary_label : boundary_label -> boundary_label -> int
val boundary_label_to_string : boundary_label -> string
val boundary_label_to_short : boundary_label -> string

(** Elements that describe alignment sequence.

    All positions are given relative the global alignment. *)
type 'sr sequence = { start : position; s : 'sr }
and gap = { gstart : position; length : int }
and boundary = { label : boundary_label; pos : position }
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

val split_sequence : string sequence -> pos:position -> string sequence * string sequence

val split_gap : gap -> pos:position -> gap * gap

(* An alignment_sequence always begins with a Boundary. *)
type 'a alignment_sequence = 'a alignment_element list

val al_seq_to_string
    : ?enclose:bool
    -> ?sep:string
    -> string alignment_sequence
    -> string

module Alteration : sig

  type per_segment = { full : bool; type_ : boundary_label; start : position; end_ : position; }
  val equal_per_segment : per_segment -> per_segment -> bool
  val compare_per_segment : per_segment -> per_segment -> int
  val per_segment_to_string : per_segment -> string
  val per_segment_list_to_string : per_segment list -> bytes
  type t = { allele : string; why : string; distance : float; positions : per_segment list; }
  val to_yojson : t -> Yojson.Safe.json
  val of_yojson : Yojson.Safe.json -> (t, string) result
  val to_string : t -> string

end (* Alteration *)

module Parser : sig

  (** [find_header_lines] parses an in channel of an alignment file up to the
   *  release info and alignment date, or returns None if they are not found.
   *)
  val find_header_lines : in_channel -> (string * string) option

  (* Alternate allele information *)
  type alt =
    { allele : string
    (** The name: ex. "B*15:148" *)

    ; seq : string alignment_sequence
    (** The sequence elements. *)

    ; alters : Alteration.t list
    (* Alterations; a result with have these empty. See Alter_MSA. *)
    }

  val sort_alts_by_nomenclature : alt list -> alt list

  type result =
    { release     : string
    (** What IMGT/HLA database release version this file was parsed from. *)

    ; align_date  : string
    (** When the sequences were aligned by IMGT. *)

    ; locus       : Nomenclature.locus
    (** Which locus does this file represent. *)

    ; reference   : string
    (** The name of the reference allele *)

    ; ref_elems   : string alignment_sequence
    (** The sequence elements of the reference. *)

    ; alt_elems   : alt list
    (** Sequence elements of alternative alleles.*)
    }

  (* Report invariant parsing violations to stdout. *)
  val report : bool ref

  val lookup_allele : result -> string -> alt

  (** Does the source contain alignment for a gDNA (UTR's, Exons & Introns)
      or from cDNA (just Exons)? *)
  type boundary_schema =
    | GDNA
    | CDNA

  (** Parse an input channel. *)
  val from_in_channel : boundary_schema
                      -> Nomenclature.locus
                      -> in_channel
                      -> result

  (** Parse an alignment file.

      @param boundary_schema if not supplied is based upon the filename
             suffix (ie. 'gen' -> GDNA,  'nuc', 'prot' -> CNDA *)
  val from_file : ?boundary_schema:boundary_schema
                -> string
                -> result

  (** The number of positions between (all) [Start]'s and [End]'s. *)
  val sequence_length : 'a alignment_sequence -> int

  val in_order_invariant
    : string alignment_sequence
    -> (unit, string alignment_element * string alignment_element) Pervasives.result

end (* Parser *)

module Boundaries : sig

  type marker =
    { label       : boundary_label                  (* What kind of Boundary? *)
    ; position    : position              (* Position of the boundary marker. *)
    ; length      : int                           (* Length to next boundary. *)
    ; seq_length  : int
    }

  val marker_to_string : marker -> string

  (* val before_start : marker -> bool *)

  val to_boundary
      : ?offset:int
      -> marker
      -> 'a alignment_element

  val all_boundaries_before_start_or_end
      : 'a alignment_sequence
      -> 'a alignment_sequence

  val first_boundary_before_start
      : 'a alignment_sequence
      -> 'a alignment_sequence

  val grouped
      : string alignment_sequence
      -> (marker * string alignment_sequence) list

  val ungrouped
      : (marker * string alignment_sequence) list
      -> string alignment_sequence

end (* Boundaries *)

(* Compute per Boundary sequence information. *)
val allele_sequences : reference:string alignment_sequence
                     -> allele:string alignment_sequence
                     -> (Boundaries.marker * string) list

(** [allele_sequence boundary_character reference allele ()] will convert the
    reference and allele alignment elements into a string representing the
    allele sequence.

    @param boundary_char if specified will be inserted into the resulting
    sequence.*)
val allele_sequence : ?boundary_char:char -> reference:string alignment_sequence ->
    allele:string alignment_sequence -> unit -> string

val reference_sequence_from_ref_alignment_elements : ?boundary_char:char ->
  string alignment_sequence -> string

(** Will return the sequence of the reference. *)
val reference_sequence : ?boundary_char:char -> Parser.result -> string

(* Segments are parts of a DNA sequence that have relevant biological
   interpretation, specifically: UTR, intron and exons. We want to measure
   things (such as string distance) on a per-segment basis. *)
module Segments : sig

  type relationship =
    | Missing
    | Partial of int                                       (* sequence length *)
    | Full of int             (* sequence length, might be > reference length *)

  type 'a t =
    { seq_length    : int
    ; mismatches    : int
    ; relationship  : 'a
    }

  val distances : reference:string alignment_sequence
                -> allele:string alignment_sequence
                -> relationship t list

  val distances_between : reference:string alignment_sequence
                        -> allele1:string alignment_sequence
                        -> allele2:string alignment_sequence
                        -> (relationship * relationship) t list

end (* Segments *)
