(** Distances

  Different logic and methods for thinking about the distance between HLA
  alleles. For constructing a uniform view of a specific HLA-loci we need to
  sometimes borrow (loosely impute) missing data from other alleles. We want
  to borrow this data from alleles that are close to another allele. This
  closeness is determined by minimizing some distance. *)

open Util

type aligned_sequence = string MSA.alignment_element list

(** The two types of distance logic that we currently support. *)
type logic =
  | Trie
  (** For a specific HLA loci construct a trie of the alleles based upon their
      (upto) 8 digit specification. The distance is always fixed at 1 and only
      the closest allele in the Trie is returned. *)

  | WeightedPerSegment
  (** Break apart alleles, as presented via their [aligned_sequence] into their
      segments (biologically relevant components: UTR, Exon & Intron).
      Distance is the number of mismatches between two alleles in a given segment
      weighted by segment sequence length. *)

val pp_logic : Format.formatter -> logic -> unit
val show_logic : logic -> string

(** Compute the distances for one specific allele. *)
val one : reference:string
        -> reference_sequence:aligned_sequence
        -> allele:aligned_sequence
        -> candidates:aligned_sequence StringMap.t
        -> logic
        -> ((string * float) list, string) result

(** Compute the distances for all the alleles [targets]. *)
val compute : reference:string
            -> reference_sequence:aligned_sequence
            -> targets:aligned_sequence StringMap.t
            -> candidates:aligned_sequence StringMap.t
            -> logic
            -> ((string * float) list StringMap.t, string) result
