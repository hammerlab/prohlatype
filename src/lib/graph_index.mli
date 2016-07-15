
open Ref_graph

(** Essentially private API that will be hidden. *)
type ksublength = [ `Part of int | `Whole ]
type 'a kmer_substring = { index : int; length : 'a; }
val whole : index:int -> [> `Whole ] kmer_substring
val part : index:int -> int -> [> `Part of int ] kmer_substring
val fold_over_kmers_in_string :
  k:int ->
  f:('a -> [> `Part of int | `Whole ] kmer_substring -> 'a) ->
  init:'a -> string -> 'a

val fold_over_kmers_in_graph :
  k:int ->
  f:('a -> alignment_position * sequence -> [> `Whole ] kmer_substring -> 'a) ->
  init:'a ->
  extend:(alignment_position * sequence -> [> `Part of int ] kmer_substring -> 'b option -> 'b) ->
  close:('a -> alignment_position * sequence -> [> `Part of int ] kmer_substring -> 'b -> 'a) ->
  Ref_graph.t -> 'a

val kmer_counts : k:int -> Ref_graph.t -> int Kmer_table.t

(** Public API *)

(** How we describe positions in the string graph. *)
type position =
  { alignment : alignment_position    (* Where, in the alignment the sequence starts. *)
  ; sequence  : sequence              (* The sequence. *)
  ; offset    : int                   (* An offset into the the sequence. *)
  }

val specific_position : Ref_graph.t -> Alleles.allele -> alignment_position ->
  (position, string) result
    
(** A graph index *)
type t = position list Kmer_table.t

(** Index a graph. *)
val create : k:int -> Ref_graph.t -> t

(** Returns where the kmer starts the given sequence. *)
val starting_with : t -> string -> (position list, string) result

(** What [k] are we using to index the graph. *)
val k : t -> int

val lookup: ?n:int -> t ->  string -> (position list, string) result
