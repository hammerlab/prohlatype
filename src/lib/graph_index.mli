
open Ref_graph

type t

type position = alignment_position * sequence * int (* offset into sequence *)

val create : k:int -> Ref_graph.G.t -> t

val starting_with : t -> string -> (position list, string) result
