(** A partition map is a data structure for a map over a partition of elements.

  Specifically, if we know (and specifically can enumerate) the elements of
  a set this data structure allows a mapping from elements to the values.
  Internally, it maintains partitions: representations of sets of the elements
  that partitions the entire universe. The most important operation is the
  {merge} of 2 (or more {merge4}) such partition maps.
*)

(* We construct partition map's in {descending} order than then convert them
   into the {ascending} order for merging. *)
type ascending
type descending

type ('o, +'a) t

(* Empty constructors. *)
val empty_d : (descending, 'a) t
val empty_a : (ascending, 'a) t

(* Initializers. These take a value and either assume a entry (the 'first' one in
   the descending case) or all of them (pass the size of the partition) in the
   ascending case. *)
val init_first_d : 'a -> (descending, 'a) t
val init_all_a : size:int -> 'a -> (ascending, 'a) t

val to_string : (_, 'a) t -> ('a -> string) -> string

(* Conversions. *)
val ascending : (descending, 'a) t -> (ascending, 'a) t

(* Observe a value for the next element. *)
val add : 'a -> (descending, 'a) t -> (descending, 'a) t

(* [get t i] returns the value associated  with the [i]'th element.

   @raise {Not_found} if [i] is outside the range [0, (size t)). *)
val get : (ascending, 'a) t -> int -> 'a

(* Merge partition maps. Technically these are {map}'s but they are purposefully
  named merge since they're only implemented for {ascending} partition maps. *)
val merge : (ascending, 'a) t
          -> (ascending, 'b) t
          -> ('a -> 'b -> 'c)
          -> (ascending, 'c) t

val merge4 : (ascending, 'a) t
            -> (ascending, 'b) t
            -> (ascending, 'c) t
            -> (ascending, 'd) t
            -> ('a -> 'b -> 'c -> 'd -> 'e)
            -> (ascending, 'e) t

(* Fold over the values. *)
val fold_values : (_, 'a) t
                -> init:'b
                -> f:('b -> 'a -> 'b)
                -> 'b

(* Map the values, the internal storage doesn't change. *)
val map : ('o, 'a) t -> f:('a -> 'b) -> ('o, 'b) t

(* Iterate over the entries and values. *)
val iter_set : ('o, 'a) t -> f:(int -> 'a -> unit) -> unit


(** Diagnostic methods. These are not strictly necessary for the operations of
    the Parametric PHMM but are exposed here for interactive use. *)

(* The size of the partition. Specifically, if [size t = n] then [get t i] will
   succeed for [0, n).  *)
val size : (_, 'a) t -> int

(* The number of unique elements in the underlying assoc . *)
val length : (_, 'a) t -> int

val descending : (ascending, 'a) t -> (descending, 'a) t

(** Set a value. *)
val set : (ascending, 'a) t -> int -> 'a -> (ascending, 'a) t
