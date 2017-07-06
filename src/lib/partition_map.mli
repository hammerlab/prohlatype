(** A partition map is a data structure for a map over a partition of a set.

  Specifically, if we know (and specifically can enumerate) the elements of
  a set this data structure allows a mapping of entries in the set to values.
  The most important, and fast, operation is the {merge} of 2
  (or more {merge4}) such partition maps.  
*)

type ascending
type descending

type ('o, +'a) t

val empty : (descending, 'a) t
val place_holder : (ascending, 'a) t

val first_d : 'a -> (descending, 'a) t

val init : int -> 'a -> (ascending, 'a) t

val to_string : (_, 'a) t -> ('a -> string) -> string

val size_a : (ascending, 'a) t -> int
val size_d : (descending, 'a) t -> int

val length : (_, 'a) t -> int

val ascending : (descending, 'a) t -> (ascending, 'a) t
val descending : (ascending, 'a) t -> (descending, 'a) t

val add : 'a -> (descending, 'a) t -> (descending, 'a) t

val get : (ascending, 'a) t -> int -> 'a
val set : (ascending, 'a) t -> int -> 'a -> (ascending, 'a) t

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

val fold_values : (_, 'a) t
          -> init:'b
          -> f:('b -> 'a -> 'b)
          -> 'b

val map : ('c, 'a) t -> f:('a -> 'b) -> ('c, 'b) t

val iter_set : ('c, 'a) t -> f:(int -> 'a -> unit) -> unit

