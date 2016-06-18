
type 'a t

val k : 'a t -> int

val make : int -> 'a -> 'a t

val init : int -> (string -> 'a) -> 'a t

val update : ('a -> 'a) -> 'a t -> string -> unit

val update_index : ('a -> 'b -> 'b) -> 'b t -> 'a -> int -> unit

val distr : int t -> int array

val lookup : 'a t -> string -> 'a

val cross_boundary : ('a * string * int) list t -> (string * ('a * string * int) list) list
