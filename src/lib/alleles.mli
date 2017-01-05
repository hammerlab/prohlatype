(** To keep an efficient representation of which allele's are present along an
    edge we split the data in two. The first pieces is a constant index, while
    the second is an efficent representation of an allele set or map that relies
    upon the index. The idea is that we have only one index (per graph) and many
    edge sets/maps. Alleles are represented by strings. *)

type allele = string

val compare : allele -> allele -> int
val equal : allele -> allele -> bool

type index

(** [index list_of_alleles] will create an [index]. *)
val index : allele list -> index

val to_alleles : index -> allele list

module Set : sig

  type t

  (** [init index] create an empty set. *)
  val init : index -> t

  (** [singleton index allele] will create an edge set with just [allele]. *)
  val singleton : index -> allele -> t

  (** [set index t allele] will make sure that [allele] is
      in [t], specifically [is_set index t allele] will be [true]. *)
  val set : index -> t -> allele -> t

  val unite : into:t -> t -> unit

  (** [clear index t allele] will make sure that [allele] is not
      in [t], specifically [is_set index t allele] will be
      [false]. *)
  val clear : index -> t -> allele -> t

  (** [is_set index t allele] is [allele] in [t]. *)
  val is_set : index -> t -> allele -> bool

  val fold : index -> f:('a -> allele -> 'a) -> init:'a -> t -> 'a

  val iter : index -> f:(allele -> unit) -> t -> unit

  (** [cardinal t] returns the number of alleles found in [t]. *)
  val cardinal : t -> int

  (** For graph construction. *)
  val empty : unit -> t
  val compare : t -> t -> int
  val equals : t -> t -> bool

  (** [union e1 e2] will return an edge set with all alleles found in
      [e1] and/or [e2]. *)
  val union : t -> t -> t

  (** [inter e1 e2] will return an edge set with alleles found in both
      [e1] and [e2]. *)
  val inter : t -> t -> t

  (** [diff e1 e2] will return an edge set with alleles found in [e1] but not
      in [e2]. *)
  val diff : t -> t -> t

  (** [complement index t] returns a set of all the alleles not in [t].*)
  val complement : index -> t -> t

  (** Construct a string of all the alleles found in the edge set. *)
  val to_string : ?compress:bool -> index -> t -> string

  (** Construct a string with all of the alleles {b not} found in the edge set.

      @param prefix will determine if the complement string is prefixed
        (ie. "Complement of ") *)
  val complement_string : ?compress:bool -> ?prefix:string -> index -> t -> string

  (** [to_human_readable] uses heuristics to pick a shorter string
      representation of the edge set.

      @param compress Use allele run compression
        (ex. ["A*01:01"; "A*01:02"; .. "A*01:05"] -> "A*01:01-05"), defaults to true.

      @param max_length Trim the string, defaults to first 500 characters.
      @param complement Whether to consider printing the complement string
        (defaults to `Yes):
        [`No]       -> Never print the complement string.
        [`Yes]      -> If the complement set has fewer alleles then print it with
                        "C. of" as a prefix
        [`Prefix s] -> If the complement set has fewer alleles then print it with
                        [s] as a prefix
      *)
  val to_human_readable : ?compress:bool -> ?max_length:int ->
    ?complement:[ `No | `Prefix of string | `Yes] -> index -> t -> string

end (* Set *)

(** Maps over indexed alleles. *)
module Map : sig

  type 'a t

  (** [make index default_value]. *)
  val make : index -> 'a -> 'a t

  (** [init index f]. *)
  val init : index -> (allele -> 'a) -> 'a t

  (** [get index map allele]. *)
  val get : index -> 'a t -> allele -> 'a

  (** [cardinal m] size of the map [m]. *)
  val cardinal : 'a t -> int

  (** [update2 source dest u] update the value of [dest] with the result of
      [u] and the values of the same allele of [source] and [dest]. *)
  val update2 : source:'a t -> dest:'b t -> ('a -> 'b -> 'b) -> unit

  (** [fold index f init map] fold over all alleles found in the [map]. *)
  val fold : index -> f:('a -> 'b -> allele -> 'a) -> init:'a -> 'b t -> 'a

  (** [iter index f map] iter over all allele assignments in the [map]. *)
  val iter : index -> f:('b -> allele -> unit) -> 'b t -> unit

  (** [map index f cm] return a new map based on applying [f] to each element
       of [cm]. *)
  val map : index -> f:('a -> allele -> 'b) -> 'a t -> 'b t

  (** [fold_wa index f init map] fold over all values found in the [map],
      without regard for the allele key (wa = without allele). *)
  val fold_wa : f:('a -> 'b -> 'a) -> init:'a -> 'b t -> 'a

  (** [iter_wa index f map] iter over all allele assignments in the [map],
      without regard for the allele key (wa = without allele). *)
  val iter_wa : f:('b -> unit) -> 'b t -> unit

  (** [map_wa f m] maps the values of the map [m] without regard to the allele
      key (wa = without allele) .*)
  val map_wa : f:('a -> 'b) -> 'a t -> 'b t

  (** [map2_wa f m1 m2] map from two maps. *)
  val map2_wa : f:('a -> 'b -> 'c) -> 'a t -> 'b t -> 'c t

  (** [values_assoc index m] compress (invert) the values found in [m] into
      an association list. *)
  val values_assoc : index -> 'a t -> ('a * Set.t) list

  (** [update_from_and_fold set f init map] updates and folds over the [set]
      elements of [map].*)
  val update_from_and_fold : Set.t -> f:('a -> 'b -> 'b * 'a) -> init:'a -> 'b t -> 'a

end (* Map *)
