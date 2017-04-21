(** To keep an efficient representation of which allele's are present along an
    edge we split the data in two. The first pieces is a constant index, while
    the second is an efficent representation of an allele set or map that relies
    upon the index. The idea is that we have only one index (per graph) and many
    edge sets/maps. Alleles are represented by strings. *)

type allele = string

val compare : allele -> allele -> int
val equal : allele -> allele -> bool

type index

val length : index -> int

(** [index list_of_alleles] will create an [index]. *)
val index : allele list -> index

val to_alleles : index -> allele list

(** Call this before any useful, non-graph construction functions in
    Set or Map function. *)
val setup : index -> unit

module Set : sig

  type set

  (** For graph construction, you can call these without the [setup] above. *)
  val empty : set
  val compare : set -> set -> int
  val equals : set -> set -> bool

  val copy : set -> set

  val init : unit -> set

  (** [set set allele] will make sure that [allele] is
      in [set], specifically [is_set set allele] will be [true].
      [set] is mutated. *)
  val set : set -> allele -> set

  (** [singleton allele] will create an edge set with just [allele]. *)
  val singleton : allele -> set

  (** [clear set allele] will make sure that [allele] is not
      in [set], specifically [is_set set allele] will be
      [false]. [set] is mutated. *)
  val clear : set -> allele -> set

  (** [is_set set allele] is [allele] in [set]. *)
  val is_set : set -> allele -> bool

  val fold : f:('a -> allele -> 'a) -> init:'a -> set -> 'a

  val iter : f:(allele -> unit) -> set -> unit

  (** [cardinal set] returns the number of alleles found in [set]. *)
  val cardinal : set -> int


  (** [union e1 e2] will return an edge set with all alleles found in
      [e1] and/or [e2]. *)
  val union : set -> set -> set

  val unite : into:set -> set -> unit

  (** [inter e1 e2] will return an edge set with alleles found in both
      [e1] and [e2]. *)
  val inter : set -> set -> set

  (** [diff e1 e2] will return an edge set with alleles found in [e1] but not
      in [e2]. *)
  val diff : set -> set -> set

  (** [complement set] returns a set of all the alleles not in [set].*)
  val complement : set -> set

  (** Construct a string of all the alleles found in the edge set. *)
  val to_string : ?compress:bool -> set -> string

  (** Construct a string with all of the alleles {b not} found in the edge set.

      @param prefix will determine if the complement string is prefixed
        (ie. "Complement of ") *)
  val complement_string : ?compress:bool -> ?prefix:string -> set -> string

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
    ?complement:[ `No | `Prefix of string | `Yes] -> set -> string

end (* Set *)

(** Maps over indexed alleles. *)
module Map : sig

  type 'a map

  val to_array : 'a map -> 'a array

  (** [make default_value]. *)
  val make : 'a -> 'a map

  (** [init f]. *)
  val init : (allele -> 'a) -> 'a map

  (** [get map allele]. *)
  val get : 'a map -> allele -> 'a

  (** [update2 source dest u] update the value of [dest] with the result of
      [u] and the values of the same allele of [source] and [dest]. *)
  val update2 : source:'a map -> dest:'b map -> ('a -> 'b -> 'b) -> unit

  (** [fold f init map] fold over all alleles found in the [map]. *)
  val fold : f:('a -> 'b -> allele -> 'a) -> init:'a -> 'b map -> 'a

  (** [iter f map] iter over all allele assignments in the [map]. *)
  val iter : f:('b -> allele -> unit) -> 'b map -> unit

  (** [map f cm] return a new map based on applying [f] to each element
       of [cm]. *)
  val map : f:('a -> allele -> 'b) -> 'a map -> 'b map

  (** [fold_wa f init map] fold over all values found in the [map],
      without regard for the allele key (wa = without allele). *)
  val fold_wa : f:('a -> 'b -> 'a) -> init:'a -> 'b map -> 'a

  (** [iter_wa f map] iter over all allele assignments in the [map],
      without regard for the allele key (wa = without allele). *)
  val iter_wa : f:('b -> unit) -> 'b map -> unit

  (** [map_wa f m] maps the values of the map [m] without regard to the allele
      key (wa = without allele) .*)
  val map_wa : f:('a -> 'b) -> 'a map -> 'b map

  (** [map2_wa f m1 m2] map from two maps. *)
  val map2_wa : f:('a -> 'b -> 'c) -> 'a map -> 'b map -> 'c map

  (** [values_assoc m] compress (invert) the values found in [m] into
      an association list. *)
  val values_assoc : 'a map -> ('a * Set.set) list

  val update_from : Set.set -> f:('a -> 'a) -> 'a map -> unit

  (** [update_from_and_fold set f init map] updates and folds over the [set]
      elements of [map].*)
  val update_from_and_fold : Set.set -> f:('a -> 'b -> 'b * 'a) -> init:'a -> 'b map -> 'a

  (** [choose m] return a single element from the map. *)
  val choose : 'a map -> 'a

end (* Map *)

module Selection :
sig
  type t =
    | Regex of string
    | Specific of string
    | Without of string
    | Number of int
  val equal : t -> t -> Ppx_deriving_runtime.bool
  val compare : t -> t -> Ppx_deriving_runtime.int
  val show : t -> string
  val pp : Format.formatter -> t -> unit
  val to_string : t -> string
  val list_to_string : t list -> string
  val apply_to_assoc : t list -> (string * 'a) list -> (string * 'a) list
end (* Selection *)

module Input :
sig
  type t =
    | AlignmentFile of (string * bool)
    | MergeFromPrefix of (string * Distances.logic * bool)
  val imputed : t -> bool
  val to_string : t -> string
  val construct :
    t ->
    (Mas_parser.result * (Util.StringMap.key * Util.StringMap.key) list,
      string)
    result
end (* Input *)
