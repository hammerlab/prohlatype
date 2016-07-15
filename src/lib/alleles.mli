(** To keep an efficient representation of which allele's are present along an
    edge we split the data in two. The first pieces is a constant index, while
    the second is an efficent representation of an allele set or map that relies
    upon the index. The idea is that we have only one index (per graph) and many
    edge sets/maps. Alleles are represented by strings. *)

type allele = string

val compare : allele -> allele -> int
val equal : allele -> allele -> bool
val to_string : allele -> string

type index

(** [index list_of_alleles] will create an [index]. *)
val index : allele list -> index

module Set : sig

  type t

  (** [init index] create an empty set. *)
  val init : index -> t

  (** [singleton index allele] will create an edge set with just [allele]. *)
  val singleton : index -> allele -> t

  (** [set index t allele] will make sure that [allele] is
      in [t], specifically [is_set index t allele] will be [true].
      *)
  val set : index -> t -> allele -> unit

  (** [clear index t allele] will make sure that [allele] is not
      in [t], specifically [is_set index t allele] will be
      [false]. *)
  val clear : index -> t -> allele -> unit

  (** [is_set index t allele] is [allele] in [t]. *)
  val is_set : index -> t -> allele -> bool

  val fold : index -> f:('a -> allele -> 'a) -> init:'a -> t -> 'a

  val iter : index -> f:(allele -> unit) -> t -> unit

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

  (** [complement index t] returns a set of all the alleles not in [t].*)
  val complement : index -> t -> t

  val is_empty : t -> bool

  (** [any t] are any of the alleles in the set? *)
  val any : t -> bool

  (** [all t] are all of the alleles in the set? *)
  val all : index -> t -> bool

  (** Construct a string of all the alleles found in the edge set. *)
  val to_string : ?compress:bool -> index -> t -> string

  (** Construct a string with all of the alleles {b not} found in the edge set.

      @param complement_prefix will determine if the complement string is
             prefixed (ie. "Complement of ") *)
  val complement_string : ?compress:bool -> ?complement_prefix:string -> index
      -> t -> string

  (** [to_human_readable] picks the shorter (fewer number of alleles) string
      representation of the edge set.

      If the algorithm chooses to represent the set with a complement "C. of"
      will be prepended unless otherwise specified. *)
  val to_human_readable : ?compress:bool -> ?max_length:int -> ?complement_prefix:string ->
    index -> t -> string

end

(** Maps over indexed alleles. *)
module Map : sig

  type 'a t

  (** [make index default_value]. *)
  val make : index -> 'a -> 'a t

  (** [init index f]. *)
  val init : index -> (allele -> 'a) -> 'a t

  (** [get index map allele]. *)
  val get : index -> 'a t -> allele -> 'a

  (** [update_all set map f] apply [f] to all alleles in [map] where whether
      they are in [set] is passed the first arg to [f]. *)
  val update_all : Set.t -> 'a t -> (bool -> 'a -> 'a) -> unit

  (** [update_from set map f] apply [f] to all alleles in [map] that are
      in [set]. *)
  val update_from : Set.t -> 'a t -> ('a -> 'a) -> unit

  val update_spec : index -> 'a t -> allele -> ('a -> 'a) -> unit

  val update2 : 'a t -> 'b t -> ('a -> 'b -> 'b) -> unit

  (** [fold index f init amap] fold over all alleles found in the [map]. *)
  val fold : index -> f:('a -> 'b -> allele -> 'a) -> init:'a -> 'b t -> 'a

end
