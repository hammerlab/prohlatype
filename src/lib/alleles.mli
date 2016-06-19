(** To keep an efficient representation of which allele's are present along an
    edge we split the data in two. The first pieces is a constant index, while
    the second is an efficent representation of an allele set or map that relies
    upon the index. The idea is that we have only one index (per graph) and many
    edge sets/maps. Alleles are represented by strings. *)

type index

(** [index list_of_alleles] will create an [index]. *)
val index : string list -> index

module Set : sig

  type t

  (** [init index] create an empty set. *)
  val init : index -> t

  (** [singleton index allele] will create an edge set with just [allele]. *)
  val singleton : index -> string -> t

  (** [set index t allele] will make sure that [allele] is
      in [t], specifically [is_set index t allele] will be [true].
      *)
  val set : index -> t -> string -> unit

  (** [clear index t allele] will make sure that [allele] is not
      in [t], specifically [is_set index t allele] will be
      [false]. *)
  val clear : index -> t -> string -> unit

  (** [is_set index t allele] is [allele] in [t]. *)
  val is_set : index -> t -> string -> bool

  val fold : index -> f:('a -> string -> 'a) -> init:'a -> t -> 'a

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

  (** Construct a string of all the alleles found in the edge set. *)
  val to_string : index -> t -> string

  (** Construct a string with all of the alleles {b not} found in the edge set.

      @param annotate will prepend "Complement of " to the string description.
        defaults to [true].  *)
  val complement_string : ?annotate:bool -> index -> t -> string

  (** [to_human_readable] picks the shorter (fewer number of alleles) string
      representation of the edge set. *)
  val to_human_readable : index -> t -> string

end

(** Maps over indexed alleles. *)
module Map : sig

  type 'a t

  (** [make index default_value] *)
  val make : index -> 'a -> 'a t

  (** [update_from set map f] apply [f] to all alleles in [map] that are
      in [set]. *)
  val update_from : Set.t -> 'a t -> ('a -> 'a) -> unit

  (** [fold index f init amap] fold over all alleles found in the [map]. *)
  val fold : index -> f:('a -> 'b -> string -> 'a) -> init:'a -> 'b t -> 'a

end
