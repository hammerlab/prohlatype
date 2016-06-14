(** Constructing a bijection between a k-mer and [0, 4^k). *)

val char_to_int : char -> int
(** [char_to_int c] converts a nucleotide [c] code to an integer. *)

val int_to_char : int -> char
(** [int_to_char i] converts an integer [i] to the appropriate nucleotide,
    reversing [char_to_int] *)

val encode : string -> int
(** [encode text] maps the nucleotides in [text] to a unique integer. *)

val decode : k:int -> int -> string
(** [decode ~k p] maps the integer pattern of [p] which encodes a kmer of
    length [k] back to a string of nuceotides. Reverses [decode]. *)

val encode_sub : string -> pos:int -> len:int -> int
(** [encode_sub text ~pos ~len] computes 
    [encode (String.sub ~pos ~len text] efficiently. *)

val encode_extend : string -> pos:int -> len:int -> int -> int

val reverse_complement : k:int -> int -> int
(** [reverse_complement k p] computes the pattern of [p]'s reverse complement,
    where [p] is a pattern for a [k]-mer. *)

val pow4 : int -> int
