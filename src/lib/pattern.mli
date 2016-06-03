(** Constructing a bijection between a k-mer and [0, 4^k). *)

val char_to_int : char -> int
(** [char_to_int c] converts a nucleotide [c] code to an integer. *)

val int_to_char : int -> char
(** [int_to_char i] converts an integer [i] to the appropriate nucleotide,
    reversing [char_to_int] *)

val pat_to_int : string -> int
(** [pat_to_int text] maps the nucleotides in [text] to a unique integer. *)

val int_to_pat : k:int -> int -> string
(** [int_to_pat ~k p] maps the integer pattern of [p] which encodes a kmer of
    length [k] back to a string of nuceotides. Reverses [int_to_pat]. *)

val pat_to_int_sub : string -> pos:int -> len:int -> int
(** [pat_to_int_sub text ~pos ~len] computes 
    [pat_to_int (String.sub ~pos ~len text] efficiently. *)

val reverse_complement : k:int -> int -> int
(** [reverse_complement k p] computes the pattern of [p]'s reverse complement,
    where [p] is a pattern for a [k]-mer. *)

val pow4 : int -> int
