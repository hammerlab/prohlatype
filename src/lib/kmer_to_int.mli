(** Constructing a bijection between a k-mer and [0, 4^k).

    Failure is via Invalid_argument exceptions. *)

val char_to_int : char -> int
(** [char_to_int c] converts a nucleotide [c] code to an integer.

    @raise Invalid_argument if char is not 'A', 'C', 'G' or 'T' *)

val int_to_char : int -> char
(** [int_to_char i] converts an integer [i] to the appropriate nucleotide,
    reversing [char_to_int].

    @raise Invalid_argument if argument is outside of [0,3]. *)

val encode : ?pos:int -> ?len:int -> ?ext:int -> string -> int

val encode_long : ?pos:int -> ?len:int -> ?ext:int64 -> string -> int64

(** [encode_N_tolerant text] converts the nucleotides in [text] to a list of integers
    that encode. The list grows for every 'N' encountered.

    @param pos  Which position in [text] to start encoding (defaults to 0).
    @param len  How many positions of [text] to encode (defaults to length of
                 [text]).
    @param exts  The current starting value (defaults to [0]). This allows
                extending an index. For example [encode_N_tolerant "ACGT"] is
                equivalent to [encode_N_tolerant "GT"
                ~exts:(encode_N_tolerant "AC")].

    @raise Invalid_argument for unsupported character sets of {char_to_int}. *)

val encode_N_tolerant : ?pos:int -> ?len:int -> ?exts:int list -> string -> int list

val decode : k:int -> int -> string
(** [decode ~k p] converts the integer pattern of [p] which encodes a kmer of
    length [k] back to a string of nuceotides. Reverses [decode]. *)

val reverse_complement : k:int -> int -> int
(** [reverse_complement k p] computes the pattern of [p]'s reverse complement,
    where [p] is a pattern for a [k]-mer. *)

(** [pow b n] returns [b ^ n]. *)
val pow : int -> int -> int

(** Neighbors that are [d] distance away from a [k] encoded value [e].

    @param Takes a list of indices in (0, k] that should be skipped. *)
val neighbors : ?skip:int list -> ?d:int -> k:int -> int -> int array

(** Encode all of the neighbors of a string that are [d] distance away.
    Automatically, "N-tolerant". *)
val encode_neighbors : ?len:int -> d:int -> string -> int array
