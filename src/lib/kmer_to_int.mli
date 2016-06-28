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
(** [encode text] converts the nucleotides in [text] to a unique integer.
 
    @param pos  Which position in [text] to start encoding (defaults to 0).
    @param len  How many positions of [text] to encode (defaults to length of
                 [text]).
    @param ext  The current starting value (defaults to 0). This allows
                extending an index. For example [encode "ACGT"] is equivalent
                to [encode "GT" ~ext:(encode "AC")].

    @raise Invalid_argument for unsupported character sets of {char_to_int}. *)

val decode : k:int -> int -> string
(** [decode ~k p] converts the integer pattern of [p] which encodes a kmer of
    length [k] back to a string of nuceotides. Reverses [decode]. *)
 
val reverse_complement : k:int -> int -> int
(** [reverse_complement k p] computes the pattern of [p]'s reverse complement,
    where [p] is a pattern for a [k]-mer. *)

(** Returns powers of 4. *)
val pow4 : int -> int

val neighbors : k:int -> int -> int array
