(** Assorted biological logic centralized. *)

open Util

module Sequence = struct

  type type_ =
    | GDNA
    | CDNA
    | Protein
  [@@deriving show { with_path = false } ]

  let is_nucleotide ?s  = function
    | 'A' | 'C' | 'G' | 'T' -> true
    (* We'll parse N as a nucleotide, but why does IMGT report this? *)
    | 'N' -> Option.iter s ~f:(eprintf "IMGT reports that %s has an 'N' as a base!");
             true
    | _   -> false

  let is_amino_acid =
    function
    | 'A'| 'C'| 'D'| 'E'| 'F'| 'G'| 'H'| 'I'| 'K'| 'L'
    | 'M'| 'N'| 'P'| 'Q'| 'R'| 'S'| 'T'| 'V'| 'W'| 'Y' -> true
    | _ -> false

  let is_valid_character ?s st c =
    match st with
    | GDNA | CDNA -> is_nucleotide ?s c
    | Protein     -> is_amino_acid c

end (* Sequence *)

module Gene_region = struct

  (* Ripe for stricter typing *)
  type t =
    | UTR5
    | UTR3
    | Exon of int
    | Intron of int
    [@@deriving eq, ord, show, yojson]

  let next_gdna_sequence = function
    | UTR5     -> Some (Exon 1)
    | UTR3     -> None
    | Exon n   -> Some (Intron n)
    | Intron n -> Some (Exon (n + 1))

  let next_cdna_sequence = function
    | Exon n   -> Some (Exon (n + 1))
    | UTR5     -> None
    | UTR3     -> None
    | Intron n -> None

  let to_string = function
    | UTR5     -> "UTR5"
    | UTR3     -> "UTR3"
    | Exon n   -> sprintf "Exon %d" n
    | Intron n -> sprintf "Intron %d" n

  let to_short = function
    | UTR5     -> "U5"
    | UTR3     -> "U3"
    | Exon n   -> sprintf "X%d" n
    | Intron n -> sprintf "I%d" n

  let next st prev =
    let opt =
      match st with
      | Sequence.GDNA    -> next_gdna_sequence prev
      | Sequence.CDNA
      | Sequence.Protein -> next_cdna_sequence prev
    in
    let msg =
      sprintf "Asked for next boundary label of %s with %s!"
        (to_string prev) (Sequence.show_type_ st)
    in
    Option.value_exn opt ~msg

end (* Gene_region *)

module Base = struct

  type t =
    | A
    | C
    | G
    | T
    [@@deriving eq,show]

  let of_char = function
    | 'A' -> A
    | 'C' -> C
    | 'G' -> G
    | 'T' -> T
    |  c  -> invalid_argf "Unsupported base: %c" c

  let to_char = function
    | A       -> 'A'
    | C       -> 'C'
    | G       -> 'G'
    | T       -> 'T'

end (* Base *)


