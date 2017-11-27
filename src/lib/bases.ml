(* What are the possible states of the DNA sequence. *)

open Util

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
