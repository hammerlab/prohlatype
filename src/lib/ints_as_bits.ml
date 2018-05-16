(** A utility module where we utilize the bits of an integer.
 *
 * This is the just the common code that is used in Bitvector.
 *)

let width = Sys.word_size - 1

(** Index from the least significant bit. We create masks to zero out the most significant
    bits that aren't used to store values. This is necessary when we are
    constructing or negating a bit vector. *)
let lsb_masks =
  let a = Array.make (width + 1) 0 in
  for i = 1 to width do
    a.(i) <- a.(i-1) lor (1 lsl (i - 1))
  done;
  a

let all_ones = lsb_masks.(width)

(* count the 1 bits in [n]. See https://en.wikipedia.org/wiki/Hamming_weight *)
let count_bits n =
  let rec recurse count n =
    if n = 0 then count else recurse (count+1) (n land (n-1))
  in
  recurse 0 n

let to_binary_string i =
  Bytes.init width (fun j ->
    let m = 1 lsl (width - j - 1) in
    if i land m = m then '1' else '0')
