
open Util

let dseq ?(alph_size=4) ~msg ?(er=0.01) =
  let n = String.length msg in
  let er_ub = 1. /. (float (alph_size * n))  in
  if er < 0.0 then
    invalid_argf "er less than 0.0 %f" er
  else if er >= er_ub then invalid_argf "er: %f greater than upper bound" er er_ub;
  let lmp = log (er /. (float (alph_size - 1))) in
  let lcp = log (1. -. er) in
  fun s ->
    let l = String.length s in
    if l <> n then invalid_argf "wrong message length: %d" l else
      let rec loop index c m =
        if index = n then
          exp ((float c) *. lcp +. (float m) *. lmp)
        else if (String.get_exn msg ~index = String.get_exn s ~index) then
          loop (index + 1) (c + 1) m
        else
          loop (index + 1) c (m + 1)
      in
      loop 0 0 0
