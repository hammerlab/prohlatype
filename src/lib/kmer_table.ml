
open Nonstd

let invalid_argf fmt = ksprintf invalid_arg fmt

let k (n,_) = n

let make k e =
  let n = Pattern.pow4 k in
  k, Array.make n e

let init k f =
  let n = Pattern.pow4 k in
  k, Array.init n ~f:(fun i -> f (Pattern.decode ~k i))

let update f ((k, t) as tb) s =
  assert (String.length s = k);
  let j = Pattern.encode s in
  t.(j) <- f t.(j);
  tb

let update_index f ((_, t) as tb) state index =
  t.(index) <- f state t.(index);
  tb

let distr (_, t) =
  let mx = Array.fold_left ~f:max ~init:0 t in
  let c = Array.make (mx + 1) 0 in
  for i = 0 to Array.length t - 1 do
    let j = t.(i) in
    c.(j) <- c.(j) + 1
  done;
  c

let lookup (k, t) s =
  let n = String.length s in
  if n <> k then
    invalid_argf "String length %d doesn't match table: %d" n k
  else
    t.(Pattern.encode s)

