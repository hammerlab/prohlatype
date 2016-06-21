
open Nonstd
open Util

type 'a t = { k     : int
            ; table : 'a array
            }

type index = int

let k {k;_} = k

let make k e =
  let n = Kmer_to_int.pow4 k in
  { k; table = Array.make n e }

let init k f =
  let n = Kmer_to_int.pow4 k in
  { k; table = Array.init n ~f:(fun i -> f (Kmer_to_int.decode ~k i)) }

let update_kmer f s { k; table} =
  assert (String.length s = k);
  let j = Kmer_to_int.encode s in
  table.(j) <- f table.(j)

let update f {table;_} index =
  table.(index) <- f table.(index)

let lookup {k; table} s =
  let n = String.length s in
  if n <> k then
    invalid_argf "String length %d doesn't match table: %d" n k
  else
    table.(Kmer_to_int.encode s)

let distr { table; _} =
  let mx = Array.fold_left ~f:max ~init:0 table in
  let c = Array.make (mx + 1) 0 in
  for i = 0 to Array.length table - 1 do
    let j = table.(i) in
    c.(j) <- c.(j) + 1
  done;
  c

let fold ~f ~init { table; _} =
  Array.fold_left ~f ~init table

let iter ~f { table; _} =
  Array.iter ~f table
