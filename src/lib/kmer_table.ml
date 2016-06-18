
open Nonstd
open Util

type 'a t = { k     : int
            ; table : 'a array
            }

let k {k;_} = k

let make k e =
  let n = Kmer_to_int.pow4 k in
  { k; table = Array.make n e }

let init k f =
  let n = Kmer_to_int.pow4 k in
  { k; table = Array.init n ~f:(fun i -> f (Kmer_to_int.decode ~k i)) }

let update f { k; table} s =
  assert (String.length s = k);
  let j = Kmer_to_int.encode s in
  table.(j) <- f table.(j)

let update_index f {table;_} state index =
  table.(index) <- f state table.(index)

let distr { table; _} =
  let mx = Array.fold_left ~f:max ~init:0 table in
  let c = Array.make (mx + 1) 0 in
  for i = 0 to Array.length table - 1 do
    let j = table.(i) in
    c.(j) <- c.(j) + 1
  done;
  c

let lookup {k; table} s =
  let n = String.length s in
  if n <> k then
    invalid_argf "String length %d doesn't match table: %d" n k
  else
    table.(Kmer_to_int.encode s)

let cross_boundary {k; table} =
  Array.fold_left table ~init:(0, []) ~f:(fun (i,acc) lst ->
      let p = Kmer_to_int.decode ~k i in
      let ni = i + 1 in
      let cross = List.filter lst ~f:(fun (_, s, o) -> o > String.length s - k) in
      match cross with
      | []   -> (ni, acc)
      | glst -> (ni, (p, glst) :: acc))
  |> snd

