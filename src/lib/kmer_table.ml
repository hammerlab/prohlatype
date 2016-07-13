
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

let update f {table;_} index =
  table.(index) <- f table.(index)

let correct_size tbl s =
  let n = String.length s in
  let k = k tbl in
  if n <> k then
    invalid_argf "String length %d doesn't match table: %d" n k
  else
    ()

let update_kmer f tbl s =
  correct_size tbl s;
  update f tbl (Kmer_to_int.encode s)

let lookup {table; _} index =
  table.(index)

let lookup_kmer tbl s =
  correct_size tbl s;
  lookup tbl (Kmer_to_int.encode s)

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

let foldi ~f ~init { table; _} =
  Array.fold_left ~init:(0, init) ~f:(fun (i, a) v -> (i + 1, f a i v)) table
  |> snd

let iteri ~f =
  foldi ~init:() ~f:(fun () p v -> f p v)

let folds ~f ~init tbl =
  let d = Kmer_to_int.decode ~k:(k tbl) in
  foldi ~f:(fun a i v -> f a (d i) v) ~init tbl

let iters ~f =
  folds ~init:() ~f:(fun () p v -> f p v)
