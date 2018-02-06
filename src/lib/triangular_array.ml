(* When we only want to store (N^2) / 2 elements. *)
open Util

module Triangular = struct

  (* Not exactly the triangular number, but T_(n-1):
      just have less subtractions this way. *)
  let number n =
    n * (n - 1) / 2

  let inverse t =
    let tf = float t in
    (truncate (sqrt (1.0 +. 8.0 *. tf)) - 1) / 2

end (* Triangular *)

let k n i j =
  Triangular.(number n - number (n - i) + j - i - 1)

let kinv n k =
  let open Triangular in
  let tn = number n in
  let i  = n - 1 - (inverse (tn - k - 1)) - 1 in
  let j  = k + i + 1 - number n + number (n - i) in
  i, j

(* Public Api *)
type 'a t =
  { n : int
  ; a : 'a array
  }

let size t =
  Array.length t.a

let make max_j z =
  let n = max_j in
  let s = Triangular.number max_j in
  { n ; a = Array.make s z }

let get t i j =
  assert (i < j);
  Array.get t.a (k t.n i j)

(*let set t i j v =
  assert (i < j);
  Array.set t.a (k t.n i j) v*)

let update t i j ~f =
  assert (i < j);
  let ki = k t.n i j in
  Array.set t.a ki (f (Array.get t.a ki))

let foldi_left t ~init ~f =
  Array.fold_left t.a ~init:(0, init)
    ~f:(fun (k, acc) v ->
        let i, j = kinv t.n k in
        (k + 1, f acc i j v))
  |> snd

let fold_left t ~init ~f =
  Array.fold_left t.a ~init ~f

