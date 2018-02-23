(* When we only want to store (N^2) / 2 elements. *)
open Util

let number n =
  n * (n + 1) / 2

let inverse t =
  let tf = float t in
  (truncate (sqrt (1.0 +. 8.0 *. tf)) - 1) / 2

(* Not exactly the triangular number, but T_(n-1):
    just have less subtractions this way. *)
let number_m1 n =
  n * (n - 1) / 2

let inverse_m1 t =
  let tf = float t in
  (truncate (sqrt (1.0 +. 8.0 *. tf)) + 1) / 2

module Indices = struct

  let upper n i j =
    number_m1 n - number_m1 (n - i) + j - i - 1

  let upper_inverse n k =
    let tn = number_m1 n in
    let i  = n - 1 - (inverse_m1 (tn - k - 1)) in
    let j  = k + i + 1 - number_m1 n + number_m1 (n - i) in
    i, j

  let full_upper n i j =
    if i >= n then
      invalid_argf "First index %d greater than or equal to max value %d."
        i n
    else if j >= n then
      invalid_argf "Second index %d greater than or equal to max value %d."
        j n
    else
      number n - number (n - i) + j - i

  let full_upper_inverse n k =
    let tn = number n in
    let i  = n - 1 - (inverse (tn - k - 1)) in
    let j  = k + i - number n + number (n - i) in
    i, j

end (* Indices *)

module Array = struct

  (* Public Api *)
  type 'a t =
    { full : bool
    ; n : int
    ; a : 'a array
    }

  let size t =
    Array.length t.a

  let make full max_j z =
    let n = max_j in
    let s = if full then number max_j else number_m1 max_j in
    { full; n ; a = Array.make s z }

  let get t i j =
    if t.full then begin
      assert (i <= j);
      Array.get t.a (Indices.full_upper t.n i j)
    end else begin
      assert (i < j);
      Array.get t.a (Indices.upper t.n i j)
    end

  let update t i j ~f =
    if t.full then begin
      assert (i <= j);
      let ki = Indices.full_upper t.n i j in
      Array.set t.a ki (f (Array.get t.a ki))
    end else begin
      assert (i < j);
      let ki = Indices.upper t.n i j in
      Array.set t.a ki (f (Array.get t.a ki))
    end

  let foldi_left t ~init ~f =
    if t.full then begin
      Array.fold_left t.a ~init:(0, init)
        ~f:(fun (k, acc) v ->
            let i, j = Indices.full_upper_inverse t.n k in
            (k + 1, f acc i j v))
      |> snd
    end else begin
      Array.fold_left t.a ~init:(0, init)
        ~f:(fun (k, acc) v ->
            let i, j = Indices.upper_inverse t.n k in
            (k + 1, f acc i j v))
      |> snd
    end

  let fold_left t ~init ~f =
    Array.fold_left t.a ~init ~f

end (* Array *)
