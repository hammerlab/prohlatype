(* Test that the triangular numbering logic gives the correct 1-1
 * indexing scheme. *)

open Prohlatype

(* Test the k <-> i, j mapping logic. *)
let test_upper n =
  let open Triangular_array in
  let kn = upper_triangular_index n in
  let ki = upper_triangular_inverse n in
  let kprev = ref (-1) in
  for i = 0 to n - 2 do
    for j = i + 1 to n - 1 do
      let k = kn i j in
      let ni, nj = ki k in
      assert (!kprev + 1 = k);
      assert (i = ni && j = nj);
      kprev := k
      (*printf "i: %d\tj: %d\t:k: %d\ ---- ni: %d\tnk: %d %b\n"
        i j k ni nj (i = ni && j = nj) *)
    done
  done

let test_full_upper n =
  let open Triangular_array in
  let kn = full_upper_triangular_index n in
  let ki = full_upper_triangular_inverse n in
  let kprev = ref (-1) in
  for i = 0 to n - 1 do
    for j = i to n - 1 do
      let k = kn i j in
      let ni, nj = ki k in
      assert (!kprev + 1 = k);
      assert (i = ni && j = nj);
      kprev := k
      (*printf "i: %d\tj: %d\t:k: %d\ ---- ni: %d\tnk: %d %b\n"
        i j k ni nj (i = ni && j = nj) *)
    done
  done


let () =
  let module Tn = Triangular_array.Triangular in
  assert Tn.(number 0 |> inverse = 0);  (* NOT true for number_m1 *)
  for i = 1 to 1000 do
    printf "%d\t%!" i;
    if  Tn.(inverse (number i) <> i) then
      failwithf "Not %d the same %d %d\n"
        i (Tn.number i) Tn.(inverse (number i));

    if  Tn.(inverse_m1 (number_m1 i) <> i) then
      failwithf "Not %d the same %d %d\n"
        i (Tn.number_m1 i) Tn.(inverse_m1 (number_m1 i))
  done;
  printf "just upper\n%!";
  test_upper 30000;
  printf "full upper\n%!";
  test_full_upper 30000
