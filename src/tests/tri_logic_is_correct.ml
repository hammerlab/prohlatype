(* Test that the triangular numbering logic gives the correct 1-1
 * indexing scheme. *)

open Prohlatype

(* Test the k <-> i, j mapping logic. *)
let test n =
  let open Triangular_array in
  let kn = k n in
  let ki = kinv n in
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

let () =
  test 30000
