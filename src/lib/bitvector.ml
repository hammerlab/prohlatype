(** Imperative Bitvectors.  *)

open Util
module Iab = Ints_as_bits

type t =
  { mask              : int
  ; final_length      : int
  ; last_array_index  : int
  ; mutable a     : int array
  }
[@@deriving eq,ord]

let empty =
  { a            = [| |]
  ; mask         = 0
  ; final_length = 0
  ; last_array_index = 0
  }

let create ~size default =
  let mask, final_length, array_length =
    let r = size mod Iab.width in
    Iab.lsb_masks.(r)
    , (if r = 0 then Iab.width else r)
    , (if r = 0 then size / Iab.width else size / Iab.width + 1)
  in
  let a = if default
    then Array.make array_length Iab.all_ones
    else Array.make array_length 0
  in
  (* adjust last bits *)
  let last_array_index = array_length - 1 in
  if default && mask <> 0 then Array.unsafe_set a last_array_index mask;
  { a; mask; final_length; last_array_index }

let copy bv = { bv with a = Array.copy bv.a }

let capacity bv = Iab.width * (bv.last_array_index + 1)

let cardinal bv =
  let n = ref 0 in
  for i = 0 to bv.last_array_index do
    n := !n + Iab.count_bits bv.a.(i)         (* MSB of last element are all 0 *)
  done;
  !n

let is_empty bv =
  try
    for i = 0 to bv.last_array_index do
      if bv.a.(i) <> 0 then raise Exit     (* MSB of last element are all 0 *)
    done;
    true
  with Exit ->
    false

(* Use the safe array access and throw an exception if out of bounds.
    This can be compiled away with -unsafe. *)
let get bv i =
  let n = i / Iab.width in
  let j = i mod Iab.width in
  bv.a.(n) land (1 lsl j) <> 0

let set bv i =
  let n = i / Iab.width in
  let j = i mod Iab.width in
  bv.a.(n) <- bv.a.(n) lor (1 lsl j)

let reset bv i =
  let n = i / Iab.width in
  let j = i mod Iab.width in
  bv.a.(n) <- bv.a.(n) land (lnot (1 lsl j))

let flip bv i =
  let n = i / Iab.width in
  let j = i mod Iab.width in
  bv.a.(n) <- bv.a.(n) lxor (1 lsl j)

let clear bv =
  Array.fill bv.a 0 (bv.last_array_index + 1) 0

let iter bv f =
  for n = 0 to bv.last_array_index - 1 do
    let j = Iab.width * n in
    for i = 0 to Iab.width - 1 do
      let () = f (j+i) (bv.a.(n) land (1 lsl i) <> 0) in
      ()
    done
  done;
  let j = Iab.width * bv.last_array_index in
  for i = 0 to bv.final_length - 1 do
    let () = f (j + i) (bv.a.(bv.last_array_index) land (1 lsl i) <> 0) in
    ()
  done

let iter_true bv f =
  iter bv (fun i b -> if b then f i else ())

let to_list bv =
  let l = ref [] in
  iter_true bv (fun i -> l := i :: !l);
  !l

let to_sorted_list bv =
  List.rev (to_list bv)

(* Interpret these as indices. *)
let of_list l =
  let bv = create ~size:(List.length l) false in
  List.iter ~f:(set bv) l;
  bv

exception FoundFirst of int

let first_exn bv =
  try
    iter_true bv (fun i -> raise (FoundFirst i));
    raise Not_found
  with FoundFirst i ->
    i

let first bv =
  try Some (first_exn bv)
  with Not_found -> None

let filter bv p =
  iter_true bv
    (fun i -> if not (p i) then reset bv i)

let negate_self b =
  for n = 0 to b.last_array_index do
    Array.unsafe_set b.a n (lnot (Array.unsafe_get b.a n))
  done;
  if b.mask <> 0 then
    Array.unsafe_set b.a b.last_array_index
      (b.mask land (Array.unsafe_get b.a b.last_array_index))

let negate b =
  let a = Array.map (lnot) b.a in
  if b.mask <> 0 then begin
    Array.unsafe_set a b.last_array_index
      (b.mask land (Array.unsafe_get a b.last_array_index))
  end;
  { b with a }

let union_into ~into bv =
  for i = 0 to bv.last_array_index do
    Array.unsafe_set into.a i
      ((Array.unsafe_get into.a i) lor (Array.unsafe_get bv.a i))
  done

(* To avoid potentially 2 passes, figure out what we need to copy. *)
let union b1 b2 =
  let into = copy b1 in
  union_into ~into b2;
  into

(* Underlying size shrinks for inter. *)
let inter_into ~into bv =
  for i = 0 to bv.last_array_index do
    Array.unsafe_set into.a i
      ((Array.unsafe_get into.a i) land (Array.unsafe_get bv.a i))
  done

let inter b1 b2 =
  let into = copy b1 in
  inter_into ~into b2;
  into

let diff_into ~into bv =
  for i = 0 to bv.last_array_index do
    Array.unsafe_set into.a i
      ((Array.unsafe_get into.a i) land (lnot (Array.unsafe_get bv.a i)))
  done

let diff in_ not_in =
  let into = copy in_ in
  diff_into ~into not_in;
  into

let inter_diff b1 b2 =
  let inter_into = copy b1 in
  let diff_into  = copy b1 in
  let same = ref true in
  let no_i = ref true in   (* no intersection? *)
  for i = 0 to b1.last_array_index do
    let b1i = Array.unsafe_get inter_into.a i in
    let b2i = Array.unsafe_get b2.a i in
    let ii  = b1i land b2i in
    let di  = b1i land (lnot b2i) in
    Array.unsafe_set inter_into.a i ii;
    Array.unsafe_set diff_into.a i di;
    same := !same && ii = b1i;
    no_i := !no_i && ii = 0;
  done;
  inter_into, diff_into, !same, !no_i

let select bv arr =
  let l = ref [] in
  begin try
      iter_true bv
        (fun i ->
          if i >= bv.last_array_index
          then raise Exit
          else l := arr.(i) :: !l)
    with Exit -> ()
  end;
  !l

let selecti bv arr =
  let l = ref [] in
  begin try
      iter_true bv
        (fun i ->
          if i > bv.last_array_index
          then raise Exit
          else l := (arr.(i), i) :: !l)
    with Exit -> ()
  end;
  !l

let print out bv =
  Format.pp_print_string out "bv {";
  iter bv
    (fun _i b ->
      Format.pp_print_char out (if b then '1' else '0')
    );
  Format.pp_print_string out "}"
