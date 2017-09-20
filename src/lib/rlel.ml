(* Run Length Encoded Lists: Encode a sequence of elements by keep track of runs:
   adjacent elements that have the same value. *)

open Util

type 'a pair =
  { value   : 'a
  ; length  : int (* > 0 *)
  }
and 'a t =
  { hd  : 'a pair
  ; tl  : 'a pair list
  }
  [@@deriving show]

let total_length { hd; tl} =
  hd.length + List.fold_left tl ~init:0 ~f:(fun s { length; _} -> s + length)

let expand { hd = {value; length}; tl} =
  let rec expand_rlp v acc = function
    | 0 -> acc
    | n -> expand_rlp v (v :: acc) (n - 1)
  in
  let rec loop acc = function
    | []                    ->
        List.rev acc
    | { value; length} :: t ->
        loop (expand_rlp value acc length) t
  in
  loop (expand_rlp value [] length) tl

let rec until_different value =
  let rec loop length = function
    | h :: t when h = value -> loop (length + 1) t
    | lst (* when h <> value
    | [] *)                 -> { value; length}, lst
  in
  loop 1

(* Run length encode a list. *)
let encode = function
  | []          -> invalid_arg "Rlel.encode: empty list"
  | value :: t  ->
      let hd, nt = until_different value t in
      let rec loop acc = function
        | []     -> List.rev acc
        | h :: t -> let rlp, nt = until_different h t in
                    loop (rlp :: acc) nt
      in
      { hd
      ; tl = loop [] nt}

let init value =
  { hd = { value; length = 1}
  ; tl = []
  }

let append_r v = function
  | { hd = { value; length} ; tl = [] } when v = value ->
      { hd = { value ; length = length + 1}; tl = []}
  | { hd ; tl = { value; length } :: t } when v = value ->
      { hd; tl = { value; length = length + 1} :: t }
  | { hd ; tl } ->
      { hd; tl = { value = v; length = 1} :: tl }

let finish_r { hd; tl} =
  {hd; tl = List.rev tl}

let fold_map l ~f ~init =
  let n_hd_value, acc = f init l.hd.value in
  let ntl, facc =
    List.fold_left l.tl ~init:([], acc)
      ~f:(fun (l, acc) e ->
            let nvalue, nacc = f acc e.value in
            { e with value = nvalue } :: l, nacc)
  in
  { hd = { l.hd with value = n_hd_value}; tl = List.rev ntl }, facc

let align c1 c2 l1 l2 =
  if c1.length < c2.length then
    c1.length
    , l1
    , { c2 with length = c2.length - c1.length } :: l2
  else if c1.length > c2.length then
    c2.length
    , { c1 with length = c1.length - c2.length } :: l1
    , l2
  else (* c1.length = c2.length *)
    c2.length
    , l1
    , l2

let fold_map2_same_length l1 l2 ~f ~init =
  let hvalue, nacc = f init l1.hd.value l2.hd.value in
  let length, nl1, nl2 = align l1.hd l2.hd l1.tl l2.tl in
  let rec loop lst acc = function
    | [], []             -> { hd = { value = hvalue; length}; tl = List.rev lst}, acc
    | h1 :: t1, h2 :: t2 -> let nvalue, nacc = f acc h1.value h2.value in
                            let length, l1, l2 = align h1 h2 t1 t2 in
                            let npair = { value = nvalue; length} in
                            loop (npair :: lst) nacc (l1, l2)
    | _,        _        -> invalid_arg "different lengths"
  in
  loop [] nacc (nl1, nl2)

let expand_into_array ~f ~update ret rl =
  let rec fill_value i length v =
    for j = i to i + length - 1 do ret.(j) <- update ret.(i) v done
  in
  let rec loop i = function
    | []                    -> ()
    | { value; length} :: t -> fill_value i length (f value);
                                loop (i + length) t
  in
  fill_value 0 rl.hd.length (f rl.hd.value);
  loop rl.hd.length rl.tl
