
module String = Sosa.Native_string
open Util

let char_to_int = function
  | 'A' -> 0
  | 'C' -> 1
  | 'G' -> 2
  | 'T' -> 3
  | x   -> invalid_argf "char_to_int: %c" x

let int_to_char = function
  | 0 -> 'A'
  | 1 -> 'C'
  | 2 -> 'G'
  | 3 -> 'T'
  | x -> invalid_argf "int_to_char: %d" x

let int_at s index =
  char_to_int (String.get_exn s ~index)

let encode ?(pos=0) ?len ?(ext=0) s =
  let r = ref ext in
  let len = Option.value len ~default:(String.length s) in
  for i = pos to pos + len - 1 do
    r := 4 * !r  + (int_at s i);
  done;
  !r

let encode_long ?(pos=0) ?len ?(ext=0L) s =
  let r = ref ext in
  let len = Option.value len ~default:(String.length s) in
  for i = pos to pos + len - 1 do
    r := Int64.mul 4L (Int64.add !r (Int64.of_int (int_at s i)));
  done;
  !r

let int_at_or_N s index =
  let c = String.get_exn s ~index in
  if c = 'N' then None else Some (char_to_int c)

let encode_N_tolerant_with_pos ?(pos=0) ?len ?(exts=[0]) s =
  let len = Option.value len ~default:(String.length s) in
  let rec loop r n_positions i =
    if i = pos + len then
      r, n_positions
    else
      match int_at_or_N s i with
      | None    ->  (* there is an 'N' !*)
          let nr = List.concat_map ~f:(fun e -> List.map [0;1;2;3] ~f:((+) (4 * e))) r in
          let nn_positions = i :: n_positions in
          loop nr nn_positions (i + 1)
      | Some c  ->
          let nr = List.map ~f:(fun e -> 4 * e + c) r in
          loop nr n_positions (i + 1)
  in
  loop exts [] pos

let encode_N_tolerant ?pos ?len ?exts s =
  fst (encode_N_tolerant_with_pos ?pos ?len ?exts s)

let decode ~k p =
  let rec loop s kp index =
    if kp = 0 || index < 0 then
      s
    else
      loop (String.set_exn s ~index ~v:(int_to_char (kp mod 4)))
        (kp / 4) (index - 1)
  in
  loop (String.make k 'A') p (k - 1)

let reverse_complement ~k i =
  let c = function
          | 0 -> 3 | 1 -> 2 | 2 -> 1 | 3 -> 0
          | x -> invalid_argf "complement:%d" x
  in
  let rec loop p index rp =
    if index <= 0 then
      rp
    else
      loop (p / 4) (index - 1) (c (p mod 4) + (rp * 4))
  in loop i k 0

let rec pow i n =
  if n = 0 then 1
  else if n = 1 then i
  else i * pow i (n - 1)

let others = function
  | 0 -> [| 1; 2; 3|]
  | 1 -> [| 0; 2; 3|]
  | 2 -> [| 0; 1; 3|]
  | 3 -> [| 0; 1; 2|]
  | x   -> invalid_argf "others: %d" x

(* The coefficients, in base 4, used to encode the index. *)
let coefficients ~k p =
  Array.init k (fun i ->
    let o = k - i - 1 in
    (p / (pow 4 o)) mod 4)

(* Coefficients but in reverse order. *)
let rev_coefficients ~k p =
  Array.init k (fun i ->
    (p / (pow 4 i)) mod 4)

let prod f t =
  let rec p a i = if i >= t then a * t else p (i * a) (i + 1) in
  p 1 f

(*
It is unfortunate but there is a different meaning to 'k' in this function:
Only here do we mean the number of ways to choose from n objects. But in the
rest of the code it refers to the size of the k in Kmer encoding (length of
string). And we will frequently choose from this k.

TODO: Would it be faster to use a cached version?
      And/or use the recursive additive form instead? *)
let n_c_k n k =
  if k = 0 then 0 else
    (prod (n - k + 1) n) / (prod 1 k)

let num_neighbors ~k ~d =
  let num_position_combs = n_c_k k d in
  let distance_per_position = pow 3 d in
  distance_per_position * num_position_combs

let combination_maccaffery n d =
  let rec maximize a b x =
    if (n_c_k a b ) <= x then a else maximize (a-1) b x
  in
  let rec iterate i x =
    if i = 0 then
      []
    else
      let m = maximize n i x in
      m :: iterate (i - 1) (x - (n_c_k m i))
  in
  Array.init (n_c_k n d) ~f:(iterate d)

let combination_maccaffery_arr n d =
  let rec maximize a b x =
    if (n_c_k a b ) <= x then a else maximize (a-1) b x
  in
  let iterate x =
    let a = Array.make d 0 in
    let rec loop i x =
      if i = 0 then a else begin
        let m = maximize n i x in
        a.(i - 1) <- m;
        loop (i - 1) (x - (n_c_k m i))
      end
    in
    loop d x
  in
  Array.init (n_c_k n d) ~f:(fun i -> iterate i)

let one_neighbors ~k e =
  let coefs = coefficients ~k e in
  Array.init (k * 3) (fun i ->
    let p = i / 3 in
    let c = coefs.(p) in
    let o = others c in
    let j = i mod 3 in
    let p4 = pow 4 (k - p - 1) in
    e + (o.(j) - c) * p4)

let threes k i =
  let r = Array.make k 0 in
  let rec loop kp i =
    if kp = 0 || i < 0 then r else begin
      r.(i) <- kp mod 3;
      loop (kp / 3) (i - 1)
    end
  in
  loop i (k - 1)

(* Neighbors that are [d] distance away from a [k] encoded value [e]. *)
let neighbors ?skip ?(d=1) ~k e =
  (* Keep them in reverse order so that as the index increases so does the
     power of the base. *)
  let coefs = rev_coefficients ~k e in
  (* comb : array of powers that is being mutated, d of them.
     diff : array of indices into other values to apply, also d of them. *)
  let apply comb diff =
    let rec loop a i =
      if i = d then a else
        let ci = comb.(i) in
        let cur_coeff = coefs.(ci) in
        let other_coeff = (others cur_coeff).(diff.(i)) in
        let p4 = pow 4 ci in
        let b = a + (other_coeff - cur_coeff) * p4 in
        (*Printf.printf "i: %d other_coeff: %d, cur_coeff: %d, p4: %d a: %d -> b: %d\n"
          i other_coeff cur_coeff p4 a b;*)
        loop b (i + 1)
    in
    loop e 0
  in
  if d < 0 then
    invalid_argf "Specified distance %d is less than zero" d
  else
    let eff_k, combs =
      match skip with
      | None   -> k, combination_maccaffery_arr k d
      | Some l ->
          let l = List.dedup ~compare l in
          match List.find l ~f:(fun i -> i >= k) with
          | Some bl -> invalid_argf "Asked to skip index %d that is larger or equal to k: %d" bl k
          | None    ->
              let len = List.length l in
              let eff_k = k - len in
              let com = combination_maccaffery_arr eff_k d in
              let s = List.init k ~f:(fun i -> i)
                      |> List.filter ~f:(fun i -> not (List.mem ~set:l i))
                      |> Array.of_list
              in
              eff_k, Array.map com ~f:(Array.map ~f:(fun i -> s.(i)))
    in
    let num_position_combs = n_c_k eff_k d in
    let distance_per_position = pow 3 d in
    let size = distance_per_position * num_position_combs in
    if size = 0 then [||] else
      let diffi = Array.init distance_per_position ~f:(threes d) in
      assert (Array.length combs = num_position_combs);
      Array.init size ~f:(fun i ->
        let k = i / distance_per_position in
        let j = i mod distance_per_position in
        apply combs.(k) diffi.(j))

let encode_neighbors ?len ~d s =
  let k = Option.value len ~default:(String.length s) in
  let ens, skip = encode_N_tolerant_with_pos ?len s in
  if d = 0 then
    Array.of_list ens
  else begin
    List.map ens ~f:(neighbors ~skip ~d ~k)
    |> Array.concat
  end
