
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

let rec pow4 n =
  if n = 0 then 1
  else if n = 1 then 4
  else 4 * pow4 (n - 1)

let others = function
  | 0 -> [| 1; 2; 3|]
  | 1 -> [| 0; 2; 3|]
  | 2 -> [| 0; 1; 3|]
  | 3 -> [| 0; 1; 2|]
  | x   -> invalid_argf "others: %d" x

let coefficients ~k p =
  Array.init k (fun i ->
    let o = k - i - 1 in
    (p / (pow4 o)) mod 4)

let neighbors ~k e =
  let coefs = coefficients ~k e in
  Array.init (k * 3) (fun i ->
    let p = i / 3 in
    let c = coefs.(p) in
    let o = others c in
    let j = i mod 3 in
    let p4 = pow4 (k - p - 1) in
    e + (o.(j) - c) * p4)

