
module String = Sosa.Native_string

let ia f = Printf.kprintf (fun s -> raise (Invalid_argument s)) f ;;

let char_to_int = function
  | 'A' -> 0
  | 'C' -> 1
  | 'G' -> 2
  | 'T' -> 3
  | x   -> ia "char_to_int: %c" x

let int_to_char = function
  | 0 -> 'A'
  | 1 -> 'C'
  | 2 -> 'G'
  | 3 -> 'T'
  | x -> ia "int_to_char: %d" x

let int_at s index = 
  char_to_int (String.get_exn s ~index)

let pat_to_int s =
  let r = ref 0 in
  for i = 0 to String.length s - 1 do
    r := 4 * !r  + (int_at s i);
  done;
  !r

let int_to_pat ~k p =
  let rec loop s kp index =
    if kp = 0 || index < 0 then
      s
    else
      loop (String.set_exn s ~index ~v:(int_to_char (kp mod 4)))
        (kp / 4) (index - 1)
  in
  loop (String.make k 'A') p (k - 1)

let pat_to_int_sub text ~pos ~len =
  let rec loop i acc =
    if i = len then acc
    else loop (i + 1) (acc * 4 + (int_at text (pos + i)))
  in
  loop 0 0

let reverse_complement ~k i =
  let c = function
          | 0 -> 3 | 1 -> 2 | 2 -> 1 | 3 -> 0
          | x -> ia "complement:%d" x
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

