
module Array = ArrayLabels

type t =
  { mutable size : int
  ; mutable data : int64 array
  }

let capacity t =
  64 * (Array.length t.data)

let empty () =
  { size = 0
  ; data = [||]
  }

let array_size n =
  n / 64 + (if n mod 64 = 0 then 0 else 1)

let create_ sfun i n =
  if n < 0 then invalid_arg ("Bitset."^sfun^": negative size");
  { size = n
  ; data = Array.make (array_size n) i
  }

let create =
  create_ "create" 0L

(* Output *)
let inner_buffer =
  Buffer.create 64

let write_int64_as_bitstring_to_buffer ?(length=64) buf c =
  let rc = ref c in
  for i = 1 to length do
    Buffer.add_char buf
      (if Int64.logand !rc 1L = 1L then '1' else '0');
    rc := Int64.shift_right_logical !rc 1
  done

let int64_to_bitstring i =
  Buffer.clear inner_buffer;
  write_int64_as_bitstring_to_buffer inner_buffer i;
  Buffer.contents inner_buffer

let to_string {size; data} =
  let b = Buffer.create size in
  for i = 0 to (size - 1) / 64 - 1 do
    write_int64_as_bitstring_to_buffer b data.(i)
  done;
  let remaining_length = size mod 64 in
  if remaining_length <> 0 then
    write_int64_as_bitstring_to_buffer
      ~length:remaining_length b data.(Array.length data - 1);
  Buffer.contents b

let copy t =
  { size = t.size
  ; data = Array.copy t.data
  }

let extend t n =
  if n < t.size then
    ()
  else if n <= capacity t then
    t.size <- n
  else begin
    let ndata = Array.make (array_size n) 0L in
    Array.blit ~src:t.data ~src_pos:0 ~dst:ndata ~dst_pos:0
      ~len:(Array.length t.data);
    t.size <- n;
    t.data <- ndata
  end

let masks =
  Array.init 64 (Int64.shift_left 1L)

type bit_op =
  | Set
  | Unset
  | Toggle

let rec valid_position_value_and_mask extendable sfun t x =
  let pos = x / 64 in
  if pos < 0 then
    invalid_arg ("Bitset."^sfun^": negative index")
  else if x < t.size then
    let cur = Array.unsafe_get t.data pos in
    let delta = x mod 64 in
    let mask = masks.(delta) in
    pos, cur, mask
  else if extendable then begin
    extend t (x+1);
    valid_position_value_and_mask false sfun t x
  end else
    invalid_arg (Printf.sprintf "Bitset.%s: index past length: %d, size %d"
                  sfun x t.size )

let apply_bit_op sfun t x =
  let pos, cur, mask = valid_position_value_and_mask true sfun t x in
  let bit = Int64.logand cur mask <> 0L in
  let set = Array.unsafe_set t.data pos in
  function
  | Set    -> if not bit then set (Int64.logor cur mask)
  | Unset  -> if bit then set (Int64.logxor cur mask)
  | Toggle -> set (Int64.logxor cur mask)

(* Mutable *)
let set t x = apply_bit_op "set" t x Set

let unset t x = apply_bit_op "unset" t x Unset

let toggle t x = apply_bit_op "toggle" t x Toggle

let mem t x =
  let _, cur, mask = valid_position_value_and_mask false "mem" t x in
  Int64.logand cur mask <> 0L

let add x t =
  let dup = copy t in
  set dup x; dup

let remove x t =
  let dup = copy t in
  unset dup x; dup

let put t =
  function
  | true -> set t
  | false -> unset t

let create_full n =
  let t = create_ "create_full" (-1L) n in
  (* Fix the tail *)
  for i = n to (capacity t) - 1 do
    unset t i
  done;
  t

let compare t1 t2 =
  let len1 = Array.length t1.data in
  let len2 = Array.length t2.data in
  let mlen = min len1 len2 in
  let rec loop i =
    if i = mlen then begin
      if len1 < len2 then
        loop_second i
      else if len1 > len2 then
        loop_first i
      else
        0
    end else
      let d = Int64.compare (Array.unsafe_get t1.data i) (Array.unsafe_get t2.data i) in
      if d <> 0 then
        d
      else
        loop (i + 1)
  and loop_first i =
    if i >= len1 then 0 else
      let d = Int64.compare (Array.unsafe_get t1.data i) 0L in
      if d <> 0 then d else loop_first (i + 1)
  and loop_second i =
    if i >= len2 then 0 else
      let d = Int64.compare 0L (Array.unsafe_get t2.data i) in
      if d <> 0 then d else loop_second (i + 1)
  in
  loop 0

let equal t1 t2 =
  compare t1 t2 = 0

let num_bits_int64 =
  let rec loop c n =
    if n = 0L then c else
      loop (c + 1) Int64.(logand n (sub n 1L))
  in
  loop 0

let count t =
  Array.fold_left t.data ~init:0 ~f:(fun s c -> s + num_bits_int64 c)

let mtable =
  [| 0;  1;  2; 53;  3;  7; 54; 27
  ;  4; 38; 41;  8; 34; 55; 48; 28
  ; 62;  5; 39; 46; 44; 42; 22;  9
  ; 24; 35; 59; 56; 49; 18; 29; 11
  ; 63; 52;  6; 26; 37; 40; 33; 47
  ; 61; 45; 43; 21; 23; 58; 17; 10
  ; 51; 25; 36; 32; 60; 20; 57; 16
  ; 50; 31; 19; 15; 30; 14; 13; 12
  |]

(* 64 Bit Debruin Sequence.*)
let magic =
  0x022fdd63cc95386dL

let clear_masks =
  Array.map masks ~f:(fun m -> Int64.mul m (-1L))

let next_set_bit64 ?starting_at i =
  let i =
    match starting_at with
    | None               -> i
    | Some p when p < 0  -> invalid_arg "Bitset.next_set_bit64 starting_at below 0"
    | Some p when p > 63 -> 0L
    | Some p             -> Int64.logand i clear_masks.(p)
  in
  if i = 0L then
    None
  else
    (*let j = Int64.(shift_right_logical (mul (logand i (mul i (-1L))) magic) 58) in *)
    let j = Int64.(shift_right_logical (mul (logand i (add (lognot i) 1L)) magic) 58) in
    Some (Array.unsafe_get mtable (Int64.to_int j))

(* Find the first set bit in the bit array *)
let rec next_set_bit t x =
  if x < 0 then
    invalid_arg "Bitset.next_set_bit"
  else
    let pos = x / 64 in
    if pos < Array.length t.data then begin
      let delta = x mod 64 in
      let cur = Array.unsafe_get t.data pos in
      match next_set_bit64 ~starting_at:delta cur with
      | None   -> next_set_bit t (64 * (pos + 1))
      | Some p -> Some (pos * 64 + p)
    end else
      None

let fold t ~f ~init =
  let rec loop acc = function
    | None   -> acc
    | Some p -> loop (f acc p) (next_set_bit t (p + 1))
  in
  loop init (next_set_bit t 0)

let iter t ~f =
  fold t ~init:() ~f:(fun () p -> f p)

type set_op =
  | Inter
  | Diff
  | Unite
  | DiffSym

let apply_set_op op t1 t2 =
  let len1 = Array.length t1.data in
  let len2 = Array.length t2.data in
  let mlen = min len1 len2 in
  let rec loop i =
    if i = mlen then begin
      if len1 < len2 && (op = Unite || op = DiffSym) then
        extend_first_with_second ()
      else
        ()
    end else
      let v1 = Array.unsafe_get t1.data i in
      let v2 = Array.unsafe_get t2.data i in
      let vc =
        match op with
        | Inter   -> Int64.logand v1 v2
        | Diff    -> Int64.(logand v1 (lognot v2))
        | Unite   -> Int64.logor v1 v2
        | DiffSym -> Int64.logxor v1 v2
      in
      Array.unsafe_set t1.data i vc;
      loop (i + 1)
  and extend_first_with_second () =
    extend t1 (len2 * 64);
    Array.blit ~src:t2.data ~src_pos:len1 ~dst:t1.data ~dst_pos:len1 ~len:(len2 - len1)
  in
  loop 0

let intersect t1 t2 = apply_set_op Inter t1 t2

let differentiate t1 t2 = apply_set_op Diff t1 t2

let unite t1 t2 = apply_set_op Unite t1 t2

let differentiate_sym t1 t2 = apply_set_op DiffSym t1 t2

let biop_with_copy f a b =
  let a' = copy a in
  f a' b;
  a'

let inter a b =
  biop_with_copy intersect a b

let union a b =
  biop_with_copy unite a b

let diff a b =
  biop_with_copy differentiate a b

let sym_diff a b =
  biop_with_copy differentiate_sym a b
