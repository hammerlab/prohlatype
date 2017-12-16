
(* Start will occupy the higher order bits. This way a comparison can test
  * for equality (the most common case) and determine the Before/After path
  * without looking at end.
  *
  * This means that the new range must start at 1 and 0 is no longer a valid
  * element in our partition.
  *)
type t = int64 [@@deriving eq,show,show]

(* Masks *)
let lower32 = 4294967295L
let upper32 = -4294967296L

type start = S of int64 [@@unboxed]
type end_ = E of int64 [@@unboxed]

let start_of_int64 i =
  S (Int64.shift_left i 32)

let int64_of_start (S i) =
  Int64.(shift_right_logical (logand i upper32) 32)

let start_p (t : t) : start =
  S (Int64.logand t upper32)

let end_of_int64 i =
  E i

let int64_of_end (E i) =
  i

let end_p t =
  E (Int64.logand t lower32)

let make64 (S s) (E e) : t =
  Int64.logor s e

let start t =
  Int64.to_int (int64_of_start (start_p t))

let end_ t =
  Int64.to_int (int64_of_end (end_p t))

let width t =
  (end_ t) - (start t) + 1

let inside i t =
  (start t) <= i && i <= (end_ t)

let to_string_interval t =
  sprintf "(%d,%d)" (start t) (end_ t)

let to_string t =
  sprintf "(%d,%d)" (start t) (end_ t)
(*
let to_string (t : t) =
  sprintf "%Ld" t
   *)

type inter_diff_res =
    | Before
    | After

    (* Different possible intersections. The way to read these is that the
       first word describes the relationship of result to the intersection
       for the first interval argument to {inter_diff} and the second word
       for the second argument.

       For example, if inter_diff i1 i2 = EmptySplit {inter; before, after}
       means that there is no resulting difference for the intersection of i1
       and i2 as it pertains to i1, it is empty, but that the intersection
       splits i2 into {before} and {after}.

       If i1 = (3,3) and i2 = (2,6) then
       inter_diff i1 i2 =
         EmptySplit {inter = (3,3); before = (2,2); after = (4,6)}.

       As a consequence the record labels are changing and may pertain to
       different elements of {i1} or {i2}. *)

    | BeforeAfter of { inter: t; before: t; after : t }
    | BeforeEmpty of { inter: t; before: t }

    | EmptyAfter of  { inter: t; after: t }
    | EmptyBefore of { inter: t; before: t }
    | EmptyEmpty of  { inter: t }
    | EmptySplit of  { inter: t; before: t; after: t}

    | AfterEmpty of  { inter: t; after: t }
    | AfterBefore of { inter: t; after: t;  before: t}

    | SplitEmpty of  { inter: t; before: t; after : t}
    [@@deriving eq, ord, show]

  let start_of_end e = 
    start_of_int64 (int64_of_end e)

  let end_of_start s =
    end_of_int64 (int64_of_start s)

  let end_pred (E e) =
    E (Int64.pred e)

  let start_1 = 4294967296L

  let start_pred (S s) =
    S (Int64.sub s start_1)

  let start_succ (S s) =
    S (Int64.add s start_1)

  let end_succ (E s) =
    E (Int64.succ s)

  let inter_diff t1 t2 =
    if t1 = t2 then
      EmptyEmpty { inter = t1 }
    else if t1 < t2 then (* s1 <= s2 *)
      begin
        let s1 = start_p t1 in
        let s2 = start_p t2 in
        let e1 = end_p t1 in
        let e2 = end_p t2 in
        if s1 = s2 then (* -> e1 < e2 ----> s1 = s2 <= e1 < e2 *)
              let after = make64 (start_of_end (end_succ e1)) e2 in
              EmptyAfter { inter = t1; after }
        else (* s1 < s2 --------> e1 ? e2  we don't know!*)
          begin  
            if e1 = e2 then (* -------> s1 < e1 = s2 = e2 *) 
                  let before = make64 s1 (end_pred (end_of_start s2)) in
                  BeforeEmpty {before; inter = t2 }
            else if e1 < e2 (* ----- s1 <= e1 ? s2 <= e2 *) then
              begin
                if e1 < end_of_start s2 then (* s1 <= e1 < s2 <= e2 *)
                  Before
                else (* e1 >= s2 --------> s1 < s2 <= e1 < e2 *)
                  let before = make64 s1 (end_pred (end_of_start s2)) in
                  let inter  = make64 s2 e1 in
                  let after  = make64 (start_of_end (end_succ e1)) e2 in
                  BeforeAfter { before; inter; after }
              end
            else (* e1 > e2    ----- s1 < s2 <= e2 < e1 *)
                  let before = make64 s1 (end_pred (end_of_start s2)) in
                  let after  = make64 (start_of_end (end_succ e2)) e1 in
                  SplitEmpty { before; after; inter = t2 }
          end
      end
    else (* t1 > t2 ------------ s1 >= s2  --- s2 <= s1 *)
      begin
        let s1 = start_p t1 in
        let s2 = start_p t2 in
        let e1 = end_p t1 in
        let e2 = end_p t2 in
        if s1 = s2 then (* --> e1 > e2 -- e2 < e1  ---- s1 = s2 <= e2 < e1*)
              let after = make64 (start_of_end (end_succ e2)) e1 in
              AfterEmpty { inter = t2; after }
        else (* s1 > s2  --- s2 < s1 --> e2 < e1 *)
          begin
            if e1 = e2 then (* -----> s2 < s1 = e1 = e2 *)
                  let before = make64 s2 (end_pred (end_of_start s1)) in
                  EmptyBefore { inter = t1; before }
            else if e1 > e2 then (* e2 < e1 ----> s2 <= e2 ? s1 <= e1 *)
              begin
                if e2 < end_of_start s1 then (* s2 <= e2 < s1 <= e2 *)
                  After
                else
                  let before = make64 s2 (end_pred (end_of_start s1)) in
                  let inter  = make64 s1 e2 in
                  let after  = make64 (start_of_end (end_succ e2)) e1 in
                  AfterBefore { before; inter; after }
              end
            else (* e1 < e2       s2 <= s1 <= e1 < e2 *)
                  let before = make64 s2 (end_pred (end_of_start s1)) in
                  let after  = make64 (start_of_end (end_succ e1)) e2 in
                  EmptySplit { after; before; inter = t1; }
          end
 
      end

(* Interfacing to ints *)
let make s e =
  if e < s then
    invalid_argf "Not in order %d < %d" e s
  else if e < 1 then
    invalid_argf "End %d less than  1" e
  else if s < 1 then
    invalid_argf "Start %d less than 1" e
  else
    let s64 = start_of_int64 (Int64.of_int s) in
    let e64 = end_of_int64 (Int64.of_int e) in
    make64 s64 e64

let extend_one t =
  Int64.succ t

let merge t1 t2 =
  let s1 = start_p t1 in 
  let s2 = start_p t2 in 
  let e1 = end_p t1 in
  let e2 = end_p t2 in
  if end_succ e1 = (end_of_start s2) then
    Some (make64 s1 e2)
  else if end_succ e2 = (end_of_start s1) then
    Some (make64 s2 e1)
  else
    None

let merge3 t1 t2 t3 =
  match merge t1 t2 with
  | Some p1d -> merge p1d t3
  | None     -> None

let merge_exn t1 t2 =
  match merge t1 t2 with
    | Some m -> m
    | None   -> invalid_argf "not merge-able %s %s"
                    (to_string_interval t1) (to_string_interval t2)


let sample_inter_diff () =
  let mx = 100 in
  let ic = 50 in
  let pair () =
    let s = 1 + Random.int mx in
    let e = s + Random.int ic in
    make s e
  in
  let p1 = pair () in
  let p2 = pair () in
  p1, p2, inter_diff p1 p2

let test_inter_diff n =
  let rec error tr p1 p2 i p =
    if p then
      `Error (show_inter_diff_res tr, to_string p1, to_string p2)
    else
      loop (i + 1)
  and loop i =
    if i = n then `Ok else
      let p1, p2, tr = sample_inter_diff () in
      match tr with
      | Before  ->
          if end_ p1 < start p2 then
            loop i
          else
            `ThereShouldBeAnIntersect (p1, p2)
      | After  ->
          if end_ p2 < start p1 then
            loop i
          else
            `ThereShouldBeAnIntersect (p1, p2)
      | BeforeAfter { inter; before; after }          ->
          error tr p1 p2 i (merge before inter <> Some p1  &&
                              merge inter after <> Some p2)
      | BeforeEmpty { inter; before }                 ->
          error tr p1 p2 i (merge before inter <> Some p1)
      | EmptyAfter { inter; after}                    ->
          error tr p1 p2 i (merge inter after <> Some p2)
      | EmptyBefore { inter; before }                 ->
          error tr p1 p2 i (merge before inter <> Some p2)
      | EmptyEmpty { inter }                          ->
          error tr p1 p2 i (inter <> p1 || inter <> p2)
      | EmptySplit { inter; before; after}            ->
          error tr p1 p2 i (merge3 before inter after <> Some p2)
      | AfterEmpty { inter; after }                   ->
          error tr p1 p2 i (merge inter after <> Some p1)
      | AfterBefore { inter; after; before} ->
          error tr p1 p2 i (merge inter after <> Some p1 &&
                            merge before inter <> Some p2)
      | SplitEmpty { inter; before; after }  ->
          error tr p1 p2 i (merge3 before inter after <> Some p1)
  in
  loop 0


