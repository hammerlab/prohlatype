(* Ordered pairs. *)

open Util

type t = int * int [@@deriving eq, ord, show]

let start (s, _) = s
let end_ (_, e) = e

let width (s, e) = e - s + 1

let make s e =
  if e < s then
    invalid_argf "Not in order %d < %d" e s
  else
    (s, e)

let inside i (s, e) =
  s <= i && i <= e

let to_string (s, e) =
  sprintf "(%d,%d)" s e

type inter_diff_res =
  (* There is no intersection, the first interval comes Before the second
      or it comes After the second. *)
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

let inter_diff t1 t2 =
  let s1, e1 = t1 in
  let s2, e2 = t2 in
  if s1 = e1 then begin
  (* Handle the cases where the first interval is 1 wide to eliminate some
      corner cases. *)

    if s1 > s2 then begin                                           (* S2_1 *)
      if e2 < s1 then                                               (* E2_1 *)
        After
      else if e2 = s1 then                                          (* E2_2 *)
        EmptyBefore { inter = t1;                               before = (s2, s1 - 1) }
      else begin (* e2 > s1 *)                                      (* E2_3 *)
        EmptySplit  { inter = t1;                               before = (s2, s1 - 1); after = (e1 + 1, e2) }
      end

    end else if s1 = s2 then begin                                  (* S2_2 *)
      if e2 < s1 then                                               (* E2_1 *)
        invalid_argf  "broken invariant (%d,%d) (%d, %d)" s1 e1 s2 e2
      else if e2 = s1 then                                          (* E2_2 *)
        EmptyEmpty { inter = t1 }
      else begin (* e2 > s1 *)                                      (* E2_3 *)
        EmptyAfter { inter = t1;                                after = (e1 + 1, e2) }
      end

    end else (*if s1 < s2 then *) begin                             (* S2_3 *)
      if e2 <= s1 then                                        (* E2_1, E2_2 *)
        invalid_argf  "broken invariant (%d,%d) (%d, %d)" s1 e1 s2 e2
      else begin (* e2 > s1 *)                                      (* E2_3 *)
        Before
      end
    end

  end else begin (* s1 < e1 *)
    if s1 > s2 then                                                 (* S2_1 *)
      if e2 < s1 then                                               (* E2_1 *)
        After
      else if (*e2 >= s1 && *) e2 < e1 then                   (* E2_2, E2_3 *)
        AfterBefore { inter = (s1, e2); after = (e2 + 1, e1);   before = (s2, s1 - 1) }
      else if e2 = e1 then                                          (* E2_4 *)
        EmptyBefore { inter = t1;                               before = (s2, s1 - 1) }
      else (* e2 > e1 *)                                            (* E2_5 *)
        EmptySplit  { inter = t1;                               before = (s2, s1 - 1); after = (e1 + 1, e2) }

    else if s1 = s2 then                                            (* S2_2 *)
      if e2 < s1 then                                               (* E2_1 *)
        invalid_argf  "broken invariant (%d,%d) (%d, %d)" s1 e1 s2 e2
      else if (*e2 >= s1 && *) e2 < e1 then                   (* E2_2, E2_3 *)
        AfterEmpty { inter = (s1, e2); after = (e2 + 1, e1); }
      else if e2 = e1 then                                          (* E2_4 *)
        EmptyEmpty { inter = t1 }
      else (* e2 > e1 *)                                            (* E2_5 *)
        EmptyAfter { inter = t1;                                after = (e1 + 1, e2) }

    else if s1 < s2 && s2 < e1 then                                 (* S2_3 *)
      if e2 <= s1 then                                        (* E2_1, E2_2 *)
        invalid_argf  "broken invariant (%d,%d) (%d, %d)" s1 e1 s2 e2
      else if (*e2 > s1 && *) e2 < e1 then                          (* E2_3 *)
        SplitEmpty { inter = t2;        before = (s1, s2-1);    after = (e2+1, e1) }
      else if e2 = e1 then                                          (* E2_4 *)
        BeforeEmpty { inter = t2;       before = (s1, s2-1) }
      else (* e2 > e1 *)                                            (* E2_5 *)
        BeforeAfter { inter = (s2, e1); before = (s1, s2-1);    after = (e1 + 1, e2) }

    else if e1 = s2 then                                            (* S2_4 *)
      if e2 < e1 then                                   (* E2_1, E2_2, E2_3 *)
        invalid_argf  "broken invariant (%d,%d) (%d, %d)" s1 e1 s2 e2
      else if e2 = e1 then                                          (* E2_4 *)
        BeforeEmpty { inter = t2;       before = (s1, s2 - 1); }
      else (* e2 > e1 *)                                            (* E2_5 *)
        BeforeAfter { inter = (s2, e1); before = (s1, s2 - 1);  after = (e1 + 1, e2) }

    else (*if e1 < s2 then *)                                       (* S2_5 *)
      if e2 <= e1 then                            (* E2_1, E2_2, E2_3, E2_4 *)
        invalid_argf  "broken invariant (%d,%d) (%d, %d)" s1 e1 s2 e2
      else (* e2 > e1 *)                                            (* E2_5 *)
        Before
  end

let extend (s, e) n =
  if e + 1 = n then
    Some (s, n)
  else
    None

let extend_one (s, e) =
  (s, e + 1)

let merge p1 p2 =
  let s1, e1 = p1 in
  let s2, e2 = p2 in
  if e1 + 1 = s2 then
    Some (s1, e2)
  else if e2 + 1 = s1 then
    Some (s2, e1)
  else
    None

let merge3 p1 p2 p3 =
  match merge p1 p2 with
  | Some p1d -> merge p1d p3
  | None     -> None

let merge_exn p1 p2 =
  match merge p1 p2 with
  | Some m -> m
  | None   -> invalid_argf "not merge-able %s %s"
                  (to_string p1) (to_string p2)

(* Test logic *)
let sample_inter_diff () =
  let mx = 100 in
  let ic = 50 in
  let pair () =
    let s = Random.int mx in
    let e = s + Random.int ic in
    (s, e)
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

(* This is a more composable version of inter_diff. We either find an
    intersection or we find something that is 'before' any possible
    intersection and should be discarded.

    The return results are
    1. Possible before of first interval.
    2. Possible before of second interval.
    3. Possible interseciton.
    4. After of first interval, prepend to tail.
    5. After of second interval, prepend to tail.  *)
let split_inter_diff i1 i2 =
  match inter_diff i1 i2 with
  | After                                 -> None,        Some i2,     None,       Some i1,    None
  | AfterBefore { after; inter; before }  -> None,        Some before, Some inter, Some after, None
  | AfterEmpty { after; inter }           -> None,        None,        Some inter, Some after, None

  | Before                                -> Some i1,     None,        None,       None,       Some i2
  | BeforeAfter { before; inter; after }  -> Some before, None,        Some inter, None,       Some after
  | BeforeEmpty { before; inter }         -> Some before, None,        Some inter, None,       None

  | EmptyBefore { inter; before }         -> None,        Some before, Some inter, None,       None
  | EmptySplit { inter; before; after }   -> None,        Some before, Some inter, None,       Some after
  | EmptyAfter { inter; after }           -> None,        None,        Some inter, None,       Some after
  | EmptyEmpty { inter }                  -> None,        None,        Some inter, None,       None

  | SplitEmpty { before; after; inter }   -> Some before, None,        Some inter, Some after, None

let prepo p1 p2 = match p1, p2 with
  | None,          None          -> None
  | Some _,        None          -> p1
  | None,          Some _        -> p2
  | Some (s1, e1), Some (s2, e2) -> Some (s1, e2)                 (* assert (e1 + 1 = s2); *)

let split_inter_diff3 i1 i2 i3 =
  let b1, b2, i, a1, a2 = split_inter_diff i1 i2 in
  match i with
  | None     -> b1, b2, None, i, a1, a2, Some i3
  | Some i12 ->
    let b12, b3, i123, a12, a3 = split_inter_diff i12 i3 in
    prepo b1 b12,
    prepo b2 b12,
    b3,
    i123,
    prepo a12 a1,
    prepo a12 a2,
    a3

let split_inter_diff4 i1 i2 i3 i4 =
  let b1, b2, b3, i, a1, a2, a3 = split_inter_diff3 i1 i2 i3 in
  match i with
  | None      -> b1, b2, b3, None, i, a1, a2, a3, Some i4
  | Some i123 ->
    let b123, b4, i1234, a123, a4 = split_inter_diff i123 i4 in
    prepo b1 b123,
    prepo b2 b123,
    prepo b3 b123,
    b4,
    i1234,
    prepo a123 a1,
    prepo a123 a2,
    prepo a123 a3,
    a4

let aligned_inter_diff3 i1 i2 i3 =
  let prep p =
    let s, _e = p in
    function
      | None         -> Some p
      | Some (_s, e) -> Some (s, e)                               (* assert (_e + 1 = _s); *)
  in
  let i, m1, m2 =
    match inter_diff i1 i2 with
    | Before
    | After
    | AfterBefore _ | BeforeAfter _
    | BeforeEmpty _ | EmptyBefore _
    | EmptySplit _ | SplitEmpty _     -> invalid_argf "Different lengths!"
    | EmptyAfter { inter; after = a } -> inter, None,   Some a
    | EmptyEmpty { inter }            -> inter, None,   None
    | AfterEmpty { after = a; inter } -> inter, Some a, None
  in
  match inter_diff i i3 with
  | Before
  | After
  | AfterBefore _ | BeforeAfter _
  | BeforeEmpty _ | EmptyBefore _
  | EmptySplit _ | SplitEmpty _     -> invalid_argf "Different lengths!"
  | EmptyAfter { inter; after = a } -> inter, m1,        m2,        Some a
  | EmptyEmpty { inter }            -> inter, m1,        m2,        None
  | AfterEmpty { after = a; inter } -> inter, prep a m1, prep a m2, None

let aligned_inter_diff4 i1 i2 i3 i4 =
  let prep p =
    let s, _e = p in
    function
      | None         -> Some p
      | Some (_s, e) -> Some (s, e)                               (* assert (_e + 1 = _s); *)
  in
  let i, m1, m2 =
    match inter_diff i1 i2 with
    | Before
    | After
    | AfterBefore _ | BeforeAfter _
    | BeforeEmpty _ | EmptyBefore _
    | EmptySplit _ | SplitEmpty _     -> invalid_argf "Different lengths!"
    | EmptyAfter { inter; after = a } -> inter, None,   Some a
    | EmptyEmpty { inter }            -> inter, None,   None
    | AfterEmpty { after = a; inter } -> inter, Some a, None
  in
  let i, m1, m2, m3 =
    match inter_diff i i3 with
    | Before
    | After
    | AfterBefore _ | BeforeAfter _
    | BeforeEmpty _ | EmptyBefore _
    | EmptySplit _ | SplitEmpty _     -> invalid_argf "Different lengths!"
    | EmptyAfter { inter; after = a } -> inter, m1,        m2,        Some a
    | EmptyEmpty { inter }            -> inter, m1,        m2,        None
    | AfterEmpty { after = a; inter } -> inter, prep a m1, prep a m2, None
  in
  match inter_diff i i4 with
  | Before
  | After
  | AfterBefore _ | BeforeAfter _
  | BeforeEmpty _ | EmptyBefore _
  | EmptySplit _ | SplitEmpty _     -> invalid_argf "Different lengths!"
  | EmptyAfter { inter; after = a } -> inter, m1,        m2,        m3,        Some a
  | EmptyEmpty { inter }            -> inter, m1,        m2,        m3,        None
  | AfterEmpty { inter; after = a } -> inter, prep a m1, prep a m2, prep a m3, None

let iter (s, e) ~f =
  for i = s to e do f i done

let fold (s, e) ~init ~f =
  let acc = ref init in
  for i = s to e do acc := f !acc i done;
  !acc

