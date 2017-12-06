(** A partition map is a data structure for a map over a partition of a set.

  The integer intervals are stored in a sorted order, so that finding
  intersections and merging can be accomplised by traversing the lists in the
  same order. *)

open Util

module Interval = struct

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

end (* Interval *)

module Set = struct

  type t = Interval.t list

  let of_interval (i : Interval.t) = [i]

  let invariant =
    let open Interval in
    let rec loop = function
      | []  -> true
      | h :: [] -> true
      | h1 :: h2 :: t ->
          match inter_diff h1 h2 with
          | Before -> loop (h2 :: t)
          | _      -> false
    in
    loop

  let to_string  =
    string_of_list ~sep:";" ~f:(fun i -> sprintf "%s" (Interval.to_string  i))

  let size = List.fold_left ~init:0 ~f:(fun a i -> a + Interval.width i)

  let inside i l =
    List.exists l ~f:(Interval.inside i)

  let universal = function
    | [i] -> true
    | _   -> false

  let compare s1 s2 =
    match s1, s2 with
    | i1 :: _ , i2 :: _ -> compare (Interval.start i1) (Interval.start i2)
    | _                 -> assert false

  (* Should I make the sets a real pair? *)
  let first_pos = function
    | []      -> assert false
    | s :: _  -> Interval.start s

  (* The elements in a set are ordered.  *)
  let split_if_in ii l =
    let open Interval in
    let before_r = ref None in
    let rec loop = function
      | []      -> l
      | h :: t  ->
          begin match inter_diff ii h with
          | Before                         -> l            (* l is ordered. *)
          | After                          -> h :: loop t
          | EmptyAfter { after; _}         -> before_r := Some [];
                                              after :: t
          | EmptyBefore { before; _ }      -> before_r := Some (of_interval before);
                                              t
          | EmptyEmpty  _                  -> before_r := Some [];
                                              t
          | EmptySplit { before; after; _} -> before_r := Some (of_interval before);
                                              after :: t
          | BeforeAfter _
          | BeforeEmpty _
          | AfterEmpty _
          | AfterBefore _
          | SplitEmpty _                   -> assert false
          end
    in
    let rest = loop l in
    !before_r, rest

  (* Zip together two non-intersecting (separate) sets. *)
  let merge_separate =
    let open Interval in
    let rec start l1 l2 = match l1, l2 with
      | _,  []            -> l1
      | [], _             -> l2
      | h1 :: t1
      , h2 :: t2 ->
          begin match inter_diff h1 h2 with
          | Before        -> loop h1 t1 l2
          | After         -> loop h2 l1 t2
          | BeforeAfter _
          | BeforeEmpty _
          | EmptyAfter _
          | EmptyBefore _
          | EmptyEmpty _
          | EmptySplit _
          | AfterEmpty _
          | AfterBefore _
          | SplitEmpty _  -> invalid_argf "Not separate!"
          end
    and loop ps l1 l2 = match l1, l2 with
      | [],  []             ->  [ps]
      | h1 :: t1, []        ->  begin match merge ps h1 with
                                | None -> ps :: l1
                                | Some m1 -> loop m1 t1 []
                                end
      | [],       h2 :: t2  ->  begin match merge ps h2 with
                                | None -> ps :: l2
                                | Some m2 -> loop m2 [] t2
                                end
      | h1 :: t1, h2 :: t2  ->
          begin match inter_diff h1 h2 with
          | Before          ->  begin match merge ps h1 with
                                | None    -> ps :: loop h1 t1 l2
                                | Some m1 -> loop m1 t1 l2
                                end
          | After           ->  begin match merge ps h2 with
                                | None    -> ps :: loop h2 l1 t2
                                | Some m2 -> loop m2 l1 t2
                                end
          | BeforeAfter _
          | BeforeEmpty _
          | EmptyAfter _
          | EmptyBefore _
          | EmptyEmpty _
          | EmptySplit _
          | AfterEmpty _
          | AfterBefore _
          | SplitEmpty _  -> invalid_argf "Not separate!"
          end
    in
    start

  let all_intersections =
    let open Interval in
    let rec loop l1 l2 = match l1, l2 with
      | _,  []                                   -> [],           l1,              l2
      | [], _                                    -> [],           l1,              l2
      | h1 :: t1
      , h2 :: t2 ->
          begin match inter_diff h1 h2 with
          | Before                               -> let i, r1, r2 = loop t1 l2 in
                                                    i,            (h1 :: r1),     r2
          | After                                -> let i, r1, r2 = loop l1 t2 in
                                                    i,            r1,             (h2 :: r2)

          | BeforeAfter { before; inter; after } -> let i, r1, r2 = loop t1 (after :: t2) in
                                                    (inter :: i), (before :: r1), r2
          | BeforeEmpty { before; inter }        -> let i, r1, r2 = loop t1 t2 in
                                                    (inter :: i), (before :: r1), r2

          | EmptyAfter { inter; after }          -> let i, r1, r2 = loop t1 (after :: t2) in
                                                    (inter :: i), r1,             r2
          | EmptyBefore { inter; before }        -> let i, r1, r2 = loop t1 t2 in
                                                    (inter :: i), r1,             (before :: r2)
          | EmptyEmpty { inter }                 -> let i, r1, r2 = loop t1 t2 in
                                                    (inter :: i), r1,             r2
          | EmptySplit { before; after; inter }  -> let i, r1, r2 = loop t1 (after :: t2) in
                                                    (inter :: i), r1,             (before :: r2)

          | AfterEmpty { after; inter }          -> let i, r1, r2 = loop (after :: t1) t2 in
                                                    (inter :: i), r1,             r2
          | AfterBefore { after; inter; before } -> let i, r1, r2 = loop (after :: t1) t2 in
                                                    (inter :: i), r1,             (before :: r2)
          | SplitEmpty { after; before; inter }  -> let i, r1, r2 = loop (after :: t1) t2 in
                                                    (inter :: i), (before :: r1), r2
          end
    in
    loop

  let must_match_at_beginning s1 s2 =
    match s1, s2 with
    | [], []  -> invalid_argf "Empty sets!"
    | [], s   -> invalid_argf "Different lengths! s2: %s" (to_string s)
    | s , []  -> invalid_argf "Different lengths! s1: %s" (to_string s)
    | h1 :: t1
    , h2 :: t2 ->
      let open Interval in
      match inter_diff h1 h2 with
      | Before        -> invalid_argf "First %s is Before second: %s"
                            (to_string h1) (to_string h2)
      | After         -> invalid_argf "Second %s is Before First: %s"
                            (to_string h2) (to_string h1)
      | BeforeAfter _
      | BeforeEmpty _
      | SplitEmpty _  -> invalid_argf "First %s is Partially Before Second : %s"
                            (to_string h1) (to_string h2)
      | AfterBefore _
      | EmptyBefore _
      | EmptySplit _  -> invalid_argf "Second %s is Partially Before First : %s"
                            (to_string h2) (to_string h1)
      | EmptyAfter { inter; after } -> let i, r1, r2 = all_intersections t1 (after :: t2) in
                                       (inter :: i), r1, r2
      | EmptyEmpty { inter }        -> let i, r1, r2 = all_intersections t1 t2 in
                                       (inter :: i), r1, r2
      | AfterEmpty { after; inter } -> let i, r1, r2 = all_intersections (after :: t1) t2 in
                                       (inter :: i), r1, r2

  let prepend_if_not_none o l =
    match o with
    | None -> l
    | Some i -> i :: l

  let all_intersections3 =
    let rec loop l1 l2 l3 = match l1, l2, l3 with
      | [],  _,  _
      |  _, [],  _
      |  _,  _, []  -> [], l1, l2, l3
      | h1 :: t1
      , h2 :: t2
      , h3 :: t3    ->
        let b1, b2, b3, i, a1, a2, a3 = Interval.split_inter_diff3 h1 h2 h3 in
        let nt1 = prepend_if_not_none a1 t1 in
        let nt2 = prepend_if_not_none a2 t2 in
        let nt3 = prepend_if_not_none a3 t3 in
        let il, r1, r2, r3 = loop nt1 nt2 nt3 in
        prepend_if_not_none i il
        , prepend_if_not_none b1 r1
        , prepend_if_not_none b2 r2
        , prepend_if_not_none b3 r3
    in
    loop

  let must_match_at_beginning3 s1 s2 s3 =
    match s1, s2, s3 with
    | [],  _,  _  -> invalid_argf "Empty 1"
    |  _, [],  _  -> invalid_argf "Empty 2"
    |  _,  _, []  -> invalid_argf "Empty 3"
    | h1 :: t1
    , h2 :: t2
    , h3 :: t3    ->
        let inter, ho1, ho2, ho3 = Interval.aligned_inter_diff3 h1 h2 h3 in
        let nt1 = prepend_if_not_none ho1 t1 in
        let nt2 = prepend_if_not_none ho2 t2 in
        let nt3 = prepend_if_not_none ho3 t3 in
        let il, r1, r2, r3 = all_intersections3 nt1 nt2 nt3 in
        inter :: il, r1, r2, r3

  let all_intersections4 =
    let rec loop l1 l2 l3 l4 = match l1, l2, l3, l4 with
      | [],  _,  _,  _
      |  _, [],  _,  _
      |  _,  _, [],  _
      |  _,  _,  _, []  -> [], l1, l2, l3, l4
      | h1 :: t1
      , h2 :: t2
      , h3 :: t3
      , h4 :: t4        ->
        let b1, b2, b3, b4, i, a1, a2, a3, a4 = Interval.split_inter_diff4 h1 h2 h3 h4 in
        let nt1 = prepend_if_not_none a1 t1 in
        let nt2 = prepend_if_not_none a2 t2 in
        let nt3 = prepend_if_not_none a3 t3 in
        let nt4 = prepend_if_not_none a4 t4 in
        let il, r1, r2, r3, r4 = loop nt1 nt2 nt3 nt4 in
        prepend_if_not_none i il
        , prepend_if_not_none b1 r1
        , prepend_if_not_none b2 r2
        , prepend_if_not_none b3 r3
        , prepend_if_not_none b4 r4
    in
    loop

  let must_match_at_beginning4 s1 s2 s3 s4 =
    match s1, s2, s3, s4 with
    | [],  _,  _,  _ -> invalid_argf "Empty 1"
    |  _, [],  _,  _ -> invalid_argf "Empty 2"
    |  _,  _, [],  _ -> invalid_argf "Empty 3"
    |  _,  _,  _, [] -> invalid_argf "Empty 4"
    | h1 :: t1
    , h2 :: t2
    , h3 :: t3
    , h4 :: t4 ->
        let inter, ho1, ho2, ho3, ho4 = Interval.aligned_inter_diff4 h1 h2 h3 h4 in
        let nt1 = prepend_if_not_none ho1 t1 in
        let nt2 = prepend_if_not_none ho2 t2 in
        let nt3 = prepend_if_not_none ho3 t3 in
        let nt4 = prepend_if_not_none ho4 t4 in
        let il, r1, r2, r3, r4 = all_intersections4 nt1 nt2 nt3 nt4 in
        inter :: il, r1, r2, r3, r4

  let fold t ~init ~f =
    List.fold_left t ~init ~f:(fun init interval ->
      Interval.fold interval ~init ~f)

  let iter t ~f =
    fold t ~init:() ~f:(fun () i -> f i)

end (* Set *)

type ascending  (*= Ascending*)
type descending (*= Descending*)

(* Things start out in descending order when we construct the partition, but
  when we 'reverse' it they are constructed into an ascending order that is
  better for merging. *)

type +'a asc =
  | E                       (* empty, merging against this will fail! *)
  | U of Set.t * 'a           (* universal, hold size - 1 = end_ *)
  | S of (Set.t * 'a) list

type +'a desc = (Interval.t * 'a) list

type (_, +'a) t =
  | Asc   : 'a asc -> (ascending, 'a) t
  | Desc  : 'a desc -> (descending, 'a) t

let ascending_invariant =
  let rec loop = function
    | []           -> false             (* fail? don't expect empty lists *)
    | (s, _) :: [] -> Set.invariant s
    | (s1, _) :: (s2, v) :: t ->
        Set.invariant s1
        && Set.first_pos s1 < Set.first_pos s2
        && loop ((s2, v) :: t)
  in
  loop

let invariant : type o a. (o, a) t -> bool =
  function
    | Asc E      -> true
    | Asc U _    -> true
    | Asc (S la) -> ascending_invariant la
    | Desc _     -> true                               (* being a bit lazy atm. *)

(* Empty Constructors *)
let empty_d = Desc []
let empty_a = Asc E

(* Initializers *)
let init_first_d v =
  let i = Interval.make 0 0 in
  Desc [i, v]

let init_all_a ~size v =
  let i = Interval.make 0 (size - 1) in
  Asc (U (Set.of_interval i, v))

(* Properties *)
let asc_to_string la to_s =
  string_of_list la ~sep:"; " ~f:(fun (s, v) ->
      sprintf "[%s]:%s" (Set.to_string s) (to_s v))

let desc_to_string ld to_s =
  string_of_list ld ~sep:";" ~f:(fun (i, v) ->
    sprintf "%s:%s" (Interval.to_string i) (to_s v))

let to_string: type o a. (o, a) t -> (a -> string) -> string =
  fun t to_s -> match t with
    | Asc E         -> "Empty!"
    | Asc (U (s,v)) -> sprintf "%s:%s" (Set.to_string s) (to_s v)
    | Asc (S l)     -> asc_to_string l to_s
    | Desc ld       -> desc_to_string ld to_s

let size_a = List.fold_left ~init:0 ~f:(fun a (s, _) -> a + Set.size s)

let size_d = function
  | []            -> 0
  | ((i, _) :: _) -> Interval.end_ i + 1

let size : type o a. (o, a) t -> int = function
  | Asc E          -> 0
  | Asc (U (s, _)) -> Set.size s
  | Asc (S l)      -> size_a l
  | Desc l         -> size_d l

let length : type o a. (o, a) t -> int = function
  | Asc E     -> 0
  | Asc (U _) -> 1
  | Asc (S l) -> List.length l
  | Desc l    -> List.length l

let assoc_remove_and_get el list =
  let rec loop acc = function
    | []                      -> None
    | (e, v) :: t when e = el -> Some (v, (List.rev acc @ t))
    | h :: t                  -> loop (h :: acc) t
  in
  loop [] list

(* Conversion *)

(* [merge_or_add_to_end eq s v l] rebuild the elements of [l] such that if
   any of the values (snd) [eq v] then merge the sets [s] and (fst). If no
   values are equal add to the end of l. *)
let merge_or_add_to_end eq s v l =
  let rec loop = function
    | []     -> [s, v]
    | h :: t ->
        let s0, v0 = h in
        if eq v v0 then
          (Set.merge_separate s0 s, v0) :: t
        else
          h :: loop t
  in
  loop l

let map_with_full_check eq l ~f =
  List.fold_left l ~init:[] ~f:(fun acc (s, v) ->
      merge_or_add_to_end eq s (f v) acc)

(* let ascending_t l =
  List.fold_left l ~init:[] ~f:(fun acc (i, v) ->
    match assoc_remove_and_get v acc with
    | None           -> (v, Set.of_interval i) :: acc
    | Some (s, nacc) -> (v, i :: s) :: nacc)      (* Set abstraction is leaky atm, this reverses. *)
  |> List.map ~f:(fun (v, s) -> Set.first_pos s, s, v)
  |> List.sort ~cmp:(fun (p1, _, _) (p2, _, _) -> compare p1 p2)
  |> List.map ~f:(fun (_, s, v) -> (s,v)) *)

let ascending_t eq l =
  List.fold_left l ~init:[] ~f:(fun acc (i, v) ->
    merge_or_add_to_end eq (Set.of_interval i) v acc)
  |> List.sort ~cmp:(fun (s1, _) (s2, _) -> Set.compare s1 s2)

let asc_sets_to_str s =
  asc_to_string s (fun _ -> "")

let ascending eq = function
  | Desc l ->
    let a = ascending_t eq l in                             (* assert (ascending_invariant l); *)
    match a with
    | []      -> invalid_arg "Empty descending!"
    | [s,v]   -> if Set.universal s then
                   Asc (U (s, v))
                 else
                  invalid_argf "Single set but not universal? %s" (Set.to_string s)
    | lst     -> Asc (S lst)

let descending = function
  | Asc E          -> invalid_arg "Can't convert empty to descending"
  | Asc (U (s, v)) -> Desc [List.hd_exn s, v]
  | Asc (S l)      ->
      List.map l ~f:(fun (s, v) -> List.map s ~f:(fun i -> i, v))
      |> List.concat
      |> List.sort ~cmp:(fun (i1, _) (i2, _) -> Interval.compare i2 i1)
      |> fun l -> Desc l

(* Getters/Setters *)
let add v = function
  | Desc []                         -> Desc ((Interval.make 0 0, v) :: [])
  | Desc ((s, ov) :: t) when v = ov -> Desc ((Interval.extend_one s, v) :: t)
  | Desc (((s, _) :: _) as l)       -> let e = 1 + Interval.end_ s in
                                       Desc ((Interval.make e e, v) :: l)


let get t i = match t with
  | Asc E          -> invalid_arg "Can't get from empty"
  | Asc (U (_, v)) -> v
  | Asc (S l)      ->
      let rec loop = function
        | []          -> raise Not_found      (* Need a better failure mode. *)
        | (s, v) :: t ->
            if Set.inside i s then
              v
            else
              loop t
      in
      loop l

let set t i v = match t with
  | Asc E          -> invalid_arg "Can't set from empty"
  | Asc (U (s, _)) -> Asc (U (s, v))
  | Asc (S l)      ->
    let open Interval in
    let ii = make i i in
    let rec loop l = match l with
      | []      -> raise Not_found
      | h :: t  ->
          let s, ov = h in
          if v = ov && Set.inside i s then              (* No use splitting *)
            l
          else
            match Set.split_if_in ii s with
            | None,    _     -> h :: loop t
            | Some [], after -> (Set.of_interval ii, v) :: (after, ov) :: t
            (* Technically this isn't scanning l again to find the
               appropriate set for {v}, we're just inserting it and maintaing
               the invariant that the sets inside are ordered.
               I'm not actually using this method in ParPHMM so I'll avoid
               a better implementation for now. *)
            | Some be, after -> (be @ after, ov) :: (Set.of_interval ii, v) :: t
    in
    Asc (S (loop l))

let insert s v l =
  let sl = Set.first_pos s in
  let rec loop l = match l with
    | []     -> [s, v]
    | h :: t -> let so, _ = h in
                if sl < Set.first_pos so then
                  (s, v) :: l
                else
                  h :: loop t
  in
  loop l

let insert_if_not_empty s v l =
  if s = [] then
    l
  else
    insert s v l

let map_with_just_last_check ~f = function
  | []  -> []
  | (s,v) :: t ->
    let rec loop ps pv = function
      | []         -> [ps, pv]
      | (s,v) :: t ->
          let nv = f v in
          if nv = pv then
            loop (Set.merge_separate ps s) pv t
          else
            (ps, pv) :: loop s nv t
    in
    loop s (f v) t

(* The reason for all this logic. *)
let rec start2 eq f l1 l2 = match l1, l2 with
  | [],     []      -> []
  | [],      s      -> invalid_argf "Different lengths! l2: %s" (asc_sets_to_str s)
  |  s,     []      -> invalid_argf "Different lengths! l1: %s" (asc_sets_to_str s)
  | (s1, v1) :: t1
  , (s2, v2) :: t2  ->
      let intersect, r1, r2 = Set.must_match_at_beginning s1 s2 in
      let nt1 = insert_if_not_empty r1 v1 t1 in
      let nt2 = insert_if_not_empty r2 v2 t2 in
      let acc = [intersect, (f v1 v2)] in
      loop2 eq f acc nt1 nt2
and loop2 eq f acc l1 l2 = match l1, l2 with
  | [],     []      -> acc
  | [],      s      -> invalid_argf "Different lengths! l2: %s" (asc_sets_to_str s)
  |  s,     []      -> invalid_argf "Different lengths! l1: %s" (asc_sets_to_str s)
  | (s1, v1) :: t1
  , (s2, v2) :: t2  ->
      let intersect, r1, r2 = Set.must_match_at_beginning s1 s2 in
      let nt1 = insert_if_not_empty r1 v1 t1 in
      let nt2 = insert_if_not_empty r2 v2 t2 in
      let nv = f v1 v2 in
      let nacc = merge_or_add_to_end eq intersect nv acc in
      loop2 eq f nacc nt1 nt2
(* TODO: There is a bug here where I'm not checking for matching ends.
  * I should functorize or somehow parameterize the construction of these
  * such that I don't worry about this. *)
and merge ~eq t1 t2 f =
  match t1, t2 with
  | Asc E          , _
  | _              , Asc E           -> invalid_argf "Can't merge empty"
  | Asc (U (s, v1)), Asc (U (_s,v2)) -> Asc (U (s, f v1 v2))
  | Asc (U (_, v1)), Asc (S l2)      -> Asc (S (map_with_just_last_check l2 ~f:(fun v2 -> f v1 v2)))
  | Asc (S l1),      Asc (U (_, v2)) -> Asc (S (map_with_just_last_check l1 ~f:(fun v1 -> f v1 v2)))
  | Asc (S l1),      Asc (S l2)      -> Asc (S (start2 eq f l1 l2))

let rec start3 eq f l1 l2 l3 =
  match l1, l2, l3 with
  | [],     [],     []  -> []
  | [],      s,      _  -> invalid_argf "Different lengths! l2: %s" (asc_sets_to_str s)
  |  _,     [],      s  -> invalid_argf "Different lengths! l3: %s" (asc_sets_to_str s)
  |  s,      _,     []  -> invalid_argf "Different lengths! l1: %s" (asc_sets_to_str s)
  | (s1, v1) :: t1
  , (s2, v2) :: t2
  , (s3, v3) :: t3      ->
      let intersect, r1, r2, r3 = Set.must_match_at_beginning3 s1 s2 s3 in
      let nt1 = insert_if_not_empty r1 v1 t1 in
      let nt2 = insert_if_not_empty r2 v2 t2 in
      let nt3 = insert_if_not_empty r3 v3 t3 in
      let acc = [intersect, (f v1 v2 v3)] in
      loop3 eq f acc nt1 nt2 nt3
and loop3 eq f acc l1 l2 l3 =
  match l1, l2, l3 with
  | [],     [],     []  -> acc     (* We insert at the end, thereby preserving order *)
  | [],      s,      _  -> invalid_argf "Different lengths! l2: %s" (asc_sets_to_str s)
  |  _,     [],      s  -> invalid_argf "Different lengths! l3: %s" (asc_sets_to_str s)
  |  s,      _,     []  -> invalid_argf "Different lengths! l1: %s" (asc_sets_to_str s)
  | (s1, v1) :: t1
  , (s2, v2) :: t2
  , (s3, v3) :: t3      ->
      let intersect, r1, r2, r3 = Set.must_match_at_beginning3 s1 s2 s3 in
      let nt1 = insert_if_not_empty r1 v1 t1 in
      let nt2 = insert_if_not_empty r2 v2 t2 in
      let nt3 = insert_if_not_empty r3 v3 t3 in
      let nv = f v1 v2 v3 in
      let nacc = merge_or_add_to_end eq intersect nv acc in
      loop3 eq f nacc nt1 nt2 nt3
and merge3 ~eq t1 t2 t3 f =
  match t1, t2, t3 with
  | Asc E          , _              , _
  | _              , Asc E          , _
  | _              , _              , Asc E           -> invalid_argf "Can't merge3 empty"

  | Asc (U (s, v1)), Asc (U (_, v2)), Asc (U (_, v3)) -> Asc (U (s, f v1 v2 v3))

  | Asc (U (_, v1)), Asc (U (_, v2)), Asc (S l3)      -> Asc (S (map_with_full_check eq l3 ~f:(fun v3 -> f v1 v2 v3)))
  | Asc (U (_, v1)), Asc (S l2),      Asc (U (_, v3)) -> Asc (S (map_with_full_check eq l2 ~f:(fun v2 -> f v1 v2 v3)))
  | Asc (S l1),      Asc (U (_, v2)), Asc (U (_, v3)) -> Asc (S (map_with_full_check eq l1 ~f:(fun v1 -> f v1 v2 v3)))

  | Asc (U (_, v1)), Asc (S l2),      Asc (S l3)      -> Asc (S (start2 eq (fun v2 v3 -> f v1 v2 v3) l2 l3))
  | Asc (S l1),      Asc (U (_, v2)), Asc (S l3)      -> Asc (S (start2 eq (fun v1 v3 -> f v1 v2 v3) l1 l3))
  | Asc (S l1),      Asc (S l2),      Asc (U (_, v3)) -> Asc (S (start2 eq (fun v1 v2 -> f v1 v2 v3) l1 l2))

  | Asc (S l1),      Asc (S l2),      Asc (S l3)      -> Asc (S (start3 eq f l1 l2 l3))

let rec start4 eq f l1 l2 l3 l4 =
  match l1, l2, l3, l4 with
  | [],     [],     [],     []      -> []
  | [],      s,      _,      _      -> invalid_argf "Different lengths! l2: %s" (asc_sets_to_str s)
  |  _,     [],      s,      _      -> invalid_argf "Different lengths! l3: %s" (asc_sets_to_str s)
  |  _,      _,     [],      s      -> invalid_argf "Different lengths! l4: %s" (asc_sets_to_str s)
  |  s,      _,      _,     []      -> invalid_argf "Different lengths! l1: %s" (asc_sets_to_str s)
  | (s1, v1) :: t1
  , (s2, v2) :: t2
  , (s3, v3) :: t3
  , (s4, v4) :: t4                  ->
      let intersect, r1, r2, r3, r4 = Set.must_match_at_beginning4 s1 s2 s3 s4 in
      let nt1 = insert_if_not_empty r1 v1 t1 in
      let nt2 = insert_if_not_empty r2 v2 t2 in
      let nt3 = insert_if_not_empty r3 v3 t3 in
      let nt4 = insert_if_not_empty r4 v4 t4 in
      let acc = [intersect, (f v1 v2 v3 v4)] in
      loop4 eq f acc nt1 nt2 nt3 nt4
and loop4 eq f acc l1 l2 l3 l4 =
  match l1, l2, l3, l4 with
  | [],     [],     [],     []      -> acc     (* We insert at the end, thereby preserving order *)
  | [],      s,      _,      _      -> invalid_argf "Different lengths! l2: %s" (asc_sets_to_str s)
  |  _,     [],      s,      _      -> invalid_argf "Different lengths! l3: %s" (asc_sets_to_str s)
  |  _,      _,     [],      s      -> invalid_argf "Different lengths! l4: %s" (asc_sets_to_str s)
  |  s,      _,      _,     []      -> invalid_argf "Different lengths! l1: %s" (asc_sets_to_str s)
  | (s1, v1) :: t1
  , (s2, v2) :: t2
  , (s3, v3) :: t3
  , (s4, v4) :: t4                  ->
      let intersect, r1, r2, r3, r4 = Set.must_match_at_beginning4 s1 s2 s3 s4 in
      let nt1 = insert_if_not_empty r1 v1 t1 in
      let nt2 = insert_if_not_empty r2 v2 t2 in
      let nt3 = insert_if_not_empty r3 v3 t3 in
      let nt4 = insert_if_not_empty r4 v4 t4 in
      let nv = f v1 v2 v3 v4 in
      let nacc = merge_or_add_to_end eq intersect nv acc in
      loop4 eq f nacc nt1 nt2 nt3 nt4
(* This method is tail recursive, and by default we pay the cost of inserting
   an element at the end, each time, Hopefully, merging, due to {eq}, instead into
   the accumulator-list will effectively constrain the size of the resulting
   list such that the cost is amortized. *)
and merge4 ~eq t1 t2 t3 t4 f =
  match t1, t2, t3, t4 with
  | Asc E,          _,              _,              _
  | _,              Asc E,          _,              _
  | _,              _,              Asc E,          _
  | _,              _,              _,              Asc E           -> invalid_argf "Can't merge empty4"

  (* 0 S's, 4 U's *)
  | Asc (U (s,v1)), Asc (U (_,v2)), Asc (U (_,v3)), Asc (U (_,v4))  -> Asc (U (s, f v1 v2 v3 v4))

  (* 1 S's, 3 U's *)
  | Asc (S l1),     Asc (U (_,v2)), Asc (U (_,v3)), Asc (U (_,v4))  -> Asc (S (map_with_full_check eq l1 ~f:(fun v1 -> f v1 v2 v3 v4)))
  | Asc (U (_,v1)), Asc (S l2),     Asc (U (_,v3)), Asc (U (_,v4))  -> Asc (S (map_with_full_check eq l2 ~f:(fun v2 -> f v1 v2 v3 v4)))
  | Asc (U (_,v1)), Asc (U (_,v2)), Asc (S l3),     Asc (U (_,v4))  -> Asc (S (map_with_full_check eq l3 ~f:(fun v3 -> f v1 v2 v3 v4)))
  | Asc (U (_,v1)), Asc (U (_,v2)), Asc (U (_,v3)), Asc (S l4)      -> Asc (S (map_with_full_check eq l4 ~f:(fun v4 -> f v1 v2 v3 v4)))

  (* 2 S's, 2 U's *)
  | Asc (S l1),     Asc (S l2),     Asc (U (_,v3)), Asc (U (_,v4))  -> Asc (S (start2 eq (fun v1 v2 -> f v1 v2 v3 v4) l1 l2))
  | Asc (S l1),     Asc (U (_,v2)), Asc (S l3),     Asc (U (_,v4))  -> Asc (S (start2 eq (fun v1 v3 -> f v1 v2 v3 v4) l1 l3))
  | Asc (S l1),     Asc (U (_,v2)), Asc (U (_,v3)), Asc (S l4)      -> Asc (S (start2 eq (fun v1 v4 -> f v1 v2 v3 v4) l1 l4))
  | Asc (U (_,v1)), Asc (S l2),     Asc (S l3),     Asc (U (_,v4))  -> Asc (S (start2 eq (fun v2 v3 -> f v1 v2 v3 v4) l2 l3))
  | Asc (U (_,v1)), Asc (S l2),     Asc (U (_,v3)), Asc (S l4)      -> Asc (S (start2 eq (fun v2 v4 -> f v1 v2 v3 v4) l2 l4))
  | Asc (U (_,v1)), Asc (U (_,v2)), Asc (S l3),     Asc (S l4)      -> Asc (S (start2 eq (fun v3 v4 -> f v1 v2 v3 v4) l3 l4))

  (* 3 S's, 1 U's *)
  | Asc (S l1),     Asc (S l2),     Asc (S l3),     Asc (U (_,v4))  -> Asc (S (start3 eq (fun v1 v2 v3 -> f v1 v2 v3 v4) l1 l2 l3))
  | Asc (S l1),     Asc (S l2),     Asc (U (_,v3)), Asc (S l4)      -> Asc (S (start3 eq (fun v1 v2 v4 -> f v1 v2 v3 v4) l1 l2 l4))
  | Asc (S l1),     Asc (U (_,v2)), Asc (S l3),     Asc (S l4)      -> Asc (S (start3 eq (fun v1 v3 v4 -> f v1 v2 v3 v4) l1 l3 l4))
  | Asc (U (_,v1)), Asc (S l2),     Asc (S l3),     Asc (S l4)      -> Asc (S (start3 eq (fun v2 v3 v4 -> f v1 v2 v3 v4) l2 l3 l4))

  (* 4 S's, 0 U's *)
  | Asc (S l1),     Asc (S l2),     Asc (S l3),     Asc (S l4)      -> Asc (S (start4 eq f l1 l2 l3 l4))

let fold_values : type o a. (o, a) t -> init:'b -> f:('b -> a -> 'b) -> 'b =
  fun l ~init ~f -> match l with
    | Desc ld         -> List.fold_left ld ~init ~f:(fun acc (_i, v) -> f acc v)
    | Asc E           -> invalid_arg "Can't fold_values on empty!"
    | Asc (U (_, v))  -> f init v
    | Asc (S la)      -> List.fold_left la ~init ~f:(fun acc (_l, v) -> f acc v)

let fold_set_sizes_and_values : type o a. (o, a) t -> init:'b -> f:('b -> int -> a -> 'b) -> 'b =
  fun l ~init ~f ->
    let ascf = List.fold_left ~init ~f:(fun acc (l, v) -> f acc (Set.size l) v) in
    match l with
    | Desc ld         -> ascf (ascending_t (fun x y -> x = y) ld)     (* TODO: Probably could be faster. *)
    | Asc E           -> invalid_arg "Can't fold_set_sizes_and_values on empty!"
    | Asc (U (s, v))  -> Set.fold s ~init ~f:(fun acc l -> f acc l v)
    | Asc (S la)      -> ascf la

let fold_indices_and_values :
  type o a. (o, a) t -> init:'b -> f:('b -> int -> a -> 'b) -> 'b =
  fun l ~init ~f -> match l with
    | Desc ld -> List.fold_left ld ~init ~f:(fun init (l, v) ->
                    Interval.fold l ~init ~f:(fun acc i -> f acc i v))
    | Asc E          -> invalid_arg "Can't fold_indices_and_values on empty!"
    | Asc (U (s, v)) -> Set.fold s ~init ~f:(fun acc i -> f acc i v)
    | Asc (S la)     -> List.fold_left la ~init ~f:(fun init (s, v) ->
                          Set.fold s ~init ~f:(fun acc i -> f acc i v))


let map : type o a. (o, a) t -> ('b -> 'b -> bool) -> f:(a -> 'b) -> (o, 'b) t =
  fun t eq ~f -> match t with
    | Desc ld         -> Desc (List.map_snd ld ~f)
    | Asc E           -> invalid_argf "Can't map empty!"
    | Asc (U (s, v))  -> Asc (U (s, f v))
    | Asc (S la)      -> Asc (S (map_with_full_check eq la ~f))

let iter_set : type o a. (o, a) t -> f:(int -> a -> unit) -> unit =
  fun t ~f -> match t with
    | Desc ld ->
        List.iter ld ~f:(fun (s, v) ->
          Interval.iter s ~f:(fun i -> f i v))
    | Asc E   -> invalid_argf "Can't iter_set empty"
    | Asc (U (s, v))  -> Set.iter s ~f:(fun i -> f i v)
    | Asc (S la)      -> List.iter la ~f:(fun (l, v) ->
                          List.iter l ~f:(Interval.iter ~f:(fun i -> f i v)))

let to_array = function
  | Asc E           -> invalid_argf "Can't to_array empty"
  | Asc (U (s, v))  -> Array.make (Set.size s) v
  | Asc (S la)      ->
    match la with
    | []  -> [||]
    | h :: t ->
        let s, v = h in
        let n = List.fold_left la ~init:0 ~f:(fun a (s, _) -> a + Set.size s) in
        let r = Array.make n v in
        let fill s v = Set.iter s ~f:(fun i -> r.(i) <- v) in
        fill s v;
        List.iter t ~f:(fun (s, v) -> fill s v);
        r
