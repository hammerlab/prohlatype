
open Util

module Interval (*: (sig

  type t = int * int

  val make : int -> int -> t

  val to_string : t -> string

  type inter_diff_res =
    | Before
    | After
    (* Different possible intersections. The way to read these is that the
       first word describes the relationship of result to the intersection
       for the first interval argument to {inter_diff} and the second word
       for the second argument:

       For example, if inter_diff i1 i2 = EmptySplit {inter; before, after}
       means that there is no resulting difference for the intersection of i1
       and i2 as it pertains to i1, it is empty, but that the intersection
       splits i2 into {before} and {after}.

       If i1 = (3,3) and i2 = (2,6) then
       inter_diff i1 i2 =
         EmptySplit {inter = (3,3); before = (2,2); after = (4,6)}.

       As a consequence the record labels are changing and may pertain
       to the the remaining element of {i1} or {i2}.

    *)
    | BeforeAfter of { inter: t; before: t; after : t }
    | BeforeEmpty of { inter: t; before: t }

    | EmptyAfter of  { inter: t; after: t }
    | EmptyBefore of { inter: t; before: t }
    | EmptyEmpty of  { inter: t }
    | EmptySplit of  { inter: t; before: t; after: t}

    | AfterEmpty of  { inter: t; after: t }
    | AfterBefore of { inter: t; after: t;  before: t}

    | SplitEmpty of  { inter: t; before: t; after : t}

  val inter_diff : t -> t -> inter_diff_res

  val extend : t -> int -> t option
  val merge : t -> t -> t option

end*) = struct

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
    | Before
    | After
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

  let inter_diff t1  t2 =
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

  let before p1 p2 =
    let _, e1 = p1 in
    let s2, _ = p2 in
    e1 < s2

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

  let aligned_inter_diff4 i1 i2 i3 i4 =
    let prep p =
      let s, _e = p in
      function
        | None         -> Some p
        | Some (_s, e) -> (*assert (_e + 1 = _s);*) Some (s, e)
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
      | AfterEmpty { inter; after = a } -> inter, Some a, None
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
      | AfterEmpty { inter; after = a } -> inter, prep a m1, prep a m2, None
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

  let inside i l =
    List.exists l ~f:(Interval.inside i)

  (* The elements in a set are ordered. *)
  let split_if_in ii l =
    assert (invariant l);
    let open Interval in
    let found = ref false in
    let rec loop = function
      | []      -> l
      | h :: t  ->
          begin match inter_diff ii h with
          | Before                         -> l            (* l is ordered. *)
          | After                          -> h :: loop t
          | EmptyAfter { after; _}         -> found := true;
                                              after :: t
          | EmptyBefore { before; _ }      -> found := true;
                                              before :: t
          | EmptyEmpty  _                  -> found := true;
                                              t
          | EmptySplit { before; after; _} -> found := true;
                                              before :: after :: t
          | BeforeAfter _
          | BeforeEmpty _
          | AfterEmpty _
          | AfterBefore _
          | SplitEmpty _                   -> assert false
          end
    in
    let s = loop l in
    !found, s

  (* Should I make the sets a real pair? *)
  let first_pos = function
    | []      -> assert false
    | s :: _  -> Interval.start s

  let merge l1 l2 =
    assert (invariant l1);
    assert (invariant l2);
    let rec loop ll1 ll2 =
      match ll1 with
      | []        -> ll2
      | h1 :: t1  ->
          begin match ll2 with
          | []        -> ll1
          | h2 :: t2  ->
              begin match Interval.merge h1 h2 with
              | Some m -> m :: loop t1 t2
              | None   -> if Interval.before h1 h2 then
                            h1 :: loop t1 ll2
                          else (* after *)
                            h2 :: loop ll1 t2
              end
          end
    in
    let ml = loop l1 l2 in
    assert (invariant ml);
    ml

  let inter_diff l1 l2 =
    assert (invariant l1);
    assert (invariant l2);
    let rec loop ll1 ll2 =
      match ll1 with
      | []        -> ll2
      | h1 :: t1  ->
          begin match ll2 with
          | []        -> ll1
          | h2 :: t2  ->
              begin match Interval.merge h1 h2 with
              | Some m -> m :: loop t1 t2
              | None   -> if Interval.before h1 h2 then
                            h1 :: loop t1 ll2
                          else (* after *)
                            h2 :: loop ll1 t2
              end
          end
    in
    let ml = loop l1 l2 in
    assert (invariant ml);
    ml

end

module Pm : sig
  type order
  type ascending
  type descending

  type ('o, +'a) t

  val empty : (descending, 'a) t
  val place_holder : (ascending, 'a) t

  val first_d : 'a -> (descending, 'a) t

  val init : int -> 'a -> (ascending, 'a) t

  val to_string : (_, 'a) t -> ('a -> string) -> string

  val size_a : (ascending, 'a) t -> int
  val size_d : (descending, 'a) t -> int

  val length : (_, 'a) t -> int

  val ascending : (descending, 'a) t -> (ascending, 'a) t
  val descending : (ascending, 'a) t -> (descending, 'a) t

  val add : 'a -> (descending, 'a) t -> (descending, 'a) t

  val get : (ascending, 'a) t -> int -> 'a
  val set : (ascending, 'a) t -> int -> 'a -> (ascending, 'a) t

  val merge : (ascending, 'a) t
            -> (ascending, 'b) t
            -> ('a -> 'b -> 'c)
            -> (ascending, 'c) t

  val merge4 : (ascending, 'a) t
             -> (ascending, 'b) t
             -> (ascending, 'c) t
             -> (ascending, 'd) t
             -> ('a -> 'b -> 'c -> 'd -> 'e)
             -> (ascending, 'e) t

  val fold : (_, 'a) t
           -> init:'b
           -> f:('b -> Interval.t -> 'a -> 'b)
           -> 'b

  val map : ('c, 'a) t -> f:('a -> 'b) -> ('c, 'b) t

  val iter_set : ('c, 'a) t -> f:(int -> 'a -> 'b) -> unit

end = struct

  type order
  type ascending = private order
  type descending = private order
  type ('o, 'a) t = (Interval.t * 'a) list

  let ascending_invariant (l : (ascending, 'a) t)  =
    let rec loop = function
      | []
      | _ :: [] -> true
      | (s1, _) :: (s2, v) :: t ->
          Interval.start s1 < Interval.start s2 &&
          Interval.end_ s1 + 1 = Interval.start s2 &&
              loop ((s2, v) :: t)
    in
    loop l

  let empty = []
  let place_holder = []

  let first_d v =
    [Interval.make 0 0, v]

  let init size v =
    [Interval.make 0 size, v]

  let to_string l to_s =
    List.map l ~f:(fun (i, v) ->
        sprintf "%s:%s" (Interval.to_string  i) (to_s v))
    |> String.concat ~sep:";"

  let size_a l =
    let rec loop a = function
      | []           -> a
      | (i, _) :: tl -> loop (a + Interval.width i) tl
    in
    loop 0 l

  let size_d l = match l with
    | []           -> 0
    | (s, _)  :: _ -> Interval.end_ s + 1

  let length = List.length

  let add v l = match l with
    | []                       -> (Interval.make 0 0, v) :: []
    | (s, ov) :: t when v = ov -> (Interval.extend_one s, v) :: t
    | (s, _) :: _              -> let e = 1 + Interval.end_ s in
                                  (Interval.make e e, v) :: l

  let ascending = List.rev
  let descending = List.rev

  let get l i =
    let rec loop = function
      | []          -> raise Not_found      (* Need a better failure mode. *)
      | (s, v) :: t ->
          if Interval.inside i s then
            v
          else
            loop t
    in
    loop l

  let set l i v =
    let open Interval in
    let ii = make i i in
    let rec loop l = match l with
      | []     -> raise Not_found
      | h :: t ->
          let s, ov = h in
          match inter_diff ii s with
          | Before                              -> raise Not_found
          | After                               -> loop2 s ov t
          | EmptyAfter { inter; after }         -> if v = ov then l else
                                                      (inter, v) :: (after, ov) :: t
          | EmptyBefore { before; inter }       -> if v = ov then l else
                                                      (before, ov) :: (inter, v) :: t
          | EmptyEmpty  _                       -> if v = ov then l else
                                                      (s, v) :: t
          | EmptySplit { before; after; inter } -> if v = ov then l else
                                                      (before, ov) :: (inter, v) :: (after, ov) :: t
          | AfterBefore _ | AfterEmpty _
          | BeforeAfter _ | BeforeEmpty _
          | SplitEmpty _                        -> assert false
    and loop2 ps pv l = match l with
      | []            -> raise Not_found
      | h :: t        ->
          let s, ov = h in
          begin match inter_diff ii s with
          | Before                              -> raise Not_found
          | After                               -> if t = [] then
                                                    raise Not_found
                                                  else
                                                    (ps, pv) :: loop2 s ov t
          | EmptyAfter { inter; after }         -> if v = ov then
                                                    (ps, pv) :: l
                                                  else if v = pv then
                                                    (merge_exn ps inter, pv) :: (after, ov) :: t
                                                  else
                                                    (ps, pv) :: (inter, v) :: (after, ov) :: t
          | EmptyBefore { before; inter }       -> if v = ov then
                                                    (ps, pv) :: l
                                                  else
                                                    (ps, pv) :: (before, ov) :: (inter, v) :: t
          | EmptyEmpty  _                       -> if v = ov then
                                                    (ps, pv) :: l
                                                  else
                                                    (ps, pv) :: (s, v) :: t
          | EmptySplit { before; after; inter } -> if v = ov then
                                                    (ps, pv) :: l
                                                  else
                                                    (ps, pv) :: (before, ov) :: (inter, v) :: (after, v) :: t
          | AfterBefore _ | AfterEmpty _
          | BeforeAfter _ | BeforeEmpty _
          | SplitEmpty _                        -> assert false
          end
    in
    loop l

  let debug_ref = ref false

  let test ~n ~m =
    let values = [| 'A'; 'B'; 'C' |] in
    let vs = Array.length values in
    let rec loop p k =
      if k = m then Ok p else
        let i = Random.int n in
        let j = Random.int vs in
        if !debug_ref then printf "%d %c\n" i values.(j);
        let np = set p i values.(j) in
        let vv = get np i in
        if vv <> values.(j) then
          Error (p, np, i, j, false)
        else if not (ascending_invariant np) then
          Error (p, np, i, j, true)
        else
          loop np (k + 1)
    in
    loop (init n values.(0)) 0

  let merge l1 l2 f =
    let rec loop l1 l2 =
      match l1, l2 with
      | [], []  -> []
      | [], s   -> invalid_argf "Different lengths! l2: %s" (to_string s (fun _ -> ""))
      | s , []  -> invalid_argf "Different lengths! l1: %s" (to_string s (fun _ -> ""))
      | h1 :: t1, h2 :: t2 ->
          let i1, v1 = h1 in
          let i2, v2 = h2 in
          let open Interval in
          match inter_diff i1 i2 with
          | Before
          | After
          | AfterBefore _
          | BeforeAfter _
          | BeforeEmpty _
          | EmptyBefore _
          | EmptySplit _
          | SplitEmpty _                -> invalid_argf "Different lengths!"
          | EmptyAfter { inter; after = a } -> (inter, f v1 v2) :: loop t1              ((a, v2) :: t2)
          | EmptyEmpty { inter }            -> (inter, f v1 v2) :: loop t1              t2
          | AfterEmpty { inter; after = a } -> (inter, f v1 v2) :: loop ((a, v1) :: t1) t2
      in
      loop l1 l2

  let merge4 l1 l2 l3 l4 f =
    let open Interval in
    let rec loop l1 l2 l3 l4 =
      match l1, l2, l3, l4 with
      | [], [], [], [] -> []
      | [],  _,  _,  _ -> invalid_argf "Different lengths!"
      |  _, [],  _,  _ -> invalid_argf "Different lengths!"
      |  _,  _, [],  _ -> invalid_argf "Different lengths!"
      |  _,  _,  _, [] -> invalid_argf "Different lengths!"
      | h1 :: t1
      , h2 :: t2
      , h3 :: t3
      , h4 :: t4 ->
          let i1, v1 = h1 in
          let i2, v2 = h2 in
          let i3, v3 = h3 in
          let i4, v4 = h4 in
          let inter, m1, m2, m3, m4 = Interval.aligned_inter_diff4 i1 i2 i3 i4 in
          (*printf "--%s %s %s %s -> %s %s %s %s %s--\n%!"
            (Interval.to_string i1)
            (Interval.to_string i2)
            (Interval.to_string i3)
            (Interval.to_string i4)
            (Interval.to_string inter)
            (match m1 with | None -> "none" | Some i -> (Interval.to_string i))
            (match m2 with | None -> "none" | Some i -> (Interval.to_string i))
            (match m3 with | None -> "none" | Some i -> (Interval.to_string i))
            (match m4 with | None -> "none" | Some i -> (Interval.to_string i)); *)
          let lt1 = match m1 with | None -> t1 | Some j1 -> (j1, v1) :: t1 in
          let lt2 = match m2 with | None -> t2 | Some j2 -> (j2, v2) :: t2 in
          let lt3 = match m3 with | None -> t3 | Some j3 -> (j3, v3) :: t3 in
          let lt4 = match m4 with | None -> t4 | Some j4 -> (j4, v4) :: t4 in
          (* TODO, add loop2 like logic to merge when possible. *)
          (inter, f v1 v2 v3 v4) :: loop lt1 lt2 lt3 lt4
      in
      loop l1 l2 l3 l4

  let fold l ~init ~f =
    List.fold_left l ~init ~f:(fun acc (s, v) -> f acc s v)

  let map = List.map_snd

  let iter_set l ~f =
    List.iter l ~f:(fun (s, v) ->
      Interval.iter s ~f:(fun i -> f i v))

end (* Pm *)
