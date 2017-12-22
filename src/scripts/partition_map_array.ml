(** A partition map is a data structure for a map over a partition of a set.

  The integer intervals are stored in a sorted order, so that finding
  intersections and merging can be accomplised by traversing the lists in the
  same order. *)

open Util

module Interval : sig

  type t = private int

  val hack : int -> t

  val compare : t -> t -> int

  val make : int -> int -> t

  val extend_one : t -> t

  val width : t -> int

  val inside : int -> t -> bool

  val split_if_inside : int -> t -> t * t * t

  val start : t -> int

  val end_ : t -> int

  val to_string : t -> string

  val is_none : t -> bool

  val none : t

  (* [strictly_before t1 t2] is true if t1's end is one less than t2's start. *)
  val strictly_before : t -> t -> bool

  val before_separate : t -> t -> bool

  val merge : t -> t -> t

  val split_inter_diff2 : t -> t ->
                          t * t * t * t * t
  val split_inter_diff3 : t -> t -> t ->
                          t * t * t * t * t * t * t
  val split_inter_diff4 : t -> t -> t -> t ->
                          t * t * t * t * t * t * t * t * t

  val aligned_inter_diff2 : t -> t ->
                            t * t * t
  val aligned_inter_diff3 : t -> t -> t ->
                            t * t * t * t
  val aligned_inter_diff4 : t -> t -> t -> t ->
                            t * t * t * t * t

  val fold : t -> init:'a -> f:('a -> int -> 'a) -> 'a

  val iter : t -> f:(int -> unit) -> unit

end = struct

  (* Start will occupy the higher order bits. This way a comparison can test
  * for equality (the most common case) and determine the Before/After path
  * without looking at end.
  *
  * This means that the new range must start at 1 and 0 is no longer a valid
  * element in our partition. We'll use 0 to encode a None
  *)
  type t = int [@@deriving ord]

  let hack x = x

  let none = 0

  let is_none t =
    t = none

  (* Masks *)
  let lower32 = 4294967295
  let upper32 = -4294967296

  let start_of_int i =
    (i lsl 32)

  let int_of_start  i =
    ((i land upper32) lsr 32)

  let start_p t =
    t land upper32

  let end_of_int i =
    i

  let int_of_end i =
    i

  let end_p t =
    t land lower32

  let make64 s e =
    s lor e

  let start_of_end e =
    start_of_int (int_of_end e)

  let end_of_start s =
    end_of_int (int_of_start s)

  (* Interfacing to ints *)
  let start t =
    (int_of_start (start_p t)) - 1

  let end_ t =
    (int_of_end (end_p t)) - 1

  let width t =
    (end_ t) - (start t) + 1

  let inside i t =
    (start t) <= i && i <= (end_ t)

  let to_string t =
    if is_none t then
      "none"
    else
      sprintf "(%d,%d)" (start t) (end_ t)

  let make s e =
    if e < s then
      invalid_argf "Not in order end: %d < start: %d" e s
    else if e < 0 then
      invalid_argf "End %d less than 0" e
    else if s < 0 then
      invalid_argf "Start %d less than 0" e
    else
      let s64 = start_of_int (s + 1) in
      let e64 = end_of_int (e + 1) in
      make64 s64 e64

  let make_o s e =
    if e < 0 || s < 0 then
      invalid_argf "Start %d or End %d less than 0" s e
    else if e < s then
      none
    else
      let s64 = start_of_int (s + 1) in
      let e64 = end_of_int (e + 1) in
      make64 s64 e64

  let extend_one t =
    succ t

  let merge t1 t2 =
    let s1 = start_p t1 in
    let s2 = start_p t2 in
    let e1 = end_p t1 in
    let e2 = end_p t2 in
    if succ e1 = (end_of_start s2) then
      make64 s1 e2
    else if succ e2 = (end_of_start s1) then
      make64 s2 e1
    else
      none

  let merge3 t1 t2 t3 =
    let p1d = merge t1 t2 in
    if is_none p1d then
      none
    else
      merge p1d t3

  let merge_exn t1 t2 =
    let m = merge t1 t2 in
    if is_none m then
      invalid_argf "not merge-able %s %s"
        (to_string t1) (to_string t2)
    else
      m

  let prep p1 p2 =
    if is_none p1 then
      p2
    else if is_none p2 then
      p1
    else
      make64 (start_p p1) (end_p p2)

(** START HERE. *)
  let split_if_inside i t =
    let s = start t in
    let e = end_ t in
    if s = e && i = s then
        none,   t,                none
    else
      make_o s (i - 1)
      , make_o (max s i) (min i e)
      , make_o (i + 1) e
(*
      if i < s then
        none,   none,             none
      else if i = s then
        none,   make i i,         make (i + 1) e
      else if i < e then
        make s (i - 1), make i i, make (i + 1) e
      else if i = e then
        make s (i - 1), none, 
      else (* i > e *)
        none,         none
        *)

  let strictly_before t1 t2 =
    t1 < t2 && (end_p t1) < (pred (end_of_start (start_p t2)))

  (* This is a very dangerous method, you better be sure that the intervals
  * are separate! *)
  let before_separate t1 t2 =
    t1 < t2

  let split_inter_diff2 t1 t2 =
    if t1 = t2 then
      (* EmptyEmpty *) none, none, t1,  none, none
      (* Only works if s1, s2 < 2^31 otherwise t1, t2 is negative. *)
    else if t1 < t2 then (* s1 <= s2 *)
      begin
        let s1 = start_p t1 in
        let s2 = start_p t2 in
        let e1 = end_p t1 in
        let e2 = end_p t2 in
        if s1 = s2 then (* -> e1 < e2 ----> s1 = s2 <= e1 < e2 *)
          let after = make64 (start_of_end (succ e1)) e2 in
          (* EmptyAfter *)none, none, t1, none, after
        else (* s1 < s2 --------> e1 ? e2  we don't know!*)
          begin
            if e1 = e2 then (* -------> s1 < e1 = s2 = e2 *)
              let before = make64 s1 (pred (end_of_start s2)) in
              (* BeforeEmpty *) before, none, t2,  none, none
            else if e1 < e2 (* ----- s1 <= e1 ? s2 <= e2 *) then
              begin
                if e1 < end_of_start s2 then (* s1 <= e1 < s2 <= e2 *)
                  (* Before *) t1, none, none, none, t2
                else (* e1 >= s2 --------> s1 < s2 <= e1 < e2 *)
                  let before = make64 s1 (pred (end_of_start s2)) in
                  let inter  = make64 s2 e1 in
                  let after  = make64 (start_of_end (succ e1)) e2 in
                  (* BeforeAfter *) before, none, inter,  none, after
              end
            else (* e1 > e2    ----- s1 < s2 <= e2 < e1 *)
              let before = make64 s1 (pred (end_of_start s2)) in
              let after  = make64 (start_of_end (succ e2)) e1 in
              (* SplitEmpty *) before, none, t2,  after,  none
          end
      end
    else (* t1 > t2 ------------ s1 >= s2  --- s2 <= s1 *)
      begin
        let s1 = start_p t1 in
        let s2 = start_p t2 in
        let e1 = end_p t1 in
        let e2 = end_p t2 in
        if s1 = s2 then (* --> e1 > e2 -- e2 < e1  ---- s1 = s2 <= e2 < e1*)
          let after = make64 (start_of_end (succ e2)) e1 in
          (* AfterEmpty *) none, none, t2,  after,  none
        else (* s1 > s2  --- s2 < s1 --> e2 < e1 *)
          begin
            if e1 = e2 then (* -----> s2 < s1 = e1 = e2 *)
              let before = make64 s2 (pred (end_of_start s1)) in
              (* EmptyBefore *) none, before, t1,  none, none
            else if e1 > e2 then (* e2 < e1 ----> s2 <= e2 ? s1 <= e1 *)
              begin
                if e2 < end_of_start s1 then (* s2 <= e2 < s1 <= e2 *)
                  (* After *)  none, t2,     none, t1,     none
                else
                  let before = make64 s2 (pred (end_of_start s1)) in
                  let inter  = make64 s1 e2 in
                  let after  = make64 (start_of_end (succ e2)) e1 in
                  (* AfterBefore *) none, before, inter,  after,  none
              end
            else (* e1 < e2       s2 <= s1 <= e1 < e2 *)
              let before = make64 s2 (pred (end_of_start s1)) in
              let after  = make64 (start_of_end (succ e1)) e2 in
              (* EmptySplit *) none, before, t1,  none, after
          end
      end

  let split_inter_diff3 i1 i2 i3 =
    let b1, b2, i, a1, a2 = split_inter_diff2 i1 i2 in
    if is_none i then
      b1, b2, none, i, a1, a2, i3
    else
      let b12, b3, i123, a12, a3 = split_inter_diff2 i i3 in
      prep b1 b12,
      prep b2 b12,
      b3,
      i123,
      prep a12 a1,
      prep a12 a2,
      a3

  let split_inter_diff4 i1 i2 i3 i4 =
    let b1, b2, b3, i, a1, a2, a3 = split_inter_diff3 i1 i2 i3 in
    if is_none i then
      b1, b2, b3, none, i, a1, a2, a3, i4
    else
      let b123, b4, i1234, a123, a4 = split_inter_diff2 i i4 in
      prep b1 b123,
      prep b2 b123,
      prep b3 b123,
      b4,
      i1234,
      prep a123 a1,
      prep a123 a2,
      prep a123 a3,
      a4

  let aligned_inter_diff2 t1 t2 =
    if t1 = t2 then
      t1, none, none    (* EmptyEmpty *)
    (* Only works if s1, s2 < 2^31 otherwise t1, t2 is negative. *)
    else if t1 < t2 then (* s1 = s2 -> e1 < e2 *)
      let e1 = end_p t1 in
      let e2 = end_p t2 in
      let after = make64 (start_of_end (succ e1)) e2 in
      t1, none, after     (* EmptyAfter *)
    else (* t1 > t2 ----    s1 = s2 -> e1 > e2 *)
      let e1 = end_p t1 in
      let e2 = end_p t2 in
      let after = make64 (start_of_end (succ e2)) e1 in
      t2, after, none     (* AfterEmpty *)

  let aligned_inter_diff3 i1 i2 i3 =
    let i, m1, m2 = aligned_inter_diff2 i1 i2 in
    let ni, ma, m3 = aligned_inter_diff2 i i3 in
    ni, prep ma m1, prep ma m2, m3

  let aligned_inter_diff4 i1 i2 i3 i4 =
    let i, m1, m2, m3 = aligned_inter_diff3 i1 i2 i3 in
    let ni, ma, m4 = aligned_inter_diff2 i i4 in
    ni, prep ma m1, prep ma m2, prep ma m3, m4

  let iter t ~f =
    for i = (start t) to (end_ t) do f i done

  let fold t ~init ~f =
    let acc = ref init in
    for i = (start t) to (end_ t) do acc := f !acc i done;
    !acc

end (* Interval *)

module Set = struct

  type t = Interval.t list

  let of_interval i =                 (* This isn't the cleanest abstraction ... *)
    if Interval.is_none i then
      []
    else
      [i]

  let invariant =
    let open Interval in
    let rec loop = function
      | []  -> true
      | h :: [] -> true
      | h1 :: h2 :: t ->
          if strictly_before h1 h2 then
            loop (h2 :: t)
          else
            false
    in
    loop

  let to_string  =
    string_of_list ~sep:";" ~f:Interval.to_string

  let size =
    List.fold_left ~init:0 ~f:(fun a i -> a + Interval.width i)

  let length =
    List.length

  let inside i l =
    List.exists l ~f:(Interval.inside i)

  let universal = function
    | [i] -> true
    | _   -> false

  let compare s1 s2 =
    match s1, s2 with
    | i1 :: _ , i2 :: _ -> Interval.compare i1 i2
    | _                 -> assert false

  (* Should I make the sets a real pair? *)
  let first_pos = function
    | []      -> assert false
    | s :: _  -> Interval.start s

  let cons_if_nnone o l =
    if Interval.is_none o then l else o :: l

  (* The elements in a set are ordered.  *)
  let split_if_in ii l =
    let open Interval in
    let rec loop = function
      | []      -> None, []
      | h :: t  ->
          let b1, b2, i, a1, a2 = split_inter_diff2 ii h in
          if Interval.is_none i then begin
            if h = a1 then (* After *)
              let o, nt = loop t in
              o, h :: nt
            else (* Before *)
              None, l
          end else
            Some (of_interval b2), (cons_if_nnone a2 t)
    in
    loop l

  (* Zip together two non-intersecting (separate) sets. *)
  let merge_separate =
    let open Interval in
    let rec start l1 l2 = match l1, l2 with
      | _,  []            -> l1
      | [], _             -> l2
      | h1 :: t1
      , h2 :: t2 ->
          if before_separate h1 h2 then
            loop h1 t1 l2
          else
            loop h2 l1 t2
    and loop ps l1 l2 = match l1, l2 with
      | [],       []        ->  [ps]
      | h1 :: t1, []        ->  let m1 = merge ps h1 in
                                if is_none m1 then
                                   ps :: l1
                                else
                                   loop m1 t1 []
      | [],       h2 :: t2  ->  let m2 = merge ps h2 in
                                if is_none m2 then
                                  ps :: l2
                                else
                                  loop m2 [] t2
      | h1 :: t1, h2 :: t2  ->
          if before_separate h1 h2 then begin
            let m1 = merge ps h1 in
            if is_none m1 then
              ps :: loop h1 t1 l2
            else
              loop m1 t1 l2
          end else begin
            let m2 = merge ps h2 in
            if is_none m2 then
              ps :: loop h2 l1 t2
            else
              loop m2 l1 t2
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
        let b1, b2, inter, a1, a2 = split_inter_diff2 h1 h2 in
        let i, r1, r2 = loop (cons_if_nnone a1 t1) (cons_if_nnone a2 t2) in
        (cons_if_nnone inter i) , (cons_if_nnone b1 r1) , (cons_if_nnone b2 r2)
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
      let inter, m1, m2 = aligned_inter_diff2 h1 h2 in
      let i, r1, r2 = all_intersections (cons_if_nnone m1 t1) (cons_if_nnone m2 t2) in
      (inter :: i), r1, r2

  let all_intersections3 =
    let rec loop l1 l2 l3 = match l1, l2, l3 with
      | [],  _,  _
      |  _, [],  _
      |  _,  _, []  -> [], l1, l2, l3
      | h1 :: t1
      , h2 :: t2
      , h3 :: t3    ->
        let b1, b2, b3, i, a1, a2, a3 = Interval.split_inter_diff3 h1 h2 h3 in
        let nt1 = cons_if_nnone a1 t1 in
        let nt2 = cons_if_nnone a2 t2 in
        let nt3 = cons_if_nnone a3 t3 in
        let il, r1, r2, r3 = loop nt1 nt2 nt3 in
        cons_if_nnone i il
        , cons_if_nnone b1 r1
        , cons_if_nnone b2 r2
        , cons_if_nnone b3 r3
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
        let nt1 = cons_if_nnone ho1 t1 in
        let nt2 = cons_if_nnone ho2 t2 in
        let nt3 = cons_if_nnone ho3 t3 in
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
        let nt1 = cons_if_nnone a1 t1 in
        let nt2 = cons_if_nnone a2 t2 in
        let nt3 = cons_if_nnone a3 t3 in
        let nt4 = cons_if_nnone a4 t4 in
        let il, r1, r2, r3, r4 = loop nt1 nt2 nt3 nt4 in
        cons_if_nnone i il
        , cons_if_nnone b1 r1
        , cons_if_nnone b2 r2
        , cons_if_nnone b3 r3
        , cons_if_nnone b4 r4
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
        let nt1 = cons_if_nnone ho1 t1 in
        let nt2 = cons_if_nnone ho2 t2 in
        let nt3 = cons_if_nnone ho3 t3 in
        let nt4 = cons_if_nnone ho4 t4 in
        let il, r1, r2, r3, r4 = all_intersections4 nt1 nt2 nt3 nt4 in
        inter :: il, r1, r2, r3, r4

  let fold_intervals t ~init ~f =
    List.fold_left t ~init ~f

  let fold t ~init ~f =
    List.fold_left t ~init ~f:(fun init interval ->
      Interval.fold interval ~init ~f)

  let iter t ~f =
    fold t ~init:() ~f:(fun () i -> f i)

end (* Set *)

(* Fix the size of the arrays used to store partition map arrays.
 * The intuition is that (1) it won't fragment the available memory if we have
 * fixed array sizes and (2) having excess capacity in arrays will avoid
 * these reallocations as we can fill in the empty space. *)
let sizes = [ 2; 4; 8; 16; 32; 64; 128; 256; 512; 1024 ]

let closest_arrays_size n =
  match List.find sizes ~f:(fun s -> n < s) with
  | Some s  -> s
  | None    -> eprintf "asked for a size greater than 1024: %d\n" n; n

let ia_transition_table = function
  | 2   -> 4
  | 4   -> 8
  | 8   -> 16
  | 16  -> 32
  | 32  -> 64
  | 64  -> 128
  | 128 -> 256
  | 256 -> 512
  | 512 -> 1024
  | n   -> eprintf "Odd ia transition: %d\n" n;
           truncate ((float n) *. 1.1)

(*
module type Equaled_type = sig

  type t

  val zero : t

  val equal : t -> t -> bool

  val to_string : t -> string

end
*)

module Pmarr = struct
  
  type 'a t =
    { interval_array  : Interval.t array
    ; end_boundaries  : int array
    (* Indices into interval_array when a given set stops. *)
    ; value_array     : 'a array
    (* What we're storing in sequential order. *)
    }

  (* recursive for loop. *)
  let from ~start ~end_ ~f ~init =
    let rec loop a i =
      if i > end_ then
        a
      else
        loop (f a i) (i + 1)
    in
    loop init start

  let fold_per_set t ~f ~init ~f_iv ~init_iv =
    let end_ = Array.length t.end_boundaries in
    let rec loop set_idx a start =
      if set_idx >= end_ then
        a
      else
        let end_ = t.end_boundaries.(set_idx) in
        let fs = from ~start ~end_ ~init:init_iv
                  ~f:(fun a i -> f_iv a t.interval_array.(i))
        in
        let na = f a fs t.value_array.(set_idx) in
        loop (set_idx + 1) na (end_ + 1)
    in
    loop 0 init 0

  let fold_per_interval t ~f ~init =
    let end_ = Array.length t.end_boundaries in
    let rec loop set_idx acc start =
      if set_idx >= end_ then
        acc
      else
        let end_ = t.end_boundaries.(set_idx) in
        let valu = t.value_array.(set_idx) in
        let nacc = from ~start ~end_ ~init:acc
                    ~f:(fun a i -> f a t.interval_array.(i) valu)
        in
        loop (set_idx + 1) nacc (end_ + 1)
    in
    loop 0 init 0

  (* 1. Within a set each interval is strictly before the next.
   * 2. Between subsequent sets, the first (sets are non-empty) interval
   *    is also before (not necessarily strictly) the next.
   *)
  let ascending_invariant et_to_string et_equal t =
    let module M = struct exception No end in
    let check_within value first previous current = 
      if Interval.strictly_before previous current then
        `PrevSet (current, first, value)
      else begin
        eprintf "within %s prev interval: %s not before: %s"
          (et_to_string value)
          (Interval.to_string previous)
          (Interval.to_string current);
        raise M.No
      end
    in
    let check_between_sets pv nv previous_first current_first = 
      if Interval.before_separate previous_first current_first then
        `NewSet (current_first, nv)
      else begin
        eprintf "from %s to %s prev first interval: %s not before: %s"
          (et_to_string pv)
          (et_to_string nv)
          (Interval.to_string previous_first)
          (Interval.to_string current_first);
        raise M.No
      end
    in 
    try
      fold_per_interval t ~init:`Start ~f:(fun s i nv ->
        match s with
        | `Start                ->
            `NewSet (i, nv)
        | `NewSet (fi, pv)    ->
            if et_equal pv nv then
              check_within nv fi fi i
            else
              check_between_sets pv nv fi i
        | `PrevSet (pi, fi, pv) ->
            if et_equal pv nv then
              check_within nv fi pi i 
            else
              check_between_sets pv nv fi i)
      |> function
          | `Start      (* Should I fail in this 'empty' case? *)
          | `NewSet _
          | `PrevSet  _ -> true
    with M.No           -> false

  let fold t ~f ~init =
    fold_per_interval t ~init ~f:(fun init interval value ->
      Interval.fold interval ~init ~f:(fun a i -> f a i value))

  let set_to_string ilst =
    string_of_list ilst ~sep:";" ~f:Interval.to_string

  let to_string t et_to_string =
    fold_per_set t
      (* Accumulate intervals into a list . *)
      ~init_iv:[] ~f_iv:(fun l i -> i :: l)
      (* Convert that list to a string and pair it with value. *) 
      ~init:[] ~f:(fun acc ilst v ->
        (sprintf "[%s]:%s" (set_to_string (List.rev ilst)) (et_to_string v)) :: acc)
    |> List.rev
    |> string_of_list ~sep:"; " ~f:id

  (* Width of the partition space. *)
  let size t =
    fold_per_interval t ~init:0 ~f:(fun s i _v -> s + Interval.width i)

  (* Number of unique values. *)
  let length t et_equal =
    let fs =
      fold_per_interval t ~init:`St ~f:(fun s _i v ->
        match s with
        | `St         -> `Pv (v, 0)
        | `Pv (pv, c) -> if et_equal v pv then
                           s
                         else
                           `Pv (v, c + 1))
    in
    match fs with
    | `St        -> 0
    | `Pv (_, c) -> c

  let of_ascending_set_list lst =
    let ilst = List.map ~f:fst lst in
    let vlst = List.map ~f:snd lst in
    let value_array = Array.of_list vlst in
    let number_of_intervals = List.fold_left ilst ~init:0 ~f:(fun i s -> i + Set.length s) in
    let intervals_size = closest_arrays_size number_of_intervals in
    let interval_array = Array.make intervals_size Interval.none in
    let end_boundaries = Array.make (Array.length value_array) 0 in
    let _last_pos_pair =
      List.fold_left ilst ~init:(0, 0) ~f:(fun (p, set_idx) s ->
        let np = Set.fold_intervals s ~init:p ~f:(fun p s ->
                  interval_array.(p) <- s;
                  p + 1)
        in
        end_boundaries.(set_idx) <- np - 1;
        np, set_idx + 1)
    in
    { interval_array
    ; end_boundaries
    ; value_array
    }

  let get (type a) (t : a t) i =
    let module M = struct exception Found of a end in
    try
      fold_per_interval t ~init:() ~f:(fun () interval v ->
        if Interval.inside i interval then
          raise (M.Found v)
        else
          ());
      raise Not_found
    with M.Found a ->
      a

  (* Modifiers and new constructors. *)
  let start_idx end_boundaries set_idx =
    if set_idx = 0 then
      0
    else
      end_boundaries.(set_idx - 1) + 1

  let for_each_index_in_set end_boundaries set_idx ~f ~init =
    let start = start_idx end_boundaries set_idx in
    let end_  = end_boundaries.(set_idx) in
    from ~start ~end_ ~f ~init

  let shift_ia interval_array starting_at =
    let len =  Array.length interval_array in
    let nia =
      if Interval.is_none (interval_array.(len - 1)) then   (* Safe to shift *)
        interval_array
      else begin
        let new_size = ia_transition_table len in
        let nia = Array.make new_size Interval.none in
        Array.blit ~src:interval_array ~src_pos:0
                   ~dst:nia ~dst_pos:0 ~len:(starting_at - 1);
        nia
      end
    in
    Array.blit ~src:interval_array ~src_pos:starting_at
               ~dst:nia ~dst_pos:(starting_at + 1)
               ~len:(len - starting_at - 1);
    nia

  let shift_boundaries end_boundaries set_idx =
    for j = set_idx to Array.length end_boundaries - 1 do
      end_boundaries.(j) <- succ end_boundaries.(j)
    done

  let shift_value_array value_array set_idx =
    let end_ = Array.length value_array in
    let rec loop i v =
      if i > end_ then
        ()
      else begin
        let nv = value_array.(i) in
        value_array.(i) <- v;
        loop (i + 1) nv
      end
    in
    loop (set_idx + 1) value_array.(set_idx)
 
  let one_bigger zero a =
    let len  = Array.length a in
    let copy = Array.make (Array.length a + 1) zero in
    Array.blit ~src:a ~src_pos:0 ~dst:copy ~dst_pos:0 ~len;
    copy

  let find_and_remove_from_singleton_set ~start ~end_ singleton interval_array =
    let si = Interval.make singleton singleton in
    if start = end_ then  
      let sii_b, sii_a = split_if_inside singleton interval_array.(start) in

      if interval_array.(start) = interval then
        true              (* found it *)
        , true            (* need to shrink *)
        , Interval.none   (* nothing to add. *)
      else
        let 


      if i = 
      if 
    else
    from ~start ~end_ ~init:(false, Interval.none) =


  (* nv doesn't equal current value:
   *  - inteval_array will be different,
   *  - end_boundaries will have different values and might be bigger!
   *  - value_array might be bigger. *)
  let set_different (type a) (t : a t) (et_equal : a -> a -> bool) i nv =
    let value_length = Array.length t.value_array in
    let insert add_me interval_array end_boundaries =
      let same_value_index_opt = Array.findi t.value_array ~f:(et_equal nv) in
      match same_value_index_opt with
      | None    ->                                                (* Add to ? *)
          let end_boundaries = one_bigger 0 end_boundaries in
          let value_array    = one_bigger nv t.value_array in
          let last_va        = Array.length value_array - 1 in
          let rec over_all_values2 set_idx =
            if set_idx = last_va then begin                     (* Last one *)
              let j = end_boundaries.(set_idx - 1) + 1 in
              interval_array.(j) <- add_me;
              end_boundaries.(set_idx) <- j;
              value_array.(last_va) <- nv;
              { interval_array ; end_boundaries ; value_array }
            end else begin
              let first_index = start_idx end_boundaries set_idx in
              let first_interval = interval_array.(first_index) in
              if Interval.before_separate add_me first_interval then begin
                let interval_array = shift_ia interval_array first_index in
                shift_boundaries end_boundaries set_idx;
                shift_value_array value_array set_idx;
                interval_array.(first_index) <- add_me;
                { interval_array ; end_boundaries ; value_array }
              end else
                over_all_values2 (set_idx + 1)
            end
          in
          over_all_values2 0
      | Some set_idx ->                                         (* Add here *)
          let remaining_interval =
            for_each_index_in_set end_boundaries set_idx
              ~init:add_me
              ~f:(fun add_me j ->
                    if Interval.is_none add_me then                   (* Done *)
                      add_me
                    else
                      let interval = interval_array.(j) in
                      if Interval.strictly_before add_me interval then begin
                        interval_array.(j) <- add_me;
                        interval
                      end else
                        let mm = Interval.merge add_me interval in
                        if Interval.is_none mm then         (* Strictly after *)
                          add_me
                        else begin
                          interval_array.(j) <- mm;
                          Interval.none                             (* Fin. *)
                        end)
          in
          if Interval.is_none remaining_interval then
            { interval_array
            ; end_boundaries
            ; value_array = t.value_array
            }
          else begin
            let interval_array = shift_ia interval_array end_boundaries.(set_idx) in
            interval_array.(end_boundaries.(set_idx)) <- remaining_interval;
            shift_boundaries end_boundaries set_idx;
            { interval_array
            ; end_boundaries
            ; value_array = t.value_array
            }
          end
    in
    let rec over_all_values interval_array end_boundaries set_idx =
      if set_idx >= value_length then
        raise Not_found
      else
        let add_me = Interval.make i i in

        let found, remaining_interval =
          for_each_index_in_set end_boundaries set_idx
            ~init:(false, Interval.none)
            ~f:(fun (found, prev_remains) j ->
                  let interval = interval_array.(j) in
                  if found then begin
                    if Interval.is_none prev_remains then
                      (found, prev_remains)
                    else begin                            (* shift to the end *)
                      interval_array.(j) <- prev_remains;
                      found, interval
                    end
                  end else begin         (* not found -> prev_remains is none *)
                    let b, a = Interval.split_if_inside i interval in
                    if Interval.is_none b then
                      if Interval.is_none a then            (* Didn't find it *)
                        false, a
                      else begin                 (* It matched interval start *)
                        interval_array.(j) <- a;
                        true, Interval.none
                      end
                    else begin 
                      interval_array.(j) <- b;
                      true, a
                    end
                  end)
        in
        let () = printf "at set_idx: %d found: %b remaining_interval: %s\n%!"
                  set_idx found (Interval.to_string remaining_interval)
        in
        if not found then
          over_all_values interval_array end_boundaries (set_idx + 1)
        else if Interval.is_none remaining_interval then begin   (* perfect fit *)
          value_array.(

            { interval_array
            ; end_boundaries
            ; value_array = t.value_array
            }
 
          insert add_me interval_array end_boundaries
        else begin (* shift and then add new set *)
          let interval_array = shift_ia interval_array end_boundaries.(set_idx) in
          interval_array.(end_boundaries.(set_idx)) <- remaining_interval;
          shift_boundaries end_boundaries set_idx;
          let add_me = Interval.make i i in
          insert add_me interval_array end_boundaries 
        end
    in
    over_all_values (Array.copy t.interval_array) t.end_boundaries 0 

  let set t et_equal i nv =
    let current_value = get t i in
    if et_equal current_value nv then
      t
    else
      set_different t et_equal i nv

      
  let map t ~f zero et_equal =
    let value_length = Array.length t.value_array in
    let new_value_array = Array.make value_length zero in
    let finish = function
      | None    ->        (* Didn't create a new copy, great! *)
          { interval_array = t.interval_array
          ; end_boundaries = t.end_boundaries
          ; value_array    = new_value_array
          }
      | Some (interval_array, end_boundaries) ->
          let eb_length = Array.length end_boundaries in
          let value_array =
            if eb_length < value_length then            (* Fewer states! *)
              Array.sub new_value_array ~pos:0 ~len:eb_length
            else
              new_value_array
          in
          { interval_array
          ; end_boundaries
          ; value_array
          }
    in
    let merge_or_extend dest_i source_i ia_eb_opt =         (* dest_i < source_i *)
      let interval_array, end_boundaries =
        match ia_eb_opt with
        | None    ->
            Array.copy t.interval_array
            , Array.copy t.end_boundaries
        | Some (interval_array, end_boundaries) ->
            interval_array
            , end_boundaries
      in
      let interval_array =
        for_each_index_in_set t.end_boundaries source_i
          ~init:interval_array ~f:(fun interval_array j ->
            let source_interval = t.interval_array.(j) in
            let merged, remaining_interval =
              for_each_index_in_set end_boundaries dest_i ~init:(false, Interval.none)
                ~f:(fun (merged, prev_remains) i ->
                      let dest_interval = interval_array.(i) in
                      if merged then begin
                        if Interval.is_none prev_remains then
                          (merged, prev_remains)
                        else begin
                          interval_array.(i) <- prev_remains;
                          merged, dest_interval
                        end
                      end else begin      (* not merged *)
                        if Interval.strictly_before source_interval dest_interval then begin
                          interval_array.(i) <- source_interval;
                          true, dest_interval
                        end else
                          let mm = Interval.merge source_interval dest_interval in
                          if Interval.is_none mm then
                            false, Interval.none
                          else begin
                            interval_array.(i) <- mm;
                            true, Interval.none
                          end
                      end)
            in
            if merged then begin
              if Interval.is_none remaining_interval then
                interval_array
              else begin
                let interval_array = shift_ia interval_array end_boundaries.(dest_i) in
                interval_array.(end_boundaries.(dest_i)) <- remaining_interval;
                shift_boundaries end_boundaries dest_i;
                interval_array
              end
            end else (* not merged*) begin   (* Comes strictly after last interval. *)
             let interval_array = shift_ia interval_array end_boundaries.(dest_i) in
             interval_array.(end_boundaries.(dest_i)) <- source_interval;
             shift_boundaries end_boundaries dest_i;
             interval_array
            end)
      in
      Some (interval_array, end_boundaries)
    in
    let find_place nv stop_index value_array =
      let rec loop i =
        if i = stop_index then
          None
        else if et_equal value_array.(i) nv then
          Some i  
        else
          loop (i + 1)
      in
      loop 0
    in
    let rec over_all_values set_idx ia_eb_opt =
      if set_idx >= value_length then
        finish ia_eb_opt
      else
        let nv = f t.value_array.(set_idx) in
        match find_place nv set_idx new_value_array with
        | None  ->                                           (* new index, or no merge. *)
            new_value_array.(set_idx) <- nv;
            over_all_values (set_idx + 1) ia_eb_opt
        | Some dest_index ->                                       (* old index, merge. *)
            let new_ia_ab_opt = merge_or_extend dest_index set_idx ia_eb_opt in
            over_all_values (set_idx + 1) new_ia_ab_opt
    in
    over_all_values 0 None


end (* Pmarr *)

(* Things start out in descending order when we construct the partition, but
 * when we 'reverse' it they are constructed into an ascending order that is
 * better for merging. 
 *)

type ascending  (*= Ascending*)
type descending (*= Descending*)

type 'a asc =
  | E                       (* empty, merging against this will fail! *)
  | U of Interval.t * 'a           (* universal, hold size - 1 = end_ *)
  | S of 'a Pmarr.t

type 'a desc = (Interval.t * 'a) list

type (_, 'a) t =
  | Asc   : 'a asc -> (ascending, 'a) t
  | Desc  : 'a desc -> (descending, 'a) t

let invariant
  : type o a. (a -> string) -> (a -> a -> bool) -> (o, a) t -> bool =
  fun et_to_string et_equal t -> match t with
    | Asc E      -> true
    | Asc U _    -> true
    | Asc (S pa) -> Pmarr.ascending_invariant et_to_string et_equal pa
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
  Asc (U (i, v))

(* Properties *)
(*let asc_to_string la to_s =
  string_of_list la ~sep:"; " ~f:(fun (s, v) ->
      sprintf "[%s]:%s" (Set.to_string s) (to_s v)) *)

let desc_to_string ld to_s =
  string_of_list ld ~sep:";" ~f:(fun (i, v) ->
    sprintf "%s:%s" (Interval.to_string i) (to_s v))

let to_string: type o a. (o, a) t -> (a -> string) -> string =
  fun t to_s -> match t with
    | Asc E         -> "Empty!"
    | Asc (U (i,v)) -> sprintf "%s:%s" (Interval.to_string i) (to_s v)
    | Asc (S pa)    -> Pmarr.to_string pa to_s
    | Desc ld       -> desc_to_string ld to_s

let size_d = function
  | []            -> 0
  | ((i, _) :: _) -> Interval.end_ i + 1

let size : type o a. (o, a) t -> int = function
  | Asc E          -> 0
  | Asc (U (i, _)) -> Interval.width i
  | Asc (S pa)     -> Pmarr.size pa
  | Desc l         -> size_d l

let length : type o a. (a -> a -> bool) -> (o, a) t -> int =
  fun et_equal t -> match t with
    | Asc E       -> 0
    | Asc (U _)   -> 1
    | Asc (S pa)  -> Pmarr.length pa et_equal 
    | Desc l      -> List.length l

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

(*
let map_with_full_check eq l ~f =
  List.fold_left l ~init:[] ~f:(fun acc (s, v) ->
      merge_or_add_to_end eq s (f v) acc)
*)

let ascending_t eq l =
  List.fold_left l ~init:[] ~f:(fun acc (i, v) ->
    merge_or_add_to_end eq (Set.of_interval i) v acc)
  |> List.sort ~cmp:(fun (s1, _) (s2, _) -> Set.compare s1 s2)

(*let asc_sets_to_str s =
  asc_to_string s (fun _ -> "")
  *)

let ascending eq = function
  | Desc l ->
    let a = ascending_t eq l in                             (* assert (ascending_invariant l); *)
    match a with
    | []      -> invalid_arg "Empty descending!"
    | [s,v]   -> if Set.universal s then
                   Asc (U (List.hd_exn s, v))
                 else
                  invalid_argf "Single set but not universal? %s" (Set.to_string s)
    | lst     -> Asc (S (Pmarr.of_ascending_set_list lst))

let descending = function
  | Asc E          -> invalid_arg "Can't convert empty to descending"
  | Asc (U (i, v)) -> Desc [i, v]
  | Asc (S pm)     ->
      Pmarr.fold_per_interval pm ~init:[] ~f:(fun acc interval value ->
        (interval, value) :: acc)
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
  | Asc (S pm)     -> Pmarr.get pm i

let set t et_equal i nv = match t with
  | Asc E           ->
      invalid_arg "Can't set from empty"
  | Asc (U (s, pv)) ->
      if et_equal pv nv then
        t
      else
        let b, a = Interval.split_if_inside i s in
        let asc_lst =
          if Interval.is_none b then
            if Interval.is_none a then
              invalid_argf "Outside %d range %s" i (Interval.to_string s)
            else
              [ [Interval.make i i], nv; [a], pv]
          else
            if Interval.is_none a then
              [ [b], pv; [Interval.make i i], nv]
            else
              [ [b; a], pv; [Interval.make i i], nv]
        in
        Asc (S (Pmarr.of_ascending_set_list asc_lst))
  | Asc (S pm)      ->
      Asc (S (Pmarr.set pm et_equal i nv))

let map t ~f zero et_equal  = match t with
  | Asc E          -> Asc E
  | Asc (U (s, v)) -> Asc (U (s, f v))
  | Asc (S pm)     -> Asc (S (Pmarr.map pm ~f zero et_equal))

let v1 = init_all_a ~size:100 "A"

let to_s v = to_string v id

let () =
  printf "v1 = %s\n" (to_s v1)

let v2 =
  set v1 (=) 20 "B"

let () =
  printf "v2 = %s\n" (to_s v2)

let v2_no_op =
  set v2 (=) 20 "B"

let () =
  printf "v2_no_op = %s\n" (to_s v2_no_op)


let v3 =
  set v2 (=) 20 "C"

let () =
  printf "v3 = %s\n" (to_s v3)
(*
let v4 =
  set v2 (=) 40 "D"

let () =
  printf "v4 = %s\n" (to_s v4)
*)



(*
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
        *)


