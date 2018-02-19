(** Test partition map *)
(*BISECT-IGNORE-BEGIN*)

open QCheck

let is_invalid f =
  try
    ignore (f ());
    false
  with (Invalid_argument _) ->
    true

let tests =
  ref []

let add_test ~count ~name a p =
  tests := (Test.make ~count ~name a p) :: !tests

let strictly_neg_int =
  make Gen.(map (fun i -> i - 1) neg_int)

let () = add_test ~count:20
  ~name:"Interval construction fail on negative values."
  (pair strictly_neg_int strictly_neg_int)
  (fun (start, end_) -> is_invalid (fun () -> Interval.make ~start ~end_))

let larger_than_max_value = int_range (Interval.max_value + 1) max_int
let () = add_test ~count:20
  ~name:"Interval construction fail on values too large."
  (pair larger_than_max_value larger_than_max_value)
  (fun (start, end_) -> is_invalid (fun () -> Interval.make ~start ~end_))

let valid_range =
  int_range 0 Interval.max_value

let () = add_test ~count:20
  ~name:"Interval construction fail on out of order values."
  (pair valid_range valid_range)
  (fun (start, end_) ->
     assume (start > end_);
     is_invalid (fun () -> Interval.make ~start ~end_))

(* Compress the range so that we can actually have a non-zero
 * chance of intersects *)

let max_range = 1000
let max_value = max_range - 1

let valid_interval_pair_gen =
  let start = Gen.int_range 0 max_value in
  Gen.(start >>= fun s ->
        map (fun e -> s, e) (int_range s max_value))

let valid_interval_pair =
  make ~print:(fun (s, e) -> sprintf "(%d,%d)" s e)
    valid_interval_pair_gen

let () = add_test ~count:100
  ~name:"Interval construction for non negative works and not none."
  valid_interval_pair
  (fun (start, end_) ->
     let i = Interval.make ~start ~end_ in
     not (Interval.is_none i))

let () = add_test ~count:100
  ~name:"Interval start property is correct."
  valid_interval_pair
  (fun (start, end_) ->
     let i = Interval.make ~start ~end_ in
     start = Interval.start i)

let () = add_test ~count:100
  ~name:"Interval end property is correct."
  valid_interval_pair
  (fun (start, end_) ->
     let i = Interval.make ~start ~end_ in
      end_ = Interval.end_ i)

let valid_interval_gen =
  Gen.(map (fun (start, end_) -> Interval.make ~start ~end_)
    valid_interval_pair_gen)

let valid_interval =
  make ~print:Interval.to_string valid_interval_gen

let () = add_test ~count:100
  ~name:"Interval width property is correct."
  valid_interval
  (fun i -> Interval.(width i = end_ i - start i + 1))

let () = add_test ~count:100
  ~name:"Interval average is inside."
  valid_interval
  (fun i -> Interval.(inside ((start i + end_ i) / 2) i))

let () = add_test ~count:100
  ~name:"Interval extending just changes the end."
  valid_interval
  (fun i ->
     let open Interval in
     let ip1 = extend_one i in
     start i = start ip1 && end_ i + 1 = end_ ip1)

let () = add_test ~count:100
  ~name:"Interval merging works."
  (pair valid_interval valid_interval)
  (fun (i1, i2) ->
     let open Interval in
     let m = merge i1 i2 in
     if end_ i1 + 1 = start i2 then
       start m = start i1 && end_ m = end_ i2
     else if end_ i2 + 1 = start i1 then
       start m = start i2 && end_ m = end_ i1
     else
       is_none m)

let k = Triangular.Indices.full_upper
let ki = Triangular.Indices.full_upper_inverse

let generate_all_pairs n =
  List.init n ~f:(fun i ->
      List.init n ~f:(fun j ->
          i, j, k n i j))
  |> List.concat
  |> List.filter ~f:(fun (i, j,_) -> (i <= j))

let manually n i1 i2 =
  generate_all_pairs n
  |> List.filter_map ~f:(fun (i,j,_k) ->
      if Interval.inside i i1 && Interval.inside j i2 then
        Some (i,j)
      else None)

let expand_ranges l =
  List.map l ~f:(Interval.to_cross_indices max_range)
  |> List.concat

let in_order = function
  | []       -> true
  | p1 :: ps -> List.fold_left ps ~init:(true, p1)
                  ~f:(fun (io, prev) i ->
                      let nio = io && Interval.strict_before prev i in
                      nio, i)
                |> fst

let () =
  add_test ~count:100
    ~name:"Intervals can recover cross indices."
    valid_interval
    (fun i ->
      let cpl = Interval.cpair max_range i i in
      let iplst = expand_ranges cpl in
      let iplst2 = manually max_range i i in
      iplst = iplst2)

let valid_interval_after rs after =
  let bound = max_range - after in
  let start = after + Random.State.int rs bound in
  if start = max_value then
    Interval.make ~start ~end_:start
  else
    let bound = max_range - start - 1 in
    let end_ = start + 1 + Random.State.int rs bound in
    Interval.make ~start ~end_

let valid_set_gen rs =
  let i = Gen.generate1 ~rand:rs valid_interval_gen in
  let rec loop acc =
    if Random.State.bool rs then
      List.rev acc
    else
      let e = Interval.end_ (List.hd_exn acc) in
      (* We can't add last element, it would be merge-able. *)
      if e >= max_value - 1 then
        List.rev acc
      else
        let ni = valid_interval_after rs (e + 2) in
        loop (ni :: acc)
  in
  loop [i]

let valid_set =
  make ~print:Set.to_string valid_set_gen

let all_sequential_pairs ilst =
  let rec loop l = match l with
    | []     -> []
    | h :: t -> List.map l ~f:(fun i -> h, i) @ loop t
  in
  loop ilst

let manually_lst n iplst =
  let inside_a_pair i j =
    List.exists iplst ~f:(fun (i1, i2) ->
        Interval.inside i i1 && Interval.inside j i2)
  in
  generate_all_pairs n
  |> List.filter_map ~f:(fun (i,j,_k) ->
      if inside_a_pair i j then
        Some (i,j)
      else None)

let () =
  add_test ~count:50
    ~name:"Cross pairing sets with itself round trip."
    valid_set
    (fun s ->
      let plst = Set.cpair max_range s in
      let t = expand_ranges plst in
      let iplst = all_sequential_pairs s in
      let m = manually_lst max_range iplst in
      in_order plst && m = t)

let complement_set set =
  let rec loop p acc = function
    | []     ->
      let nacc =
        if p > max_value then
          acc
        else
          (Interval.make ~start:p ~end_:max_value) :: acc
      in
      List.rev nacc
    | h :: t ->
        let s = Interval.start h in
        let np = Interval.end_ h + 1 in
        if p < s then
          let nacc = Interval.make ~start:p ~end_:(s - 1) :: acc in
          loop np nacc t
        else
          loop np acc t
  in
  loop 0 [] set

let non_empty_subset_gen lst rs =
  match lst with
  | []      -> invalid_arg "subsets of non empty sets only!"
  | o :: [] -> None
  | lst     -> Some (List.filter lst ~f:(fun _ ->
                  if Random.bool rs then true else false))

let paired =
  Gen.(valid_set_gen >>= fun s ->
        let c = complement_set s in
        if Set.compare s c < 0 then
          return (s, c)
        else
          return (c, s))

(* Return a set and it's complement, by design the first set will have the
 * lower element, so that when we compute cross pairs it should be the first
 * argument. *)
let non_intersecting_pair_opt =
  let print (s, c) =
    sprintf "%s - %s" (Set.to_string s) (Set.to_string c)
  in
  make ~print paired

let ordered_sequential_pairs lst1 lst2 =
  let rec loop acc l1 l2 = match l1, l2 with
    | [],       _
    | _,        []       -> List.rev acc
    | h1 :: t1, h2 :: t2 ->
      if Interval.before_separate h1 h2 then
        let nacc = List.rev_map ~f:(fun i2 -> h1, i2) l2 @ acc in
        loop nacc t1 l2
      else
        let nacc = List.rev_map ~f:(fun i1 -> h2, i1) l1 @ acc in
        loop nacc l1 t2
  in
  loop [] lst1 lst2

let () =
  add_test ~count:50
    ~name:"Cross pairing separate sets round trip."
    non_intersecting_pair_opt
    (fun (s, c) ->
      let plst = Set.cpair_separate max_range s c in
      let t = expand_ranges plst in
      let orp = ordered_sequential_pairs s c in
      let m = manually_lst max_range orp in
      in_order plst && m = t)

let () =
  let u = Triangular.number max_range - 1 in
  let t = Set.of_interval (Interval.make 0 u) in
  add_test ~count:10
    ~name:"Cross pairing separate and against self should yield universal."
    non_intersecting_pair_opt
    (fun (s, c) ->
       let sm = Set.cpair max_range s in
       let cm = Set.cpair max_range c in
       let sc_m = Set.cpair_separate max_range s c in
       let mgd =
         (Set.merge_separate (Set.merge_separate sm cm)
            sc_m)
       in
       mgd = t)

let without l n =
  List.filteri l ~f:(fun i _ -> i <> n)

let reorder s1 s2 =
  if Set.compare s1 s2 <= 0 then
    s1, s2
  else
    s2, s1

let subsampled_pair_gen =
  Gen.(paired >>= fun (s, c) ->
        let l_s = List.length s in
        let l_c = List.length c in
        let params = triple bool (int_bound l_s) (int_bound l_c) in
        map (fun (first, i_s, i_c) ->
              if first && l_s > 1 then
                reorder (without s i_s) c
              else if l_c > 1 then
                reorder s (without c i_c)
              else
                (* We'd have to change the generating code to make sure our
                 * elements have more than 1 set. *)
                (s, c))
          params)

let subsampled_non_intersecting_pairs =
  let print (s, c) =
    sprintf "%s - %s" (Set.to_string s) (Set.to_string c)
  in
  make ~print subsampled_pair_gen

(* Similar to "Cross pairing separate sets round trip." but now we make sure
 * that the two sets do not partition the entire universe, ie. there are
 * elements missing. *)
let () =
  add_test ~count:50
    ~name:"Cross pairing separate (not full) sets round trip."
    subsampled_non_intersecting_pairs
    (fun (s, c) ->
      let plst = Set.cpair_separate max_range s c in
      let t = expand_ranges plst in
      let orp = ordered_sequential_pairs s c in
      let m = manually_lst max_range orp in
      in_order plst && m = t)

let () =
  QCheck_runner.run_tests_main (List.rev !tests)

(*BISECT-IGNORE-END*)
