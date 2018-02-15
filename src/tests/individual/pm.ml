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
     else
       is_none m)

let k = Triangular_array.full_upper_triangular_index
let ki = Triangular_array.full_upper_triangular_inverse
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

let expand_ranges_back_to_pairs n plst =
  List.map plst ~f:(fun i ->
      let s = Interval.start i in
      let w = Interval.width i in
      List.init w ~f:(fun j -> ki n (j + s)))
  |> List.concat

let () =
  add_test ~count:50
  ~name:"Pairing round trip."
  (pair valid_interval valid_interval)
  (fun (i1, i2) ->
     let m = manually max_range i1 i2 in
     let plst = Interval.cpair max_range i1 i2 in
     let t = expand_ranges_back_to_pairs max_range plst in
     match plst with
     | []       -> m = []
     | p1 :: ps ->
        let in_order, _last =
          List.fold_left ps ~init:(true, p1) ~f:(fun (io, prev) i ->
            let nio = io && Interval.strict_before prev i in
            nio, i)
        in
        in_order && t = m)

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

let manually_lst n ilst =
  let iplst = all_sequential_pairs ilst in
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
  add_test ~count:1000
    ~name:"Pairing sets round trip."
    valid_set
    (fun s ->
      let plst = Set.cpair max_range s in
      let t = expand_ranges_back_to_pairs max_range plst in
      let m = manually_lst max_range s in
      match plst with
      | []       -> m = []
      | p1 :: ps ->
          let in_order, _last =
            List.fold_left ps ~init:(true, p1) ~f:(fun (io, prev) i ->
              let nio = io && Interval.strict_before prev i in
              nio, i)
          in
          in_order && t = m)

let () =
  QCheck_runner.run_tests_main (List.rev !tests)

(*BISECT-IGNORE-END*)
