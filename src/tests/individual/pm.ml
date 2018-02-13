(** Test partition map *)
(*BISECT-IGNORE-BEGIN*)

open QCheck

let is_invalid f =
  try
    ignore (f ());
    false
  with (Invalid_argument _) ->
    true

let count = 20

let strictly_neg_int =
  make Gen.(map (fun i -> i - 1) neg_int)

let t1 = Test.make ~count
  ~name:"Interval construction fail on negative values."
  (pair strictly_neg_int strictly_neg_int)
  (fun (start, end_) -> is_invalid (fun () -> Interval.make ~start ~end_))

let t2 =
  let larger_than_max_value = int_range (Interval.max_value + 1) max_int in
  Test.make ~count
    ~name:"Interval construction fail on values too large."
    (pair larger_than_max_value larger_than_max_value)
    (fun (start, end_) -> is_invalid (fun () -> Interval.make ~start ~end_))

let valid_range =
  int_range 0 Interval.max_value

let t3 = Test.make ~count
  ~name:"Interval construction fail on out of order values."
  (pair valid_range valid_range)
  (fun (start, end_) ->
     assume (start > end_);
     is_invalid (fun () -> Interval.make ~start ~end_))

(* Compress the range so that we can actually have a non-zero
 * chance of intersects *)

let count = 100

let valid_interval_gen =
  let start_range = Gen.int_range 0 1000 in
    Gen.(map (fun (s, i) -> s, s + i) (pair start_range nat))

let valid_interval_pair =
  make ~print:(fun (s, e) -> sprintf "(%d,%d)" s e)
    valid_interval_gen

let t4 = Test.make ~count
  ~name:"Interval construction for non negative works and not none."
  valid_interval_pair
  (fun (start, end_) ->
     let i = Interval.make ~start ~end_ in
     not (Interval.is_none i))

let t5 = Test.make ~count
  ~name:"Interval start property is correct."
  valid_interval_pair
  (fun (start, end_) ->
     let i = Interval.make ~start ~end_ in
     start = Interval.start i)

let t6 = Test.make ~count
  ~name:"Interval end property is correct."
  valid_interval_pair
  (fun (start, end_) ->
     let i = Interval.make ~start ~end_ in
      end_ = Interval.end_ i)

let valid_interval =
  make ~print:Interval.to_string
    Gen.(map (fun (start, end_) -> Interval.make ~start ~end_)
           valid_interval_gen)

let t7 = Test.make ~count
  ~name:"Interval width property is correct."
  valid_interval
  (fun i -> Interval.(width i = end_ i - start i + 1))

let t8 = Test.make ~count
  ~name:"Interval average is inside."
  valid_interval
  (fun i -> Interval.(inside ((start i + end_ i) / 2) i))

let t9 = Test.make ~count
  ~name:"Interval extending just changes the end."
  valid_interval
  (fun i ->
     let open Interval in
     let ip1 = extend_one i in
     start i = start ip1 && end_ i + 1 = end_ ip1)

let t10 = Test.make ~count
  ~name:"Interval merging works."
  (pair valid_interval valid_interval)
  (fun (i1, i2) ->
     let open Interval in
     let m = merge i1 i2 in
     if end_ i1 + 1 = start i2 then
       start m = start i1 && end_ m = end_ i2
     else
       is_none m)

let tests =
  [ t1
  ; t2
  ; t3
  ; t4
  ; t5
  ; t6
  ; t7
  ; t8
  ; t9
  ; t10
  ]

let () =
  QCheck_runner.run_tests_main tests

(*BISECT-IGNORE-END*)
