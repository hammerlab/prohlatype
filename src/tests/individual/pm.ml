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

let t4 = Test.make ~count
  ~name:"Interval construction for non negative works and not none."
  (pair valid_range valid_range)
  (fun (start, end_) ->
     assume (start <= end_);
     let i = Interval.make ~start ~end_ in
     not (Interval.is_none i))
 
let tests =
  [ t1
  ; t2
  ; t3
  ; t4
  ]

let () =
  QCheck_runner.run_tests_main tests

(*BISECT-IGNORE-END*)
