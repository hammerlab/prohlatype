(* From Improving SNP discovery by base alignment quality *)
open Util

type probabilities =
  { gap_open      : float     (* alpha *)
  ; gap_extension : float     (* beta *)
  ; ending        : float     (* gamma *)
  }

let default_probabilities ?(read_length=100) () =
  { gap_open      = 0.001
  ; gap_extension = 0.1
  ; ending        = 1. /. (2.0 *. float read_length)
  }

let match_idx     = 0
let insertion_idx = 1
let deletion_idx  = 2
let start_end_idx = 3 (* Compress the two nodes ala the Samtools note. *)

module Float = struct
  let ( + ) x y = x +. y
  let ( - ) x y = x -. y
  let ( * ) x y = x *. y
  let ( / ) x y = x /. y
end

let to_transition_matrix ?(model_probs=`Default) read_length ref_length =
  let mp =
    match model_probs with
    | `Default                 -> default_probabilities ~read_length ()
    | `WithAverageReadLength l -> default_probabilities ~read_length:l ()
    | `Spec model              -> model
  in
  (* Rename the probabilities to follow the convention in the paper.
     Still uncertain if it actually makes it easier to understand transition
     probability calculations. *)
  let alpha = mp.gap_open in
  let beta  = mp.gap_extension in
  let gamma = mp.ending in
  let tm = Array.make_matrix ~dimx:4 ~dimy:4 0. in
  let ll = float ref_length in
  let open Float in
  tm.(match_idx).(match_idx)            <- (1. - 2. * alpha) * (1. - gamma);
  tm.(match_idx).(insertion_idx)        <- alpha * (1. - gamma);
  tm.(match_idx).(deletion_idx)         <- alpha * (1. - gamma);
  tm.(match_idx).(start_end_idx)        <- gamma;

  tm.(insertion_idx).(match_idx)        <- (1. - beta) * (1. - gamma);
  tm.(insertion_idx).(insertion_idx)    <- beta * (1. - gamma);
  (*tm.(insertion_idx).(deletion_idx)   <- 0.; *)
  tm.(insertion_idx).(start_end_idx)    <- gamma;

  tm.(deletion_idx).(match_idx)         <- (1. - beta);
  (*tm.(deletion_idx).(insertion_idx)   <- 0.; *)
  tm.(deletion_idx).(deletion_idx)      <- beta;
  (*tm.(deletion_idx).(end_idx)         <- 0.; *)

  tm.(start_end_idx).(match_idx)        <- (1. - alpha) / ll;
  tm.(start_end_idx).(insertion_idx)    <- alpha / ll;
  (*tm.(start_end_idx).(deletion_idx)   <- 0.; *)
  (*tm.(start_end_idx).(start_end_idx)  <- 0.; *)
  tm

(* TODO:
    - Make N tolerant for sequences/reads. This would require changing p_c_m to
      return a uniform error probability (ie 1) ?
    - Transition to computing everything in logs.
*)
let forward ?(normalize=true) ?model_probs ~refs ~read read_probs =
  let read_length = String.length read in   (* l *)
  let ref_length = String.length refs in    (* L *)
  let tm = to_transition_matrix ?model_probs read_length ref_length in
  let insert_prob = 0.25 in (* There are 4 characters, assume equality. *)
  let p_c_m i k =
    if read.[i] = refs.[k] then
      1. -. read_probs.(i)
    else
      read_probs.(i) /. 3.
  in
  let fm = Array.make_matrix ~dimx:read_length ~dimy:(ref_length * 3) 0. in
  let mi i = 3 * i in
  let ii i = 3 * i + 1 in
  let di i = 3 * i + 2 in
  let over_reference ~init ~g col f =
    let rec loop s k =
      if k = ref_length then
        s
      else begin
        let mv, iv, dv = f k in
        fm.(col).(mi k) <- mv;
        fm.(col).(ii k) <- iv;
        fm.(col).(di k) <- dv;
        loop (g s mv iv dv) (k + 1)
      end
    in
    loop init 0
  in
  let sum_over_reference = over_reference ~init:0. ~g:(fun s mv iv dv -> s +. mv +. iv +. dv) in
  let iter_over_reference = over_reference ~init:() ~g:(fun u _ _ _ -> u) in
  let s1 = sum_over_reference 0
            (fun k -> tm.(start_end_idx).(match_idx) *. p_c_m 0 k
                    , tm.(start_end_idx).(insertion_idx) *. insert_prob
                    , 0.)
  in
  if normalize then
    iter_over_reference 0 (fun k -> fm.(0).(mi k) /. s1, fm.(0).(ii k) /. s1, fm.(0).(di k) /. s1);
  let rl1 = read_length - 1 in
  let f_m i = function
    | 0 -> 0.
    | k -> (p_c_m i k) *.
              (  tm.(match_idx).(match_idx)      *. fm.(i-1).(mi (k - 1))
              +. tm.(insertion_idx).(match_idx)  *. fm.(i-1).(ii (k - 1))
              +. tm.(deletion_idx).(match_idx)   *. fm.(i-1).(di (k - 1)))
  in
  let i_m i k =
    (* Notice that we do not transition between successive inserts: we're just
       extending the current (same k) insert. *)
    insert_prob *.
      ( tm.(match_idx).(insertion_idx)      *. fm.(i-1).(mi k)
      +.tm.(insertion_idx).(insertion_idx)  *. fm.(i-1).(ii k))
  in
  let d_m i = function
    | 0 -> 0.
    | k -> ( tm.(match_idx).(deletion_idx)     *. fm.(i).(mi (k - 1))
           +.tm.(deletion_idx).(deletion_idx)  *. fm.(i).(di (k - 1)))
  in
  for i = 1 to rl1 do
    let s = sum_over_reference i (fun k -> f_m i k, i_m i k, d_m i k) in
    if normalize then
      iter_over_reference i (fun k -> fm.(i).(mi k) /. s, fm.(i).(ii k) /. s, fm.(i).(di k) /. s)
  done;
  let fl = sum_over_reference rl1
            (fun k ->
              tm.(match_idx).(start_end_idx) *. fm.(rl1).(mi k)
              , tm.(insertion_idx).(start_end_idx) *. fm.(rl1).(ii k)
              , 0. (* since tm.(deletion_idx).(start_end_idx) = 0. *))
  in
  fm, fl

(* From Accurately Computing log(1 − exp( −| a |))
  Assessed by the Rmpfr package by Martin Maechler *)
let log1mexp =
  let l2 = log 2. in
  fun x ->
    if x <= l2 then
      log (-1. *. expm1 (-.x))
    else
      log1p (-1. *. exp (-.x))

let log1pexp x =
  if x <= -37. then
    exp x
  else if x <= 18. then
    log1p (exp x)
  else if x <= 33. then
    x +. exp (-.x)
  else
    x

