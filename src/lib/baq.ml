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


type emission_prob = int -> int -> float
type transition_matrix = float array array

type recurrences =
  { start   : emission_prob -> insert_prob:float -> transition_matrix -> int -> float * float * float
  ; match_  : emission_prob -> transition_matrix -> float array array -> i:int -> int -> float
  ; insert  : insert_prob:float -> transition_matrix -> float array array -> i:int -> int -> float
  ; delete  : transition_matrix -> float array array -> i:int -> int -> float
  ; end_    : transition_matrix -> float array array -> i:int -> int -> float * float * float
  }

(* We are storing match/insert/delete states as sequential rows in the same
   matrix. These methods are helper indicies. *)
let mi i = 3 * i
let ii i = 3 * i + 1
let di i = 3 * i + 2

let forward_recurrences =
  { start   = begin fun emission_prob ~insert_prob tm k ->
                tm.(start_end_idx).(match_idx) *. (emission_prob 0 k)
                , tm.(start_end_idx).(insertion_idx) *. insert_prob
                , 0.
              end
  ; match_  = begin fun emission_prob tm dm ~i -> function
                | 0 -> 0.
                | k -> (emission_prob i k) *.
                          (  tm.(match_idx).(match_idx)      *. dm.(i-1).(mi (k - 1))
                          +. tm.(insertion_idx).(match_idx)  *. dm.(i-1).(ii (k - 1))
                          +. tm.(deletion_idx).(match_idx)   *. dm.(i-1).(di (k - 1)))
              end
   (* Notice that we do not transition between successive inserts: we're just
      extending the current (same k) insert. *)
  ; insert  = begin fun ~insert_prob tm dm ~i k ->
                insert_prob *.
                  (   tm.(match_idx).(insertion_idx)      *. dm.(i-1).(mi k)
                  +.  tm.(insertion_idx).(insertion_idx)  *. dm.(i-1).(ii k))
              end
  ; delete  = begin fun tm dm ~i -> function
                | 0 -> 0.
                | k -> ( tm.(match_idx).(deletion_idx)    *. dm.(i).(mi (k - 1))
                       +.tm.(deletion_idx).(deletion_idx) *. dm.(i).(di (k - 1)))
              end
  ; end_    = begin fun tm dm ~i k ->
                tm.(match_idx).(start_end_idx) *. dm.(i).(mi k)
                , tm.(insertion_idx).(start_end_idx) *. dm.(i).(ii k)
                , 0. (* since tm.(deletion_idx).(start_end_idx) = 0. *)
              end
  }

(* TODO:
    - Make N tolerant for sequences/reads. This would require changing p_c_m to
      return a uniform error probability (ie 1) ?
    - Transition to computing everything in logs. *)
let forward_gen recurrences ~normalize ?model_probs ~refs ~read read_probs =
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
  let m = Array.make_matrix ~dimx:(read_length + 1) ~dimy:(ref_length * 3) 0. in
  let over_reference ~init ~g row f =
    let rec loop s k =
      if k = ref_length then
        s
      else begin
        let mv, iv, dv = f k in
        m.(row).(mi k) <- mv;
        m.(row).(ii k) <- iv;
        m.(row).(di k) <- dv;
        loop (g s mv iv dv) (k + 1)
      end
    in
    loop init 0
  in
  let sum_over_reference = over_reference ~init:0. ~g:(fun s mv iv dv -> s +. mv +. iv +. dv) in
  let iter_over_reference = over_reference ~init:() ~g:(fun u _ _ _ -> u) in
  let normalize row by =
    if normalize then
      iter_over_reference row
        (fun k -> m.(row).(mi k) /. by
                , m.(row).(ii k) /. by
                , m.(row).(di k) /. by)
    else
      ()
  in
  let rl1 = read_length - 1 in
  let s1 = sum_over_reference 0 (recurrences.start p_c_m ~insert_prob tm) in
  normalize 0 s1;
  let f_m = recurrences.match_ p_c_m tm m in
  let i_m = recurrences.insert ~insert_prob tm m in
  let d_m = recurrences.delete tm m in
  for i = 1 to rl1 do
    let s = sum_over_reference i (fun k -> f_m ~i k, i_m ~i k, d_m ~i k) in
    normalize i s;
  done;
  let fl = sum_over_reference read_length (recurrences.end_ tm m ~i:rl1) in
  normalize read_length fl;
  m, fl

let forward ?(normalize=true) =
  forward_gen forward_recurrences ~normalize

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

let max3 x y z = max x (max y z)

let viterbi_recurrences =
  { forward_recurrences with
    match_  = begin fun emission_prob tm dm ~i -> function
                | 0 -> 0.
                | k -> (emission_prob i k) *.
                          (max3
                            (tm.(match_idx).(match_idx)      *. dm.(i-1).(mi (k - 1)))
                            (tm.(insertion_idx).(match_idx)  *. dm.(i-1).(ii (k - 1)))
                            (tm.(deletion_idx).(match_idx)   *. dm.(i-1).(di (k - 1))))
              end
  ; insert  = begin fun ~insert_prob tm dm ~i k ->
                insert_prob *.
                  (max (tm.(match_idx).(insertion_idx)      *. dm.(i-1).(mi k))
                       (tm.(insertion_idx).(insertion_idx)  *. dm.(i-1).(ii k)))
              end
  ; delete  = begin fun tm dm ~i -> function
                | 0 -> 0.
                | k -> max (tm.(match_idx).(deletion_idx)     *. dm.(i).(mi (k - 1)))
                           (tm.(deletion_idx).(deletion_idx)  *. dm.(i).(di (k - 1)))
              end
  }

let max_pos ?upto = function
  | [||]  -> invalid_argf "Empty array"
  | a     -> let upto = Option.value upto ~default:(Array.length a) in
             let rec loop bv bi i =
               if i >= upto then bi else
                 let v = a.(i) in
                 if v >= bv then
                   loop v i (i + 1)
                 else
                   loop bv bi (i + 1)
             in
             loop a.(0) 0 0

let viterbi ?(normalize=true) ?model_probs ~refs ~read read_probs =
  let vm, _ =
    forward_gen viterbi_recurrences ~normalize ?model_probs
      ~refs ~read read_probs
  in
  let rec prev_match i pos acc =
    let nacc = `M (i, pos / 3) :: acc in
    if i = 0 then
      nacc
    else
      let pm = vm.(i-1).(pos - 3) in
      let pi = vm.(i-1).(pos - 2) in
      let pd = vm.(i-1).(pos - 1) in
      if pm >= pi then begin
        if pm >= pd then
          prev_match (i-1) (pos - 3) nacc
        else
          prev_delete (i-1) (pos - 1) nacc
      end else if pi >= pd then
        prev_insert (i-1) (pos - 1) nacc
      else
        prev_delete (i-1) (pos - 1) nacc
  and prev_insert i pos acc =
    let nacc = `I i :: acc in
    if i = 0 then
      nacc
    else
      let pm = vm.(i-1).(pos - 1) in
      let pi = vm.(i-1).(pos) in
      if pm >= pi then
        prev_match (i-1) (pos - 1) nacc
      else
        prev_insert (i-1) pos nacc
  and prev_delete i pos acc =       (* Can't start on a delete! *)
    let nacc = `D (pos / 3) :: acc in
    let pm = vm.(i).(pos - 5) in
    let pd = vm.(i).(pos - 3) in
    if pm >= pd then
      prev_match i (pos - 5) nacc
    else
      prev_delete i (pos - 3) nacc
  in
  let last = Array.length vm - 2 in
  let pos = max_pos vm.(last) in
  let pth =
    match pos mod 3 with
    | 0 -> prev_match last pos []
    | 1 -> prev_insert last pos []
    | 2 -> prev_delete last pos []
    | x -> assert false
  in
  vm, pth

let expand_path ?(insert_char ='-') ~refs ~read =
  let rec loop sref sread = function
    | []            -> String.of_character_list (List.rev sref),
                       String.of_character_list (List.rev sread)
    | `M (i,p) :: t -> loop (String.get_exn refs ~index:p :: sref)
                            (String.get_exn read ~index:i :: sread) t
    | `I i :: t     -> loop (insert_char :: sref)
                            (String.get_exn read ~index:i :: sread) t
    | `D _ :: t     -> loop sref sread t
  in
  loop [] []
