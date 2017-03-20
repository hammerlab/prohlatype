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

module Float = struct
  let ( + ) x y = x +. y
  let ( - ) x y = x -. y
  let ( * ) x y = x *. y
  let ( / ) x y = x /. y
end

module TransitionMatrix = struct

  let match_idx     = 0
  let insertion_idx = 1
  let deletion_idx  = 2
  let start_end_idx = 3 (* Compress the two nodes ala the Samtools note. *)

  let init ?(model_probs=`Default) ~ref_length read_length =
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
    let to_index = function
      | `Match      -> match_idx
      | `Insert     -> insertion_idx
      | `Delete     -> deletion_idx
      | `StartOrEnd -> start_end_idx
    in
    fun f t -> tm.(to_index f).(to_index t)

  type state = [ `Match | `Insert | `Delete | `StartOrEnd ]
  type t = state -> state -> float
end (* TransitionMatrix *)

type emission_prob = int -> int -> float

type fwd_recurrences =
  { start   : emission_prob -> insert_prob:float -> TransitionMatrix.t -> int -> float * float * float
  ; match_  : emission_prob -> TransitionMatrix.t -> float array array -> i:int -> int -> float
  ; insert  : insert_prob:float -> TransitionMatrix.t -> float array array -> i:int -> int -> float
  ; delete  : TransitionMatrix.t -> float array array -> i:int -> int -> float
  ; end_    : TransitionMatrix.t -> float array array -> i:int -> int -> float * float * float
  }

(* We are storing match/insert/delete states as sequential rows in the same
   matrix. These methods are helper indicies. *)
let mi i = 3 * i
let ii i = 3 * i + 1
let di i = 3 * i + 2

let forward_recurrences =
  { start   = begin fun emission_prob ~insert_prob tm k ->
                (tm `StartOrEnd `Match) *. (emission_prob 0 k)
                , (tm `StartOrEnd `Insert) *. insert_prob
                , 0.
              end
  ; match_  = begin fun emission_prob tm dm ~i -> function
                | 0 -> 0.
                | k -> (emission_prob i k) *.
                          (  (tm `Match `Match)   *. dm.(i-1).(mi (k - 1))
                          +. (tm `Insert `Match)  *. dm.(i-1).(ii (k - 1))
                          +. (tm `Delete `Match)  *. dm.(i-1).(di (k - 1)))
              end
   (* Notice that we do not transition between successive inserts: we're just
      extending the current (same k) insert. *)
  ; insert  = begin fun ~insert_prob tm dm ~i k ->
                insert_prob *.
                  (   (tm `Match `Insert)   *. dm.(i-1).(mi k)
                  +.  (tm `Insert `Insert)  *. dm.(i-1).(ii k))
              end
  ; delete  = begin fun tm dm ~i -> function
                | 0 -> 0.
                | k -> ( (tm `Match `Delete)  *. dm.(i).(mi (k - 1))
                       +.(tm `Delete `Delete) *. dm.(i).(di (k - 1)))
              end
  ; end_    = begin fun tm dm ~i k ->
                (tm `Match `StartOrEnd) *. dm.(i).(mi k)
                , (tm `Insert `StartOrEnd) *. dm.(i).(ii k)
                , 0. (* (since (tm `Delete `StartOrEnd) = 0. *)
              end
  }

let create_workspace_matrix ~read_length ~ref_length =
  Array.make_matrix ~dimx:(read_length + 1) ~dimy:(ref_length * 3) 0.0

(* TODO:
    - Make N tolerant for sequences/reads. This would require changing p_c_m to
      return a uniform error probability (ie 1) ?
    - Transition to computing everything in logs. *)
let forward_gen recurrences ?m ~normalize ?model_probs ~refs ~read read_probs =
  let read_length = String.length read in   (* l *)
  let ref_length = String.length refs in    (* L *)
  let tm = TransitionMatrix.init ?model_probs ~ref_length read_length in
  let insert_prob = 0.25 in (* There are 4 characters, assume equality. *)
  let p_c_m i k =
    if read.[i] = refs.[k] then
      1. -. read_probs.(i)
    else
      read_probs.(i) /. 3.
  in
  let m =
    match m with
    | Some m -> m
    | None   -> create_workspace_matrix ~read_length ~ref_length
  in
  let over_row ~init ~g ~row f =
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
  let sum_over = over_row ~init:0. ~g:(fun s mv iv dv -> s +. mv +. iv +. dv) in
  let iter_over = over_row ~init:() ~g:(fun u _ _ _ -> u) in
  let normalize row by =
    if normalize then
      iter_over row
        (fun k -> m.(row).(mi k) /. by
                , m.(row).(ii k) /. by
                , m.(row).(di k) /. by)
    else
      ()
  in
  let rl1 = read_length - 1 in
  let s1 = sum_over ~row:0 (recurrences.start p_c_m ~insert_prob tm) in
  normalize 0 s1;
  let f_m = recurrences.match_ p_c_m tm m in
  let i_m = recurrences.insert ~insert_prob tm m in
  let d_m = recurrences.delete tm m in
  for i = 1 to rl1 do
    let s = sum_over ~row:i (fun k -> f_m ~i k, i_m ~i k, d_m ~i k) in
    normalize i s;
  done;
  let fl = sum_over ~row:read_length (recurrences.end_ tm m ~i:rl1) in
  (*normalize read_length fl;*)
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
                            ((tm `Match `Match)      *. dm.(i-1).(mi (k - 1)))
                            ((tm `Insert `Match)  *. dm.(i-1).(ii (k - 1)))
                            ((tm `Delete `Match)   *. dm.(i-1).(di (k - 1))))
              end
  ; insert  = begin fun ~insert_prob tm dm ~i k ->
                insert_prob *.
                  (max ((tm `Match `Insert)      *. dm.(i-1).(mi k))
                       ((tm `Insert `Insert)  *. dm.(i-1).(ii k)))
              end
  ; delete  = begin fun tm dm ~i -> function
                | 0 -> 0.
                | k -> max ((tm `Match `Delete)     *. dm.(i).(mi (k - 1)))
                           ((tm `Delete `Delete)  *. dm.(i).(di (k - 1)))
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

let recover_path vm =
  let rec prev_match i pos acc =
    (*printf "prev_match i: %d pos: %d: v: %f\n" i pos vm.(i).(pos);*)
    let nacc = `M (i, pos / 3) :: acc in
    if i = 0 then
      nacc
    else
      let pm = pos - 3 in
      let pi = pos - 2 in
      let pd = pos - 1 in
      let im1 = i - 1 in
      (*printf "prev_match M: i: %d pos: %d: v: %f\n" im1 pm vm.(im1).(pm);
      printf "prev_match I: i: %d pos: %d: v: %f\n" im1 pi vm.(im1).(pi);
      printf "prev_match D: i: %d pos: %d: v: %f\n" im1 pd vm.(im1).(pd);*)
      if vm.(im1).(pm) >= vm.(im1).(pi) then begin
        if vm.(im1).(pm) >= vm.(im1).(pd) then
          prev_match im1 pm nacc
        else
          prev_delete im1 pd nacc
      end else if vm.(im1).(pi) >= vm.(im1).(pd) then
        prev_insert im1 pi nacc
      else
        prev_delete im1 pd nacc
  and prev_insert i pos acc =
    (*printf "prev_insert i: %d pos: %d v: %f\n" i pos vm.(i).(pos);*)
    let nacc = `I i :: acc in
    if i = 0 then
      nacc
    else
      let im1 = i - 1 in
      let pm = pos - 1 in
      let pi = pos in
      if vm.(im1).(pm) >= vm.(im1).(pi) then
        prev_match im1 pm nacc
      else
        prev_insert im1 pi nacc
  and prev_delete i pos acc =       (* Can't start on a delete! *)
    (*printf "prev_delete i: %d pos: %d v: %f\n" i pos vm.(i).(pos);*)
    let nacc = `D (pos / 3) :: acc in
    let pm = pos - 5 in
    let pd  =pos - 3 in
    if vm.(i).(pm) >= vm.(i).(pd) then
      prev_match i pm nacc
    else
      prev_delete i pd nacc
  in
  let last = Array.length vm - 2 in
  let pos = max_pos vm.(last) in
  match pos mod 3 with
  | 0 -> prev_match last pos []
  | 1 -> prev_insert last pos []
  | 2 -> prev_delete last pos []
  | x -> assert false

let viterbi ?(normalize=true) ?model_probs ~refs ~read read_probs =
  let vm, _ =
    forward_gen viterbi_recurrences ~normalize ?model_probs
      ~refs ~read read_probs
  in
  vm, recover_path vm

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

type bwd_recurrences =
  { last    : TransitionMatrix.t -> int -> float * float * float
  ; first   : emission_prob -> insert_prob:float -> TransitionMatrix.t ->
                float array array -> int -> float * float * float
  ; match_  : emission_prob -> insert_prob:float -> TransitionMatrix.t ->
                float array array -> i:int -> int -> float
  ; insert  : emission_prob -> insert_prob:float -> TransitionMatrix.t ->
                float array array -> i:int -> int -> float
  ; delete  : emission_prob -> TransitionMatrix.t -> float array array ->
                i:int -> int -> float
  }

let backward_recurrences ref_size =
  { last    = begin fun tm _ ->
                (tm `Match `StartOrEnd )
                , (tm `Insert `StartOrEnd)
                , 0.
              end
  ; first   = begin fun emission_prob ~insert_prob tm dm k ->
                  (tm `StartOrEnd `Match) *. (emission_prob 0 k) *. dm.(1).(mi k)
                , (tm `StartOrEnd `Insert) *. insert_prob *. dm.(1).(ii k)
                , 0.
              end
  ; match_  = begin fun emission_prob ~insert_prob tm dm ~i k ->
                let dmm, dmd =
                  if k + 1 = ref_size then
                    0., 0.
                  else
                    dm.(i+1).(mi (k+1)), dm.(i).(di (k+1))
                in
                let emission_p = emission_prob i (k + 1) in
                   emission_p   *. (tm `Match `Match)   *. dmm
                +. insert_prob  *. (tm `Match `Insert)  *. dm.(i+1).(ii k)
                +.                 (tm `Match `Delete ) *. dmd
              end
  ; insert  = begin fun emission_prob ~insert_prob tm dm ~i k ->
                let dmm =
                  if k + 1 = ref_size then
                    0.
                  else
                    dm.(i+1).(mi (k+1))
                in
                let emission_p = emission_prob i (k + 1) in
                   emission_p   *. (tm `Insert `Match)  *. dmm
                +. insert_prob  *. (tm `Insert `Insert) *. dm.(i+1).(ii k)
              end
  ; delete  = begin fun emission_prob tm dm ~i k->
                let dmm, dmd =
                    if k + 1 = ref_size then
                      0., 0.
                    else
                      dm.(i+1).(mi (k+1)), dm.(i).(di (k+1))
                in
                let emission_p = emission_prob i (k + 1) in
                emission_p *. (tm `Delete `Match)  *. dmm
                +.            (tm `Delete `Delete) *. dmd
              end
 }

let backward_gen recurrences ~normalize ?model_probs ~refs ~read read_probs =
  let read_length = String.length read in   (* l *)
  let ref_length = String.length refs in    (* L *)
  let tm = TransitionMatrix.init ?model_probs ~ref_length read_length in
  let insert_prob = 0.25 in
  let p_c_m i k =
    if k = ref_length || i = ref_length then
      0.0
    else if read.[i] = refs.[k] then
      1. -. read_probs.(i)
    else
      read_probs.(i) /. 3.
  in
  let m = Array.make_matrix ~dimx:(read_length + 1) ~dimy:(ref_length * 3) 1. in
  let over_row ~init ~g row f =
    let rec loop s k =
      if k < 0 then
        s
      else begin
        let mv, iv, dv = f k in
        m.(row).(mi k) <- mv;
        m.(row).(ii k) <- iv;
        m.(row).(di k) <- dv;
        loop (g s mv iv dv) (k - 1)
      end
    in
    loop init (ref_length - 1)
  in
  let sum_over = over_row ~init:0. ~g:(fun s mv iv dv -> s +. mv +. iv +. dv) in
  let iter_over = over_row ~init:() ~g:(fun u _ _ _ -> u) in
  let normalize row by =
    if normalize then
      iter_over row
        (fun k -> m.(row).(mi k) /. by
                , m.(row).(ii k) /. by
                , m.(row).(di k) /. by)
    else
      ()
  in
  let rl1 = read_length - 1 in
  let recurrences = recurrences ref_length in
  let s1 = sum_over read_length (recurrences.last tm) in
  normalize read_length s1;
  let f_m = recurrences.match_ p_c_m ~insert_prob tm m in
  let i_m = recurrences.insert p_c_m ~insert_prob tm m in
  let d_m = recurrences.delete p_c_m tm m in
  for i = rl1 downto 1 do
    let s = sum_over i (fun k -> f_m ~i k, i_m ~i k, d_m ~i k) in
    normalize i s;
  done;
  let bl = sum_over 0 (recurrences.first p_c_m ~insert_prob tm m) in
  (*normalize 0 bl;*)
  m, bl

let backward ?(normalize=true) =
  backward_gen ~normalize backward_recurrences

let posterior ?(normalize=true) ?model_probs ~refs ~read read_probs =
  let fm, lh = forward ~normalize ?model_probs ~refs ~read read_probs in
  let bm, _ = backward ~normalize ?model_probs ~refs ~read read_probs in
  let mm =
    Array.mapi fm ~f:(fun i fmr ->
      Array.mapi fmr ~f:(fun j v -> v *. bm.(i).(j) /. lh))
  in
  (* Could do the path recover without constructing mm. *)
  let path =
    Array.init (Array.length mm - 1) ~f:(fun i ->
      let j = max_pos mm.(i) in
      match j mod 3 with
      | 0 -> `M (i, j / 3)
      | 1 -> `I i
      | 2 -> `D (j / 3)
      | _ -> assert false)
    |> Array.to_list
  in
  mm, path
