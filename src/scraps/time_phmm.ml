
open StdLabels
open Core_bench.Std

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

  let init ?(model_probs=`Default) read_length ref_length =
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
  let tm = TransitionMatrix.init ?model_probs read_length ref_length in
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
  let over_row ~init ~g row f =
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
  let s1 = sum_over 0 (recurrences.start p_c_m ~insert_prob tm) in
  normalize 0 s1;
  let f_m = recurrences.match_ p_c_m tm m in
  let i_m = recurrences.insert ~insert_prob tm m in
  let d_m = recurrences.delete tm m in
  for i = 1 to rl1 do
    let s = sum_over i (fun k -> f_m ~i k, i_m ~i k, d_m ~i k) in
    normalize i s;
  done;
  let fl = sum_over read_length (recurrences.end_ tm m ~i:rl1) in
  (*normalize read_length fl;*)
  m, fl

let forward ?(normalize=true) =
  forward_gen forward_recurrences ~normalize


let refs =
  "CCTCGTTCAGGGCGATGTAATCCTTGCCGTCGTAGGCGGACTGGTCATGCCCGCGGAGGAGGCGCCCGTCGGGCCCCAGGTCGCAGCCATACATCCTCTG"

let full_refs =
  "CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGGCGTGGCTCTCAGGGTCTCAGGCCCCGAAGGCGGTGTATGGATTGGGGAGTCCCAGCCTTGGGGATTCCCCAACTCCGCAGTTTCTTTTCTCCCTCTCCCAACCTACGTAGGGTCCTTCATCCTGGATACTCACGACGCGGACCCAGTTCTCACTCCCATTGGGTGTCGGGTTTCCAGAGAAGCCAATCAGTGTCGTCGCGGTCGCTGTTCTAAAGTCCGCACGCACCCACCGGGACTCAGATTCTCCCCAGACGCCGAGGATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGTGAGTGCGGGGTCGGGAGGGAAACCGCCTCTGCGGGGAGAAGCAAGGGGCCCTCCTGGCGGGGGCGCAGGACCGGGGGAGCCGCGCCGGGAGGAGGGTCGGGCAGGTCTCAGCCACTGCTCGCCCCCAGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTGAGTGACCCCGGCCCGGGGCGCAGGTCACGACCCCTCATCCCCCACGGACGGGCCAGGTCGCCCACAGTCTCCGGGTCCGAGATCCACCCCGAAGCCGCGGGACTCCGAGACCCTTGTCCCGGGAGAGGCCCAGGCGCCTTTACCCGGTTTCATTTTCAGTTTAGGCCAAAAATCCCCCCGGGTTGGTCGGGGCGGGGCGGGGCTCGGGGGACTGGGCTGACCGCGGGGTCGGGGCCAGGTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCCGCGGGTACCGGCAGGACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGATCACCAAGCGCAAGTGGGAGGCGGTCCATGCGGCGGAGCAGCGGAGAGTCTACCTGGAGGGCCGGTGCGTGGACGGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGGGTACCAGGGGCCACGGGGCGCCTCCCTGATCGCCTATAGATCTCCCGGGCTGGCCTCCCACAAGGAGGGGAGACAATTGGGACCAACACTAGAATATCACCCTCCCTCTGGTCCTGAGGGAGAGGAATCCTCCTGGGTTTCCAGATCCTGTACCAGAGAGTGACTCTGAGGTTCCGCCCTGCTCTCTGACACAATTAAGGGATAAAATCTCTGAAGGAGTGACGGGAAGACGATCCCTCGAATACTGATGAGTGGTTCCCTTTGACACCGGCAGCAGCCTTGGGCCCGTGACTTTTCCTCTCAGGCCTTGTTCTCTGCTTCACACTCAATGTGTGTGGGGGTCTGAGTCCAGCACTTCTGAGTCTCTCAGCCTCCACTCAGGTCAGGACCAGAAGTCGCTGTTCCCTTCTCAGGGAATAGAAGATTATCCCAGGTGCCTGTGTCCAGGCTGGTGTCTGGGTTCTGTGCTCTCTTCCCCATCCCGGGTGTCCTGTCCATTCTCAAGATGGCCACATGCGTGCTGGTGGAGTGTCCCATGACAGATGCAAAATGCCTGAATTTTCTGACTCTTCCCGTCAGACCCCCCCAAGACACATATGACCCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGAGAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTCTGCCCAAGCCCCTCACCCTGAGATGGGGTAAGGAGGGAGATGGGGGTGTCATGTCTCTTAGGGAAAGCAGGAGCCTCTCTGGAGACCTTTAGCAGGGTCAGGGCCCCTCACCTTCCCCTCTTTTCCCAGAGCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCCGTGATGTGGAGGAGGAAGAGCTCAGGTGGAGAAGGGGTGAAGGGTGGGGTCTGAGATTTCTTGTCTCACTGAGGGTTCCAAGCCCCAGCTAGAAATGTGCCCTGTCTCATTACTGGGAAGCACCTTCCACAATCATGGGCCGACCCAGCCTGGGCCCTGTGTGCCAGCACTTACTCTTTTGTAAAGCACCTGTTAAAATGAAGGACAGATTTATCACCTTGATTACGGCGGTGATGGGACCTGATCCCAGCAGTCACAAGTCACAGGGGAAGGTCCCTGAGGACAGACCTCAGGAGGGCTATTGGTCCAGGACCCACACCTGCTTTCTTCATGTTTCCTGATCCCGCCCTGGGTCTGCAGTCACACATTTCTGGAAACTTCTCTGGGGTCCAAGACTAGGAGGTTCCTCTAGGACCTTAAGGCCCTGGCTCCTTTCTGGTATCTCACAGGACATTTTCTTCCCACAGATAGAAAAGGAGGGAGTTACACTCAGGCTGCAAGTAAGTATGAAGGAGGCTGATGCCTGAGGTCCTTGGGATATTGTGTTTGGGAGCCCATGGGGGAGCTCACCCACCCCACAATTCCTCCTCTAGCCACATCTTCTGTGGGATCTGACCAGGTTCTGTTTTTGTTCTACCCCAGGCAGTGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGTAAAGGTGAGAGCTTGGAGGGCCTGATGTGTGTTGGGTGTTGGGTGGAACAGTGGACACAGCTGTGCTATGGGGTTTCTTTGCGTTGGATGTATTGAGCATGCGATGGGCTGTTTAAGGTGTGACCCCTCACTGTGATGGATATGAATTTGTTCATGAATATTTTTTTCTATAGTGTGAGACAGCTGCCTTGTGTGGGACTGAGAGGCAAGAGTTGTTCCTGCCCTTCCCTTTGTGACTTGAAGAACCCTGACTTTGTTTCTGCAAAGGCACCTGCATGTGTCTGTGTTCGTGTAGGCATAATGTGAGGAGGTGGGGAGAGCACCCCACCCCCATGTCCACCATGACCCTCTTCCCACGCTGACCTGTGCTCCCTCTCCAATCATCTTTCCTGTTCCAGAGAGGTGGGGCTGAGGTGTCTCCATCTCTGTCTCAACTTCATGGTGCACTGAGCTGTAACTTCTTCCTTCCCTATTAAAA" ;;

let read = 
  "GGGCTCGGGGGACGGGGCTGACCGCGGGGCCGGGGCCAGGGTCTCACATCATCCAGAGGATGTATGGCTGCGACCTGGGGCCCGACGGGCGCCTCCTCCG"

let read_prob =
  [|0.000398107170553497352; 0.000398107170553497352; 0.000398107170553497352;
    0.000199526231496887876; 0.000199526231496887876; 0.000199526231496887876;
    0.000199526231496887876; 0.000199526231496887876; 0.000158489319246111419;
    0.000125892541179416744; 0.000125892541179416744; 0.000158489319246111419;
    0.000125892541179416744; 7.94328234724282208e-05; 7.94328234724282208e-05;
    7.94328234724282208e-05; 7.94328234724282208e-05; 7.94328234724282208e-05;
    7.94328234724282208e-05; 7.94328234724282208e-05; 7.94328234724282208e-05;
    7.94328234724282208e-05; 0.000125892541179416744; 0.000125892541179416744;
    0.000199526231496887876; 0.000199526231496887876; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000398107170553497352;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000501187233627272528;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000251188643150957955; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000398107170553497352;
    0.000316227766016837939; 0.000251188643150957955; 0.000251188643150957955;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000501187233627272528;
    0.000501187233627272528; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000316227766016837939; 0.000316227766016837939;
    0.000316227766016837939; 0.000501187233627272528; 0.000316227766016837939;
    0.000501187233627272528; 0.000316227766016837939; 0.000316227766016837939;
    0.000398107170553497352; 0.000398107170553497352; 0.000316227766016837939;
    0.000316227766016837939|]

let () =
  let read_length = String.length read in
  let ref_length  = String.length refs in
  let workspace = create_workspace_matrix ~read_length ~ref_length in
  let f () =
    ignore (forward ~m:workspace ~refs ~read read_prob)
  in

  let full_ref_length  = String.length full_refs in
  let full_workspace = create_workspace_matrix ~read_length ~ref_length:full_ref_length in
  let f2 () =
    ignore (forward ~m:full_workspace ~refs:full_refs ~read read_prob)
  in
  Core.Command.run (Bench.make_command
    [ Bench.Test.create ~name:"forward pass" f
    ; Bench.Test.create ~name:"forward pass against full reference" f2
  ])
