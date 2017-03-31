(* Parametric Profile Hidden Markov Model.
   "Parameterize" the match/insert/delete states by different alleles.

  TODO: Remove type annotations and transition to an interface.
*)

open Util

let array_findi v a =
  let n = Array.length a in
  let rec loop i =
    if i >= n then raise Not_found
    else if a.(i) = v then i
    else loop (i + 1)
  in
  loop 0

(* What are the possible states of the alleles. *)
module BaseState = struct

  type t =
    | Unknown         (* Treat like an 'N', anything goes. *)
    | Gap_            (* Nothing goes! *)
    | A
    | C
    | G
    | T
    [@@deriving show]

  let of_char = function
    | 'A' -> A
    | 'C' -> C
    | 'G' -> G
    | 'T' -> T
    |  c  -> invalid_argf "Unsupported base: %c" c

  let to_char = function
    | Unknown -> 'N'
    | Gap_    -> '_'
    | A       -> 'A'
    | C       -> 'C'
    | G       -> 'G'
    | T       -> 'T'

end (* BaseState *)

(* Run Length Encoded Lists: Encode a sequence of elements by keep track of runs:
   adjacent elements that have the same value. *)
module Rlel = struct

  type 'a pair =
    { value   : 'a
    ; length  : int (* > 0 *)
    }
  and 'a t =
    { hd  : 'a pair
    ; tl  : 'a pair list
    }
    [@@deriving show]

  let total_length { hd; tl} =
    hd.length + List.fold_left tl ~init:0 ~f:(fun s { length; _} -> s + length)

  let expand { hd = {value; length}; tl} =
    let rec expand_rlp v acc = function
      | 0 -> acc
      | n -> expand_rlp v (v :: acc) (n - 1)
    in
    let rec loop acc = function
      | []                    ->
          List.rev acc
      | { value; length} :: t ->
          loop (expand_rlp value acc length) t
    in
    loop (expand_rlp value [] length) tl

  let rec until_different value =
    let rec loop length = function
      | h :: t when h = value -> loop (length + 1) t
      | lst (* when h <> value
      | [] *)                 -> { value; length}, lst
    in
    loop 1

  (* Run length encode a list. *)
  let encode = function
    | []          -> invalid_arg "Empty"
    | value :: t  ->
        let hd, nt = until_different value t in
        let rec loop acc = function
          | []     -> List.rev acc
          | h :: t -> let rlp, nt = until_different h t in
                      loop (rlp :: acc) nt
        in
        { hd
        ; tl = loop [] nt}

  let init value =
    { hd = { value; length = 1}
    ; tl = []
    }

  let append_r v = function
    | { hd = { value; length} ; tl = [] } when v = value ->
        { hd = { value ; length = length + 1}; tl = []}
    | { hd ; tl = { value; length } :: t } when v = value ->
        { hd; tl = { value; length = length + 1} :: t }
    | { hd ; tl } ->
        { hd; tl = { value = v; length = 1} :: tl }

  let finish_r { hd; tl} =
    {hd; tl = List.rev tl}

  let fold_map l ~f ~init =
    let n_hd_value, acc = f init l.hd.value in
    let ntl, facc =
      List.fold_left l.tl ~init:([], acc)
        ~f:(fun (l, acc) e ->
              let nvalue, nacc = f acc e.value in
              { e with value = nvalue } :: l, nacc)
    in
    { hd = { l.hd with value = n_hd_value}; tl = List.rev ntl }, facc

  let align c1 c2 l1 l2 =
    if c1.length < c2.length then
      c1.length
      , l1
      , { c2 with length = c2.length - c1.length } :: l2
    else if c1.length > c2.length then
      c2.length
      , { c1 with length = c1.length - c2.length } :: l1
      , l2
    else (* c1.length = c2.length *)
      c2.length
      , l1
      , l2

  let fold_map2_same_length l1 l2 ~f ~init =
    let hvalue, nacc = f init l1.hd.value l2.hd.value in
    let length, nl1, nl2 = align l1.hd l2.hd l1.tl l2.tl in
    let rec loop lst acc = function
      | [], []             -> { hd = { value = hvalue; length}; tl = List.rev lst}, acc
      | h1 :: t1, h2 :: t2 -> let nvalue, nacc = f acc h1.value h2.value in
                              let length, l1, l2 = align h1 h2 t1 t2 in
                              let npair = { value = nvalue; length} in
                              loop (npair :: lst) nacc (l1, l2)
      | _,        _        -> invalid_arg "different lengths"
    in
    loop [] nacc (nl1, nl2)

  let fill_array ~f ret rl =
    let rec fill_value v i length =
      for j = i to i + length - 1 do ret.(j) <- v done
    in
    let rec loop i = function
      | []                    -> ()
      | { value; length} :: t -> fill_value (f value) i length;
                                 loop (i + length) t
    in
    fill_value (f rl.hd.value) 0 rl.hd.length;
    loop rl.hd.length rl.tl



end (* Rlel *)

type position_map = (Mas_parser.position * int) list

(*** Construction
  1. From Mas_parser.result -> Run-Length-Encoded-Lists Array of BaseState.t 's.
    a. Figure out the BaseState.t of reference sequence and a position map
       (position in Mas_parser.result to index into final array)
       [initialize_base_array_and_position_map].
    b. Start Run-Length encoded lists and extend them with each alternate
       allele.


 ***)

(* Figure out the BaseState.t of the reference and aggregate a position map:
   (position in Mas_parser.result * index into base_state array.) list.

  The mapping (position into base_state array) is then recovered by iterating
  over this position map as we move along Mas_parser.alignment_element's for
  the alleles. *)
let initialize_base_array_and_position_map mp =
  let sequence_to_base_states_array s =
    Array.init (String.length s) ~f:(fun index ->
      BaseState.of_char (String.get_exn s ~index))
  in
  let maybe_add_unknown last_known_mas_pos p pacc acc pos =
    let d = pos - last_known_mas_pos in
    if d > 0 then
      last_known_mas_pos + d
      , p + d
      , (last_known_mas_pos, p) :: pacc
      , Array.make d BaseState.Unknown :: acc
    else
      last_known_mas_pos, p, pacc, acc
  in
  let open Mas_parser in
  let rec loop lkmp p pacc acc = function
    | []                          -> Array.concat (List.rev acc),
                                     List.rev ((lkmp, p) :: pacc)
    (* We're not explicitly tracking what is {un}known through this sequence and
       assuming that it is well formed. So we're only adding unknowns after we
       transition (maybe back) to known state, ie, after Start *)
    | Start pos :: t              -> let nlkmp, np, npacc, nacc = maybe_add_unknown lkmp p pacc acc pos in
                                     loop nlkmp np npacc nacc t
    | End pos :: t                -> loop pos (p + pos - lkmp) ((lkmp, p) :: pacc) acc t
    | Boundary { pos; _ } :: t    -> loop (pos + 1) (p + pos - lkmp) ((lkmp, p) :: pacc) acc t
    | Sequence { s; start } :: t  -> let l = String.length s in
                                     loop (start + l) (p + l) ((lkmp, p) :: pacc)
                                        (sequence_to_base_states_array s :: acc) t
    | Gap { length; gstart } :: t -> loop (gstart + length) (p + length) ((lkmp, p) :: pacc)
                                        (Array.make length BaseState.Gap_ :: acc) t
  in
  match mp.Mas_parser.ref_elems with
  | Start s :: t -> loop s 0 [] [] t
  | e :: _       -> invalid_argf "Reference not at Start : %s" (al_el_to_string e)
  | []           -> invalid_argf "Empty reference sequence!"

(* Remove redundant (difference between the two doesn't change) positions.
   This step is not strictly necessary. *)
let reduce_position_map : position_map -> position_map = function
  | []          -> invalid_arg "empty"
  | (p, a) :: t ->
      List.fold_left t ~init:(a - p, [p,a])
        ~f:(fun (d, acc) (p,a) ->
            let d1 = a - p in
            if d <> d1 then (d1, (p,a) :: acc) else (d, acc))
      |> snd
      |> List.rev

(* Helper method to create an actual function for computing the index into
   the base state array. This is useful for debugging between the Mas_parser
   positions and the index into BaseState array. Assumes a 'reduced' (via
   [reduce_position_map]) position map. *)
let to_position_map : position_map -> (int -> int option) = function
  | [] -> invalid_arg "empty"
  | (p, a) :: t ->
      (* m in reverse *)
      let _, lp, m =
        List.fold_left t ~init:(p - a, min_int, [p, a])
          ~f:(fun (d, _, acc) (p2, a2) ->
              let d2 = p2 - a2 in
              if d2 <> d then d2, p2, (p2, a2) :: acc else d, p2, acc)
      in
      fun x ->
        if x >= lp then None else
          List.find_map m ~f:(fun (p,o) ->
            if p <= x then Some (x + o - p) else None)

(* The position map is just a list of the correct sequential offsets.
   They change as we encounter new boundaries (for UTR/exon/intron breaks).
   When we need to know the correct position we walk (recurse down) this list
   to find the most recent difference and compute the position there.
   We can discard previous differences because we use this method as we merge
   the elements of the alternate alleles in [add_alternate_allele]. *)
let rec position_and_advance sp (pos_map : position_map) =
  match pos_map with
  | []                                                -> invalid_argf "reached end of sequence map: %d" sp
  | (p1, o1) :: (p2, _) :: _ when p1 <= sp && sp < p2 -> sp - p1 + o1, pos_map
  | (p1, o1) :: []           when p1 <= sp            -> sp - p1 + o1, pos_map
  | h :: t                                            -> position_and_advance sp t

(* Add an allele's Mas_parser.sequence to the current run-length encoded state. *)
let add_alternate_allele ~position_map reference_states allele_instr rl_arr =
  let last_pos = Array.length rl_arr in
  let extend_to start lp v =
    for i = start to lp - 1 do
      rl_arr.(i) <- Rlel.append_r v rl_arr.(i)
    done
  in
  let extend_from_base start lp =
    for i = start to lp - 1 do
      rl_arr.(i) <- Rlel.append_r reference_states.(i) rl_arr.(i)
    done
  in
  let rec loop position_map lp =
    let open Mas_parser in function
    | []                          -> ()
    | Start p :: t                ->
        let ap, position_map = position_and_advance p position_map in
        extend_to lp ap BaseState.Unknown;
        loop position_map ap t
    | End p :: t                  ->
        let ap, position_map = position_and_advance p position_map in
        extend_from_base lp ap;
        extend_to ap last_pos BaseState.Unknown;    (* Preemptively set the rest to Unknown. *)
        loop position_map ap t
    | Boundary { pos; _ } :: t    ->
        let ap, position_map = position_and_advance pos position_map in
        extend_from_base lp ap;
        loop position_map ap t
    | Sequence { start; s } :: t  ->
        let ap, position_map = position_and_advance start position_map in
        extend_from_base lp ap;
        let n = String.length s in
        for index = 0 to n - 1 do
          let st = BaseState.of_char (String.get_exn s ~index) in
          rl_arr.(index + ap) <- Rlel.append_r st rl_arr.(index + ap)
        done;
        loop position_map (ap + n) t
    | Gap { gstart; length } :: t ->
        let ap, position_map = position_and_advance gstart position_map in
        extend_from_base lp ap;
        for index = 0 to length - 1 do
          rl_arr.(index + ap) <- Rlel.append_r BaseState.Gap_ rl_arr.(index + ap)
        done;
        loop position_map (ap + length) t
  in
  loop position_map 0 allele_instr

let build_allele_and_rls ?spec ?n mp =
  let rg, pmap = initialize_base_array_and_position_map mp in
  let position_map = reduce_position_map pmap in
  (* reference goes first *)
  let rl_arr = Array.map rg ~f:(Rlel.init) in
  let init = mp.Mas_parser.reference :: [] in

  let alleles = List.sort ~cmp:(fun (a1,_) (a2,_) -> compare a1 a2) mp.Mas_parser.alt_elems in
  let add_these =
    match spec with
    | Some set -> List.filter alleles ~f:(fun (a, _) -> List.mem a ~set)
    | None ->
      match n with
      | None -> alleles
      | Some n -> List.take alleles n
  in
  let add_allele = add_alternate_allele ~position_map rg in
  List.fold_left add_these ~init ~f:(fun acc (name, allele_els) ->
    add_allele allele_els rl_arr;
    (name :: acc))
  |> fun alleles ->
      List.rev alleles
      (* Remember that the encodings in rl_arr are backward right now! *)
      , Array.map rl_arr ~f:Rlel.finish_r

let recover_allele_map al rl =
  List.map2 ~f:(fun a v -> a, v) al (Rlel.expand rl)

(* At the end of the day we need to translate various states (whether they be
   actual BaseState.t's (as when calculating an emission probability) or a
   transition (between two positions in the variant graph or k of the PHMM), to
   an index into an array where we store the value (probability of said event).
   To determine this index, we track the different states as we encounter them
   (growing) and then [finalize] the states into an order (index) so that we
   can quickly [lookup] the index.
*)
module IndexMap = struct

  type 'a growing =
    { max : int
    ; ord : 'a list     (* assignment order, reversed. *)
    }

  let empty = { max = 0; ord = []}

  let is_empty g =
    g.ord = []

  let add gm p =
    let rec loop i = function
      | []                  -> gm.max, { max = gm.max + 1; ord = p :: gm.ord }
      | hp :: _ when hp = p -> i, gm
      | _ :: ht             -> loop (i - 1) ht
    in
    loop (gm.max - 1) gm.ord

  type 'a t = 'a array

  let finalize gm =
    Array.of_list (List.rev gm.ord)

  let lookup v im =
    array_findi v im

end (* IndexMap. *)

(* If we consider the state of an allele position by position, we want to model
   two aspects:
   1. There is a normal flow of base's that describe the allele sequence. In
      this case the allele path should flow with the match/insert/delete markov
      hidden states.
   2. There isa gap in the allele, with respect to some given global alignment.
      In this case we want the to only allow a transition from the
      match/insert/delete states at position k to position k + n where n is the
      length of the gap. *)
type transition_state =
  | Regular of int      (* index into the emission/forward prob *)
  | Gapped of int * int (* Gap length, previous valid index *)

let extend_gap ~length index =
  Gapped (length + 1, index)

(* Given a run-length encoded list of base states return:
  1. an array of possible emissions: the unique non-Gap states.
  2. transition state run-length encoded. *)
let init_emissions_and_transitions rlst =
  let to_first_transition gm = function
    | BaseState.Gap_  -> invalid_arg "Can't have a gap on the first position!"
    | st              -> let index, gm = IndexMap.add gm st in
                         Regular index, gm
  in
  Rlel.fold_map rlst ~init:IndexMap.empty ~f:to_first_transition
  |> fun (trlst, gm) -> trlst, IndexMap.finalize gm

(* To model the transition, between different _reference_ states, the 'k' of
   the PHMM.*)
type previous_base =
  { state : int       (* How far back, usually -1. *)
  ; index : int       (* Into index previous output, [0,size_of_previous_transition) *)
  } [@@deriving show]
type transition_index =
  { previous  : previous_base option    (* Optional because the first row doesn't have it. *)
  ; current   : int
  } [@@deriving show]

let to_f regular =
  let open BaseState in
  fun im_pair p b ->
    match (p, b) with
    | (Regular index),          Gap_    -> extend_gap ~length:1 index, im_pair
    | (Gapped (length, index)), Gap_    -> extend_gap ~length index, im_pair
    | (Regular index),          not_gap -> regular not_gap im_pair { state = -1; index }
    | (Gapped (length, index)), not_gap -> regular not_gap im_pair { state = -length; index }

(* Given a transition state (from, the source) run-length encoded list and a
   base-state (to, destination) run-length encoded list, return:
   1. Emission array/IndexMap
   2. Transition states run-length encoded
   3. Transition IndexMap

Throw invalid_arg if t1 and t2 do not encode sequences of the same length! *)
let merge_start_positions ptsrl bsrl =
  let regular base_state (gm, trans_gm) previous =
    let current, ngm = IndexMap.add gm base_state in
    let trans_index, ntrans_gm =
      IndexMap.add trans_gm { previous = Some previous; current }
    in
    (Regular trans_index), (ngm, ntrans_gm)
  in
  let f = to_f regular in
  Rlel.fold_map2_same_length ptsrl bsrl ~init:(IndexMap.empty, IndexMap.empty) ~f
  |> fun (mgd_rlst, (gm, tgm)) ->
      if IndexMap.is_empty tgm then
        invalid_argf "Didn't identify any non-Gapped transitions!"
      else
        IndexMap.finalize gm
        , mgd_rlst
        , IndexMap.finalize tgm

 (* Given a transition state (from, the source) run-length encoded list and a
   base-state (to, destination) run-length encoded list and an IndexMap of the
   base state run-length list (previously computed): return

   1. Transition state run-length encoded.
   2. Transition indexing scheme, IndexMap.  *)
let merge ptsrl bsrl im  =
  let regular base_state trans_gm previous =
    let trans_index, ntrans_gm =
      IndexMap.add trans_gm
        { previous = Some previous
        ; current  = IndexMap.lookup base_state im  (* Indicies are known! *)
        }
    in
    (Regular trans_index), ntrans_gm
  in
  Rlel.fold_map2_same_length ptsrl bsrl ~init:IndexMap.empty ~f:(to_f regular)
  |> fun (rlst, trans_gm) -> rlst, (IndexMap.finalize trans_gm)

(* In the following "transitions" is used in 2 ways. The first, and the one in
   this type describes an indexing scheme between successives match, insert,
   and delete nodes. The second, as is more common in the TransitionMatrix of
   Phmm, are the state _transition_ probabilities underlying the Markov model.
*)

type emissions = BaseState.t array [@@deriving show]
type transitions = transition_index array [@@deriving show]
type conf =
  { read_size       : int
  ; emissions_a     : emissions array
  ; transitions_m   : transitions array array
  ; final_run_len   : transition_state Rlel.t array
  ; max_transition  : int
  }

(* TODO: For the cells where i > k, the transitions should be identical (same
   states). One potential optimization would be to not compute them. *)
let build_matrix ?(read_size=100) rlarr =
  if read_size <= 0 then
    invalid_argf "read size must be greater than 0: %d" read_size;
  let bigK = Array.length rlarr in
  let transitions_m = Array.make_matrix ~dimx:bigK ~dimy:(read_size - 1) [||] in
  let etrl = Rlel.init (Regular 1) in
  let first_run_len = Array.make bigK etrl in
  let final_run_len = Array.make bigK etrl in
  let trans_rtl_a = Array.make (read_size - 1) etrl in
  let max_transitions_ref = ref 0 in
  let update_max_transitions trans_arr =
    max_transitions_ref := max !max_transitions_ref (Array.length trans_arr);
    (*printf "max_transitions: %d\n" !max_transitions_ref *)
  in
  let emissions_a =
    match Array.to_list rlarr with
    | []     -> invalid_arg "Empty array"
    | h :: t ->
        let trl, em = init_emissions_and_transitions h in
        update_max_transitions em;  (* for first state # transitions = # emissions. *)
        let id_start i = { previous = None; current = i} in
        transitions_m.(0).(0) <- Array.init (Array.length em) ~f:id_start;
        (*printf "k: %d start  rtl total length: %d\n" k (total_length trl); *)
        first_run_len.(0) <- trl;
        let _k, _ftrl, emlst =
          List.fold_left t ~init:(1, trl, [em]) ~f:(fun (k, trl, acc) brl ->
            let em, ntrl, ti = merge_start_positions trl brl in
            first_run_len.(k) <- ntrl;
            update_max_transitions ti;
            transitions_m.(k).(0) <- ti;
            (k + 1, ntrl, em ::acc))
        in
        Array.of_list (List.rev emlst)
  in
  (* The first transitions row is special, since it has no input. *)
  final_run_len.(0) <- first_run_len.(0) ;
  for i = 1 to read_size - 2 do
    transitions_m.(0).(i) <- transitions_m.(0).(0);
    trans_rtl_a.(i) <- first_run_len.(0);
  done;
  let trans_rtl_b = Array.copy trans_rtl_a in
  let merge_fill k =
    let base_rls = rlarr.(k) in
    let ems = emissions_a.(k) in
    let source, dest =
      if k mod 2 = 1 then
        trans_rtl_a, trans_rtl_b
      else
        trans_rtl_b, trans_rtl_a
    in
    source.(0) <- first_run_len.(k-1);
    let rec loop i =
      if i = read_size - 1 then begin
        final_run_len.(k) <- dest.(i-1);
        ()
      end else begin
        let trl = source.(i-1) in
        (*printf "k: %d i: %d base length: %d rtl total length: %d\n" k i
          (total_length base_rls) (total_length trl); *)
        let ntrl, ti = merge trl base_rls ems in
        (*printf "k: %d i: %d after: %d\n" k i (total_length ntrl); *)
        update_max_transitions ti;
        transitions_m.(k).(i) <- ti;
        dest.(i) <- ntrl;
        loop (i + 1)
      end
    in
    (*printf "merge filling:k:%d \n%!" k; *)
    loop 1
  in
  (* Do not fill the first row (k = 0) since there are _no_ previous,
     all id_transitions.
     The forward pass should special case this code to return 0's or the
     rescaled insert values. *)
  for k = 1 to bigK - 1 do merge_fill k done;
  { read_size
  ; emissions_a
  ; transitions_m
  ; final_run_len
  ; max_transition = !max_transitions_ref
  }

(* Fill a configuration with the actual paths that we encounter along the
   transitions. This is a debug method to help understand the actual paths
   that the PHMM follows. *)
let fill_possibilities conf =
  let bigK = Array.length conf.transitions_m in
  let fm = Array.make_matrix ~dimx:bigK ~dimy:conf.read_size [||] in
  for k = 0 to bigK - 1 do
    fm.(k).(0) <- Array.map conf.transitions_m.(k).(0)
      ~f:(fun { current; _ } ->
            [ BaseState.to_char conf.emissions_a.(k).(current) ]);
    for i = 0 to conf.read_size - 2 do
      fm.(k).(i + 1) <-
        Array.map conf.transitions_m.(k).(i)
          ~f:(fun { previous ; current } ->
                (BaseState.to_char conf.emissions_a.(k).(current)) ::
                (match previous with
                  | None                  ->  []
                  | Some { index ; state} -> fm.(k + state).(i).(index)))
    done;
  done;
  Array.map fm ~f:(fun r ->
    Array.map r ~f:(fun r2 ->
      Array.map r2 ~f:(fun l ->
        String.of_character_list (List.rev l))))


(***** Forward Pass ***)

(* For every k there are 3 possible states. The cell stores values for each of
   those states based upon the transitions into this k. *)
type cell =
  { match_  : float array
  ; insert  : float array
  ; delete  : float array
  }

type workspace = cell array array

(* The actual workspace is large. Repeatedly allocating a new one for each read
   is wasteful. It would be better to recycle it between reads. Therefore,
   these methods operate by mutating the last cell. *)
type fwd_recurrences =
  { start   : emissions -> transitions -> char -> float -> cell -> unit
  ; middle  : workspace -> emissions -> transitions -> char -> float
                          -> i:int -> int -> cell  -> unit
  ; end_    : workspace -> int -> cell -> unit  (* Doesn't use the delete section. *)
  }

let to_match_prob reference_state_arr base base_error =
  let compare_against c = if base = c then 1. -. base_error else base_error /. 3. in
  let open BaseState in
  Array.map reference_state_arr ~f:(function
    | Unknown -> 1.   (* Treat like an 'N', anything goes.
                         TODO: But need to figure out a way to not reward this.
                               Specifically the probability shouldn't be higher than actual matches,
                               so perhaps 1-. base_error is better? .*)
    | A       -> compare_against 'A'
    | C       -> compare_against 'C'
    | G       -> compare_against 'G'
    | T       -> compare_against 'T'
    | Gap_    -> invalid_argf "Asked to match against a Gap_ state!")

(** Need a better name and destination for this module. *)
module L = struct

  let is_nan (x : float) = x <> x

  let log_z = nan

  let exp10 x = if is_nan x then 0. else 10. ** x
  let log10 = function | 0. -> nan | x -> log10 x

  (* in_lsp -> in log space *)
  let log10sum_in_lsp lx ly =
    if is_nan lx then
      ly
    else if is_nan ly then
      lx
    else if lx > ly then
      lx +. log10 (1. +. exp10 (ly -. lx))
    else
      ly +. log10 (1. +. exp10 (lx -. ly))

  let log10prod_in_lsp lx ly =
    if is_nan lx || is_nan ly then
      nan
    else
      lx +. ly

  let log10_one_minus_l lq =
    let p = 10. ** (-.lq) in
    log10 (1. -. p)

  (* The base error (qualities) are generally know. To avoid repeating the manual
    calculation (as described above) of the log quality to log (1. -. base error)
    we precompute these values.

  Weird. this seems to be slower! TODO: Why? )
  let log10_one_minus_l = function
    | -0.0 -> log_z
    | -0.1 -> -0.686825324380115454
    | -0.2 -> -0.432923433336248276
    | -0.3 -> -0.302062439928300397
    | -0.4 -> -0.220480830541908562
    | -0.5 -> -0.165088538626769726
    | -0.6 -> -0.125627577491815079
    | -0.7 -> -0.0966528953262047186
    | -0.8 -> -0.07494036743261491
    | -0.9 -> -0.058435173882679825
    | -1.0 -> -0.0457574905606751153
    | -1.1 -> -0.035944514242268806
    | -1.2 -> -0.0283047837831196247
    | -1.3 -> -0.0223306727357915694
    | -1.4 -> -0.0176431456736382448
    | -1.5 -> -0.0139554338820558448
    | -1.6 -> -0.0110483332892353306
    | -1.7 -> -0.00875292940206854816
    | -1.8 -> -0.0069382318574496421
    | -1.9 -> -0.00550215071190342936
    | -2.0 -> -0.00436480540245008826
    | -2.1 -> -0.00346349774554599588
    | -2.2 -> -0.00274889425384098815
    | -2.3 -> -0.00218210128532180967
    | -2.4 -> -0.0017324081870171721
    | -2.5 -> -0.00137553579921727916
    | -2.6 -> -0.00109227082153636758
    | -2.7 -> -0.000867397043708781281
    | -2.8 -> -0.000688856394105166097
    | -2.9 -> -0.000547088803770739681
    | -3.0 -> -0.000434511774017691684
    | -3.1 -> -0.000345109452404739883
    | -3.2 -> -0.000274107777278245887
    | -3.3 -> -0.000217717413117151276
    | -3.4 -> -0.000172930172032690271
    | -3.5 -> -0.000137357693108739246
    | -3.6 -> -0.000109103544996691523
    | -3.7 -> -8.66617872714958135e-05
    | -3.8 -> -6.88364918576594357e-05
    | -3.9 -> -5.46778777877239756e-05
    | -4.0 -> -4.34316198075056039e-05
    | -4.1 -> -3.44986070951026853e-05
    | -4.2 -> -2.7402993817532012e-0
    |    x -> (*eprintf "asked to compute log10_one_minus_l of %f\n" x; *)
              log10_one_minus_l_manual x
 *)

end (* L *)

let log10_1_3 =
  L.log10 (1. /. 3.)

(* TODO. Avoid the `float_of_int (Phred_score.to_int c) /. -10.` round trip
         and just use char's instead. *)
let to_match_prob_log reference_state_arr base base_error =
  let compare_against c =
    if base = c then
      L.log10_one_minus_l base_error
    else
      L.log10prod_in_lsp base_error log10_1_3
  in
  let open BaseState in
  Array.map reference_state_arr ~f:(function
    | Unknown -> 0.   (* Treat like an 'N', anything goes.
                         TODO: But need to figure out a way to not reward this.
                               Specifically the probability shouldn't be higher than actual matches,
                               so perhaps 1-. base_error is better? .*)
    | A       -> compare_against 'A'
    | C       -> compare_against 'C'
    | G       -> compare_against 'G'
    | T       -> compare_against 'T'
    | Gap_    -> invalid_argf "Asked to match against a Gap_ state!")


(*
let make_constants size v =
  let arr = Array.init size ~f:(fun i -> lazy (Array.make (i + 1) v)) in
  fun i -> Lazy.force arr.(i-1)

let wrap_array f arr = f (Array.length arr) *)

let generate_workspace conf read_size =
  let just_zeros n = Array.make n 0. in
  Array.mapi conf.transitions_m ~f:(fun k trans_row ->
    Array.init read_size ~f:(fun i ->
      let n =
        if i = 0 then
          Array.length trans_row.(0)
        else
          Array.length trans_row.(i-1)
      in
      { match_ = just_zeros n
      ; insert = just_zeros n
      ; delete = just_zeros n
      })),
  (* The final emissions should have the same size as the last transition
     column. Don't bother allocating the delete *)
  Array.map conf.transitions_m ~f:(fun trans_row ->
    let n = Array.length trans_row.(read_size-2) in
    { match_ = just_zeros n
    ; insert = just_zeros n
    ; delete = [||]
    })

let forward_recurrences tm ~insert_prob read_size =

  let t_s_m = tm `StartOrEnd `Match in
  let t_s_i = tm `StartOrEnd `Insert in
  let start_i = t_s_i *. insert_prob in
  let t_m_m = tm `Match `Match in
  let t_i_m = tm `Insert `Match in
  let t_d_m = tm `Delete `Match in

  let t_m_i = tm `Match `Insert in
  let t_i_i = tm `Insert `Insert in

  let t_m_d = tm `Match `Delete in
  let t_d_d = tm `Delete `Delete in

  let t_m_s = tm `Match `StartOrEnd in
  let t_i_s = tm `Insert `StartOrEnd in

  { start   = begin fun emissions transitions base base_error cell ->
                let ce = to_match_prob emissions base base_error in
                Array.iteri transitions ~f:(fun j { current; _} ->
                  cell.match_.(j) <- ce.(current) *. t_s_m;
                  cell.insert.(j) <- start_i;
                  cell.delete.(j) <- 0.);
              end
  ; middle  = begin fun fm emissions transitions base base_error ~i k cell ->
                let ce = to_match_prob emissions base base_error in
                Array.iteri transitions
                    ~f:(fun j { previous ; current } ->
                          let pm_i = fm.(k).(i-1).match_.(current) in
                          let pi_i = fm.(k).(i-1).insert.(current) in
                          cell.insert.(j) <- insert_prob *. (t_m_i *. pm_i +. t_i_i *. pi_i);
                          match previous with
                          | None                  -> (* at k = 0 -> mirror emissions. *)
                              cell.match_.(j) <- ce.(current) *. t_s_m;
                              cell.delete.(j) <- 0.
                          | Some { state; index } ->
                              let ks = k + state in
                              let pm_m = fm.(ks).(i-1).match_.(index) in
                              let pi_m = fm.(ks).(i-1).insert.(index) in
                              let pd_m = fm.(ks).(i-1).delete.(index) in
                              cell.match_.(j) <- ce.(current) *. (t_m_m *. pm_m +. t_i_m *. pi_m +. t_d_m *. pd_m);
                              let pm_d = fm.(ks).(i).match_.(index) in
                              let pd_d = fm.(ks).(i).insert.(index) in
                              cell.delete.(j) <- (t_m_d *. pm_d +. t_d_d *. pd_d));
              end
  ; end_    = begin fun fm k cell ->
                let fc = fm.(k).(read_size-2) in
                Array.iteri fc.match_ ~f:(fun j m ->
                  cell.match_.(j) <- m *. t_m_s;
                  cell.insert.(j) <- fc.insert.(j) *. t_i_s;
                  (* delete is empty!*))
              end
  }

let forward_recurrences_log tm ~insert_prob read_size =

  let open L in
  let t_s_m = log10 (tm `StartOrEnd `Match) in
  let t_s_i = log10 (tm `StartOrEnd `Insert) in
  let t_m_m = log10 (tm `Match `Match) in
  let t_i_m = log10 (tm `Insert `Match) in
  let t_d_m = log10 (tm `Delete `Match) in

  let t_m_i = log10 (tm `Match `Insert) in
  let t_i_i = log10 (tm `Insert `Insert) in

  let t_m_d = log10 (tm `Match `Delete) in
  let t_d_d = log10 (tm `Delete `Delete) in

  let t_m_s = log10 (tm `Match `StartOrEnd) in
  let t_i_s = log10 (tm `Insert `StartOrEnd) in

  let linsert_prob = log10 insert_prob in
  let start_i = log10prod_in_lsp t_s_i linsert_prob in

  { start   = begin fun emissions transitions base base_error cell ->
                let ce = to_match_prob_log emissions base base_error in
                Array.iteri transitions ~f:(fun j { current; _} ->
                  cell.match_.(j) <- log10prod_in_lsp ce.(current) t_s_m;
                  cell.insert.(j) <- start_i;
                  cell.delete.(j) <- log_z);
              end
  ; middle  = begin fun fm emissions transitions base base_error ~i k cell ->
                let ce = to_match_prob_log emissions base base_error in
                Array.iteri transitions
                    ~f:(fun j { previous ; current } ->
                          let pm_i = fm.(k).(i-1).match_.(current) in
                          let pi_i = fm.(k).(i-1).insert.(current) in
                          cell.insert.(j) <-
                            log10prod_in_lsp linsert_prob
                              (log10sum_in_lsp
                                (log10prod_in_lsp t_m_i pm_i)
                                (log10prod_in_lsp t_i_i pi_i));
                          match previous with
                          | None                  -> (* at k = 0 -> mirror emissions. *)
                              cell.match_.(j) <- log10prod_in_lsp ce.(current) t_s_m;
                              cell.delete.(j) <- log_z
                          | Some { state; index } ->
                              let ks = k + state in
                              let pm_m = fm.(ks).(i-1).match_.(index) in
                              let pi_m = fm.(ks).(i-1).insert.(index) in
                              let pd_m = fm.(ks).(i-1).delete.(index) in
                              cell.match_.(j) <-
                                log10prod_in_lsp ce.(current)
                                  (log10sum_in_lsp
                                    (log10sum_in_lsp
                                      (log10prod_in_lsp t_m_m pm_m)
                                      (log10prod_in_lsp t_i_m pi_m))
                                    (log10prod_in_lsp t_d_m pd_m));
                              let pm_d = fm.(ks).(i).match_.(index) in
                              let pd_d = fm.(ks).(i).insert.(index) in
                              cell.delete.(j) <-
                                log10sum_in_lsp
                                  (log10prod_in_lsp t_m_d pm_d)
                                  (log10prod_in_lsp t_d_d pd_d))
              end
  ; end_    = begin fun fm k cell ->
                let fc = fm.(k).(read_size-2) in
                Array.iteri fc.match_ ~f:(fun j m ->
                  cell.match_.(j) <- log10prod_in_lsp m t_m_s;
                  cell.insert.(j) <- log10prod_in_lsp fc.insert.(j) t_i_s;
                  (* delete is empty!*))
              end
  }


type t =
  { conf        : conf
  ; align_date  : string
  ; alleles     : string list
  ; merge_map   : (string * string) list
  }

let construct input selectors read_size =
  Alleles.Input.construct input >>= fun (mp, merge_map) ->
    let nalt_elems =
      mp.Mas_parser.alt_elems
      |> List.sort ~cmp:(fun (n1, _) (n2, _) -> Alleles.compare n1 n2)
      |> Alleles.Selection.apply_to_assoc selectors
    in
    let alleles, rlarr =
      build_allele_and_rls { mp with Mas_parser.alt_elems = nalt_elems}
    in
    let conf = build_matrix ~read_size rlarr in
    Ok { conf ; align_date = mp.Mas_parser.align_date ; alleles ; merge_map }

let forward_pass ?(logspace=true) conf read_size =
  let read_size =
    if read_size > conf.read_size then begin
      eprintf "requested read size %d greater than configuration %d, will use smaller\n"
        read_size conf.read_size;
        conf.read_size
    end else
      read_size
  in
  let sm, e = generate_workspace conf read_size in
  let bigK = Array.length conf.transitions_m in
  let tm = Phmm.TransitionMatrix.init ~ref_length:bigK read_size in
  let insert_prob = 0.25 in
  let recurrences =
    (if logspace then forward_recurrences_log else forward_recurrences)
      tm ~insert_prob read_size
  in
  fun read read_prob ->
    (* Fill in start. *)
    for k = 0 to bigK - 1 do
      recurrences.start conf.emissions_a.(k) conf.transitions_m.(k).(0)
        (String.get_exn read 0) read_prob.(0) sm.(k).(0);
    done;
    (* Fill in middle  *)
    for i = 1 to read_size - 1 do
      let base = String.get_exn read i in
      let base_prob = read_prob.(i) in
      for k = 0 to bigK - 1 do
        recurrences.middle sm conf.emissions_a.(k) conf.transitions_m.(k).(i-1)
          base base_prob ~i k sm.(k).(i)
      done;
    done;
    for k = 0 to bigK - 1 do
      recurrences.end_ sm k e.(k)
    done;
    e

let to_allele_arr ~logspace fe rtl =
  let len = Rlel.total_length rtl.(0) in
  let ret = Array.make len 0. in
  Array.iteri rtl ~f:(fun k rtl ->
      (* value is an index into match/insert probs. *)
      Rlel.fill_array ret rtl ~f:(function
        | Regular i ->
            if logspace then
              L.log10sum_in_lsp fe.(k).match_.(i) fe.(k).insert.(i)
            else
              fe.(k).match_.(i) +. fe.(k).insert.(i)
        | Gapped _  -> 0. (* an allele can end in a gap after read_length. *)));
  ret
