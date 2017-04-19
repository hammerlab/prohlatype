
open Common

type base_state =
  | Unknown         (* Treat like an 'N', anything goes. *)
  | Gap_            (* Nothing goes! *)
  | A
  | C
  | G
  | T
  [@@deriving show]

let char_to_base_state = function
  | 'A' -> A
  | 'C' -> C
  | 'G' -> G
  | 'T' -> T
  |  c  -> invalid_argf "Unsupported base: %c" c

let base_state_to_char = function
  | Unknown -> 'N'
  | Gap_    -> '_'
  | A       -> 'A'
  | C       -> 'C'
  | G       -> 'G'
  | T       -> 'T'

let init_from_masp_ref =
  let sequence_to_base_states_array s =
    Array.init (String.length s) ~f:(fun index ->
      char_to_base_state (String.get_exn s ~index))
  in
  let maybe_add_unknown last_known_mas_pos p pacc acc pos =
    let d = pos - last_known_mas_pos in
    if d > 0 then
      last_known_mas_pos + d
      , p + d
      , (last_known_mas_pos, p) :: pacc
      , Array.make d Unknown :: acc
    else
      last_known_mas_pos, p, pacc, acc
  in
  let open Mas_parser in
  let rec loop lkmp p pacc acc = function
    | []                          -> List.rev acc, List.rev ((lkmp, p) :: pacc)
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
                                        (Array.make length Gap_ :: acc) t
  in
  function
  | Start s :: t -> loop s 0 [] [] t
  | e :: _       -> invalid_argf "Not at Start : %s" (al_el_to_string e)
  | []           -> invalid_argf "Empty"

type 'a run_length_pair =
  { value   : 'a
  ; length  : int (* > 0 *)
  }
and 'a rl_list =
  { hd  : 'a run_length_pair
  ; tl  : 'a run_length_pair list
  }
  [@@deriving show]

let total_rl_length { hd; tl} =
  hd.length + List.fold_left tl ~init:0 ~f:(fun s { length; _} -> s + length)

let expand_rl { hd = {value; length}; tl} =
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

let reduce_position_map = function
  | []  -> invalid_arg "empty"
  | (p, a)  :: t ->
      List.fold_left t ~init:(a - p, [p,a])
        ~f:(fun (d, acc) (p,a) ->
            let d1 = a - p in
            if d <> d1 then (d1, (p,a) ::acc) else (d, acc))
      |> snd
      |> List.rev

let to_position_map = function
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

let init_rl value =
  { hd = { value; length = 1}
  ; tl = []
  }

let append_to_rl v = function
  | { hd = { value; length} ; tl = [] } when v = value ->
      { hd = { value ; length = length + 1}; tl = []}
  | { hd ; tl = { value; length } :: t } when v = value ->
      { hd; tl = { value; length = length + 1} :: t }
  | { hd ; tl } ->
      { hd; tl = { value = v; length = 1} :: tl }

let rec to_array_pos sp splst =
  match splst with
  | []                                                -> invalid_argf "reached end of sequence map: %d" sp
  | (p1, o1) :: (p2, _) :: _ when p1 <= sp && sp < p2 -> sp - p1 + o1, splst
  | (p1, o1) :: []           when p1 <= sp            -> sp - p1 + o1, splst
  | h :: t                                            -> to_array_pos sp t

let extend_run_length_encoded rmap allele_mas_instr ref_pos_arr rl_arr =
  let last_pos = Array.length rl_arr in
  let extend_to start lp v =
    for i = start to lp - 1 do
      rl_arr.(i) <- append_to_rl v rl_arr.(i)
    done
  in
  let extend_from_base start lp =
    for i = start to lp - 1 do
      rl_arr.(i) <- append_to_rl ref_pos_arr.(i) rl_arr.(i)
    done
  in
  let rec loop rmap lp =
    let open Mas_parser in function
    | []                          -> ()
    | Start p :: t                ->
        let ap, rmap = to_array_pos p rmap in
        extend_to lp ap Unknown;
        loop rmap ap t
    | End p :: t                  ->
        let ap, rmap = to_array_pos p rmap in
        extend_from_base lp ap;
        extend_to ap last_pos Unknown;    (* Preemptively set the rest to Unknown. *)
        loop rmap ap t
    | Boundary { pos; _ } :: t    ->
        let ap, rmap = to_array_pos pos rmap in
        extend_from_base lp ap;
        loop rmap ap t
    | Sequence { start; s } :: t  ->
        let ap, rmap = to_array_pos start rmap in
        extend_from_base lp ap;
        let n = String.length s in
        for index = 0 to n - 1 do
          let st = char_to_base_state (String.get_exn s ~index) in
          rl_arr.(index + ap) <- append_to_rl st rl_arr.(index + ap)
        done;
        loop rmap (ap + n) t
    | Gap { gstart; length } :: t ->
        let ap, rmap = to_array_pos gstart rmap in
        extend_from_base lp ap;
        for index = 0 to length - 1 do
          rl_arr.(index + ap) <- append_to_rl Gap_ rl_arr.(index + ap)
        done;
        loop rmap (ap + length) t
  in
  loop rmap 0 allele_mas_instr


let build_allele_and_rls ?spec ?n mp =
  let rg_lst, pmap = init_from_masp_ref mp.Mas_parser.ref_elems in
  let rpmap = reduce_position_map pmap in
  let rg = Array.concat rg_lst in
  (* reference goes first *)
  let rl_arr = Array.map rg ~f:(init_rl) in
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
  List.fold_left add_these ~init ~f:(fun acc (name, allele_els) ->
    extend_run_length_encoded rpmap allele_els rg rl_arr;
    (name :: acc))
  |> fun alleles ->
      List.rev alleles
      (* Remember that the encodings in rl_arr are backward right now! *)
      , Array.map rl_arr ~f:(fun {hd; tl} -> {hd ; tl = List.rev tl })

let recover_allele_map al rl =
  List.map2 ~f:(fun a v -> a, v) al (expand_rl rl)

type 'a growing_map =
  { max : int
  ; ord : 'a list     (* assignment order, reversed. *)
  }

let empty_growing_map = { max = 0; ord = []}

let add_to_growing_map gm p =
  let rec loop i = function
    | []                  -> gm.max, { max = gm.max + 1; ord = p :: gm.ord }
    | hp :: _ when hp = p -> i, gm
    | _ :: ht             -> loop (i - 1) ht
  in
  loop (gm.max - 1) gm.ord

type 'a index_map = 'a array

let growing_to_index gm =
  Array.of_list (List.rev gm.ord)

let array_findi v a =
  let n = Array.length a in
  let rec loop i =
    if i >= n then raise Not_found
    else if a.(i) = v then i
    else loop (i + 1)
  in
  loop 0

type transition_state =
  | Regular of int      (* index into the emission/forward prob *)
  | Gapped of int * int (* Gap length, previous valid index *)

let to_emission_im gm =
  Array.of_list (List.rev gm.ord)

(* Given a run-length encoded list of base states return:
  1. an array of possible emissions. (The emissions are just the unique
   non-Gap states.
  2. transition state run-length encoded list *)
let init_emission_mp_and_rl {hd ; tl} =
  let to_value gm = function
    | Gap_   -> invalid_arg "Can't have a gap on the first position!"
    | st     -> let index, gm = add_to_growing_map gm st in
                gm, Regular index
  in
  let gm, nvalue = to_value empty_growing_map hd.value in
  let nhd = { hd with value = nvalue } in
  let final_gm, ntl =
    List.fold_left tl ~init:(gm, []) ~f:(fun (gm, acc) rl ->
      let gm, nvalue = to_value gm rl.value in
      gm, { rl with value = nvalue } :: acc)
  in
  let nrlst = { hd = nhd; tl = List.rev ntl} in
  let emission_im = to_emission_im final_gm in
  emission_im, nrlst

let align_rtl c1 c2 l1 l2 =
  if c1.length < c2.length then
    { value = (c1.value, c2.value); length = c1.length }
    , l1
    , { c2 with length = c2.length - c1.length } :: l2
  else if c1.length > c2.length then
    { value = (c1.value, c2.value); length = c2.length }
    , { c1 with length = c1.length - c2.length } :: l1
    , l2
  else (* c1.length = c2.length *)
    { value = (c1.value, c2.value); length = c2.length }
    , l1
    , l2

type previous_base =
  { state : int       (* How far back. *)
  ; index : int       (* Into index previous output! *)
  } [@@deriving show]
type transition_index =
  { previous  : previous_base option
  ; current   : int
  } [@@deriving show]

let id_transition i = { previous = None; current = i}

(* Given a transition state run-length encoded list and a base-state
   run-length encoded list, return:
   1. Emission array
   2. Transition state run-length encoded
   3. Transition indexing scheme.

Throw invalid_arg if t1 and t2 do not encode sequences of the same length! *)
let merge_start_positions ptsrl { hd; tl } =
  let merge_one gm2 trans_gm merged =
    let v1, v2 = merged.value in
    match v2 with
    | Gap_  ->  (* Entering a gap *)
        let new_state =
          match v1 with              (* since 1 is just the previous row. *)
          | Regular index          -> Gapped (2, index)
          | Gapped (length, index) -> Gapped (length + 1, index)
        in
        { value = new_state; length = merged.length}, gm2, trans_gm
    | b2    ->
        let index2, nim2 = add_to_growing_map gm2 b2 in
        let previous =
          match v1 with
          | Regular index          -> { state = -1; index }
          | Gapped (length, index) -> { state = -length; index }
        in
        let trans_index, ntrans_gm =
          add_to_growing_map trans_gm
            { previous = Some previous
            ; current  = index2
            }
        in
        { value = Regular trans_index; length = merged.length}, nim2, ntrans_gm
  in
  let merge_me, l1, l2 = align_rtl ptsrl.hd hd ptsrl.tl tl in
  let mhd, gm, tgm = merge_one empty_growing_map empty_growing_map merge_me in
  let rec loop gm tgm acc = function
    | [], []             -> if tgm.ord = [] then
                              invalid_argf "Didn't identify any non-Gapped transitions!"
                            else
                              (to_emission_im gm)
                              , { hd = mhd; tl = List.rev acc }
                              , (to_emission_im tgm)
    | h1 :: t1, h2 :: t2 -> let merge_me, l1, l2 = align_rtl h1 h2 t1 t2 in
                            let code, gm2, tgm2 = merge_one gm tgm merge_me in
                            loop gm2 tgm2 (code :: acc) (l1, l2)
    | _,        _        -> invalid_arg "different lengths"
  in
  loop gm tgm [] (l1, l2)

(*
let rl1 = encode [A; Unknown; Unknown; A; A;    C] ;;
let rl2 = encode [C; Unknown; Unknown; C; Gap_; C] ;;
let rl3 = encode [G; Unknown; G;       G; Gap_; C] ;;
let rl4 = encode [T; Unknown; Unknown; T; A;    A] ;;

let em1, rlp1 = init_emission_mp_and_rl rl1 ;;
let em2, rlp2, tm2 = merge_old rlp1 rl2 ;;
let em3, rlp3, tm3 = merge_old rlp2 rl3 ;;
let em4, rlp4, tm4 = merge_old rlp3 rl4 ;;
*)


let merge ptsrl { hd; tl } im =
  let merge_one trans_gm merged =
    let v1, v2 = merged.value in
    match v2 with
    | Gap_  ->  (* Entering a gap *)
        let new_state =
          match v1 with
          | Regular index          -> Gapped (2, index)
          | Gapped (length, index) -> Gapped (length + 1, index)
        in
        trans_gm, { value = new_state; length = merged.length}
    | b2    ->
        let index2 = array_findi b2 im in
        let previous =
          match v1 with
          | Regular index          -> { state = -1; index }
          | Gapped (length, index) -> { state = -length; index }
        in
        let trans_index, ntrans_gm =
          add_to_growing_map trans_gm
            { previous = Some previous
            ; current  = index2
            }
        in
        ntrans_gm, { value = Regular trans_index; length = merged.length}
  in
  let merge_me, l1, l2 = align_rtl ptsrl.hd hd ptsrl.tl tl in
  let tgm, mhd = merge_one empty_growing_map merge_me in
  let rec loop tgm acc = function
    | [], []             -> { hd = mhd; tl = List.rev acc }
                            , (to_emission_im tgm)
    | h1 :: t1, h2 :: t2 -> let merge_me, l1, l2 = align_rtl h1 h2 t1 t2 in
                            let tgm2, code = merge_one tgm merge_me in
                            loop tgm2 (code :: acc) (l1, l2)
    | _,        _        -> invalid_arg "different lengths"
  in
  loop tgm [] (l1, l2)

(* There "transitions" is used in 2 ways. The first, and the one in this type
   describes an indexing scheme between successives match,insert,and delete
   nodes. The second, as is more common in the TransitionMatrix of Phmm, are
   the state _transition_ probabilities underlying the Markov model. *)
type transitions = transition_index array [@@deriving show]
type conf =
  { read_size       : int
  ; emissions       : base_state array array
  ; transitions     : transitions array array
  ; final_run_len   : transition_state rl_list array
  ; max_transition  : int
  }

let build_matrix ?(read_size=100) rlarr =
  if read_size <= 0 then
    invalid_argf "read size must be greater than 0: %d" read_size;
  let bigK = Array.length rlarr in
  let transitions = Array.make_matrix ~dimx:bigK ~dimy:(read_size - 1) [||] in
  let etrl = init_rl (Regular 1) in
  let rtl_trans_matrix = Array.make_matrix ~dimx:bigK ~dimy:read_size etrl in
  let max_transitions_ref = ref 0 in
  let update_max_transitions trans_arr =
    max_transitions_ref := max !max_transitions_ref (Array.length trans_arr);
    (*printf "max_transitions: %d\n" !max_transitions_ref *)
  in
  let emissions =
    match Array.to_list rlarr with
    | []     -> invalid_arg "Empty array"
    | h :: t ->
        let em, trl = init_emission_mp_and_rl h in
        update_max_transitions em;  (* for first state # transitions = # emissions. *)
        transitions.(0).(0) <- Array.init (Array.length em) ~f:(id_transition);
        (*printf "k: %d start  rtl total length: %d\n" k (total_rl_length trl); *)
        rtl_trans_matrix.(0).(0) <- trl;
        let _k, _ftrl, emlst =
          List.fold_left t ~init:(1, trl, [em]) ~f:(fun (k, trl, acc) brl ->
            let em, ntrl, ti = merge_start_positions trl brl in
            rtl_trans_matrix.(k).(0) <- ntrl;
            update_max_transitions ti;
            transitions.(k).(0) <- ti;
            (k + 1, ntrl, em ::acc))
        in
        Array.of_list (List.rev emlst)
  in
  for i = 1 to read_size - 2 do
    transitions.(0).(i) <- transitions.(0).(0);
    rtl_trans_matrix.(0).(i) <- rtl_trans_matrix.(0).(0)
  done;
  let merge_fill k =
    let base_rls = rlarr.(k) in
    let ems = emissions.(k) in
    let rec loop i =
      if i = read_size - 1 then
        ()
      (*else if i > k then begin
        transitions.(k).(i) <- transitions.(k).(i-1);
        loop (i + 1)
      end *)else begin
        (* The rtl_trans_matrix is off by one wrt transitions and we're only
           storing them up to k, again off by one :) *)
        let trl = rtl_trans_matrix.(k-1).(i-1) in
        (*printf "k: %d i: %d base length: %d rtl total length: %d\n" k i
          (total_rl_length base_rls) (total_rl_length trl); *)
        let ntrl, ti = merge trl base_rls ems in
        (*printf "k: %d i: %d after: %d\n" k i (total_rl_length ntrl); *)
        update_max_transitions ti;
        transitions.(k).(i) <- ti;
        rtl_trans_matrix.(k).(i) <- ntrl;
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
  let final_run_len =
    Array.mapi rtl_trans_matrix
      ~f:(fun k rtl_row -> rtl_row.(read_size - 1))
  in
  { read_size
  ; emissions
  ; transitions
  ; final_run_len
  ; max_transition = !max_transitions_ref
  }, rtl_trans_matrix

let fill_possibilities conf =
  let bigK = Array.length conf.transitions in
  let fm = Array.make_matrix ~dimx:bigK ~dimy:conf.read_size [||] in
  for k = 0 to bigK - 1 do
    fm.(k).(0) <- Array.map conf.transitions.(k).(0)
      ~f:(fun { current; _ } ->
            [ base_state_to_char conf.emissions.(k).(current) ]);
    for i = 0 to conf.read_size - 2 do
      fm.(k).(i + 1) <-
        Array.map conf.transitions.(k).(i)
          ~f:(fun { previous ; current } ->
                (base_state_to_char conf.emissions.(k).(current)) ::
                (match previous with
                  | None                  ->  []
                  | Some { index ; state} -> fm.(k + state).(i).(index)))
    done;
  done;
  Array.map fm ~f:(fun r ->
    Array.map r ~f:(fun r2 ->
      Array.map r2 ~f:(fun l -> String.of_character_list (List.rev l))))

let mp = Mas_parser.from_file (to_alignment_file "A_gen") ;;

let spec =
  [ "A*01:01:01:01"
  ; "A*01:01:01:02N"
  ; "A*01:01:38L"
  ; "A*11:01:01:01"
  ; "A*02:01:01:01"
  ; "A*02:264"
  ]

(*
let alleles, rlarr = build_allele_and_rls mp ;;
let rlarr_sub = rlarr (*Array.sub rlarr ~pos:0 ~len:1000 *) ;;
let conf, tlarr = build_matrix ~read_size:10 rlarr_sub ;;
let ffa = fill_possibilities  conf;;
*)

(***** Forward Pass ***)

(* For every k there are 3 possible states. The cell stores values for each of
   those states based upon the transitions into this k. *)
type cell =
  { match_  : float array
  ; insert  : float array
  ; delete  : float array
  }

type workspace = cell array array

type emissions = base_state array [@@deriving show]

type fwd_recurrences =
  { start   : emissions -> transitions -> char -> float -> cell
  ; middle  : workspace -> emissions -> transitions -> char -> float -> i:int -> int -> cell
(*; end_    : TransitionMatrix.t -> float array array -> i:int -> int -> float * float * float
  *)
  }

let to_match_prob reference_state_arr base base_error =
  let of_char c = if base = c then 1. -. base_error else base_error /. 3. in
  Array.map reference_state_arr ~f:(function
    | Unknown -> 1.   (* Treat like an 'N', anything goes.
                         TODO: But need to figure out a way to not reward this.
                               Specifically the probability shouldn't be higher than actual matches,
                               so perhaps 1-. base_error is better? .*)
    | A       -> of_char 'A'
    | C       -> of_char 'C'
    | G       -> of_char 'G'
    | T       -> of_char 'T'
    | Gap_    -> assert false (*0. NEVER match against a gap! *))

(* let to_constant_prob reference_state_arr v =
  Array.map reference_state_arr ~f:(fun _ -> v)

let scale v =
  Array.map ~f:(fun x -> x *. v) *)

let make_constants size v =
  let arr = Array.init size ~f:(fun i -> Array.make (i + 1) v) in
  fun i ->
    if i < 1 || i > size then
      invalid_argf "More emission states: %d than possible! [1,%d]" i size
    else
      arr.(i-1)

let wrap_array f arr = f (Array.length arr)

let forward_recurrences tm ~insert_prob max_transition_size =

  let t_s_m = tm `StartOrEnd `Match in
  let t_s_i = tm `StartOrEnd `Insert in

  let t_m_m = tm `Match `Match in
  let t_i_m = tm `Insert `Match in
  let t_d_m = tm `Delete `Match in

  let t_m_i = tm `Match `Insert in
  let t_i_i = tm `Insert `Insert in

  let t_m_d = tm `Match `Delete in
  let t_d_d = tm `Delete `Delete in

  let just_zeros a = wrap_array (make_constants max_transition_size 0.) a in
  let init_deletes a = just_zeros a in
  let init_inserts a = wrap_array (make_constants max_transition_size (t_s_i *. insert_prob)) a in
  { start   = begin fun emissions transitions base base_error ->
                let ce = to_match_prob emissions base base_error in
                { match_ = Array.map transitions ~f:(fun { current; _} -> ce.(current) *. t_s_m)
                ; insert = init_inserts transitions
                ; delete = init_deletes transitions
                }
              end
  ; middle  = begin fun fm emissions transitions base base_error ~i k ->
                let ce = to_match_prob emissions base base_error in
                let match_ =
                  Array.map transitions
                      ~f:(fun { previous ; current } ->
                            match previous with
                            | None                  -> (* at k = 0 -> mirror emissions. *)
                                ce.(current) *. t_s_m
                            | Some { state; index } ->
                                let pm = fm.(k + state).(i-1).match_.(index) in
                                let pi = fm.(k + state).(i-1).insert.(index) in
                                let pd = fm.(k + state).(i-1).delete.(index) in
                                ce.(current) *. (t_m_m *. pm +. t_i_m *. pi +. t_d_m *. pd))
                in
                let insert =
                  Array.map transitions
                      ~f:(fun { current; _ } ->
                            let pm = fm.(k).(i-1).match_.(current) in
                            let pi = fm.(k).(i-1).insert.(current) in
                            insert_prob *. (t_m_i *. pm +. t_i_i *. pi  ))
                in
                let delete =
                  Array.map transitions
                      ~f:(function
                            | { previous = None; _ } -> (* A gap start -> mirror emissions *)
                                0.
                            | { previous = Some { state; index } ; _} ->
                                let pm = fm.(k + state).(i).match_.(index) in
                                let pd = fm.(k + state).(i).insert.(index) in
                                (t_m_d *. pm +. t_d_d *. pd))
                in
                { match_; insert; delete }
              end
  }

let forward_pass conf read read_prob =
  if conf.read_size <> String.length read then
    invalid_argf "read length %d doesn't match configuration read_size %d"
      (String.length read) conf.read_size
  else if conf.read_size <> Array.length read_prob then
    invalid_argf "read probability length %d doesn't match configuration read_size %d"
      (Array.length read_prob) conf.read_size
  else
    let bigK = Array.length conf.transitions in
    let ecell = { match_ = [||]; insert = [||]; delete = [||] } in
    let fm = Array.make_matrix ~dimx:bigK  ~dimy:conf.read_size ecell in
    let tm = Phmm.TransitionMatrix.init ~ref_length:bigK conf.read_size in
    let insert_prob = 0.25 in
    let recurrences = forward_recurrences tm ~insert_prob conf.max_transition in
    (* Fill in start. *)
    for k = 0 to bigK - 1 do
      fm.(k).(0) <-
        recurrences.start conf.emissions.(k) conf.transitions.(k).(0)
          (String.get_exn read 0) read_prob.(0)
    done;
    (* Fill in middle  *)
    for i = 1 to conf.read_size - 1 do
      let base = String.get_exn read i in
      let base_prob = read_prob.(i) in
      for k = 0 to bigK - 1 do
        fm.(k).(i) <-
          recurrences.middle fm conf.emissions.(k) conf.transitions.(k).(i-1)
            base base_prob ~i k
      done;
    done;
    fm

let build ?read_size ?len ?spec mp =
  let alleles, rlarr = build_allele_and_rls ?spec mp in
  let rlarr = match len with | None -> rlarr | Some len -> Array.sub rlarr ~pos:0 ~len in
  let conf, _tlr = build_matrix ?read_size rlarr in
  alleles,
  rlarr,
  conf,
  forward_pass conf

let array_same arr =
  let f = arr.(0) in
  Array.fold_left arr ~init:true ~f:(fun b v -> b && v = f)

