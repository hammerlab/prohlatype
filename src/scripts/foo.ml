
open Common

type base_state =
  | Unknown         (* Treat like an 'N', anything goes. *)
  | Gap_            (* Nothing goes! *)
  | A
  | C
  | G
  | T

let char_to_base_state = function
  | 'A' -> A
  | 'C' -> C
  | 'G' -> G
  | 'T' -> T
  |  c  -> invalid_argf "Unsupported base: %c" c

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

let total_rl_length { hd; tl} =
  hd.length + List.fold_left tl ~init:0 ~f:(fun s { length; _} -> s + length)

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

let mp =
  Mas_parser.from_file (to_alignment_file "A_gen") ;;

let build ?n mp =
  let rg_lst, pmap = init_from_masp_ref mp.Mas_parser.ref_elems in
  let rpmap = reduce_position_map pmap in
  let rg = Array.concat rg_lst in
  (* reference goes first *)
  let rl_arr = Array.map rg ~f:(init_rl) in
  let init = mp.Mas_parser.reference :: [] in

  let alleles = List.sort ~cmp:(fun (a1,_) (a2,_) -> compare a1 a2) mp.Mas_parser.alt_elems in
  let add_these =
    match n with
    | None -> alleles
    | Some n -> List.take alleles n
  in
  (* Remember that the encodings in rl_arr are backward right now! *)
  List.fold_left add_these ~init ~f:(fun acc (name, allele_els) ->
    extend_run_length_encoded rpmap allele_els rg rl_arr;
    (name :: acc))
  |> fun alleles -> List.rev alleles, rl_arr

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

type 'a index_map = 
  'a array

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
  match hd.value with
  | Gap_   -> invalid_arg "How can the first position be in a gap!"
  | st     -> let index, gm = add_to_growing_map empty_growing_map st in
              let nhd = { hd with value = Regular index } in
              let ntl, final_gm =
                List.fold_left tl ~init:([], gm) ~f:(fun (acc, gm) rl ->
                  let index, gm = add_to_growing_map gm rl.value in
                  { rl with value = Regular index } :: acc, gm)
              in
              let nrlst = { hd = nhd; tl = List.rev ntl} in
              let emission_im = to_emission_im final_gm in
              emission_im, nrlst

type transition_indeex =
  { prev_em_index : int      (* Into index previous emission. *)
  ; prev_em_state : int      (* How far back. *)
  ; cur_em_index  : int
  }
(* Given a transition state run-length encoded list and a base-state
   run-length encoded list, return:
   1. Emission array
   2. Transition state run-length encoded
   3. Transition indexing scheme.

Throw invalid_arg if t1 and t2 do not encode sequences of the same length! *)
let merge ptsrl { hd; tl } =
  let align_lengths c1 c2 l1 l2 =
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
  in
  let merge_one gm2 trans_gm merged =
    let v1, v2 = merged.value in
    match v2 with
    | Gap_  ->  (* Entering a gap *)
        let new_state =
          match v1 with
          | Regular index       -> Gapped (1, index)
          | Gapped (length, pi) -> Gapped (length + 1, pi)
        in
        { value = new_state; length = merged.length}, gm2, trans_gm
    | b2    ->
        let index2, nim2 = add_to_growing_map gm2 b2 in
        let prev_index, gap_length =
          match v1 with
          | Regular index       -> index, -1
          | Gapped (length, pi) -> pi, -length
        in
        let trans_index, ntrans_gm =
          add_to_growing_map trans_gm
            { prev_em_index = prev_index
            ; prev_em_state = gap_length
            ; cur_em_index  = index2
            }
        in
        { value = Regular trans_index; length = merged.length}, nim2, ntrans_gm
  in 
  let merge_me, l1, l2 = align_lengths ptsrl.hd hd ptsrl.tl tl in
  let mhd, gm, tgm = merge_one empty_growing_map empty_growing_map merge_me in
  let rec loop gm tgm acc = function
    | [], []             -> (to_emission_im gm)
                            , { hd = mhd; tl = List.rev acc }
                            , (to_emission_im tgm)
    | h1 :: t1, h2 :: t2 -> let merge_me, l1, l2 = align_lengths h1 h2 t1 t2 in
                            let code, gm2, tgm2 = merge_one gm tgm merge_me in
                            loop gm2 tgm2 (code :: acc) (l1, l2)
    | _,        _        -> invalid_arg "different lengths"
  in
  loop gm tgm [] (l1, l2)

let rl1 = encode [A; Unknown; Unknown; A; A;    C] ;;
let rl2 = encode [C; Unknown; Unknown; C; Gap_; C] ;;
let rl3 = encode [G; Unknown; G;       G; Gap_; C] ;;
let rl4 = encode [T; Unknown; Unknown; T; A;    A] ;;

let em1, rlp1 = init_emission_mp_and_rl rl1 ;;
let em2, rlp2, tm2 = merge rlp1 rl2 ;;
let em3, rlp3, tm3 = merge rlp2 rl3 ;;
let em4, rlp4, tm4 = merge rlp3 rl4 ;;

