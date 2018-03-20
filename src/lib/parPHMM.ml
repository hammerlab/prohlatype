(* Parametric Profile Hidden Markov Model.

   "Parameterize" the match/insert/delete states by different alleles.
*)

open Util
open Biology

(* Add a short hand for partition maps and some helper functions for to
   make the code easier. *)
module Pm = Partition_map

(* In all cases size is really just the number of alleles. *)
let pm_init_all ~number_alleles v =
  Pm.init_all_a ~size:number_alleles v

(* Construction
  1. From MSA.Parser.result -> Array of Base.t option 's.
      None if allele is in a gap.
    a. Figure out the Base.t of reference sequence and a position map
       (position in MSA.Parser.result to index into final array)
       [initialize_base_array_and_position_map].
    b. Start Partition_map descending lists so that we can extend them with
       each alternate allele.
  2. Extend the array with each new allele.

The mapping (position into base_state array) is then recovered by iterating
over this position map as we move along MSA.Parser.alignment_element's for
the alleles.
*)

type gapped_bases = Base.t option          (* None -> the allele is in a gap. *)
type base_emissions = (Pm.ascending, gapped_bases) Pm.t

type position_map = (MSA.position * int) list

let some x = Some x

let initialize_base_array_and_position_map ref_elems =
  let wrap_char c = Base.of_char c |> some in
  let sequence_to_base_states_array s =
    String.to_character_list s
    |> List.map ~f:wrap_char
    |> Array.of_list
  in
  let gap_to_base_states_array len = Array.make len None in
  let open MSA in
  let prefix = "initialize_base_array_and_position_map: " in
  let ia fmt = invalid_argf ~prefix fmt in
  let rec loop lkmp p pp pacc acc = function
    | End _pos :: []              -> Array.concat (List.rev acc),
                                     List.rev ((lkmp, p) :: pacc)
    | Start pos :: t              -> ia "second start: %d" pos
    | End pos :: t                -> ia "end with more %d" pos
    | []                          -> ia "list before End"
    | Boundary { pos; _ } :: t    -> loop pos (p + pos - lkmp) pp ((lkmp, p) :: pacc) acc t
    | Sequence { s; start } :: t  -> let l = String.length s in
                                     loop (start + l) (p + l) (-1) ((lkmp, p) :: pacc)
                                        (sequence_to_base_states_array s :: acc) t
    | Gap { length; gstart } :: t -> loop (gstart + length) (p + length) (-length-1) ((lkmp, p) :: pacc)
                                        (gap_to_base_states_array length :: acc) t
  in
  match ref_elems with
  | Start s :: t -> loop s 0 min_int [] [] t
  | e :: _       -> ia "reference not at Start %s" (al_el_to_string e)
  | []           -> ia "Empty reference sequence!"

(* Remove redundant (difference between the two doesn't change) positions.
   This step is not strictly necessary. *)
let reduce_position_map : position_map -> position_map = function
  | []          -> invalid_arg "reduce_position_map: empty"
  | (p, a) :: t ->
      List.fold_left t ~init:(a - p, [p,a])
        ~f:(fun (d, acc) (p,a) ->
            let d1 = a - p in
            if d <> d1 then (d1, (p,a) :: acc) else (d, acc))
      |> snd
      |> List.rev

(* Helper method to create an actual function for computing the index into
   the base state array. This is useful for debugging between the MSA.Parser
   positions and the index into Base array. Assumes a 'reduced' (via
   [reduce_position_map]) position map. *)
let to_position_map : position_map -> (int -> int option) = function
  | [] -> invalid_arg "to_position_map: empty"
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

let init_state =
  Array.map ~f:Pm.init_first_d

(* Add an allele's MSA.Parser.sequence to the state. *)
let add_alternate_allele ~position_map allele allele_instr reference_arr arr =
  let add_reference_value start end_ =
    for i = start to end_ - 1 do
      arr.(i) <- Pm.add reference_arr.(i) arr.(i)
    done
  in
  let add_gap start end_ =
    for i = start to end_ - 1 do
      arr.(i) <- Pm.add None arr.(i)
    done
  in
  let add_sequence start s =
    String.iteri s ~f:(fun i c ->
      let j = i + start in
      arr.(j) <- Pm.add (some (Base.of_char c)) arr.(j))
  in
  let prefix = sprintf "add_alternate_allele (%s): " allele in
  let ia fmt = invalid_argf ~prefix fmt in
  let open MSA in
  let rec loop pmap lp = function
    | End p :: []                 ->
        let ap, _position_map = position_and_advance p pmap in
        add_reference_value lp ap
    | Start p :: _                -> ia "second start: %d" p
    | []                          -> ia "didn't End"
    | End p :: t                  -> ia "end before end: %d." p
    | Boundary { pos; _ } :: t    ->
        let ap, pmap = position_and_advance pos pmap in
        add_reference_value lp ap;
        loop pmap ap t
    | Sequence { start; s } :: t  ->
        let ap, pmap = position_and_advance start pmap in
        add_reference_value lp ap;
        add_sequence ap s;
        let fap = ap + String.length s in
        loop pmap fap t
    | Gap { gstart; length } :: t ->
        let ap, pmap = position_and_advance gstart pmap in
        add_reference_value lp ap;
        let fap = ap + length in
        add_gap ap fap;
        loop pmap fap t
  in
  match allele_instr with
  | Start s :: t -> loop position_map 0 t
  | e :: _       -> ia "not at Start : %s" (al_el_to_string e)
  | []           -> ia "Empty sequence!"

let pos_to_string pm =
  Pm.to_string pm
    (function | None   -> "none"
              | Some c -> sprintf "some %c" (Base.to_char c))

(***** Forward Pass ***)

(* For every k there are 3 possible states. *)
type 'a cell =
  { match_  : 'a
  ; insert  : 'a
  ; delete  : 'a
  }

let cell_to_string f c =
  sprintf "{m: %s; i: %s; d: %s}"
    (f c.match_) (f c.insert) (f c.delete)

let float_cell_to_string = cell_to_string (sprintf "%0.3f")

(* Note that we'll embed the missing/None/gap logic inside of 'a.
   They will all come from some Ring and we'll check via {is_gap}. *)
type ('a, 'e) cell_recurrences =
  { start     : ?base_p:'a -> 'e -> 'a cell -> 'a cell
  ; fst_col   : 'e -> 'a cell -> 'a cell
  ; middle    : 'e -> match_c:('a cell)
                   -> insert_c:('a cell)
                   -> delete_c:('a cell)
                   -> 'a cell
  ; end_      : 'a cell -> 'a
  }

exception FilterIsPastThreshold of string

let filter_is_past_thresholdf fmt =
  ksprintf (fun s -> raise (FilterIsPastThreshold s)) fmt

module Filter = struct

  (* Filter's operate in two stages, we first observe a {'base} (an actual base
      and quality pair) and then as we scan the reference we see all of the
      resulting {entry}'s that we're going to insert into the matrix storing the
      state of the pass. Either method can throw the {FilterIsPastThreshold}
      exception, because we may only want to act after we've observed the best
      entry for a given base.

      [init] is used as the first observation in the pass and may be set
      approriately before we start a pass.
  *)
  type ('a, 'base, 'entry) t =
    { init  : 'a
    ; base  : 'a -> 'base -> 'a
    ; entry : 'a -> 'entry -> 'a
    }

  let join f1 f2 =
    let init  = f1.init, f2.init in
    let base (s1, s2) b = f1.base s1 b, f2.base s2 b in
    let entry (s1, s2) b = f1.entry s1 b, f2.entry s2 b in
    { init ; base ; entry }

  let empty =
    { init  = ()
    ; base  = begin fun () _ -> () end
    ; entry = begin fun () _ -> () end
    }

  let update t n_init =
    { t with init = n_init }

end (* Filter *)

(** What are the values and equations that determine how probabilities are
    calculated in the forward pass. *)
module Forward_calcations_over_cells (R : Probability.Ring) = struct

  let mismatch_p =
    R.times_one_third

  let match_p =
    R.complement_probability

  (* TODO. Avoid the `float_of_int (Phred_score.to_int c) /. -10.` round trip
      between converting to log10p and then back to log10, and just use char's
      instead for the quality calc. *)
  let to_match_prob (base, base_error) =
    let compare_against c =
      if base = 'N' then begin
        (* This is commented out because this function needs to be fast.
         * eprintf "Spotted an N with %f\n%!" base_error; *)
        R.one
      end else if base = c then
        match_p base_error
      else
        mismatch_p base_error
    in
    let open Base in
    function
    | A -> compare_against 'A'
    | C -> compare_against 'C'
    | G -> compare_against 'G'
    | T -> compare_against 'T'

  let debug_ref = ref false

  let max3 m i d =
    R.max (R.max m i) d

  let max2_by_fst a b =
    let av, al = a in
    let bv, bl = b in
    if R.(bv <= av) then
      a
    else
      b

  let max3_by_fst m i d =
    max2_by_fst (max2_by_fst m i) d

  let g ?(insert_p=Phmm.default_insert_probability) tm read_length =

    let open R in                       (* Opening R shadows '+' and '*' below*)
    let open Phmm.TransitionMatrix in
    let t_s_m = constant (tm StartOrEnd Match) in
    let t_s_i = constant (tm StartOrEnd Insert) in
    let t_m_m = constant (tm Match Match) in
    let t_i_m = constant (tm Insert Match) in
    let t_d_m = constant (tm Delete Match) in

    let t_m_i = constant (tm Match Insert) in
    let t_i_i = constant (tm Insert Insert) in

    let t_m_d = constant (tm Match Delete) in
    let t_d_d = constant (tm Delete Delete) in

    let t_m_s = constant (tm Match StartOrEnd) in
    let t_i_s = constant (tm Insert StartOrEnd) in

    (*printf "t_s_m : %s before %f\n" (R.to_string t_s_m) (tm StartOrEnd Match);
    printf "t_s_i : %s before %f\n" (R.to_string t_s_i) (tm StartOrEnd Insert);
    printf "t_m_m : %s before %f\n" (R.to_string t_m_m) (tm Match Match);
    printf "t_i_m : %s before %f\n" (R.to_string t_i_m) (tm Insert Match);
    printf "t_d_m : %s before %f\n" (R.to_string t_d_m) (tm Delete Match);

    printf "t_m_i : %s before %f\n" (R.to_string t_m_i) (tm Match Insert);
    printf "t_i_i : %s before %f\n" (R.to_string t_i_i) (tm Insert Insert);

    printf "t_m_d : %s before %f\n" (R.to_string t_m_d) (tm Match Delete);
    printf "t_d_d : %s before %f\n" (R.to_string t_d_d) (tm Delete Delete);

    printf "t_m_s : %s before %f\n" (R.to_string t_m_s) (tm Match StartOrEnd);
    printf "t_i_s : %s before %f\n" (R.to_string t_i_s) (tm Insert StartOrEnd); *)

    let insert_p = constant insert_p in
    let start_i = insert_p * t_s_i in
    let f = (* forward *)
      let start ?(base_p=one) emission_p prev_c =
        if is_gap emission_p then
          prev_c
        else
          { match_ = emission_p * t_s_m * base_p
          ; insert = start_i * base_p
          ; delete = zero (* * base_p *)
          }
      in
      let fst_col emission_p insert_c =
        if is_gap emission_p then
          (* This is a weird case, for an allele there is a gap, in the first
             position. We'll 'fail' by propogating nan. *)
          { match_ = gap
          ; insert = start_i
          ; delete = gap
          }
        else
          { match_ = emission_p * ( t_m_m * zero
                                  + t_i_m * zero
                                  + t_d_m * zero)
          ; insert = insert_p   * ( t_m_i * insert_c.match_
                                  + t_i_i * insert_c.insert)
          ; delete = (* one *  *) ( t_m_d * zero
                                  + t_d_d * zero)
          }
      in
      let middle emission_p ~match_c ~insert_c ~delete_c =
        let r =
          if is_gap emission_p then
            delete_c
          else
            { match_ = emission_p * ( t_m_m * match_c.match_
                                    + t_i_m * match_c.insert
                                    + t_d_m * match_c.delete)
            ; insert = insert_p   * ( t_m_i * insert_c.match_
                                    + t_i_i * insert_c.insert)
            ; delete = (* one *)    ( t_m_d * delete_c.match_
                                    + t_d_d * delete_c.delete)
            }
        in
        r
      in
      let end_ cell = cell.match_ * t_m_s + cell.insert * t_i_s in
      { start
      ; fst_col
      ; middle
      ; end_
      }
    in
    let pfst v (t, l) = t * v, l in               (* Multiply first *)
    let v = (* viterbi *)
      let start ?(base_p=(one,[])) emission_p prev_c =
        if is_gap emission_p then
          prev_c
        else
          { match_ = pfst emission_p (pfst t_s_m base_p)
          ; insert = pfst start_i base_p
          ; delete = zero (* * base_p *),         []
          }
      in
      let fst_col emission_p insert_c =
        if is_gap emission_p then
          (* This is a weird case, for an allele there is a gap, in the first
            position. We'll 'fail' by propogating nan. *)
          { match_ = gap,     []
          ; insert = start_i, []
          ; delete = gap,     []
          }
        else
          { match_ = emission_p * (max3 (t_m_m * zero)
                                        (t_i_m * zero)
                                        (t_d_m * zero))
                   , []
          ; insert = pfst insert_p (max2_by_fst (pfst t_m_i insert_c.match_)
                                                (pfst t_i_i insert_c.insert))
          ; delete = (* one *  *) (max (t_m_d * zero)
                                       (t_d_d * zero))
                   , []
          }
      in
      let middle emission_p ~match_c ~insert_c ~delete_c =
        if is_gap emission_p then
          delete_c
        else
          { match_ = pfst emission_p (max3_by_fst (pfst t_m_m match_c.match_)
                                                  (pfst t_i_m match_c.insert)
                                                  (pfst t_d_m match_c.delete))
          ; insert = pfst insert_p   (max2_by_fst (pfst t_m_i insert_c.match_)
                                                  (pfst t_i_i insert_c.insert))
          ; delete = (* pfst one *)  (max2_by_fst (pfst t_m_d delete_c.match_)
                                                  (pfst t_d_d delete_c.delete))
          }
      in
      let end_ cell = max2_by_fst (pfst t_m_s cell.match_)
                                  (pfst t_i_s cell.insert)
      in
      { start
      ; fst_col
      ; middle
      ; end_
      }
    in
    f, v

  let zero_cell =
    { match_ = R.zero
    ; insert = R.zero
    ; delete = R.zero
    }

  let cells_close_enough c1 c2 =
    R.close_enough c1.match_ c2.match_
    && R.close_enough c1.insert c2.insert
    && R.close_enough c1.delete c2.delete

  let cells_close_enough_cmp c1 c2 =
    let mc = R.close_enough_cmp c1.match_ c2.match_ in
    if mc <> 0 then mc else
      let ic = R.close_enough_cmp c1.insert c2.insert in
      if ic <> 0 then ic else
        R.close_enough_cmp c1.delete c2.delete

  let max_cell_by_match c1 c2 =
    if R.(c1.match_ < c2.match_) then c2 else c1

end (* Forward_calcations_over_cells *)

module ForwardFilters (R : Probability.Ring) = struct

  open Filter
  module Fc = Forward_calcations_over_cells(R)

  let debug_ref = ref false

  module Past_threshold = struct

    type t =
      { threshold       : R.t
      ; best_seen_entry : R.t
      }

    let to_string { threshold; best_seen_entry } =
      sprintf "{threshold: %s; best_seen_entry:%s }"
        (R.to_string threshold)
        (R.to_string best_seen_entry)

    let create_t threshold =
      { threshold; best_seen_entry = R.zero}

    let create threshold of_entry =
      let base {threshold; best_seen_entry} _base_doesnt_matter =
        if best_seen_entry <> R.zero && best_seen_entry < threshold then
          let msg =
            sprintf "threshold %s breached: %s"
              (R.to_string threshold)
              (R.to_string best_seen_entry)
          in
          filter_is_past_thresholdf "%s" msg
        else
          create_t threshold
      in
      let entry t e =
        { t with best_seen_entry = R.max t.best_seen_entry (of_entry e)}
      in
      let init = create_t threshold in
      { init; base; entry }

    let update_init_value last =
      create_t R.(last.threshold / last.best_seen_entry)

  end (* Past_threshold *)

  module Max_number_of_mismatches = struct

    (* Keep track of the encountered read errors, assume that {mismatches}
      of the most likely (lowest base error) bases are wrong. While the the rest
      (total number of bases - mismatches) are correct. This way we have a well
      defined value for the threshold. If all of the path probablities are below
      the threshold, stop executing! *)

    (* Assume that {errors} are in ascending order; therefore the lowest base
      error are at the front of the list; this way when we say that there are
      'n' mismatches, those n have the lowest error (highest probability of
      being right). In this model/calculation there are only transitions between
      match states, we need to scale by the transition probabilities
      (t_s_m, t_m_m) in order to make it comparable with the match value in the
      PHMM cell. *)
    let ungapped_path_prob number_mismatches ~t_s_m ~t_m_m errors =
      let to_emission_p n e =
        if n < number_mismatches then Fc.mismatch_p e else Fc.match_p e
      in
      let rec loop n p l =
        match l with
        | []     -> p
        | e :: t -> let emission_p = to_emission_p n e in
                    loop (n + 1) R.(p * emission_p * t_m_m) t
      in
      match errors with
      | []     -> R.one
      | e :: t -> let emission_p = to_emission_p 0 e in
                  loop 1 R.(t_s_m * emission_p) t

    (* Insert by maintaining the ascending order. *)
    let rec insert error l = match l with
      | []     -> error :: []
      | e :: t -> if error <= e then
                    error :: l
                  else
                    e :: (insert error t)

    let create tm ~number_mismatches of_entry =
      let open Phmm.TransitionMatrix in
      let t_s_m = R.constant (tm StartOrEnd Match) in
      let t_m_m = R.constant (tm Match Match) in
      let path_prob = ungapped_path_prob number_mismatches ~t_s_m ~t_m_m in
      (* Called at the start of each base *)
      let base (errors, max_value) (_base, base_error) =
        let threshold = path_prob errors in
        let msg = lazy
            (sprintf "max_value: %s, threshold: %s, rows: %d"
              (R.to_string max_value)
              (R.to_string threshold)
              (List.length errors))
        in
        if !debug_ref then printf "filter: %s\n%!" (Lazy.force msg);
        if max_value < threshold then
          filter_is_past_thresholdf "%s" (Lazy.force msg)
        else
          (insert base_error errors, R.zero)
      in
      let entry (es, mv) c = (es, R.max mv (of_entry c)) in
      { init  = ([], R.one)
      ; base
      ; entry
      }

  end (* Max_number_of_mismatches *)

end (* ForwardFilters *)

(* Layout logic:
   The dynamic array consists of values for each base in the read in a row
   and for each value (base, or list of bases) of the reference in a column.

   'i' is used to index into the read, hence row.
   'k' is used to index into the reference, hence column

  The workspaces are reused for the 2nd read and the current logic is configured
  so that the 2 reads have the same length [rows].
*)

module type Workspace_intf = sig

  type t

  val generate : ref_length:int -> read_length:int -> t

  val columns : t -> int
  val rows : t -> int

  type entry
  type final_entry
  type final_emission

  val get : t -> i:int -> k:int -> entry
  val set : t -> i:int -> k:int -> entry -> unit

  val get_final : t -> int -> final_entry
  val set_final : t -> int -> final_entry -> unit

  val get_emission : t -> final_emission
  val set_emission : t -> final_emission -> unit

  (* This is exported to
     1. initialize the bands.
     2. grab the maximum match. *)
  val fold_over_row
      : ?range:(int * int)
      -> t
      -> int
      -> init:'a
      -> f:('a -> entry -> 'a)
      -> 'a

  (* Necessary for bands *)
  val foldi_over_row
      : ?range:(int * int)
      -> t
      -> int
      -> init:'a
      -> f:('a -> int -> entry -> 'a)
      -> 'a

  (* Fold over the final entries. These are the probabilities of finishing at
     a given reference position. These are then averaged to compute the final
     emission value. *)
  val fold_over_final
      : ?range:(int * int)
      -> t
      -> init:'a
      -> f:('a -> final_entry -> 'a)
      -> 'a

  val foldi_over_final
      : ?range:(int * int)
      -> t
      -> init:'a
      -> f:('a -> int -> final_entry -> 'a)
      -> 'a

  (* Reset for calculation. *)
  val clear : t -> unit

end (* Workspace_intf *)

module CommonWorkspace = struct

  (* The recurrence logic (see {middle} _below_) relies only upon the previous
     state of the read, or the row (the Markov in pphMm). Therefore, for the
     bulk of the forward pass (ignoring final emission) we can conserve space
     by storing only 2 rows and alternating. Externally, the workspace allows
     access to {read_length} amount of rows; although this isn't enforced and
     is only a product of Perform_forward_calculation calling {rows} correctly.

     Internally, the workspace maps (via mod) to the appriate offset into the 2
     rows. In practice this has a small impact on run time performance: ~1-3,
     as we essentially just shift when GC work gets done. But this has a huge
     impact on total memory usage and is probably an important change as it
     will allow more instances of prohlatype to be run in parallel. *)
  type ('entry, 'final_entry, 'final_emission) w =
    { mutable forward   : 'entry array array
    ; mutable final     : 'final_entry array
    ; mutable emission  : 'final_emission
    ; read_length       : int
    }

  let last_array_index arr =
    Array.length arr - 1

  let columns ws =
    last_array_index ws.forward.(0)

  let rows ws =
    ws.read_length - 1

  let get ws ~i ~k = ws.forward.(i mod 2).(k)

  let set ws ~i ~k e = ws.forward.(i mod 2).(k) <- e

  let get_final ws k = ws.final.(k)
  let set_final ws k e = ws.final.(k) <- e

  let get_emission ws = ws.emission
  let set_emission ws e = ws.emission <- e

  let to_se range ws =
    Option.value ~default:(0, columns ws) range

  let fold_over_row ?range ws i ~init ~f =
    let start, end_ = to_se range ws  in
    let rec loop k acc =
      if k >= end_ then acc else loop (k + 1) (f acc (get ws ~i ~k))
    in
    loop start init

  let foldi_over_row ?range ws i ~init ~f =
    let start, end_ = to_se range ws  in
    let rec loop k acc =
      if k >= end_ then acc else loop (k + 1) (f acc k (get ws ~i ~k))
    in
    loop start init

  let fold_over_final ?range ws ~init ~f =
    let start, end_ = to_se range ws  in
    let rec loop k acc =
      if k >= end_ then acc else loop (k + 1) (f acc (get_final ws k))
    in
    loop start init

  let foldi_over_final ?range ws ~init ~f =
    let start, end_ = to_se range ws  in
    let rec loop k acc =
      if k >= end_ then acc else loop (k + 1) (f acc k (get_final ws k))
    in
    loop start init

end (* CommonWorkspace *)

(* Create the workspace for the forward calculation for a single,
   aka reference, allele. *)
module SingleWorkspace (R : Probability.Ring) :
  (Workspace_intf with type entry = R.t cell
                   and type final_entry = R.t
                   and type final_emission = (R.t * int)) = struct

  module Fc = Forward_calcations_over_cells(R)

  type entry = R.t cell
  type final_entry = R.t
  type final_emission = R.t * int             (* Position of max final entry. *)

  include CommonWorkspace
  type t = (entry, final_entry, final_emission) w

  let init_emission = R.zero, -1
  let generate ~ref_length ~read_length =
    { forward  = Array.init 2 (*read_length*) ~f:(fun _ -> Array.make ref_length Fc.zero_cell)
    ; final    = Array.make ref_length R.zero
    ; emission = init_emission
    ; read_length
    }

  let clear ws =
    let ref_length  = Array.length ws.final in
    let number_rows = Array.length ws.forward in
    ws.forward  <- Array.init number_rows ~f:(fun _ -> Array.make ref_length Fc.zero_cell);
    ws.final    <- Array.make ref_length R.zero;
    ws.emission <- init_emission

  let save ws =
    let fname = Filename.temp_file ~temp_dir:"." "forward_workspace" "" in
    let oc = open_out fname in
    Marshal.to_channel oc ws [];
    close_out oc;
    printf "Saved workspace to %s\n" fname

  let load fname =
    let ic = open_in fname in
    let ws : t = Marshal.from_channel ic in
    close_in ic;
    ws

end (* SingleWorkspace *)

(* Describe a path, a decoding through the PHMM. *)
module Path = struct

  type t =
    (* Start and emission values *)
    | S of int                                   (* Index into reference index. *)
    | E of int
    | M of char * Base. t (* Store the read pair, that _may_ have 'N', and Base *)
    | I of int * char                     (* index into read and read base pair *)
    | D of int * Base.t

  let to_string = function
    | S i      -> sprintf "S%d" i
    | E i      -> sprintf "E%d" i
    | M (b, r) -> sprintf "M%c%c" b (Base.to_char r)
    | I (p, b) -> sprintf "I%d%c" p b
    | D (p, r) -> sprintf "D%d%c" p (Base.to_char r)

  type summary =
    { reference : string
    ; read      : string
    ; start     : int
    ; end_      : int
    }

  let to_strings =
    let g = '_' in
    let to_pair = function
      | M (b, r) -> (b, Base.to_char r)
      | I (_, b) -> (b, g)
      | D (_, r) -> (g, Base.to_char r)
      | _        -> assert false
    in
    let rec first start reada refa = function
      | E end_ :: [] ->
          Sp.Single { reference = List.rev refa |> String.of_character_list
                    ; read = List.rev reada |> String.of_character_list
                    ; start ; end_ }
      | E end_ :: S nstart :: tl ->
          let r1 =
            { reference = List.rev refa |> String.of_character_list
            ; read = List.rev reada |> String.of_character_list
            ; start ; end_ }
          in
          let r2 = second nstart [] [] tl in
          Sp.Paired (r1, r2)
      |  _ :: []
      | [] -> invalid_argf "Didn't end with an end"

      | p :: tl -> let a, b = to_pair p in
                   first start (a :: reada) (b :: refa) tl
    and second start reada refa = function
      | E end_ :: [] ->
          { reference = List.rev refa |> String.of_character_list
          ; read = List.rev reada |> String.of_character_list
          ; start ; end_ }
      | E _ :: _t    -> invalid_arg "Things after end!"
      | S _ :: _     -> invalid_arg "More than 2 starts in path!"
      |  _ :: []
      | []           -> invalid_argf "Didn't end with an end"
      | p :: tl      -> let a, b = to_pair p in
                        second start (a :: reada) (b :: refa) tl
    in
    function
      | S start :: tl -> first start [] [] tl
      | h  :: _       -> invalid_argf "Path.to_strings didn't start with S but: %s"
                           (to_string h)

      | []            -> invalid_argf "Empty path list"

end (* Path *)

module SingleViterbiWorkspace (R : Probability.Ring) :
  (Workspace_intf with type entry = (R.t * Path.t list) cell
                   and type final_entry = R.t * Path.t list
                   and type final_emission = R.t * Path.t list
  ) = struct

  module Fc = Forward_calcations_over_cells(R)

  type entry = (R.t * Path.t list) cell
  type final_entry = R.t * Path.t list
  type final_emission = R.t * Path.t list

  include CommonWorkspace
  type t = (entry, final_entry, final_emission) w

  let empty_entry =
    { match_ = R.zero, []
    ; insert = R.zero, []
    ; delete = R.zero, []
    }
  let empty_final_entry = R.zero, []
  let empty_final_emission = R.zero, []

  let generate ~ref_length ~read_length =
    { forward  = Array.init 2 (*read_length*) ~f:(fun _ -> Array.make ref_length empty_entry)
    ; final    = Array.make ref_length empty_final_entry
    ; emission = empty_final_emission
    ; read_length
    }

  let clear ws =
    let ref_length  = Array.length ws.final in
    let number_rows = Array.length ws.forward in
    ws.forward  <- Array.init number_rows ~f:(fun _ -> Array.make ref_length empty_entry);
    ws.final    <- Array.make ref_length empty_final_entry;
    ws.emission <- empty_final_emission

  let save ws =
    let fname = Filename.temp_file ~temp_dir:"." "forward_workspace" "" in
    let oc = open_out fname in
    Marshal.to_channel oc ws [];
    close_out oc;
    printf "Saved workspace to %s\n" fname

  let load fname =
    let ic = open_in fname in
    let ws : t = Marshal.from_channel ic in
    close_in ic;
    ws

end (* SingleViterbiWorkspace *)

type 'a mt = (Pm.ascending, 'a) Pm.t

(* Create and manage the workspace for the forward pass for multiple alleles.*)
module MakeMultipleWorkspace (R : Probability.Ring) :
  (Workspace_intf with type entry = R.t cell mt
                   and type final_entry = R.t mt
                   and type final_emission = (R.t * int) mt) = struct

  type entry = R.t cell mt
  type final_entry = R.t mt
  type final_emission = (R.t * int) mt

  include CommonWorkspace
  type t = (entry, final_entry, final_emission) w

  let pa = Pm.empty_a  (* This is just a place holder that is mutated. *)

  let generate ~ref_length ~read_length =
    { forward   = Array.init 2 (*read_length*) ~f:(fun _ -> Array.make ref_length pa)
    ; final     = Array.make ref_length pa
    ; emission  = pa
    ; read_length
    }

  let clear ws =
    let ref_length     = Array.length ws.final in
    let number_rows    = Array.length ws.forward in
    ws.forward  <- Array.init number_rows ~f:(fun _ -> Array.make ref_length pa);
    ws.final    <- Array.make ref_length pa;
    ws.emission <- pa

  let save ws =
    let fname = Filename.temp_file ~temp_dir:"." "forward_workspace" "" in
    let oc = open_out fname in
    Marshal.to_channel oc ws [];
    close_out oc;
    printf "Saved workspace to %s\n" fname

  let load fname =
    let ic = open_in fname in
    let ws : t = Marshal.from_channel ic in
    close_in ic;
    ws

end (* MakeMultipleWorkspace *)

(* Ring and position functions. *)
module R_p (R : Probability.Ring) = struct

  let equal (l1, p1) (l2, p2) =
    p1 = p2 && R.(close_enough l1 l2)

  (* R.t, int, R.t when reducing/merging a pm. *)
  let equal3 (l1, p1, bl1) (l2, p2, bl2) =
    p1 = p2 && R.(close_enough l1 l2 && close_enough bl1 bl2)

  (* We're tacitly taking the _first_ best final entry. Perhaps we should warn
   * if there is ever more than one? Ex. The read is small enough that we can't
   * tell? *)

  (* Add values and track best position . *)
  let avatbp acc pos new_value =
    let sum, best_pos, value_at_best_pos = acc in
    if R.is_gap new_value then
      acc
    else if R.(value_at_best_pos < new_value) then
      R.(sum + new_value), pos, new_value
    else
      R.(sum + new_value), best_pos, value_at_best_pos

end (* R_p *)

module Lp = Probability.Log10
module Lpr = R_p(Lp)

(* Pairing the observation makes it easier to abstract the regular vs
   reverse complement access pattern. Leaving this as a pair to avoid
   redundant pairing/unpairing.

   I'll use obsp (observation pair) as the variable name. *)
type obs = char * float

type read_accessor = int -> obs

(* I'm somewhat purposefully shadowing the cell_recurrences field names. *)
type ('workspace, 'entry, 'final_entry, 'final_emission, 'base) recurrences =
  (* final_emission should have the same type as likelihood *)
  { start   : ?base_p:'final_emission -> 'workspace -> obs -> 'base -> k:int -> 'entry
  ; fst_col : 'workspace -> obs -> 'base -> i:int -> k:int -> 'entry
  ; middle  : 'workspace -> obs -> 'base -> i:int -> k:int -> 'entry
  ; end_    : 'workspace -> int -> 'final_entry
  ; final_e : range:(int * int) -> 'workspace -> 'final_emission
  }

type column_spec =
  { cs_start    : int
  ; ref_length  : int
  }

let column_spec_to_string {cs_start; ref_length; _} =
  sprintf "{cs_start: %d; ref_length: %d}"
    cs_start ref_length

module type Forward_pass_intf = sig

  type workspace
  type entry
  type final_entry
  type final_emission

  val full
      : ?rows:int
      -> ?columns:column_spec
      -> ?base_p:final_emission
      -> filter:('a, obs, entry) Filter.t
      -> workspace
      -> (workspace, entry, final_entry, final_emission, 'base) recurrences
      -> read:read_accessor
      -> reference:(int -> 'base)
      -> 'a

  val paired
      : ?rows:int
      -> ?columns:column_spec
      -> filter:('a, obs, entry) Filter.t
      -> workspace
      -> (workspace, entry, final_entry, final_emission, 'base) recurrences
      -> read1:read_accessor
      -> read2:read_accessor
      -> reference:(int -> 'base)
      -> 'a

end (* Forward_pass_intf *)

(* Given:
    - a workspace
    - its dimensions
    - recurrence functions
    - read and references acccessors

   Produce functions that actually traverse and fill in the forward pass,
   and emission values. These functions also take an optional [filter] that
   may raise an exception ({FilterIsPastThreshold}) to signal a premature
   stop to the computation. *)
module Perform_forward_calculation (W : Workspace_intf) : sig
  include Forward_pass_intf
    with type workspace = W.t
     and type entry = W.entry
     and type final_entry = W.final_entry
     and type final_emission = W.final_emission
  end = struct

  type workspace = W.t
  type entry = W.entry
  type final_entry = W.final_entry
  type final_emission = W.final_emission

  let to_se columns ws =
    Option.value_map columns ~default:(0, W.columns ws)
    ~f:(fun {cs_start; ref_length } ->
          (* last indexable position *)
          let end_ = cs_start + ref_length - 1 in
          cs_start, end_)

  (* compute the full forward pass
     @param rows how many elements of the read to fill, defaults to the
        configured size of the workspace; which should be at most
        the length of the read.
     @param base_p base probability of a path; allows us to weight the 2nd read
        based upon the likelihood (final_emission) of the first. *)
  let pass ?rows ?columns ?base_p ~filter ws recurrences ~read ~reference =
    let open Filter in
    let start, end_ = to_se columns ws in
    let rows = Option.value rows ~default:(W.rows ws) in
    let a_0 = read 0 in
    let ft = ref (filter.base filter.init a_0) in
    for k = start to end_ do
      let ne = recurrences.start ?base_p ws a_0 (reference k) ~k in
      ft := filter.entry !ft ne;
      W.set ws ~i:0 ~k ne;
    done;
    for i = 1 to rows do
      let a_i = read i in
      ft := filter.base !ft a_i;
      let ne = recurrences.fst_col ws a_i ~i ~k:start (reference start) in
      ft := filter.entry !ft ne;
      W.set ws ~i ~k:start ne;
      for k = start + 1 to end_ do
        let ne = recurrences.middle ws a_i ~i ~k (reference k) in
        ft := filter.entry !ft ne;
        W.set ws ~i ~k ne
      done
    done;
    !ft

  let final ?columns ws recurrences =
    let start, end_ = to_se columns ws in
    for k = start to end_ do
      W.set_final ws k (recurrences.end_ ws k)
    done;
    start, end_

  (* After filling in both parts of the workspace,
     compute the final emission value. *)
  let full ?rows ?columns ?base_p ~filter ws recurrences ~read ~reference =
    let ft = pass ?rows ?columns ?base_p ~filter ws recurrences ~read ~reference in
    let range = final ?columns ws recurrences in
    W.set_emission ws (recurrences.final_e ~range ws);
    ft

  (* paired pass. *)
  let paired ?rows ?columns ~filter ws recurrences ~read1 ~read2 ~reference =
    let nft = full ?rows ?columns ~filter ws recurrences ~read:read1 ~reference in
    let base_p = W.get_emission ws in
    let filter = Filter.update filter nft in
    full ?rows ?columns ~base_p ~filter ws recurrences ~read:read2 ~reference

end (* Perform_forward_calculation *)

(* We introduce a custom result type because it is more informative of what
 * happened: a pass could end because it was Filtered as opposed to
 * Error'ing. *)
module Pass_result = struct

  type 'a t =
    | Completed of 'a
    | Filtered of string
    [@@deriving show,yojson]

  let to_string c_to_s = function
    | Completed c -> sprintf "Completed: %s" (c_to_s c)
    | Filtered m  -> sprintf "Filtered: %s" m

  let map ~f = function
    | Completed c -> Completed (f c)
    | Filtered m  -> Filtered m

  let completed = function
    | Completed _ -> true
    | Filtered  _ -> false

  let wrap f =
    try
      let r = f () in
      Completed r
    with FilterIsPastThreshold msg ->
      Filtered msg

  let unwrap = function
    | Completed c -> c
    | Filtered f  -> invalid_argf "Asked to unwrap a filtered %s result." f

end (* Pass_result *)

(* A helper record to expose the functions generated by these modules. *)
type 'a passes =
  { full    : read:(int -> obs) -> 'a
  ; paired  : read1:(int -> obs) -> read2:(int -> obs) -> 'a
  }

module ForwardSingleGen (R: Probability.Ring) = struct

  module W = SingleWorkspace(R)
  module Fc = Forward_calcations_over_cells(R)
  module Ff = ForwardFilters(R)
  module Rp = R_p(R)

  let debug_ref = ref false

  module MakePasses (Ws : Workspace_intf) = struct

    module Pfc = Perform_forward_calculation(Ws)

    let p ?max_number_mismatches ?transition_ref_length
      ~prealigned_transition_model ~read_length ~ref_length ws allele_a
       to_recurrences of_entry =
      let tm, columns =
        let rltm = Option.value transition_ref_length ~default:ref_length in
        if prealigned_transition_model then
          let column_spec = { cs_start = 0; ref_length; } in
          Phmm.TransitionMatrix.prealigned ~ref_length:rltm read_length
          , Some column_spec
        else
          Phmm.TransitionMatrix.init ~ref_length:rltm read_length
          , None
      in
      let reference k = allele_a.(k) in
      let recs = to_recurrences tm in
      let full filter ~read =
        Pass_result.wrap (fun () ->
          let _ = Pfc.full ?columns ~filter ws recs ~reference ~read in
          ())
      in
      let paired filter ~read1 ~read2 =
        Pass_result.wrap (fun () ->
          let _ = Pfc.paired ?columns ~filter ws recs ~reference ~read1 ~read2 in
          ())
      in
      match max_number_mismatches with
      | None ->
          { full   = full Filter.empty
          ; paired = paired Filter.empty
          }
      | Some number_mismatches ->
          let filter = Ff.Max_number_of_mismatches.create tm
                          ~number_mismatches of_entry in
          { full   = full filter
          ; paired = paired filter
          }

  end (* MakePasses *)

  let single_recurrences ?insert_p tm read_length =
    let r, _ = Fc.g ?insert_p tm read_length in
    let start ?base_p ws obsp base ~k =
      let prev_c = if k = 0 then Fc.zero_cell else (W.get ws ~i:0 ~k:(k-1)) in
      let base_p = Option.map base_p ~f:fst in      (* Ignore best position. *)
      r.start ?base_p (Fc.to_match_prob obsp base) prev_c
    in
    let fst_col ws obsp base ~i ~k =
      r.fst_col (Fc.to_match_prob obsp base) (W.get ws ~i:(i-1) ~k)
    in
    let middle ws obsp base ~i ~k =
      let emp = Fc.to_match_prob obsp base in
      let ks = k-1 in
      let insert_c = W.get ws ~i:(i-1) ~k in
      let match_c = W.get ws ~i:(i-1) ~k:ks in
      let delete_c = W.get ws ~i ~k:ks in
      r.middle emp ~insert_c ~delete_c ~match_c
    in
    let end_ ws k = r.end_ (W.get ws ~i:(read_length-1) ~k) in
    let final_e ~range ws =
      let sum_likelihoods, best_pos, _entry_at_best_pos =
        W.foldi_over_final ~range ws ~init:(R.zero, -1, R.zero)
          ~f:Rp.avatbp
      in
      sum_likelihoods, best_pos
    in
    { start ; fst_col ; middle ; end_ ; final_e }

  module SinglePasses = MakePasses(W)

  let passes_s ?max_number_mismatches ?insert_p ?transition_ref_length
    ~prealigned_transition_model
    ~read_length ~ref_length ws allele_a =
    SinglePasses.p ?max_number_mismatches ?transition_ref_length
      ~prealigned_transition_model ~read_length ~ref_length
      ws allele_a
      (fun tm -> single_recurrences ?insert_p tm read_length)
      (fun c -> c.match_)


  module V = SingleViterbiWorkspace(R)

  let v_zero_cell =
    { match_ = R.zero, []
    ; insert = R.zero, []
    ; delete = R.zero, []
    }

  let map_snd f (a, b) = (a, f b)

  (* Viterbi *)
  let viterbi_recurrences ?insert_p tm read_length =
    let open Path in
    let _, r = Fc.g ?insert_p tm read_length in
    let annotate_cell i b r c =
      { match_  = map_snd (fun l -> M (b, r) :: l) c.match_
      ; insert  = map_snd (fun l -> I (i, b) :: l) c.insert
      ; delete  = map_snd (fun l -> D (i, r) :: l) c.delete
      }
    in
    let start ?base_p ws obsp base ~k =
      let prev_c =
        if k = 0 then
          v_zero_cell
        else
          V.get ws ~i:0 ~k:(k-1)
      in
      let base_p =
        Option.value_map base_p
          ~default:(R.one, [ S k])
          ~f:(fun (bp, pp) -> (bp, S k :: pp))
      in
      let nc = r.start ~base_p (Fc.to_match_prob obsp base) prev_c in
      annotate_cell 0 (fst obsp) base nc
    in
    let fst_col ws obsp base ~i ~k =
      let pc = V.get ws ~i:(i-1) ~k in
      let nc = r.fst_col (Fc.to_match_prob obsp base) pc in
      annotate_cell i (fst obsp) base nc
    in
    let middle ws obsp base ~i ~k =
      let emp = Fc.to_match_prob obsp base in
      let ks = k-1 in
      let insert_c = V.get ws ~i:(i-1) ~k in
      let match_c  = V.get ws ~i:(i-1) ~k:ks in
      let delete_c = V.get ws ~i ~k:ks in
      let nc = r.middle emp ~insert_c ~delete_c ~match_c in
      let pt = annotate_cell i (fst obsp) base nc in
      if !debug_ref then begin
        let fst_two_pths = function
          | []          -> "         "
          | h :: []     -> sprintf "%-9s" (Path.to_string h)
          | h :: i :: _ -> sprintf "%s - %s" (Path.to_string h) (Path.to_string i)
        in
        printf "at i: %d k: %d obs: %c -%f, base: %c pt: %s\n"
          i k
          (fst obsp)
          (snd obsp)
          (Base.to_char base)
          (cell_to_string (fun (p, l) ->
              sprintf "%s - %s" (R.to_string p) (fst_two_pths l)) pt)
      end;
      pt
    in
    let end_ ws k =
      (*let en, lst = V.get ws ~i:(read_length-1) ~k in
      r.end_ en, (E k :: lst) (* The list contains best because we chose at middle. *) *)
      let en, lst = r.end_ (V.get ws ~i:(read_length-1) ~k) in
      en, (E k :: lst)
    in
    let final_e ~range ws =
      let open R in
      V.fold_over_final ~range ws ~init:(zero, [])
        ~f:(fun ((be, _) as b) ((fe, _) as e)->
              if R.(be < fe) then e else b)
    in
    { start ; fst_col ; middle ; end_ ; final_e }

  module ViterbiPasses = MakePasses(V)

  let passes_v ?insert_p ?transition_ref_length ~prealigned_transition_model
    ~read_length ~ref_length ws allele_a =
    ViterbiPasses.p ?transition_ref_length ~prealigned_transition_model
      ~read_length ~ref_length ws allele_a
      (fun tm -> viterbi_recurrences ?insert_p tm read_length)
      (fun _ -> failwith "No filters during viterbi passes")

end (* ForwardSingleGen *)

module PosMap = Map.Make(struct type t = int [@@deriving ord] end)

let largest k a i lst = topn ~greater:(>) k a i lst
let smallest k a i lst = topn ~greater:(<) k a i lst

(* Band config *)
type band_config =
  { warmup  : int     (* How many rows of the forward space to fill before
                         using a banded pass. *)
  ; number  : int     (* How many bands to calculate. This logic seems
                         inelegant in the sense that if we reduce the
                         calculation from the full forward space to just this
                         number of bands, why not reduce it further when some
                         bands, inevitably, will have far less probability
                         mass. *)
  ; radius   : int    (* How many columns of the band to calculate. *)
  }

let band_default =
  { warmup  = 10
  ; number  = 5
  ; radius  = 3
  }

module ForwardMultipleGen (R : Probability.Ring) = struct

  let debug_ref = ref false

  (* Eh... not the best place for it. *)
  let pm_max_cell_by_match =
    Pm.fold_values ~init:R.zero ~f:(fun m c -> R.max m c.match_)

  module W = MakeMultipleWorkspace(R)

  let max_of_ws_row ?range ws r =
    W.fold_over_row ?range ws r ~init:R.zero
      ~f:(fun v fpm -> R.max v (pm_max_cell_by_match fpm))

  let max_pm_of_ws_row ?range ws r =
    W.fold_over_row ?range ws r ~init:(R.zero, Pm.empty_a)
      ~f:(fun (bv,bpm) fpm ->
            let nv = pm_max_cell_by_match fpm in
            if R.(bv < nv) then
              (nv, fpm)
            else
              (bv, bpm))
    |> snd

  let median_of_pm apm =
    let values_with_size =
      Pm.fold_set_and_values apm ~init:[]
        ~f:(fun acc st v ->
              let nv = Partition_map.Set.size st in
              insert_sorted ~greater:R.( < ) v nv acc)
    in
    let n = Pm.size apm in
    let m = n / 2 in
    let even = n mod 2 = 0 in
    let rec loop number_seen = function
      | []           -> invalid_argf "Didn't find a median after %d values" n
      | (v, nv) :: t ->
          let new_number_seen = number_seen + nv in
          if new_number_seen > m then
            v
          else if new_number_seen = m then begin
            if even then
              match t with
              | []          -> invalid_argf "Didn't find a median after odd midpoint %d values" n
              | (v1,_) :: _ -> R.((v + v) / (constant 2.))
            else (* not even *)
              v
          end else loop new_number_seen t
    in
    loop 0 values_with_size

  let maximum_positions_median_match ?range ws r =
    max_pm_of_ws_row ?range ws r
    |> fun p -> Pm.map p R.close_enough ~f:(fun c -> c.match_)
    |> median_of_pm

  (* offset and emission probabilties
  type 'a oeps = (int * 'a) Cm.t
  type 'a emission_map = 'a oeps PosMap.t
  type 'a banded_recs =
    (* This isn't the greatest design, but we need to separate (to cache) the
       emission calculation that is then passed back into the banded middle
       call. *)
    { middle_emissions : obs
                          -> base_emissions
                          -> 'a oeps
    ; banded           : W.t
                          -> 'a oeps
                          -> ?prev_row:('a cell)              (* Needed for missing/edge cases.*)
                          -> ?cur_row:(W.entry)               (* Needed for missing/edge cases.*)
                          -> i:int
                          -> k:int
                          -> W.entry
    ; spec_final_e     : int list list -> W.t -> R.t array
    }*)

  let per_allele_emission_arr len =
    Array.make len R.one

  module Fc = Forward_calcations_over_cells(R)
  module Ff = ForwardFilters(R)
  module Rp = R_p(R)

  let pm_max_cell =
    Pm.fold_values ~init:Fc.zero_cell ~f:Fc.max_cell_by_match

  let max_cell_of_ws_row ?range ws r =
    W.fold_over_row ?range ws r ~init:Fc.zero_cell
    ~f:(fun c fpm -> Fc.max_cell_by_match c (pm_max_cell fpm))

  let recurrences ?insert_p tm read_length number_alleles =
    let r, _ = Fc.g ?insert_p tm read_length in

   (* TODO: I could imagine some scenario's where it makes sense to cache,
       precompute or memoize this calculation. The # of base errors isn't
       that large (<100) and there are only 4 bases. So we could be performing
       the same lookup. *)
    let to_em_set obsp emissions =
      Pm.map emissions R.close_enough
        ~f:(Option.value_map ~default:R.gap ~f:(Fc.to_match_prob obsp))
    in
    let zero_cell_pm = pm_init_all ~number_alleles Fc.zero_cell in
    let eq = Fc.cells_close_enough in
    let start ?base_p ws obsp base ~k =
      let ems = to_em_set obsp base in
      let prev_pm = if k = 0 then zero_cell_pm else W.get ws ~i:0 ~k:(k-1) in
      let m2 =
        match base_p with
        | None    -> Pm.merge ~eq  ems prev_pm r.start         (* Not weighing alleles. *)
        | Some l  -> let nl = Pm.map l R.close_enough ~f:fst in       (* Drop best pos. *)
                     Pm.merge3 ~eq nl ems prev_pm (fun base_p emission prev_c ->
                        r.start ~base_p emission prev_c)
      in
      m2
    in
    let fst_col ws obsp emissions ~i ~k =
      Pm.merge ~eq  (to_em_set obsp emissions) (W.get ws ~i:(i-1) ~k) r.fst_col
    in
    let middle ws obsp emissions ~i ~k =
      let matches = W.get ws ~i:(i-1) ~k:(k-1) in
      let inserts = W.get ws ~i:(i-1) ~k       in
      let deletes = W.get ws ~i       ~k:(k-1) in
      let ems = to_em_set obsp emissions in
      let r   = Pm.merge4 ~eq ems matches inserts deletes
        (fun emission_p match_c insert_c delete_c ->
          r.middle emission_p ~insert_c ~delete_c ~match_c)
      in
      if !debug_ref then begin
        printf "at i: %d k: %d: o:%c,%f e: %s, ems: %s, m: %d, i: %d, d: %d, r: %d \n%s\n%!"
          i k (fst obsp) (snd obsp)
              (Pm.to_string emissions
                  (function
                    | None -> "gap"
                    | Some b -> sprintf "%c" (Base.to_char b)))
              (Pm.to_string ems (R.to_string ~precision:10))
              (Pm.length matches)
              (Pm.length inserts)
              (Pm.length deletes)
              (Pm.length r)
              (Pm.to_string r (cell_to_string R.to_string))
      end;
      r
    in
    let end_ ws k = Pm.map (W.get ws ~i:(read_length-1) ~k) R.close_enough ~f:r.end_ in
    let final_e ~range ws =
      (* CAN'T use empty_a since we're merging! *)
      let init = pm_init_all ~number_alleles (R.zero, -1, R.zero) in
      let triple_pm =
        W.foldi_over_final ~range ws ~init
          ~f:(fun e1 p e2 -> Pm.merge ~eq:Rp.equal3 e1 e2
                (fun acc l -> Rp.avatbp acc p l))
      in
      Pm.map triple_pm Rp.equal ~f:(fun (ls, bp, _bl) -> ls, bp)
    in
    (*
    let banded ws allele_ems ?prev_row ?cur_row ~i ~k =
      let with_insert inters (offset, emission_p) insert_c =
        let ks = Pervasives.(+) k offset in
        let calc insert_c delete_c match_c =
          r.middle emission_p ~insert_c ~delete_c ~match_c
        in
        let matches = W.get ws ~i:(i-1) ~k:ks in
        let deletes = W.get ws ~i ~k:ks in
        let insertsi = Cm.singleton inters insert_c in
        Cm.map3_partial ~eq insertsi
          ~by1:deletes
          ~missing1:(fun missing_deletes _insert_c ->
              let default = Cm.singleton missing_deletes Fc.zero_cell in
              Option.value_map ~default cur_row ~f:(fun as_ ->
                Option.value (Cm.get missing_deletes as_) ~default))
          ~by2:matches
          ~missing2:(fun missing_matches _insert_c _delete_c ->
            Cm.singleton missing_matches
              (Option.value prev_row ~default:Fc.zero_cell))
          ~f:calc
      in
      let inserts = W.get ws ~i:(i-1) ~k in
        Cm.concat_map2_partial ~eq allele_ems ~by:inserts
          ~missing:(fun missing_inserts ep_pair ->
              match prev_row with
              | None -> invalid_argf "At %d %d looking for inserts still missing %s"
                          k i (Cm.allele_set_to_string missing_inserts)
              | Some v -> with_insert missing_inserts ep_pair v)
          ~f:with_insert
    in
    let spec_final_e spec_cols ws =
      let em = Array.make number_alleles R.zero in
      List.iter spec_cols ~f:(fun cols ->
        List.iter cols ~f:(fun k ->
          update_emission_from_cam em (W.get_final ws k)));
      em
    in
    *)
    { start; fst_col; middle; end_; final_e}
    (*, { middle_emissions = to_em_set ; banded ; spec_final_e} *)

  module Pfc = Perform_forward_calculation(W)

  (* Banded Pass logic *
  module Bands = struct

    let debug_ref = ref false

    (* 1. Identify bands. *)
    let select_specific_band_indices ws c =
      let lg = largest c.number in        (* compares by match_, insert, delete *)
      W.fold_over_row ws c.warmup ~init:(0, Cm.init_everything [])
        ~f:(fun (k, acc) by ->
              k + 1,
              Cm.map2_partial acc ~by ~f:(fun lst c -> lg c k lst)
                ~missing:(fun s v -> Cm.singleton s v))
      |> snd

    let expand_allele_set_map l =
      Cm.to_list l
      |> List.map ~f:(fun (alleles, l) -> List.map l ~f:(fun c -> alleles, c))
      |> List.concat

    let group_by_allele_value lst =
      let rec loop as1 v1 acc = function
        | []              ->
            List.rev ((as1, v1) :: acc)
        | (as2, v2) :: tl ->
            if v1 = v2 then
              loop (Aset.union as1 as2) v1 acc tl
            else
              loop as2 v2 ((as1, v1) :: acc) tl
      in
      match List.sort ~cmp:(fun (_, v1) (_, v2) -> compare v1 v2) lst with
      | []              -> []
      | (as1, v1) :: tl -> loop as1 v1 [] tl

    (* TODO: keeping the bv = best_value through this transform is awkward, but
      seems like the most straightforward. *)
    let find_indices_above emissions inds =
      Cm.concat_map inds ~f:(fun s (bv, ilst) ->
        let i = List.hd_exn ilst in
        Cm.get_exn s emissions.(i)
        |> Cm.map ~bijective:true ~f:(fun (_bs, o) ->
            if o = min_int then
              (bv, ilst)
            else
              (bv, (i + o) :: ilst)))

    let find_indices_below increments inds =
      let bigK = Array.length increments in
      Cm.concat_map inds ~f:(fun s (bv, ilst) ->
        let i = List.hd_exn ilst in
        if i = bigK then
          Cm.singleton s (bv, ilst)
        else
          Cm.get_exn s increments.(i)
          |> Cm.map ~bijective:true ~f:(fun o ->
            bv, o :: ilst))

    let n_times n f s =
      let rec loop i a =
        if i = n then a else
          loop (i + 1) (f a)
      in
      loop 0 s

    let find_indices_above_n n emissions inds =
      n_times n (find_indices_above emissions) inds

    let find_indices_below_n n increments inds =
      n_times n (find_indices_below increments) inds

    let to_bands emissions_a increment_a c ~to_index inds =
      let lnds = Cm.map inds ~bijective:true ~f:(fun st -> st, [to_index st]) in
      let ai = find_indices_above_n c.radius emissions_a lnds in
      let bi = find_indices_below_n c.radius increment_a lnds in
      (* tl_exn -> drop the center band, so we don't duplicate it. *)
      Cm.map2 ai bi ~f:(fun (st, a) (_st2, b) -> st, a @ (List.tl_exn (List.rev b)))

    (* This step (see get_exn) partitions the current bands based on the
      last cell's. Since these value are already away from the optimal point
      in the band (that would be radius above), I'm not certain that this
      is necessary. We're using this value as just an approximation in lieu
      of filling in the forward matrix. Specifically, we can have a less strict
      cell equality test, so that the allele sets are joined together.

      This is a general notion to consider; when are cells close enough
      (difference just due to numerical rounding) that it isn't worth the
      split. *)
    let lookup_previous_values ws row bands =
      Cm.concat_map bands ~f:(fun s (bv, cols) ->
        let end_col = List.last cols (* doesn't compile! ~msg:"empty cols!"*) in
        Cm.get_exn s (W.get ws ~i:row ~k:end_col)
        |> Cm.map ~bijective:true ~f:(fun lv -> (cols, bv, lv)))

    type t =
      { cols        : int list
      ; best_value  : R.t cell
      ; last_value  : R.t cell         (* Doesn't have to occur at end_col *)
      ; alleles     : set
      } (*The order determines comparison. *)

    let to_string t =
      sprintf "cols: %s\tbv: %s \tlv: %s\n\t\ts: %s"
        (string_of_list ~sep:";" ~f:(sprintf "%d") t.cols)
        (cell_to_string R.to_string t.best_value)
        (cell_to_string R.to_string t.last_value)
        (Aset.to_human_readable t.alleles)

    (* TODO: Should we shift these down 1 ala next_band ? *)
    let setup emissions_a increment_a c ws =
      select_specific_band_indices ws c
      |> expand_allele_set_map
      |> group_by_allele_value
      (* We have to keep the Allele.Set bands separate, not in an
        Cm.t to avoid overlaps. *)
      |> List.map ~f:(fun p ->
          Cm.of_list [p]
          |> to_bands emissions_a increment_a c ~to_index:(fun (bv, i) -> i)
          |> lookup_previous_values ws c.warmup
          |> Cm.to_list
          |> List.map ~f:(fun (alleles, (cols, (best_value, _), last_value)) ->
              { cols; best_value; last_value; alleles}))
      |>  List.flatten

    (* As we fill a band we keep track of a little bit of state to determine how
      we orient the next band parameters. In particular we need to
      1. Find the highest likelihood value in the pass: best_c.match_. This helps to
          orient the next two functions.
      2. We need to know when to stop filling the band. We could use a fixed
          radius and do something like just move down 1 position per column but:
            - This doesn't account for gaps in the alleles.
            - This doesn't account for inserts/deletes that will shift the center
              of the band. In area's of ambiguity we could have 2 (or more)
              likelihood values that are close in value so we may err in
              determining the center.
          Therefore we need an adaptive strategy. We count the number of match
          values that are worse; where the likelihood is less than the best.
          Once this counter reaches the band config's radius we stop considering
          those alleles.
        3. Keep track of the calculated cols for an allele. This allows us to
          adaptively figure out the cols of the next pass by moving radius away
          from the best_col. See [find_next_row_from_fill_state].
      *)
    type fill_state =
      { best_col : int          (* where is the best, lowest match likelihood row *)
      ; best_c   : R.t cell
      ; worse    : int          (* Number of likelihoods < than best_c.match_ *)
      ; last_c   : R.t cell
      ; ncols    : int list     (* Where we're calculating. Since there might be
                                  gaps, the radius needs to look inside this list
                                  for the next start/end_col *)
      }

    let init_fill_state col cell =
      { best_col  = col
      ; best_c    = cell
      ; worse     = 0
      ; last_c    = cell
      ; ncols     = [col]
      }

    let update_fill_state col fs cell =
      if cell.match_ > fs.best_c.match_ then
        { best_col = col
        ; best_c   = cell
        ; worse    = 0
        ; last_c   = cell
        ; ncols    = col :: fs.ncols
        }
      else
        { fs with worse  = fs.worse + 1
                ; last_c = cell
                ; ncols  = col :: fs.ncols
        }

    (* cols are in reverse, descending, order! *)
    let to_next_cols radius best_col cols =
      let rec find_best acc = function
        | []     -> invalid_argf "Didn't find best row."
        | h :: t ->
            if h = best_col then
              (* TODO: These can silently take less than radius. *)
              let s = List.take t radius in
              let e = List.take acc radius in
              (List.rev s) @ (h :: e)
            else
              find_best (h :: acc) t
      in
      find_best [] cols

    let find_next_row_from_fill_state c fs =
      to_next_cols c.radius fs.best_col fs.ncols

    let next_band emissions_a increment_a c fs_map =
      Cm.concat_map fs_map ~f:(fun alleles fs ->
        (* Shift the band, by adjusting around best_col, for next column *)
        Cm.get_exn alleles increment_a.(fs.best_col)
        (* Now fill in the radius. *)
        |> to_bands emissions_a increment_a c ~to_index:(fun nr -> nr)
        |> Cm.map ~bijective:true ~f:(fun (_br,cols) -> (cols, fs.best_c, fs.last_c)))
      |> Cm.to_list
      |> List.map ~f:(fun (alleles, (cols, best_value, last_value)) ->
          { cols ; alleles ; best_value ; last_value })

    let fill_next emissions_a increment_a c (rrecs, brecs) em_map ws obsp i b col_values =
      if !debug_ref then begin
        let base, base_prob = obsp in
        printf "current bands for %c %f cols:%s at %d \n\t%s\n"
          base base_prob
          (string_of_list ~sep:";" ~f:(sprintf "%d") b.cols)
            i (to_string b)
      end;
      let cur_row = Cm.get b.alleles col_values in
      let update ?cur_row emp k alleles =
        let em_values, nem_map =
          try
            let emv = PosMap.find k emp in
            emv, emp
          with Not_found ->
            let es = brecs.middle_emissions obsp emissions_a.(k) in
            let nemp = PosMap.add ~key:k ~data:es emp in
            es, nemp
        in
        let entry =
          if k = 0 then
            (* Even though f_start returns the transitions for _all_  alleles
               the rest of the banded logic assumes that this entry is for
               only _alleles. *)
            Cm.get_exn alleles (rrecs.start obsp emissions_a.(0))
          else
            let allele_emissions = Cm.get_exn alleles em_values in
            brecs.banded ws allele_emissions ~i ~k
              (* Poor design: No harm in adding prev_row and cur_row as banded
                will only use this value in the missing case. So we're not going
                to track that we're at the right row. *)
              ~prev_row:b.last_value ?cur_row
        in
        W.set ws ~i ~k (Cm.join entry (W.get ws ~i ~k));
        nem_map, entry
      in
      match b.cols with
      | []                -> invalid_argf "empty cols"
      | start_row :: trows -> begin
          let nem_map, first_entry = update ?cur_row em_map start_row b.alleles in
          let state =
            Cm.map first_entry ~bijective:true
              ~f:(init_fill_state start_row)
          in
          let update_fill_state prev nk cur =
            Cm.map2_partial prev ~by:cur
              ~f:(update_fill_state nk)
              ~missing:(fun s p -> Cm.singleton s p)
          in
          let rec loop em_map acc fill_state not_full_alleles cur_row = function
            | []        -> invalid_argf "empty row, was there only one row?"
            | k :: cols ->
                let nem_map, entry = update ~cur_row em_map k not_full_alleles in
                let new_fill_state = update_fill_state fill_state k entry in
                if cols <> [] then                  (* Still have remaining cols to fill. *)
                  loop nem_map acc new_fill_state not_full_alleles entry cols
                else begin                          (* Filled all that we're supposed to. *)
                  let full, not_full_state =
                    Cm.partition_map new_fill_state ~f:(fun _s fs ->
                      if fs.worse >= c.radius then `Fst fs else `Snd fs)
                  in
                  let full_bands = next_band emissions_a increment_a c full in
                  if !debug_ref then begin
                    printf "new bands for k:%d at %d: %d\n" k i (List.length full_bands);
                    List.iter full_bands ~f:(fun b ->
                      printf "\t%s\n" (to_string b))
                  end;
                  let nacc = full_bands @ acc in
                  if Cm.length not_full_state = 0 ||     (* Nothing left to fill -> Done *)
                    k = Array.length increment_a then                     (* Reached end! *)
                    nacc, nem_map
                  else begin
                    Cm.fold not_full_state ~init:(nacc, nem_map)
                      ~f:(fun init alleles state ->
                            (*printf "not_full_recurse %d %s in %s\n%!"
                              k (Alleles.Set.to_human_readable alleles)
                                (Alleles.Set.to_human_readable (Cm.domain increment_a.(k))); *)
                            Cm.get_exn alleles increment_a.(k)
                            |> Cm.fold ~init ~f:(fun (acc, em_map) alleles2 next_row ->
                                loop em_map acc (Cm.singleton alleles2 state) alleles2
                                entry [next_row]))
                  end
                end
          in
          loop nem_map [] state b.alleles first_entry trows
      end

    let fill_end rrecs ws b =
      List.iter b.cols ~f:(fun k ->
        W.set_final ws k (Cm.join (rrecs.end_ ws k) (W.get_final ws k)))

    let pass c emissions_a increment_a ws recurrences last_read_index a =
      (* order matters for passing along last_col *)
      let first_bands = setup emissions_a increment_a c ws |> List.sort ~cmp:compare in
      if !debug_ref then begin
        printf "first_bands %d \n" (List.length first_bands);
        List.iter first_bands ~f:(fun t ->
          printf "\t%s\n" (to_string t))
      end;
      let banded_middle first_banded_column =
        let rec loop bands i =
          let new_bands_to_flatten, _last_em_map, _last_col_values =
            List.fold_left bands ~init:([], PosMap.empty, Cm.empty_d)
              ~f:(fun (acc, em_map, col_values) b ->
                    let nb, nem_map =
                      fill_next emissions_a increment_a c recurrences em_map ws (a i) i b col_values
                    in
                    let ncol_values =
                      List.map nb ~f:(fun t -> t.alleles, t.last_value)
                      |> Cm.of_list
                    in
                    nb :: acc, nem_map, ncol_values)
          in
          if i = last_read_index then
            bands (* We need these bands for end_ *)
          else begin
            let new_bands =
              List.flatten new_bands_to_flatten
              (* The default comparator will sort first by cols (first field),
                and within the int lists, the comparison is by the values,
                with smaller length lists taking precedence. *)
              |> List.sort ~cmp:compare
            in
            if !debug_ref then begin
              printf "bands at %d %d \n" i (List.length new_bands);
              List.iter new_bands ~f:(fun t ->
                printf "\t%d%s\n" (Hashtbl.hash t) (to_string t))
            end;
            loop new_bands (i + 1)
          end
        in
        loop first_bands first_banded_column
      in
      let rrecs, brecs = recurrences in
      let banded_end bands =
        List.iter bands ~f:(fill_end rrecs ws);
        let spec_cols = List.map bands ~f:(fun b -> b.cols) in
        W.set_emission ws (brecs.spec_final_e spec_cols ws)
      in
      banded_end (banded_middle (c.warmup + 1))

  end (* Bands *)
*)

end (* ForwardMultipleGen *)

(* Instantiate the actual passes over the 2 implemented probability rings.
   But we'll only use the Log10 in production as one can see from the
   drivers. *)
module ForwardS = ForwardSingleGen(Probability.RegularFloats)
module ForwardSLogSpace = ForwardSingleGen(Lp)

module ForwardM = ForwardMultipleGen(Probability.RegularFloats)
module ForwardMLogSpace = ForwardMultipleGen (Lp)

type t =
  { locus           : Nomenclature.locus
  ; release         : string
  ; align_date      : string
  ; alleles         : (string * MSA.Alteration.t list) array        (* Canonical order. *)
  ; number_alleles  : int
  ; emissions_a     : base_emissions array
  }

let construct input =
  if not (Alleles.Input.is_imputed input) then
    invalid_argf "Allele input MUST be imputed!"
  else begin
    let open MSA.Parser in
    Alleles.Input.construct input >>= fun mp ->
      let base_arr, pmap = initialize_base_array_and_position_map mp.ref_elems in
      let state_a = init_state base_arr in
      List.iter mp.alt_elems ~f:(fun a ->
        add_alternate_allele pmap a.allele a.seq base_arr state_a);
      let eq x y =
        match x, y with
        | None,   None      -> true
        | None, _ | _, None -> false
        | Some x, Some y    -> Base.equal x y
      in
      let emissions_a = Array.map (Pm.ascending eq) state_a in
      let alleles =
        (mp.reference, []) :: List.map mp.alt_elems ~f:(fun a -> a.allele, a.alters)
        |> Array.of_list
      in
      let number_alleles = Array.length alleles in
      Ok { release    = mp.release
         ; align_date = mp.align_date
         ; locus      = mp.locus
         ; alleles
         ; number_alleles
         ; emissions_a
         }
  end

let save_pphmm t =
  let fname = Filename.temp_file ~temp_dir:"." "pphmm" "" in
  let oc = open_out fname in
  Marshal.to_channel oc t [];
  close_out oc;
  printf "Saved ParPHMM.t to %s\n" fname

let load_pphmm fname =
  let ic = open_in fname in
  let t : t = Marshal.from_channel ic in
  close_in ic;
  t

let float_arr_to_str a =
  Array.to_list a
  |> string_of_list ~sep:";" ~f:(sprintf "%f")
  |> sprintf "[|%s|]"

(* Abstracts, behind a function the regular or reverse complement access
   pattern of a read. This is to avoid manually reversing and converting
   the read. *)
let access ?(o=0) rc read read_prob =
  let m = Array.length read_prob - 1 in
  if rc then
    fun i ->
      let i = i + o in
      complement (String.get_exn read (m - i))
      , Array.get read_prob (m - i)
  else
    fun i ->
      let i = i + o in
      String.get_exn read i
      , Array.get read_prob i

(*** Full Forward Pass *)

type viterbi_result =
  { reverse_complement : bool             (* Was the reverse complement best? *)
  ; emission           : Lp.t                    (* Final emission likelihood *)
  ; path_list          : Path.t list
  }

let lookup_allele_or_fail t ~allele =
  (* Recover allele's base sequence, aka, emission. *)
  match Array.findi t.alleles ~f:(fun (s, _) -> s = allele) with
  | None              ->
      invalid_argf "%s not found among the list of alleles!" allele
  | Some allele_index ->
      allele_index

type viterbi_proc =
  { single  : read:bytes -> read_errors:float array -> viterbi_result
  ; paired  : read1:bytes -> read_errors1:float array
            -> read2:bytes -> read_errors2:float array
            -> viterbi_result
  }

let setup_single_allele_viterbi_pass ?insert_p ~prealigned_transition_model
  read_length ~allele t =
  let allele_index = lookup_allele_or_fail t ~allele in
  let allele_a =
    Array.fold_left t.emissions_a ~init:[] ~f:(fun acc pm ->
      match Pm.get pm allele_index with
      | None   -> acc            (* Gap *)
      | Some b -> b :: acc)
    |> List.rev
    |> Array.of_list
  in
  let ref_length = Array.length allele_a in
  let ws = ForwardSLogSpace.V.generate ~ref_length ~read_length in
  let transition_ref_length = Array.length t.emissions_a in
  let passes =
    ForwardSLogSpace.passes_v ?insert_p ~transition_ref_length
      ~prealigned_transition_model ~ref_length ~read_length ws allele_a
  in
  let result reverse_complement ws =
    let emission, pl = ForwardSLogSpace.V.get_emission ws in
    { emission; reverse_complement; path_list = List.rev pl }
  in
  let most_likely_viterbi vr1 vr2 =
    if vr1.emission > vr2.emission then
      vr1
    else
      vr2
  in
  let single ~read ~read_errors =
    let _ = passes.full (access false read read_errors) in
    let r = result false ws in
    let _ = passes.full (access true read read_errors) in
    let c = result true ws in
    most_likely_viterbi r c
  in
  let paired ~read1 ~read_errors1 ~read2 ~read_errors2 =
    let _ = passes.paired
              (access false read1 read_errors1)
              (access true read2 read_errors2)
    in
    let r = result false ws in
    let _ = passes.paired
              (access true read1 read_errors1)
              (access false read2 read_errors2)
    in
    let c = result true ws in
    most_likely_viterbi r c
  in
  { single; paired }


(* The forward pass return type requires a bit more subtlety as we want
   diagnostic ability; to debug and to pass along values for downstream filters.
   Therefore these passes return a structure that can (1) execute the pass
   (ex [single]) and (2) interrogate the result
   (ex. [maximum_positions_median_match]). *)

type proc =
  { name              : string
  (* Give it a name of what we're processing. *)

  ; init_global_state : unit -> Lp.t array
  (* Allocate the right size for global state of the per allele likelihoods. *)

  ; single            : ?prev_threshold:Lp.t
                      -> ?base_p:(Lp.t * int) mt
                      (* NOTE: Be careful about passing in a base_p. It could
                         slow down the calculation dramnatically, so have think
                         carefully about whether these values couldn't be
                         factored out. *)
                      -> read:string
                      -> read_errors:float array    (* also enforce that these are log p ?*)
                      -> bool
                      (* reverse_complement *)
                      -> unit Pass_result.t
  (* Perform a forward pass. We purposefully only signal whether we fully
     completed the pass or were filtered, because we may have different uses
     for the result of a pass. The accessors below allow the user to extract
     more purposeful information. *)

  (*; best_allele_pos   : int -> per_allele_datum list
     After we perform a forward pass we might be interested in either some
     diagnostic information such as which were the best alleles or where
     was the best alignment in the loci? These support a mode where we
     want to see how the reads "map" for different loci. *)

  ; per_allele_llhd_and_pos   : unit -> (Lp.t * int) mt
  (* If we're not interested in diagnostics, but just the end result, we want
     the per allele likelihood. *)

  ; maximum_positions_median_match : unit -> Lp.t
  (* This value is used as a short-circuiting threshold filter for subsequent
   * passes. For the position with the maximum emission, chose the median emission
   * probability. The maximum position allows us to situate the read as best as we
   * can, but choosing the median (as opposed to again the max) is a bit more robust
   * since there are alleles between _loci_ that can always match against a read
   * better than another loci. This provides a slightly more robust method so that we
   * don't short-circuit the evaluation too early. *)
  }

(* Single, for one allele, forward pass *)
let setup_single_allele_forward_pass ?insert_p ?max_number_mismatches
  ~prealigned_transition_model read_length allele t =
  let allele_index = lookup_allele_or_fail t ~allele in
  (* Recover allele's base sequence, aka, emission. *)
  let allele_a =
    Array.fold_left t.emissions_a ~init:[] ~f:(fun acc pm ->
      match Pm.get pm allele_index with
      | None   -> acc            (* Gap *)
      | Some b -> b :: acc)
    |> List.rev
    |> Array.of_list
  in
  let ref_length = Array.length allele_a in
  let ws = ForwardSLogSpace.W.generate ~ref_length ~read_length in
  let maximum_positions_median_match () =
    ForwardSLogSpace.W.fold_over_row ws (read_length - 1)
      ~init:Lp.zero ~f:(fun v c -> Lp.max v c.match_)
  in
  (* For comparison against all of the alleles we want to have the same
    transition probabilities, that depends upon the reference size.
    Therefore we'll use the reference size of all the alleles; size of
    emission_a.  Perhaps we want to expose a parameter to switch to the
    'local' reference size; size of allele_a. *)
  let transition_ref_length = Array.length t.emissions_a in
  let pass =
    ForwardSLogSpace.passes_s ?insert_p ?max_number_mismatches
      ~transition_ref_length
      ~prealigned_transition_model
      ~ref_length ~read_length ws allele_a
  in
  let single ?prev_threshold ?base_p ~read ~read_errors reverse_complement =
    (* Ignore base_p for the moment, as I can't think of a good reason to
        implement this logic in the single case. *)
    let read = access reverse_complement read read_errors in
    pass.full read
  in
  let per_allele_llhd_and_pos () =
    pm_init_all ~number_alleles:1 (ForwardSLogSpace.W.get_emission ws)
  in
  let init_global_state () = [| Lp.one |] in
  { name = sprintf "Single allele %s" allele
  ; single
  ; per_allele_llhd_and_pos
  ; init_global_state
  ; maximum_positions_median_match
  }

(* Return
  1. a function to process one read
  2. the workspace that it uses
  3. a function to combine results *)
let setup_single_pass ?band ?insert_p ?max_number_mismatches
  ~prealigned_transition_model read_length t =
  let { number_alleles; emissions_a; alleles; _ } = t in
  let ref_length = Array.length emissions_a in
  let module F = ForwardMLogSpace in
  let tm =
    if prealigned_transition_model then
      Phmm.TransitionMatrix.prealigned ~ref_length read_length
    else
      Phmm.TransitionMatrix.init ~ref_length read_length
  in
  let r(*, br*) = F.recurrences ?insert_p tm read_length number_alleles in
  let ws = F.W.generate ref_length read_length in
  let last_read_index = read_length - 1 in
  let per_allele_llhd_and_pos () = F.W.get_emission ws in
  let maximum_positions_median_match () =
    F.maximum_positions_median_match ws (read_length - 1)
  in
  let reference i = emissions_a.(i) in
  let init_global_state () = F.per_allele_emission_arr t.number_alleles in
  let normal () =
    (* TODO: Refactor this to be a bit more elegant, though it isn't obvious
       how to preserve the type variability of the filters. *)
    let single ?prev_threshold ?base_p ~read ~read_errors reverse_complement =
      (* We do not have to clear the workspace, since a full pass will
         overwrite all elements of the workspace.
        F.Workspace.clear ws;*)
      let read = access reverse_complement read read_errors in
      let wrap filter =
        Pass_result.wrap (fun () ->
          ignore (F.Pfc.full ?base_p ~filter ws r ~reference ~read))
      in
      match max_number_mismatches, prev_threshold with
      | None,                   None            ->
          wrap Filter.empty
      | None,                   Some threshold  ->
          let of_entry = F.pm_max_cell_by_match in
          wrap (F.Ff.Past_threshold.create threshold of_entry)
      | Some number_mismatches, None            ->
          let of_entry = F.pm_max_cell_by_match in
          wrap (F.Ff.Max_number_of_mismatches.create tm ~number_mismatches of_entry)
      | Some number_mismatches, Some threshold  ->
          let of_entry = F.pm_max_cell_by_match in
          let filter1 = F.Ff.Past_threshold.create threshold of_entry in
          let filter2 = F.Ff.Max_number_of_mismatches.create tm ~number_mismatches of_entry in
          let filter = Filter.join filter1 filter2 in
          wrap filter
    in
    { name = sprintf "Locus: %s" (Nomenclature.show_locus t.locus)
    ; single
    ; per_allele_llhd_and_pos
    ; maximum_positions_median_match
    ; init_global_state
    }
  in
  let banded c =
    failwith "NI"
    (*
    let single ?prev_threshold rc rd rd_errors =
      (* Clear the forward/final array since the banded pass algorithm relies
         on unfilled cells to indicate boundaries (places where we use
         heuristics). *)
      F.W.clear ws;
      let read = access rc rd rd_errors in
      F.Regrular.pass ws r ~reference ~read ~rows:c.warmup;
      F.Bands.pass c emissions_a increment_a ws (r, br) last_read_index read;
      Completed ()
    in
    { single
    ; best_alleles
    ; per_allele_llhd
    ; init_global_state
    ; maximum_positions_median_match (* TODO: Not quite right for bands. *)
    } *)
  in
  match band with
  | None                                    -> normal ()
  | Some c when c.warmup >= last_read_index -> normal ()
  | Some c                                  -> banded c

(* Create an array that maps indices of the references to positions that do
   not have gaps! This prevents us from starting looking for the likelihoood
   at an index where an allele has a gap and consequently has no starting
   likelihood. *)
let has_gap =
  Pm.fold_values ~init:false ~f:(fun b gbo ->
      match gbo with
      | Some _ -> b
      | None   -> true)

let non_gapped_emissions_arr arr =
  let _index, _latest_non_gap_index, non_gap_indices =
    Array.fold_left arr ~init:(0,0,[]) ~f:(fun (i, ng, acc) pm ->
      if has_gap pm then
        (i + 1, ng, ng :: acc)
      else
        (i + 1, i, i :: acc))
  in
  Array.of_list (List.rev non_gap_indices)

(* Because of the disconnect between performing the pass (single) and
   the accessors used to get information (per_allele_lhood ... etc)
   it is helpful to keep track of some state without passing
   something along recursively. *)

module Splitting_state = struct

  type t =
    { mutable emission_pm   : (Lp.t * int) mt
    ; mutable max_med_match : Lp.t
    (* The splitting passes incur an etra start/stop emission probability
      that makes them difficult to compare against the normal version.
      We're going to keep track of these values and then remove the
      probabilities from the final emission. *)
    ; mutable scaling       : Lp.t list
    }

  let init number_alleles =
    { emission_pm = pm_init_all ~number_alleles (Lp.one, -1)
    ; max_med_match = Lp.one
    ; scaling = []
    }

  let reset ss =
    ss.emission_pm <-
      pm_init_all ~number_alleles:(Pm.size ss.emission_pm) (Lp.one, -1);
    ss.max_med_match <- Lp.one;
    ss.scaling <- []

  let add_emission ss em mm =
    (* Choose second position, that is the one farther in the match. *)
    let m (l1,_p1) (l2,p2) = Lp.(l1 * l2), p2 in
    ss.emission_pm <- Pm.merge ~eq:Lpr.equal ss.emission_pm em m;
    ss.max_med_match <- Lp.(ss.max_med_match * mm)

  let add_scaling ss v =
    ss.scaling <- v :: ss.scaling

  let finish ss =
    let pr = List.fold_left ss.scaling ~init:Lp.one
        ~f:Lp.( * ) in
    ss.emission_pm <- Pm.map ss.emission_pm Lpr.equal
      ~f:(fun (e,p) -> Lp.(e / pr),p)

end (* Splitting_state *)

(* One of the big bottlenecks of the full forward pass is that as we get
   farther into the read the total variation across alleles grows (the length
   of the Partition_map.t). It is constrained by numerical checks for equality
   (dx), but this is insufficient. But PHMM's are (somewhat) linear so that
   P(R_i|A) ~ P(R_-i|A) P(R_+i|R_-i, A) where `R_-i` is the first half and
   `R_+i` is the second half. Since R_-i is half as long as R_i the total
   variation across alleles should be less. If we calculate the 2 halfs
   individually we can reduce the growth by combining only their final
   emissions.

   if number_of_splits doesn't equally divide read_length throw an error.
*)
let setup_splitting_pass ?band ?insert_p ?max_number_mismatches
  ~prealigned_transition_model read_length number_of_splits t =
  let { number_alleles; emissions_a; alleles; _ } = t in
  let ref_length = Array.length emissions_a in
  if read_length mod number_of_splits <> 0 then
    invalid_argf "number of splits %d doesn't equally divide read_length: %d"
      number_of_splits read_length
  else
    let eff_read_length = read_length / number_of_splits in
    let module F = ForwardMLogSpace in
    (* We'll be using different parts of the full reference length space during
       the different splits, this is big enough to contain everything. Note it
       is important to realize that the recurrences will override the appropriate
       parts of the reference-columns of the array, such that as we fill it in,
       we're only seeing the new values. This avoids a costly clearing of the
       space. *)
    let ws = F.W.generate ref_length eff_read_length in
    (* Set the emission and maximum match as references that we update after
       each split run. Only {single} defined below is different between the
       prealigned_transition_model case. *)
    let ss = Splitting_state.init number_alleles in
    let merge_after_pass ?range () =
      let em = F.W.get_emission ws in
      (* Close enough ? *)
      let max_med_em = F.maximum_positions_median_match ?range ws (read_length - 1) in
      Splitting_state.add_emission ss em max_med_em;
    in
    let per_allele_llhd_and_pos () = ss.Splitting_state.emission_pm in
    let maximum_positions_median_match () = ss.Splitting_state.max_med_match in
    let reference i = emissions_a.(i) in
    let init_global_state () = F.per_allele_emission_arr t.number_alleles in
    if not prealigned_transition_model then begin
      (* Easier to understand this version first:
         The tm and r are fixed. *)
      let tm = Phmm.TransitionMatrix.init ~ref_length eff_read_length in
      let r = F.recurrences ?insert_p tm eff_read_length number_alleles in
      let single ?prev_threshold ?base_p ~read ~read_errors reverse_complement =
        let wrap filter update =
          Splitting_state.reset ss;
          let rec loop n filter =
            if n = number_of_splits then
              Splitting_state.finish ss
            else begin
              let o = eff_read_length * n in
              let read = access ~o reverse_complement read read_errors in
              let nf_init = F.Pfc.full ~filter ws r ~reference ~read in
              let nfilter = update nf_init in
              if n > 0 then begin
                let open Phmm.TransitionMatrix in
                Splitting_state.add_scaling ss (Lp.constant (tm StartOrEnd Match));
                Splitting_state.add_scaling ss (Lp.constant (tm Match StartOrEnd))
              end;
              merge_after_pass ();
              loop (n + 1) nfilter
            end
          in
          Pass_result.wrap (fun () -> loop 0 filter)
        in
        match max_number_mismatches, prev_threshold with
        | None,                   None            ->
            wrap Filter.empty (Filter.update Filter.empty)
        | None,                   Some threshold  ->
            let of_entry = F.pm_max_cell_by_match in
            let filter = F.Ff.Past_threshold.create threshold of_entry in
            let update last_value =
              Filter.update filter (F.Ff.Past_threshold.update_init_value last_value)
            in
            wrap filter update
        | Some number_mismatches, None            ->
            let of_entry = F.pm_max_cell_by_match in
            let filter = F.Ff.Max_number_of_mismatches.create tm ~number_mismatches of_entry in
            let update = Filter.update filter in
            wrap filter update
        | Some number_mismatches, Some threshold  ->
            let of_entry = F.pm_max_cell_by_match in
            let filter1 = F.Ff.Past_threshold.create threshold of_entry in
            let filter2 = F.Ff.Max_number_of_mismatches.create tm ~number_mismatches of_entry in
            let filter = Filter.join filter1 filter2 in
            let update (pt_last_value, mnm_last_value) =
              let pt_new_init = F.Ff.Past_threshold.update_init_value pt_last_value in
              Filter.update filter (pt_new_init, mnm_last_value)
            in
            wrap filter update
      in
      let name =
        sprintf "Locus %s splitting %d (not prealigned)"
            (Nomenclature.show_locus t.locus) number_of_splits
      in
      { name
      ; single
      ; per_allele_llhd_and_pos
      ; maximum_positions_median_match
      ; init_global_state
      }
    end else begin (* prealigned_transition_model *)
      let update ~cs_start ~ref_length =
        (*printf "update cs_start: %d ref_length: %d eff_read_length: %d\n"
          cs_start ref_length eff_read_length; *)
        let tm = Phmm.TransitionMatrix.prealigned ~ref_length eff_read_length in
        let r = F.recurrences ?insert_p tm eff_read_length number_alleles in
        tm, r, Some { cs_start; ref_length }
      in
      let non_gap_emissions_a = non_gapped_emissions_arr emissions_a in
      let lookup_best_emissions n range =
        let cs_start =
          if n = 0 then
            0
          else begin
            let _highest_emission, best_pos =
              Pm.fold_values (per_allele_llhd_and_pos ()) ~init:(Lp.zero, -1)
                ~f:(fun (bem, bpos) (em, pos) ->
                      if Lp.is_gap em || not Lp.(bem < em) then
                        (bem, bpos)
                      else (* Lp.(bem < em) *)
                        (em, pos))
            in
            if best_pos = -1 then begin
              invalid_argf "range: %s\n"
                (Option.value_map range ~default:"None"
                    ~f:(fun (s,e) -> sprintf "Some (%d,%d)" s e))
            end else
              best_pos
          end
        in
        let non_gapped_cs_start = non_gap_emissions_a.(cs_start) in
        let ref_length = ref_length - non_gapped_cs_start in
        non_gapped_cs_start, ref_length
      in
      (* In order for this logic to work it needs to properly incorpoate gaps.
         Otherwise a gap of 15, not uncommon, seems to easily trigger this
         check on reads of 100 with split of 4 -> eff read_length 25.
      let big_jump cs_start =
        Option.iter ~f:(fun (previous_start, _) ->
          if previous_start <> 0  (* Not the first segment that locates the read! *)
          && cs_start - previous_start > 2 * eff_read_length then
            filter_is_past_thresholdf
              "new highest emissions %d, from %d separated by more than 2x \
              effective read length: %d."
                cs_start previous_start eff_read_length
          else
            ())
      in *)
      let single ?prev_threshold ?base_p ~read ~read_errors reverse_complement =
        (* We still do not have to clear the workspace, since a during the full
           we always start wide enough that we're clearing downstream cells
           that overwrite previous values.
          F.Workspace.clear ws;*)
        let wrap to_filter =
          Splitting_state.reset ss;
          let rec loop n prev_filter_value_opt range =
            if n = number_of_splits then
              Splitting_state.finish ss
            else
              let cs_start, ref_length = lookup_best_emissions n range in
              if ref_length <= 2 then
                filter_is_past_thresholdf "Best emission is at reference end before we \
                            matched entire read."
              else if ref_length <= eff_read_length then
                filter_is_past_thresholdf
                  "Remaining reference %d is less than or equal to remaining read length %d"
                    ref_length eff_read_length
              else begin
                (* big_jump cs_start range; See above *)
                let ntm, nr, columns = update ~cs_start ~ref_length in
                let o = eff_read_length * n in
                let read = access ~o reverse_complement read read_errors in
                let filter = to_filter ntm prev_filter_value_opt in
                let nf_init = F.Pfc.full ?columns ~filter ws nr ~reference ~read in
                if n > 0 then begin
                  let open Phmm.TransitionMatrix in
                  Splitting_state.add_scaling ss (Lp.constant (ntm StartOrEnd Match));
                  Splitting_state.add_scaling ss (Lp.constant (ntm Match StartOrEnd))
                end;
                let new_range = cs_start, cs_start + ref_length in
                merge_after_pass ~range:new_range ();
                loop (n + 1) (Some nf_init) (Some new_range)
              end
          in
          Pass_result.wrap (fun () -> loop 0 None None)
        in
        match max_number_mismatches, prev_threshold with
        | None,                   None            ->
            wrap (fun _ _ -> Filter.empty)
        | None,                   Some threshold  ->
            let of_entry = F.pm_max_cell_by_match in
            let to_filter _tm =
              let filter = F.Ff.Past_threshold.create threshold of_entry in
              function
              | None    -> filter
              | Some lv -> Filter.update filter (F.Ff.Past_threshold.update_init_value lv)
            in
            wrap to_filter
        | Some number_mismatches, None            ->
            let of_entry = F.pm_max_cell_by_match in
            let to_filter tm lv_opt =
              let filter = F.Ff.Max_number_of_mismatches.create tm ~number_mismatches of_entry in
              Option.value_map ~default:filter lv_opt ~f:(Filter.update filter)
            in
            wrap to_filter
        | Some number_mismatches, Some threshold  ->
            let of_entry = F.pm_max_cell_by_match in
            let to_filter tm lv_opt =
              let filter1 = F.Ff.Past_threshold.create threshold of_entry in
              let filter2 = F.Ff.Max_number_of_mismatches.create tm ~number_mismatches of_entry in
              let filter = Filter.join filter1 filter2 in
              Option.value_map ~default:filter lv_opt
                ~f:(fun (pt_last_value, mnm_last_value) ->
                      let pt_new_init = F.Ff.Past_threshold.update_init_value pt_last_value in
                      Filter.update filter (pt_new_init, mnm_last_value))
            in
            wrap to_filter
      in
      let name =
        sprintf "Locus %s splitting %d (prealigned)"
            (Nomenclature.show_locus t.locus) number_of_splits
      in
      { name
      ; single
      ; per_allele_llhd_and_pos
      ; maximum_positions_median_match
      ; init_global_state
      }
    end
