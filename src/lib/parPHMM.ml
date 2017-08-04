(* Parametric Profile Hidden Markov Model.

   "Parameterize" the match/insert/delete states by different alleles.
*)

open Util
module Pm = Partition_map

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

(*let dx = 2.22044e-16 *)
let dx = ref 1.e-6

let close_enough x y =
  abs_float (x -. y) < !dx

(* Probability Ring where we perform the forward pass calculation. *)
module type Ring = sig

  type t
  val to_string : t -> string
  val zero : t
  val one  : t

  val gap  : t
  val is_gap : t -> bool

  val ( + ) : t -> t -> t
  val ( * ) : t -> t -> t
  val max   : t -> t -> t

  val close_enough : t -> t -> bool

  (* Special constructs necessary for the probabilistic logic. *)
  (* Convert constant probabilities. *)
  val constant : float -> t

  (* Scale a probability be a third. *)
  val times_one_third : float -> t

  (* Complement probability. *)
  val complement_probability : float -> t

end (* Ring *)

module MultiplicativeProbability = struct
  type t = float
  let zero  = 0.
  let one   = 1.

  let gap   = nan
  let is_gap x = x <> x

  let ( + ) = ( +. )
  let ( * ) = ( *. )
  let max   = max

  let close_enough x y =
    close_enough x y

  let constant x = x

  let complement_probability p =
    1. -. p

  let times_one_third p =
    p /. 3.

  let to_string = sprintf "%f"

end (* MultiplicativeProbability *)

module LogProbabilities = struct

  type t = float

  let to_string = sprintf "%f"

  let zero  = neg_infinity
  let one   = 0.  (* log10 1. *)

  let gap   = nan
  let is_gap x = x <> x

  let exp10 x = 10. ** x

  let ( * ) lx ly = lx +. ly

  let ( + ) lx ly =
         if lx = neg_infinity then ly
    else if ly = neg_infinity then lx
    else if lx > ly           then lx +. log10 (1. +. exp10 (ly -. lx))
    else (* lx < ly *)             ly +. log10 (1. +. exp10 (lx -. ly))

  let max   = max

  let close_enough x y =
    close_enough x y

  let constant = log10

  let l13 = constant (1. /. 3.)

  let times_one_third = ( * ) l13

  let complement_probability lq =
    log10 (1. -. (exp10 lq))

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

end (* LogProbabilities *)

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
type 'a cell_recurrences =
  { start     : ?base_p:'a -> 'a -> 'a cell -> 'a cell
  ; fst_col   : 'a -> 'a cell -> 'a cell
  ; middle    : 'a -> match_c:('a cell)
                 -> insert_c:('a cell)
                 -> delete_c:('a cell)
                 -> 'a cell
  ; end_      : 'a cell -> 'a
  }

exception PastThreshold of string

type ('a, 'base, 'entry) filter =
  { init  : 'a
  ; base  : 'a -> 'base -> 'a
  ; entry : 'a -> 'entry -> 'a
  }

let join_filter f1 f2 =
  let init  = f1.init, f2.init in
  let base (s1, s2) b = f1.base s1 b, f2.base s2 b in
  let entry (s1, s2) b = f1.entry s1 b, f2.entry s2 b in
  { init ; base ; entry }

(** What are the values and equations that determine how probabilities are
    calculated in the forward pass. *)
module ForwardCalcs (R : Ring) = struct

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
        (*printf "Spotted an N with %f\n%!" base_error; *)
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

    let insert_p = constant insert_p in
    let start_i = insert_p * t_s_i in
    let start ?(base_p=R.one) emission_p prev_c =
      let r =
        if R.is_gap emission_p then
          prev_c
        else
          { match_ = emission_p * t_s_m * base_p
          ; insert = start_i * base_p
          ; delete = zero (* * base_p *)
          }
      in
      let () =
        if !debug_ref then
          printf "--------start: emission_p:%s\t insert_p%s:\n\tafter: %s\n"
            (R.to_string emission_p)
            (R.to_string insert_p)
            (cell_to_string R.to_string r)
      in
      r
    in
    let f = (* forward *)
      let fst_col emission_p insert_c =
        if R.is_gap emission_p then
          (* This is a weird case, for an allele there is a gap, in the first
              position. We'll 'fail' by propogating nan. *)
          { match_ = R.gap
          ; insert = start_i
          ; delete = R.gap
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
          if R.is_gap emission_p then
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
        let () =
          if !debug_ref then
            printf "--------middle: emission:%s \
                  \n\tmatch_: %s\n\tinsert: %s\n\tdelete: %s\n\tafter : %s\n"
              (R.to_string emission_p)
              (cell_to_string R.to_string match_c)
              (cell_to_string R.to_string insert_c)
              (cell_to_string R.to_string delete_c)
              (cell_to_string R.to_string r)
        in
        r
      in
      let end_ cell = cell.match_ * t_m_s + cell.insert * t_i_s in
      { start; fst_col; middle; end_ }
    in
    let v = (* viterbi *)
      let fst_col emission_p insert_c =
        if R.is_gap emission_p then
          (* This is a weird case, for an allele there is a gap, in the first
            position. We'll 'fail' by propogating nan. *)
          { match_ = R.gap
          ; insert = start_i
          ; delete = R.gap
          }
        else
          { match_ = emission_p * (max3 (t_m_m * zero)
                                        (t_i_m * zero)
                                        (t_d_m * zero))
          ; insert = insert_p   * (R.max (t_m_i * insert_c.match_)
                                         (t_i_i * insert_c.insert))
          ; delete = (* one *  *) (R.max (t_m_d * zero)
                                         (t_d_d * zero))
          }
      in
      let middle emission_p ~match_c ~insert_c ~delete_c =
        let r =
          if R.is_gap emission_p then
            delete_c
          else
            { match_ = emission_p * (max3 (t_m_m * match_c.match_)
                                          (t_i_m * match_c.insert)
                                          (t_d_m * match_c.delete))
            ; insert = insert_p   * (R.max (t_m_i * insert_c.match_)
                                           (t_i_i * insert_c.insert))
            ; delete = (* one *)    (R.max (t_m_d * delete_c.match_)
                                           (t_d_d * delete_c.delete))
            }
        in
        let () =
          if !debug_ref then
            printf "--------middle: emission:%s \
                  \n\tmatch_: %s\n\tinsert: %s\n\tdelete: %s\n\tafter : %s\n"
              (R.to_string emission_p)
              (cell_to_string R.to_string match_c)
              (cell_to_string R.to_string insert_c)
              (cell_to_string R.to_string delete_c)
              (cell_to_string R.to_string r)
        in
        r
      in
      let end_ cell = R.max (cell.match_ * t_m_s) (cell.insert * t_i_s) in
      { start; fst_col; middle; end_ }
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

end (* ForwardCalcs *)

module ForwardFilters  (R : Ring) = struct

  module Fc = ForwardCalcs(R)

  let debug_ref = ref false

  let past_threshold_filter threshold of_entry =
    let base (threshold, best_seen_entry) _base_doesnt_matter =
      if best_seen_entry <> R.zero && best_seen_entry < threshold then
        let msg =
          sprintf "threshold %s breached: %s"
            (R.to_string threshold)
            (R.to_string best_seen_entry)
        in
        raise (PastThreshold msg)
      else
        threshold, R.zero
    in
    let entry (threshold, best_entry_value) entry =
      (threshold, R.max best_entry_value (of_entry entry))
    in
    { init = (threshold, R.zero); base ; entry }

  (* Keep track of the encountered read errors, assume that {mismatches}
     of the most likely (lowest base error) bases are wrong. While the the rest
     (total number of bases - mismatches) are correct. This way we have a well
     defined value for the threshold. If all of the path probablities are below
     the threshold, stop executing! *)

  (* Assume that {errors} are in ascending order; therefore the lowest base
     error are at the front of the list. In this model/calculation there are
     only transitions between match states, so we don't need to scale by
     transition probabilities such as t_m_m. *)
  let ungapped_path_prob number_mismatches errors =
    let to_emission_p n e =
      if n < number_mismatches then Fc.mismatch_p e else Fc.match_p e
    in
    let rec loop n p l =
      match l with
      | []     -> p
      | e :: t -> let emission_p = to_emission_p n e in
                  loop (n + 1) R.(p * emission_p) t
    in
    match errors with
    | []     -> R.one
    | e :: t -> let emission_p = to_emission_p 0 e in
                loop 1 R.((* one *) emission_p) t

  (* Insert by maintaining the ascending order. *)
  let rec insert error l = match l with
    | []     -> error :: []
    | e :: t -> if error <= e then
                  error :: l
                else
                  e :: (insert error t)

  let match_filter ~number_mismatches of_entry =
    let path_prob = ungapped_path_prob number_mismatches in
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
        raise (PastThreshold (Lazy.force msg))
      else
        (insert base_error errors, R.zero)
    in
    let entry (es, mv) c = (es, R.max mv (of_entry c)) in
    { init  = ([], R.one)
    ; base
    ; entry
    }

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

  (* TODO: Consider renaming these to be more in line with 'i', 'k' *)
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
  val fold_over_row : t -> int -> init:'a -> f:('a -> entry -> 'a) -> 'a

  val fold_over_final : t -> init:'a -> f:('a -> final_entry -> 'a) -> 'a

  (* Reset for calculation. *)
  val clear : t -> unit

end (* Workspace_intf *)

module CommonWorkspace = struct

  (* The recurrence logic (see {middle} _below_) relies only upon the previous
     state of the read, or the row (the Markov in pphMm). Therefore for the
     bulk of the forward pass (ignoring final emission) we can conserve space
     by storing only 2 rows and alternating. Externally, the workspace allows
     access to {read_length} amount of rows; although this isn't enforced and
     is only a product of ForwardPass calling {rows} correctly.

     Internally, the workspace maps (mod) to the appriate offset into the 2
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

  let fold_over_row ws i ~init ~f =
    let fk = columns ws in
    let rec loop k acc =
      if k >= fk then acc else loop (k + 1) (f acc (get ws ~i ~k))
    in
    loop 0 init

  let fold_over_final ws ~init ~f =
    let fk = columns ws in
    let rec loop k acc =
      if k >= fk then acc else loop (k + 1) (f acc (get_final ws k))
    in
    loop 0 init

end (* CommonWorkspace *)

(* Create the workspace for the forward calculation for a single,
   aka reference, allele. *)
module SingleWorkspace (R : Ring) :
  (Workspace_intf with type entry = R.t cell
                   and type final_entry = R.t
                   and type final_emission = R.t) = struct

  module Fc = ForwardCalcs(R)

  type entry = R.t cell
  type final_entry = R.t
  type final_emission = R.t

  include CommonWorkspace
  type t = (entry, final_entry, final_emission) w

  let generate ~ref_length ~read_length =
    { forward  = Array.init 2 (*read_length*) ~f:(fun _ -> Array.make ref_length Fc.zero_cell)
    ; final    = Array.make ref_length R.zero
    ; emission = R.zero
    ; read_length
    }

  let clear ws =
    let ref_length  = Array.length ws.final in
    let number_rows = Array.length ws.forward in
    ws.forward  <- Array.init number_rows ~f:(fun _ -> Array.make ref_length Fc.zero_cell);
    ws.final    <- Array.make ref_length R.zero;
    ws.emission <- R.zero

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
          Single { reference = List.rev refa |> String.of_character_list
                 ; read = List.rev reada |> String.of_character_list
                 ; start ; end_ }
      | E end_ :: S nstart :: tl ->
          let r1 =
            { reference = List.rev refa |> String.of_character_list
            ; read = List.rev reada |> String.of_character_list
            ; start ; end_ }
          in
          let r2 = second nstart [] [] tl in
          Paired (r1, r2)
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

module SingleViterbiWorkspace (R : Ring) :
  (Workspace_intf with type entry = R.t cell * Path.t list
                   and type final_entry = R.t * Path.t list
                   and type final_emission = R.t * Path.t list) = struct

  module Fc = ForwardCalcs(R)

  type entry = R.t cell * Path.t list
  type final_entry = R.t * Path.t list
  type final_emission = R.t * Path.t list

  include CommonWorkspace
  type t = (entry, final_entry, final_emission) w

  let empty_entry = Fc.zero_cell, []
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
module MakeMultipleWorkspace (R : Ring) :
  (Workspace_intf with type entry = R.t cell mt
                   and type final_entry = R.t mt
                   and type final_emission = R.t mt) = struct

  type entry = R.t cell mt
  type final_entry = R.t mt
  type final_emission = R.t mt

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
  ; fst_col : 'workspace -> obs -> 'base -> i:int -> 'entry
  ; middle  : 'workspace -> obs -> 'base -> i:int -> k:int -> 'entry
  ; end_    : 'workspace -> int -> 'final_entry
  ; final_e : 'workspace -> 'final_emission
  }

(* Given:
    - a workspace
    - its dimensions
    - recurrence functions
    - read and references acccessors

   Produce functions that actually traverse and fill in the forward pass,
   and emission values. These functions also take an optional [filter] that
   may raise an exception (such as {PastThreshold}) to signal a stop of the
   computation. *)
module ForwardPass (W : Workspace_intf) = struct

  let empty_filter =
    { init  = ()
    ; base  = begin fun () _ -> () end
    ; entry = begin fun () _ -> () end
    }

  (* compute the full forward pass
     @param rows how many elements of the read to fill, defaults to the
        configured size of the workspace; which should be at most
        the length of the read.
     @param base_p base probability of a path; allows us to weight the 2nd read
        based upon the likelihood (final_emission) of the first. *)
  let pass_f ?rows ?base_p ~filter ws recurrences ~read ~reference =
    let columns = W.columns ws in
    let rows = Option.value rows ~default:(W.rows ws) in
    let a_0 = read 0 in
    let ft = ref (filter.base filter.init a_0) in
    for k = 0 to columns do
      let ne = recurrences.start ?base_p ws a_0 (reference k) ~k in
      ft := filter.entry !ft ne;
      W.set ws ~i:0 ~k ne
    done;
    for i = 1 to rows do
      let a_i = read i in
      ft := filter.base !ft a_i;
      let ne = recurrences.fst_col ws a_i ~i (reference 0) in
      ft := filter.entry !ft ne;
      W.set ws ~i ~k:0 ne;
      for k = 1 to columns do
        let ne = recurrences.middle ws a_i ~i ~k (reference k) in
        ft := filter.entry !ft ne;
        W.set ws ~i ~k ne
      done
    done;
    !ft

  let final ws recurrences =
    let columns = W.columns ws in
    for k = 0 to columns do
      W.set_final ws k (recurrences.end_ ws k)
    done

  (* Fill in both parts of the workspace. *)
  let both_f ?rows ?base_p ~filter ws recurrences ~read ~reference =
    let ft = pass_f ?rows ?base_p ~filter ws recurrences ~read ~reference in
    final ws recurrences;
    ft

  (* After filling in both parts of the workspace,
     compute the final emission value. *)
  let full_f ?rows ?base_p ~filter ws recurrences ~read ~reference =
    let ft = both_f ?rows ?base_p ~filter ws recurrences ~read ~reference in
    W.set_emission ws (recurrences.final_e ws);
    ft

  (* paired pass.  *)
  let paired_f ?rows ~filter ws recurrences ~read1 ~read2 ~reference =
      let nft = full_f ?rows ~filter ws recurrences ~read:read1 ~reference in
      let base_p = W.get_emission ws in
      let filter = { filter with init = nft } in
      full_f ?rows ~base_p ~filter ws recurrences ~read:read2 ~reference

  (* Without the filter. *)
  let pass ?base_p ?rows ws recurrences ~read ~reference =
    ignore (pass_f ?base_p ~filter:empty_filter ?rows ws recurrences ~read ~reference)

  let full ?base_p ?rows ws recurrences ~read ~reference =
    ignore (full_f ?base_p ~filter:empty_filter ?rows ws recurrences ~read ~reference)

  let paired ?rows ws recurrences ~read1 ~read2 ~reference =
    ignore (paired_f ?rows ~filter:empty_filter ws recurrences ~read1 ~read2 ~reference)

end (* ForwardPass *)

type 'a pass_result =
  | Completed of 'a
  | Filtered of string

let pass_result_to_string c_to_s = function
  | Completed c -> sprintf "Completed: %s" (c_to_s c)
  | Filtered m  -> sprintf "Filtered: %s" m

let map_completed ~f = function
  | Completed c -> Completed (f c)
  | Filtered m  -> Filtered m

let completed = function
  | Completed _ -> true
  | Filtered  _ -> false

(* A helper record to expose the functions generated by these modules. *)
type 'a passes =
  { full    : read:(int -> obs) -> 'a
  ; paired  : read1:(int -> obs) -> read2:(int -> obs) -> 'a
  }

module ForwardSingleGen (R: Ring) = struct

  module W = SingleWorkspace(R)
  module V = SingleViterbiWorkspace(R)
  module Fc = ForwardCalcs(R)
  module Ff = ForwardFilters(R)

  let debug_ref = ref false

  type base = Base.t

  let recurrences ?insert_p tm read_length =
    let r, _ = Fc.g ?insert_p tm read_length in
    let start ?base_p ws obsp base ~k =
      let prev_c = if k = 0 then Fc.zero_cell else (W.get ws ~i:0 ~k:(k-1)) in
      r.start ?base_p (Fc.to_match_prob obsp base) prev_c
    in
    let fst_col ws obsp base ~i =
      r.fst_col (Fc.to_match_prob obsp base) (W.get ws ~i:(i-1) ~k:0)
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
    let final_e ws =
      let open R in
      W.fold_over_final ws ~init:zero ~f:(+)
    in
    { start ; fst_col ; middle ; end_ ; final_e }

  module Fp = ForwardPass(W)

  let passes ?insert_p ?transition_ref_length ?max_number_mismatches ~read_length
    ws allele_a =
    let tm =
      let ref_length = Option.value transition_ref_length
        ~default:(Array.length allele_a)
      in
      Phmm.TransitionMatrix.init ~ref_length read_length in
    let recurrences = recurrences ?insert_p tm read_length in
    let reference k = allele_a.(k) in
    match max_number_mismatches with
    | None ->
        let full ~read =
          Fp.full ws recurrences ~reference ~read; Completed ()
        in
        let paired ~read1 ~read2 =
          Fp.paired ws recurrences ~reference ~read1 ~read2; Completed ()
        in
        { full; paired }
    | Some number_mismatches ->
        let full ~read =
          let of_entry c = c.match_ in
          let filter = Ff.match_filter ~number_mismatches of_entry in
          try
            let _filter_s = Fp.full_f ~filter ws recurrences ~reference ~read in
            Completed ()
          with PastThreshold msg ->
            Filtered msg
        in
        let paired ~read1 ~read2 =
          let of_entry c = c.match_ in
          let filter = Ff.match_filter ~number_mismatches of_entry in
          try
            let _filter_s = Fp.paired_f ~filter ws recurrences ~reference ~read1 ~read2 in
            Completed ()
          with PastThreshold msg ->
            Filtered msg
        in
        { full; paired }

  (* Viterbi *)
  let viterbi_recurrences ?insert_p tm read_length =
    let open Path in
    let _, r = Fc.g ?insert_p tm read_length in
    let cell_to_path i b r c =
      if c.match_ >= c.insert then begin
        if c.match_ >= c.delete then
          M (b, r)
        else
          D (i, r)
      end else (* c.match_ < c.insert *) begin
        if c.insert >= c.delete then
          I (i,b)
        else
          D (i, r)
      end
    in
    let start ?base_p ws obsp base ~k =
      let prev_c =
        if k = 0 then
          Fc.zero_cell
        else
          fst (V.get ws ~i:0 ~k:(k-1)) (* Ignore the path of the previous cell;
                                          not an actual transitin. *)
      in
      let base_pe, prev_path =
        Option.value_map base_p
          ~default:(R.one, [ S k])
          ~f:(fun (bp, pp) -> (bp, S k :: pp))
      in
      let nc = r.start ~base_p:base_pe (Fc.to_match_prob obsp base) prev_c in
      nc, (cell_to_path 0 (fst obsp) base nc :: prev_path)
    in
    let fst_col ws obsp base ~i =
      let pc, pl = V.get ws ~i:(i-1) ~k:0 in
      let nc = r.fst_col (Fc.to_match_prob obsp base) pc in
      nc, (cell_to_path i (fst obsp) base nc :: pl)
    in
    let middle ws obsp base ~i ~k =
      let emp = Fc.to_match_prob obsp base in
      let ks = k-1 in
      let insert_c, ilst = V.get ws ~i:(i-1) ~k in
      let match_c, mlst = V.get ws ~i:(i-1) ~k:ks in
      let delete_c, dlst = V.get ws ~i ~k:ks in
      let nc = r.middle emp ~insert_c ~delete_c ~match_c in
      let pt = cell_to_path i (fst obsp) base nc in
      if !debug_ref then
        printf "at i: %d k: %d nc: %s pt: %s (%d,%d,%d)\n"
          i k (cell_to_string R.to_string nc)
            (Path.to_string pt)
            (List.length mlst)
            (List.length ilst)
            (List.length dlst);
      let pl =
        match pt with         (* Don't track _all_ paths, choose one here *)
        | M _ -> pt :: mlst
        | I _ -> pt :: ilst
        | D _ -> pt :: dlst
        | S _
        | E _ -> assert false
      in
      nc, pl
    in
    let end_ ws k =
      let en, lst = V.get ws ~i:(read_length-1) ~k in
      r.end_ en, (E k :: lst) (* The list contains best because we chose at middle. *)
    in
    let final_e ws =
      let open R in
      V.fold_over_final ws ~init:(zero, [])
        ~f:(fun ((be, _) as b) ((fe, _) as e)->
              if fe > be then e else b)
    in
    { start ; fst_col ; middle ; end_ ; final_e }

  module Fpv = ForwardPass(V)

  let passes_v ?insert_p ?transition_ref_length ~read_length ws allele_a =
    let tm =
      let ref_length = Option.value transition_ref_length
        ~default:(Array.length allele_a)
      in
      Phmm.TransitionMatrix.init ~ref_length read_length in
    let recurrences = viterbi_recurrences ?insert_p tm read_length in
    let reference k = allele_a.(k) in
    let full ~read =
      Fpv.full ws recurrences ~reference ~read
    in
    let paired ~read1 ~read2 =
      Fpv.paired ws recurrences ~reference ~read1 ~read2
    in
    { full; paired }

end (* ForwardSingleGen *)

module PosMap = Map.Make(struct type t = int [@@deriving ord] end)

let topn p k a i lst =
  let rec loop added n lst =
    if n >= k then
      []
    else
      match lst with
      | []         -> if added then [] else [a,i]
      | (u,j) :: t -> if p a u && not added then
                        (a,i) :: loop true (n + 1) lst
                      else
                        (u,j) :: loop added (n + 1)  t
  in
  loop false 0 lst

let largest k a i lst = topn (>) k a i lst
let smallest k a i lst = topn (<) k a i lst

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
  ; width   : int     (* How many columns of the band to calculate. *)
  }

let band_default =
  { warmup  = 10
  ; number  = 5
  ; width   = 3
  }

module ForwardMultipleGen (R : Ring) = struct

  let debug_ref = ref false

  (* Eh... not the best place for it. *)
  let cam_max =
    Pm.fold_values ~init:R.zero ~f:(fun m v -> R.max m v)

  let cam_max_cell =
    Pm.fold_values ~init:R.zero ~f:(fun m c -> R.max m c.match_)

  module W = MakeMultipleWorkspace(R)

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

  module Fc = ForwardCalcs(R)
  module Ff = ForwardFilters(R)

  let recurrences ?insert_p tm read_length number_alleles =
    let r, _ = Fc.g ?insert_p tm read_length in

   (* TODO: I could imagine some scenario's where it makes sense to cache,
       precompute or memoize this calculation. The # of base errors isn't
       that large (<100) and there are only 4 bases. So we could be performing
       the same lookup. *)
    let to_em_set obsp emissions =
      Pm.map emissions ~f:(Option.value_map ~default:R.gap ~f:(Fc.to_match_prob obsp))
    in
    let zero_cell_pm = Pm.init_all_a ~size:number_alleles Fc.zero_cell in
    let eq = Fc.cells_close_enough in
    let start ?base_p ws obsp base ~k =
      let ems = to_em_set obsp base in
      let prev_pm = if k = 0 then zero_cell_pm else W.get ws ~i:0 ~k:(k-1) in
      let m2 =
        match base_p with
        | None    -> Pm.merge ems prev_pm r.start       (* Not weighing alleles *)
        | Some l  -> Pm.merge3 ~eq l ems prev_pm (fun base_p emission prev_c ->
                        r.start ~base_p emission prev_c)
      in
      (*printf "start:k: %d pm_length: %d: %s\n%!" k (Pm.length m2)
        (Pm.to_string m2 (cell_to_string R.to_string)); *)
      m2
    in
    let fst_col ws obsp emissions ~i =
      Pm.merge (to_em_set obsp emissions) (W.get ws ~i:(i-1) ~k:0)
        r.fst_col
    in
    let middle ws obsp emissions ~i ~k =
      let matches = W.get ws ~i:(i-1) ~k:(k-1) in
      let inserts = W.get ws ~i:(i-1) ~k       in
      let deletes = W.get ws ~i       ~k:(k-1) in
      let ems = to_em_set obsp emissions in
      (*printf "at i: %d k: %d: e: %d, m: %d, i: %d, d: %d \n%!"
        i k (Pm.length ems) (Pm.length matches) (Pm.length inserts)
            (Pm.length deletes); *)
      Pm.merge4 ~eq ems matches inserts deletes
        (fun emission_p match_c insert_c delete_c ->
          r.middle emission_p ~insert_c ~delete_c ~match_c)
    in
    let end_ ws k =
      Pm.map (W.get ws ~i:(read_length-1) ~k) ~f:r.end_
    in
    let final_e ws =
      (* CAN'T use empty_a since we're merging! *)
      let init = Pm.init_all_a ~size:number_alleles R.zero in
      W.fold_over_final ws ~init ~f:(fun e1 e2 -> Pm.merge e1 e2 R.(+))
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

  module Regular = ForwardPass(W)

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
      let ai = find_indices_above_n c.width emissions_a lnds in
      let bi = find_indices_below_n c.width increment_a lnds in
      (* tl_exn -> drop the center band, so we don't duplicate it. *)
      Cm.map2 ai bi ~f:(fun (st, a) (_st2, b) -> st, a @ (List.tl_exn (List.rev b)))

    (* This step (see get_exn) partitions the current bands based on the
      last cell's. Since these value are already away from the optimal point
      in the band (that would be width above), I'm not certain that this
      is necessary. We're using this value as just an approximation in lieu
      of filling in the forward matrix. Specifically, we can have a less strict
      cell equality test, so that the allele sets are joined together.

      This is a general notion to consider; when are cells close enough
      (difference just due to numerical rounding) that it isn't worth the
      split. *)
    let lookup_previous_values ws row bands =
      Cm.concat_map bands ~f:(fun s (bv, cols) ->
        let end_col = List.last cols |> Option.value_exn ~msg:"empty cols!" in
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
        (String.concat ~sep:";" (List.map t.cols ~f:(sprintf "%d")))
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
          width and do something like just move down 1 position per column but:
            - This doesn't account for gaps in the alleles.
            - This doesn't account for inserts/deletes that will shift the center
              of the band. In area's of ambiguity we could have 2 (or more)
              likelihood values that are close in value so we may err in
              determining the center.
          Therefore we need an adaptive strategy. We count the number of match
          values that are worse; where the likelihood is less than the best.
          Once this counter reaches the band config's width we stop considering
          those alleles.
        3. Keep track of the calculated cols for an allele. This allows us to
          adaptively figure out the cols of the next pass by moving width away
          from the best_col. See [find_next_row_from_fill_state].
      *)
    type fill_state =
      { best_col : int          (* where is the best, lowest match likelihood row *)
      ; best_c   : R.t cell
      ; worse    : int          (* Number of likelihoods < than best_c.match_ *)
      ; last_c   : R.t cell
      ; ncols    : int list     (* Where we're calculating. Since there might be
                                  gaps, the width needs to look inside this list
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
    let to_next_cols width best_col cols =
      let rec find_best acc = function
        | []     -> invalid_argf "Didn't find best row."
        | h :: t ->
            if h = best_col then
              (* TODO: These can silently take less than width. *)
              let s = List.take t width in
              let e = List.take acc width in
              (List.rev s) @ (h :: e)
            else
              find_best (h :: acc) t
      in
      find_best [] cols

    let find_next_row_from_fill_state c fs =
      to_next_cols c.width fs.best_col fs.ncols

    let next_band emissions_a increment_a c fs_map =
      Cm.concat_map fs_map ~f:(fun alleles fs ->
        (* Shift the band, by adjusting around best_col, for next column *)
        Cm.get_exn alleles increment_a.(fs.best_col)
        (* Now fill in the width. *)
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
          (String.concat ~sep:";" (List.map b.cols ~f:(sprintf "%d"))) i
            (to_string b)
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
                      if fs.worse >= c.width then `Fst fs else `Snd fs)
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

(* Instantiate the actual passes over the 2 implemented probability rings. *)
module ForwardS = ForwardSingleGen(MultiplicativeProbability)
module ForwardSLogSpace = ForwardSingleGen(LogProbabilities)

module ForwardM = ForwardMultipleGen(MultiplicativeProbability)
module ForwardMLogSpace = ForwardMultipleGen (LogProbabilities)

type t =
  { align_date      : string
  ; alleles         : (string * Alter_MSA.info) array     (* Canonical order. *)
  ; number_alleles  : int
  ; emissions_a     : base_emissions array
  }

let construct input =
  if not (Alleles.Input.is_imputed input) then
    invalid_argf "Allele input MUST be imputed!"
  else begin
    let open MSA.Parser in
    Alleles.Input.construct input >>= fun (mp, merge_map) ->
      let { reference; ref_elems; alt_elems; align_date} = mp in
      let base_arr, pmap = initialize_base_array_and_position_map ref_elems in
      let state_a = init_state base_arr in
      List.iter alt_elems ~f:(fun (allele, allele_seq) ->
        add_alternate_allele pmap allele allele_seq
          base_arr state_a);
      let emissions_a = Array.map Pm.ascending state_a in
      let alleles =
        (reference, Alter_MSA.FullSequence) ::
         (List.map alt_elems ~f:(fun (a,_) -> a,
           List.Assoc.get a merge_map
           |> Option.value_exn ~msg:(sprintf "couldn't find %s in merge_map" a)))
        |> Array.of_list
      in
      let number_alleles = Array.length alleles in
      Ok { align_date; alleles; number_alleles; emissions_a }
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
  |> List.map ~f:(sprintf "%f")
  |> String.concat ~sep:"; "
  |> sprintf "[|%s|]"

(* Abstracts, behind a function the regular or reverse complement access
   pattern of a read. This is to avoid manually reversing and converting
   the read. *)
let access rc read read_prob =
  let m = Array.length read_prob - 1 in
  if rc then
    fun i ->
      complement (String.get_exn read (m - i))
      , Array.get read_prob (m - i)
  else
    fun i -> String.get_exn read i , Array.get read_prob i

(*** Full Forward Pass *)

type viterbi_result =
  { reverse_complement : bool       (* Was the reverse complement best? *)
  ; emission           : float      (* Final emission likelihood *)
  ; path_list          : Path.t list
  }

let setup_single_allele_viterbi_pass ?insert_p read_length ~allele t =
  let { emissions_a; alleles; _ } = t in
  (* Recover allele's base sequence, aka, emission. *)
  match Array.findi alleles ~f:(fun (s, _) -> s = allele) with
  | None              ->
      invalid_argf "%s not found among the list of alleles!" allele
  | Some allele_index ->
      let allele_a =
        Array.fold_left emissions_a ~init:[] ~f:(fun acc pm ->
          match Pm.get pm allele_index with
          | None   -> acc            (* Gap *)
          | Some b -> b :: acc)
        |> List.rev
        |> Array.of_list
      in
      let ref_length = Array.length allele_a in
      let transition_ref_length = Array.length emissions_a in
      let ws = ForwardSLogSpace.V.generate ~ref_length ~read_length in
      let passes =
        ForwardSLogSpace.passes_v ?insert_p ~transition_ref_length ~read_length
          ws allele_a
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
        passes.full (access false read read_errors);
        let r = result false ws in
        passes.full (access true read read_errors);
        let c = result true ws in
        most_likely_viterbi r c
      in
      let paired ~read1 ~read_errors1 ~read2 ~read_errors2 =
        passes.paired (access false read1 read_errors1) (access true read2 read_errors2);
        let r = result false ws in
        passes.paired (access true read1 read_errors1) (access false read2 read_errors2);
        let c = result true ws in
        most_likely_viterbi r c
      in
      single, paired

(* The forward pass return type requires a bit more subtlety as we want
   diagnostic ability; to debug and to pass along values for downstream filters.
   Therefore these passes return a structure that can (1) execute the pass
   (ex [single]) and (2) interrogate the result (ex [best_alleles]). *)

type proc =
  { init_global_state : unit -> float array
  (* Allocate the right size for global state of the per allele likelihoods. *)

  ; single            : ?prev_threshold:float
                      -> ?base_p:float mt             (* ascending partition map. *)
                      (* NOTE: Be careful about passing in a base_p. It could
                         slow down the calculation dramnatically, so have think
                         carefully about whether these values couldn't be
                         factored out. *)
                      -> reverse_complement:bool
                      -> read:string
                      -> read_errors:float array
                      -> unit pass_result
  (* Perform a forward pass. We purposefully only signal whether we fully
     completed the pass or were filtered, because we may have different uses
     for the result of a pass. The accessors below allow the user to extract
     more purposeful information. The method is called "single" to contrast
     with [paired] which performs the forward pass for paired reads. *)

  ; paired            : ?prev_threshold:float
                      -> reverse_complement:bool
                      -> read1:string
                      -> read_errors1:float array
                      -> read2:string
                      -> read_errors2:float array
                      -> unit pass_result
  (* Perform a forward pass for paired reads. Takes care of flipping either
     the first or second read information depending on reverse complement. *)

  ; best_alleles    : int -> (float * string) list
  ; best_positions  : int -> (float * int) list
  ; best_allele_pos : int -> (float * string * int) list
  (* After we perform a forward pass we might be interested in either some
     diagnostic information such as which were the best alleles or where
     was the best alignment in the loci? These support a mode where we
     want to see how the reads "map" for different loci. *)

  ; per_allele_llhd : unit -> float mt
  (* If we're not interested in diagnostics, but just the end result, we want
     the per allele likelihood: *)

  (* For paired reads when we want to know the next base. *)
  ; maximum_match   : unit -> float
  (* We might also want to get a possible threshold value for future passes. *)
  }

(* Single, for one allele, forward pass *)
let setup_single_allele_forward_pass ?insert_p ?max_number_mismatches
  read_length allele t =
  let { emissions_a; alleles; _ } = t in
  (* Recover allele's base sequence, aka, emission. *)
  match Array.findi alleles ~f:(fun (s, _) -> s = allele) with
  | None              ->
      invalid_argf "%s not found among the list of alleles!" allele
  | Some allele_index ->
      let allele_a =
        Array.fold_left emissions_a ~init:[] ~f:(fun acc pm ->
          match Pm.get pm allele_index with
          | None   -> acc            (* Gap *)
          | Some b -> b :: acc)
        |> List.rev
        |> Array.of_list
      in
      let ref_length = Array.length allele_a in
      let ws = ForwardSLogSpace.W.generate ~ref_length ~read_length in
      let best_positions n =               (* 'n' useless when single allele. *)
        let lgn = largest n in
        ForwardSLogSpace.W.fold_over_final ws ~init:(0, [])
          ~f:(fun (p, acc) l -> (p + 1, lgn l p acc))
        |> snd
      in
      let maximum_match () =
        ForwardSLogSpace.W.fold_over_row ws (read_length - 1) ~init:neg_infinity
          ~f:(fun v c -> max v c.match_)
      in
      (* For comparison against all of the alleles we want to have the same
        transition probabilities, that depends upon the reference size.
        Therefore we'll use the reference size of all the alleles; size of
        emission_a.  Perhaps we want to expose a parameter to switch to the
        'local' reference size; size of allele_a. *)
      let transition_ref_length = Array.length emissions_a in
      let pass =
        ForwardSLogSpace.passes ?insert_p ?max_number_mismatches
          ~transition_ref_length ~read_length ws allele_a
      in
      let single ?prev_threshold ?base_p ~reverse_complement ~read ~read_errors =
        (* Ignore base_p for the moment, as I can't think of a good reason to 
           implement this logic in the single case. *)
        let read = access reverse_complement read read_errors in
        pass.full read
      in
      let paired ?prev_threshold ~reverse_complement
          ~read1 ~read_errors1 ~read2 ~read_errors2 =
        let read1 = access reverse_complement read1 read_errors1 in
        let read2 = access (not reverse_complement) read2 read_errors2 in
        pass.paired read1 read2
      in
      let best_alleles _n = [ ForwardSLogSpace.W.get_emission ws, allele] in
      let best_allele_pos n =
        best_positions n |> List.map ~f:(fun (l, i) -> (l, allele, i))
      in
      let per_allele_llhd () = Pm.init_all_a ~size:1 (ForwardSLogSpace.W.get_emission ws) in
      let init_global_state () = [| LogProbabilities.one |] in
      { single; paired
      ; best_alleles; best_positions; best_allele_pos
      ; per_allele_llhd
      ; init_global_state
      ; maximum_match
      }

(* Return
  1. a function to process one read
  2. the workspace that it uses
  3. a function to combine results *)
let setup_single_pass ?band ?insert_p ?max_number_mismatches read_length t =
  let { number_alleles; emissions_a; alleles; _ } = t in
  let ref_length = Array.length emissions_a in
  let tm = Phmm.TransitionMatrix.init ~ref_length read_length in
  let module F = ForwardMLogSpace in
  let r(*, br*) = F.recurrences ?insert_p tm read_length number_alleles in
  let ws = F.W.generate ref_length read_length in
  let last_read_index = read_length - 1 in
  let best_alleles n =
    let alleles = Array.map ~f:fst alleles in
    let emission_pm = F.W.get_emission ws in
    let lgn = largest n in
    Pm.fold_indices_and_values emission_pm ~init:[]
      ~f:(fun acc i emission -> lgn emission alleles.(i) acc)
  in
  let best_allele_pos n =
    let alleles = Array.map ~f:fst alleles in
    let emission_pm = F.W.get_emission ws in
    let lgn a p = largest n a p in
    Pm.fold_indices_and_values emission_pm ~init:[]
      ~f:(fun acc i emission ->
            let _, bplst =
              F.W.fold_over_final ws ~init:(0, [])
                ~f:(fun (p, acc) fcam -> (p + 1, lgn (Pm.get fcam i) p acc))
            in
            let _best_prob, best_pos = List.hd_exn bplst in
            lgn emission (alleles.(i), best_pos) acc)
    |> List.map ~f:(fun (e, (a, p)) -> (e, a, p))
  in
  let best_positions n =
    let lgn = largest n in
    F.W.fold_over_final ws ~init:(0, [])
      ~f:(fun (p, acc) fcam ->
        (p + 1, lgn (F.cam_max fcam) p acc))
    |> snd
  in
  let maximum_match () =
    F.W.fold_over_row ws (read_length - 1)
      ~init:neg_infinity ~f:(fun v fcam -> max v (F.cam_max_cell fcam))
  in
  let reference i = emissions_a.(i) in
  let per_allele_llhd () = F.W.get_emission ws in
  let init_global_state () = F.per_allele_emission_arr t.number_alleles in
  let normal () =
    (* TODO: Refactor this to be a bit more elegant, though it isn't obvious
       how to preserve the type variability of the filters. *)
    let single ?prev_threshold ?base_p ~reverse_complement ~read ~read_errors =
      (* We do not have to clear the workspace, since a full pass will
         overwrite all elements of the workspace.
        F.Workspace.clear ws;*)
      let read = access reverse_complement read read_errors in
      let unfiltered () =
        F.Regular.full ?base_p ws r ~reference ~read;
        Completed ()
      in
      let filtered ~filter =
        try
          let _final_filter = F.Regular.full_f ?base_p ~filter ws r ~reference ~read in
          Completed ()
        with PastThreshold msg ->
          Filtered msg
      in
      match max_number_mismatches, prev_threshold with
      | None,                   None            ->
          unfiltered ()
      | None,                   Some threshold  ->
          let of_entry = F.cam_max_cell in
          let filter = F.Ff.past_threshold_filter threshold of_entry in
          filtered ~filter
      | Some number_mismatches, None            ->
          let of_entry = F.cam_max_cell in
          let filter = F.Ff.match_filter ~number_mismatches of_entry in
          filtered ~filter
      | Some number_mismatches, Some threshold  ->
          let of_entry = F.cam_max_cell in
          let filter1 = F.Ff.past_threshold_filter threshold of_entry in
          let filter2 = F.Ff.match_filter ~number_mismatches of_entry in
          let filter = join_filter filter1 filter2 in
          filtered ~filter
    in
    let paired ?prev_threshold ~reverse_complement ~read1 ~read_errors1
      ~read2 ~read_errors2 =
      (* We do not have to clear the workspace, since a full pass will
         overwrite all elements of the workspace.
        F.Workspace.clear ws;*)
      let read1 = access reverse_complement read1 read_errors1 in
      let read2 = access (not reverse_complement) read2 read_errors2 in
      let unfiltered () =
        F.Regular.paired ws r ~reference ~read1 ~read2;
        Completed ()
      in
      let filtered ~filter =
        try
          let _final_filter = F.Regular.paired_f ~filter ws r ~reference ~read1 ~read2 in
          Completed ()
        with PastThreshold msg ->
          Filtered msg
      in
      match max_number_mismatches, prev_threshold with
      | None,                   None            ->
          unfiltered ()
      | None,                   Some threshold  ->
          let of_entry = F.cam_max_cell in
          let filter = F.Ff.past_threshold_filter threshold of_entry in
          filtered ~filter
      | Some number_mismatches, None            ->
          let of_entry = F.cam_max_cell in
          let filter = F.Ff.match_filter ~number_mismatches of_entry in
          filtered ~filter
      | Some number_mismatches, Some threshold  ->
          let of_entry = F.cam_max_cell in
          let filter1 = F.Ff.past_threshold_filter threshold of_entry in
          let filter2 = F.Ff.match_filter ~number_mismatches of_entry in
          let filter = join_filter filter1 filter2 in
          filtered ~filter
    in
    { single; paired
    ; best_alleles ; best_positions ; best_allele_pos
    ; per_allele_llhd
    ; maximum_match
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
      F.Regular.pass ws r ~reference ~read ~rows:c.warmup;
      F.Bands.pass c emissions_a increment_a ws (r, br) last_read_index read;
      Completed ()
    in
    { single
    ; best_alleles
    ; best_positions
    ; per_allele_llhd
    ; init_global_state
    ; maximum_match (* TODO: Not quite right for bands. *)
    } *)
  in
  match band with
  | None                                    -> normal ()
  | Some c when c.warmup >= last_read_index -> normal ()
  | Some c                                  -> banded c
