(* Parametric Profile Hidden Markov Model.
   "Parameterize" the match/insert/delete states by different alleles.

  TODO: Remove type annotations and transition to an interface.
*)

open Util

let debug_ref = ref false

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

  (* Assume that we're working on fully imputed sequences, so no 'Unknown'
     states for an allele. *)
  type t =
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
    | A       -> 'A'
    | C       -> 'C'
    | G       -> 'G'
    | T       -> 'T'

end (* BaseState *)

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
let initialize_base_array_and_position_map aset reference ref_elems =
  let module AS = (val aset : Alleles.Set) in
  let open Mas_parser in
  let ref_set () = AS.singleton reference in
  let sequence_to_base_states_array prev_state s =
    String.to_character_list s
    |> List.mapi ~f:(fun i c ->
        let p = if i = 0 then prev_state else -1 in
        let b = BaseState.of_char c in
        [(b, p), ref_set ()])
    |> Array.of_list
  in
  let gap_to_base_states_array len = Array.make len [] in
  let open Mas_parser in
  let rec loop lkmp p pp pacc acc = function
    | End _pos :: []              -> Array.concat (List.rev acc),
                                     List.rev ((lkmp, p) :: pacc)
    | Start pos :: t              -> invalid_argf "initialize_base_array_and_position_map: second start: %d"
                                        pos
    | End pos :: t                -> invalid_argf "initialize_base_array_and_position_map: end with more %d"
                                        pos
    | []                          -> invalid_argf "initialize_base_array_and_position_map: list before End"

    | Boundary { pos; _ } :: t    -> loop (pos + 1) (p + pos - lkmp) pp ((lkmp, p) :: pacc) acc t
    | Sequence { s; start } :: t  -> let l = String.length s in
                                     loop (start + l) (p + l) (-1) ((lkmp, p) :: pacc)
                                        (sequence_to_base_states_array pp s :: acc) t
    | Gap { length; gstart } :: t -> loop (gstart + length) (p + length) (-length-1) ((lkmp, p) :: pacc)
                                        (gap_to_base_states_array length :: acc) t
  in
  match ref_elems with
  | Start s :: t -> loop s 0 min_int [] [] t
  | e :: _       -> invalid_argf "Reference not at Start : %s" (al_el_to_string e)
  | []           -> invalid_argf "Empty reference sequence!"

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
   the base state array. This is useful for debugging between the Mas_parser
   positions and the index into BaseState array. Assumes a 'reduced' (via
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

(* Add an allele's Mas_parser.sequence to the current run-length encoded state. *)
let add_alternate_allele aset reference ~position_map allele allele_instr arr =
  let module AS = (val aset : Alleles.Set) in
  let base_and_offset b o ((bp,bo), _) = b = bp && o = bo in
  let add_to_base_state i b o =
    (*printf "adding base at %d %c %d\n" i (BaseState.to_char b) o; *)
    match List.find arr.(i) ~f:(base_and_offset b o) with
    | None              -> let s = AS.singleton allele in
                           (*printf "single cell at %d for %s \n"  i allele; *)
                           arr.(i) <- ((b, o), s) :: arr.(i)
    | Some ((ab,ro), s) -> (*printf "At: %d %s to %s\n"  i allele (AS.to_string s); *)
                           ignore (AS.set s allele)
  in
  let has_reference_set (_, s) = AS.is_set s reference in
  let add_to_reference_set offset start end_ =
    (*printf "adding reference set at %d %d %d\n" offset start end_; *)
    let rec loop i offset =
      if i = end_ then offset else begin
        match List.find arr.(i) ~f:has_reference_set with
        | None              -> (* Reference in gap -> so are we! *)
                               loop (i + 1) (-1)  (* Offset becomes -1 after 1st use! *)
        | Some ((rb,ro), s) ->
            if ro = offset then begin
              ignore (AS.set s allele);
              loop (i + 1) (-1)
            end else begin
              add_to_base_state i rb offset;
              loop (i + 1) (-1)
            end
      end
    in
    loop start offset
  in
  let open Mas_parser in
  let rec loop position_map lp ~offset = function
    | End p :: []                 ->
        let ap, position_map = position_and_advance p position_map in
        let _final_offset = add_to_reference_set offset lp ap in
        (* Check that ap = last_pos ? *)
        ()
    | Start p :: _                -> invalid_argf "add_alternate_allele: second start: %d" p
    | []                          -> invalid_argf "add_alternate_allele: didn't End"
    | End p :: t                  -> invalid_argf "add_alternate_allele: end before end: %d." p

    | Boundary { pos; _ } :: t    ->
        let ap, position_map = position_and_advance pos position_map in
        loop position_map ap ~offset:(add_to_reference_set offset lp ap) t
    | Sequence { start; s } :: t  ->
        let ap, position_map = position_and_advance start position_map in
        let noffset = add_to_reference_set offset lp ap in
        let fap, foffset =
          String.fold s ~init:(ap, noffset) ~f:(fun (p, o) c ->
            add_to_base_state p (BaseState.of_char c) o;
            (p + 1, -1))
        in
        loop position_map fap ~offset:foffset t
    | Gap { gstart; length } :: t ->
        let ap, position_map = position_and_advance gstart position_map in
        let _noffset = add_to_reference_set offset lp ap in
        (* The Gap determines new offsets! *)
        loop position_map (ap + length) ~offset:(-length - 1) t
  in
  match allele_instr with
  | Start s :: t -> loop position_map 0 ~offset:min_int t
  | e :: _       -> invalid_argf "add_alternate_allele: Allele %s not at Start : %s" allele (al_el_to_string e)
  | []           -> invalid_argf "add_alternate_allele: Empty allele %s sequence!" allele

(***** Forward Pass ***)

type set = Alleles.set

(* Probability Ring where we perform the forward pass calculation. *)
module type Ring = sig

  type t
  val to_string : t -> string
  val zero : t
  val one  : t

  val ( + ) : t -> t -> t
  val ( * ) : t -> t -> t
  val max   : t -> t -> t

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
  let ( + ) = ( +. )
  let ( * ) = ( *. )
  let max   = max

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

  let exp10 x = 10. ** x

  let ( * ) lx ly = lx +. ly

  let ( + ) lx ly =
         if lx = neg_infinity then ly
    else if ly = neg_infinity then lx
    else if lx > ly           then lx +. log10 (1. +. exp10 (ly -. lx))
    else (* lx < ly *)             ly +. log10 (1. +. exp10 (lx -. ly))

  let max   = max

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

type 'a cell_recurrences =
  { start   : 'a -> 'a cell
  ; fst_col : 'a -> 'a cell -> 'a cell
  ; middle  : 'a -> insert_c:('a cell)
                 -> delete_c:('a cell)
                 -> match_c:('a cell)
                 -> 'a cell
  ; end_    : 'a cell -> 'a
  }

exception PastThreshold of string

type ('a, 'base, 'entry) filter =
  { init  : 'a
  ; base  : 'a -> 'base -> 'a
  ; entry : 'a -> 'entry -> 'a
  }

(** What are the values and equations that determine how probabilities are
    calculated in the forward pass. *)
module ForwardCalcs  (R : Ring) = struct

  let mismatch_p =
    R.times_one_third

  let match_p =
    R.complement_probability

  (* TODO. Avoid the `float_of_int (Phred_score.to_int c) /. -10.` round trip
      between converting to log10p and then back to log10, and just use char's
      instead for the quality calc. *)
  let to_match_prob (base, base_error) =
    let compare_against c =
      if base = c then
        match_p base_error
      else
        mismatch_p base_error
    in
    let open BaseState in
    function
    | A -> compare_against 'A'
    | C -> compare_against 'C'
    | G -> compare_against 'G'
    | T -> compare_against 'T'

  let debug_ref = ref false

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
    { start   = begin fun emission_p ->
                  let r =
                    { match_ = emission_p * t_s_m
                    ; insert = start_i
                    ; delete = zero
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
                end
    ; fst_col = begin fun emission_p insert_c ->
                  { match_ = emission_p * ( t_m_m * zero
                                          + t_i_m * zero
                                          + t_d_m * zero)
                  ; insert = insert_p   * ( t_m_i * insert_c.match_
                                          + t_i_i * insert_c.insert)
                  ; delete = (* one *  *) ( t_m_d * zero
                                          + t_d_d * zero)
                  }
                end
    ; middle  = begin fun emission_p ~insert_c ~delete_c ~match_c ->
                  let r =
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
                 end
   ; end_    = begin fun cell ->
                  cell.match_ * t_m_s + cell.insert * t_i_s
               end
   }

  let zero_cell =
    { match_ = R.zero
    ; insert = R.zero
    ; delete = R.zero
    }

  (* Keep track of the encountered read errors, assume that {mismatches}
     of the most likely (lowest base error) bases are wrong. While the the rest
     (total number of bases - mismatches) are correct. This way we have a well
     defined value for the threshold. If all of the path probablities are below
     the threshold, stop executing! *)

  (* Assume that {errors} are in ascending order; therefore the lowest base
     error are at the front of the list. In this model/calculation there are
     only transitions between match states, so we don't need to scale by
     transition probabilities such as t_m_m. *)
  let path_prob number_mismatches errors =
    let to_emission_p n e =
      if n < number_mismatches then mismatch_p e else match_p e
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
    let path_prob = path_prob number_mismatches in
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

end (* ForwardCalcs *)

type ('entry, 'final_entry, 'final_emission) workspace =
  { mutable forward   : 'entry array array
  ; mutable final     : 'final_entry array
  ; mutable emission  : 'final_emission
  }

(* Layout logic:
   The dynamic array consists of values for each base in the read in a row
   and for each value (base, or list of bases) of the reference in a column.

   'i' is used to index into the read, hence row.
   'k' is used to index into the reference, hence column
*)

(* Create the workspace for the forward calculation for a single,
   aka reference, allele. *)
module SingleWorkspace (R : Ring) = struct

  module Fc = ForwardCalcs(R)

  type entry = R.t cell
  type final_entry = R.t

  type t = (entry, final_entry, R.t) workspace

  let generate ~ref_length ~read_length =
    { forward  = Array.init read_length ~f:(fun _ -> Array.make ref_length Fc.zero_cell)
    ; final    = Array.make ref_length R.zero
    ; emission = R.zero
    }

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

(* Create and manage the workspace for the forward pass for multiple alleles.*)
module MakeMultipleWorkspace (R : Ring)(Cm : CAM.M) = struct

  type entry = R.t cell CAM.t
  type final_entry = R.t CAM.t

  type t = (entry, final_entry, R.t array) workspace

  let generate number_alleles ~ref_length ~read_length =
    { forward   = Array.init read_length ~f:(fun _ -> Array.make ref_length Cm.empty)
    ; final     = Array.make ref_length Cm.empty
    ; emission  = Array.make number_alleles R.zero
    }

  let clear ws =
    let ref_length       = Array.length ws.final in
    let read_length      = Array.length ws.forward in
    let number_alleles = Array.length ws.emission in
    ws.forward <- Array.init read_length ~f:(fun _ -> Array.make ref_length Cm.empty);
    ws.final   <- Array.make ref_length Cm.empty;
    Array.fill ws.emission ~pos:0 ~len:number_alleles R.zero

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
  { start   :  obs -> 'base -> 'entry
  ; fst_col : 'workspace -> obs -> 'base -> i:int -> 'entry
  ; middle  : 'workspace -> obs -> 'base -> i:int -> k:int -> 'entry
  ; end_    : 'workspace -> int -> 'final_entry
  ; final_e : 'workspace -> 'final_emission
  }

(* Given:
    - a workspace
    - it's dimensions
    - recurrence functions
    - read and references acccessors

   Produce functions that actually traverse and fill in the forward pass,
   and emission values. These functions also take an optional [filter] that
   may raise an exception (such as {PastThreshold}) to signal a stop of the
   computation. *)
module ForwardPass = struct

  let last_array_index arr =
    Array.length arr - 1

  let empty_filter =
    { init  = ()
    ; base  = begin fun () _ -> () end
    ; entry = begin fun () _ -> () end
    }

  let pass_f ?rows ~filter ws recurrences ~read ~reference =
    let columns = last_array_index ws.forward.(0) in
    let rows = Option.value rows ~default:(last_array_index ws.forward) in
    let a_0 = read 0 in
    let ft = ref (filter.base filter.init a_0) in
    for k = 0 to columns do
      let ne = recurrences.start a_0 (reference k) in
      ft := filter.entry !ft ne;
      ws.forward.(0).(k) <- ne
    done;
    for i = 1 to rows do
      let a_i = read i in
      ft := filter.base !ft a_i;
      let ne = recurrences.fst_col ws a_i ~i (reference 0) in
      ft := filter.entry !ft ne;
      ws.forward.(i).(0) <- ne;
      for k = 1 to columns do
        let ne = recurrences.middle ws a_i ~i ~k (reference k) in
        ft := filter.entry !ft ne;
        ws.forward.(i).(k) <- ne
      done
    done

  let final ws recurrences =
    let columns = last_array_index ws.final in
    for k = 0 to columns do
      ws.final.(k) <- recurrences.end_ ws k
    done

  (* Fill in both parts of the workspace. *)
  let both_f ?rows ~filter ws recurrences ~read ~reference =
    pass_f ?rows ~filter ws recurrences ~read ~reference;
    final ws recurrences

  (* After filling in both parts of the workspace,
     compute the final emission value. *)
  let full_f ?rows ~filter ws recurrences ~read ~reference =
    both_f ?rows ~filter ws recurrences ~read ~reference;
    ws.emission <- recurrences.final_e ws

  (* Without the filter. *)
  let pass = pass_f ~filter:empty_filter

  let full = full_f ~filter:empty_filter

end (* ForwardPass *)

type 'a pass_result =
  | Completed of 'a
  | Filtered of string

let pass_result_to_string c_to_s = function
  | Completed c -> sprintf "Completed: %s" (c_to_s c)
  | Filtered m  -> sprintf "Filtered: %s" m

module ForwardSingleGen (R: Ring) = struct

  module Workspace = SingleWorkspace(R)
  module Fc = ForwardCalcs(R)

  type base = BaseState.t

  let recurrences ?insert_p tm read_length =
    let open Workspace in
    let r = Fc.g ?insert_p tm read_length in
    let start obsp base = r.start (Fc.to_match_prob obsp base) in
    let fst_col ws obsp base ~i =
      r.fst_col (Fc.to_match_prob obsp base) (ws.forward.(i-1).(0))
    in
    let middle ws obsp base ~i ~k =
      let emp = Fc.to_match_prob obsp base in
      let ks = k-1 in
      let insert_c = ws.forward.(i-1).(k) in
      let match_c = ws.forward.(i-1).(ks) in
      let delete_c = ws.forward.(i).(ks) in
      r.middle emp ~insert_c ~delete_c ~match_c
    in
    let end_ ws k = r.end_ ws.forward.(read_length-1).(k) in
    let final_e ws =
      Array.fold_left ~init:R.zero ~f:(R.(+)) ws.final
      (* Array.iter ws.final ~f:(fun f -> ws.emission <- R.(ws.emission + f)) *)
    in
    { start ; fst_col ; middle ; end_ ; final_e }

  let full ?insert_p ?transition_ref_length ?max_number_mismatches ~read_length
    ws allele_a =
    let tm =
      let ref_length = Option.value transition_ref_length
        ~default:(Array.length allele_a)
      in
      Phmm.TransitionMatrix.init ~ref_length read_length in
    let recurrences = recurrences ?insert_p tm read_length in
    let reference k = allele_a.(k) in
    let open Workspace in
    match max_number_mismatches with
    | None ->
        begin fun ~read ->
          ForwardPass.full ws recurrences ~reference ~read;
          Completed ()
        end
    | Some number_mismatches ->
        begin fun ~read ->
          let of_entry c = c.match_ in
          let filter = Fc.match_filter ~number_mismatches of_entry in
          try
            ForwardPass.full_f ~filter ws recurrences ~reference ~read;
            Completed ()
          with PastThreshold msg ->
            Filtered msg
        end

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

module ForwardMultipleGen (R : Ring)(Aset: Alleles.Set) = struct

  module Cm = CAM.Make(Aset)

  (* Eh... not the best place for it. *)
  let cam_max =
    Cm.fold ~init:(R.constant neg_infinity) ~f:(fun m _s v -> R.max m v)

  let cam_max_cell =
    Cm.fold ~init:(R.constant neg_infinity) ~f:(fun m _s c -> R.max m c.match_)

  module Workspace = MakeMultipleWorkspace(R)(Cm)

  type base_emissions = (BaseState.t * int) CAM.t

  (* offset and emission probabilties *)
  type 'a oeps = (int * 'a) CAM.t
  type 'a emission_map = 'a oeps PosMap.t
  type 'a banded_recs =
    (* This isn't the greatest design, but we need to separate (to cache) the
       emission calculation that is then passed back into the banded middle
       call. *)
    { middle_emissions : obs
                          -> base_emissions
                          -> 'a oeps
    ; banded           : Workspace.t
                          -> 'a oeps
                          -> ?prev_row:('a cell)              (* Needed for missing/edge cases.*)
                          -> ?cur_row:(Workspace.entry)       (* Needed for missing/edge cases.*)
                          -> i:int
                          -> k:int
                          -> Workspace.entry
    ; spec_final_e     : int list list -> Workspace.t -> R.t array
    }

  let per_allele_emission_arr len =
    Array.make len R.one

  module Fc = ForwardCalcs(R)

  let recurrences ?insert_p tm read_length number_alleles =
    let open Workspace in
    let r = Fc.g ?insert_p tm read_length in

   (* TODO: I could imagine some scenario's where it makes sense to cache,
       precompute or memoize this calculation. The # of base errors isn't
       that large (<100) and there are only 4 bases. So we could be performing
       the same lookup. *)
    let to_em_set obsp emissions =
      Cm.map emissions ~f:(fun (b, offset) ->
        offset, Fc.to_match_prob obsp b)
    in
    let start obsp emissions =
      to_em_set obsp emissions
      |> Cm.map ~bijective:true
        ~f:(fun (_offset, emissionp) -> r.start emissionp)
    in
    let fst_col ws obsp emissions ~i =
      to_em_set obsp emissions
      |> Cm.map2 (ws.forward.(i-1).(0))
          ~f:(fun insert_c (_offset, emission_p) ->
                r.fst_col emission_p insert_c)
    in
    let middle ws obsp emissions ~i ~k =
      let inserts = ws.forward.(i-1).(k) in
      let ems = to_em_set obsp emissions in
      Cm.concat_map2 inserts ~by:ems   (* ORDER matters for performance! *)
          ~f:(fun inters insert_c (offset, emission_p) ->
                let ks = Pervasives.(+) k offset in
                let matches = ws.forward.(i-1).(ks) in
                let deletes = ws.forward.(i).(ks) in
                let insertsi = Cm.singleton inters insert_c in
                (* inserti should come before other 2 for performance. *)
                Cm.map3 insertsi deletes matches
                  ~f:(fun insert_c delete_c match_c ->
                        r.middle emission_p ~insert_c ~delete_c ~match_c))
    in
    let end_ ws k =
      Cm.map ~bijective:true ws.forward.(read_length-1).(k) ~f:r.end_
    in
    let update_emission_from_cam em l =
      let open R in
      Cm.iter_values l ~f:(fun i v -> em.(i) <- em.(i) + v)
    in
    let final_e ws =
      let em = Array.make number_alleles R.zero in
      Array.iter ws.final ~f:(update_emission_from_cam em);
      em
      (*Array.iter ws.final ~f:(update_emission_from_cam ws) *)
    in
    let banded ws allele_ems ?prev_row ?cur_row ~i ~k =
      let with_insert inters (offset, emission_p) insert_c =
        let ks = Pervasives.(+) k offset in
        let calc insert_c delete_c match_c =
          r.middle emission_p ~insert_c ~delete_c ~match_c
        in
        let matches = ws.forward.(i-1).(ks) in
        let deletes = ws.forward.(i).(ks) in
        let insertsi = Cm.singleton inters insert_c in
        Cm.map3_partial insertsi
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
      let inserts = ws.forward.(i-1).(k) in
        Cm.concat_map2_partial allele_ems ~by:inserts
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
          update_emission_from_cam em ws.final.(k)));
      em
    in
    { start; fst_col; middle; end_; final_e}
    , { middle_emissions = to_em_set ; banded ; spec_final_e}

  module Regular = ForwardPass

  (* Banded Pass logic **)
  module Bands = struct

    (* 1. Identify bands. *)
    let select_specific_band_indices ws c =
      let open Workspace in
      let lg = largest c.number in        (* compares by match_, insert, delete *)
      let ev = Cm.init_everything [] in
      Array.fold_left ws.forward.(c.warmup) ~init:(0, ev) ~f:(fun (k, acc) by ->
        let nacc =
          Cm.map2_partial acc ~by ~f:(fun lst c -> lg c k lst)
            ~missing:(fun s v -> Cm.singleton s v)
        in
        k + 1, nacc)
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
      let open Workspace in
      Cm.concat_map bands ~f:(fun s (bv, cols) ->
        let end_col = List.last cols |> Option.value_exn ~msg:"empty cols!" in
        Cm.get_exn s ws.forward.(row).(end_col)
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
      let open Workspace in
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
        ws.forward.(i).(k) <- Cm.join entry ws.forward.(i).(k);
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
      let open Workspace in
      List.iter b.cols ~f:(fun k ->
        ws.final.(k) <- Cm.join (rrecs.end_ ws k) ws.final.(k))

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
            List.fold_left bands ~init:([], PosMap.empty, Cm.empty)
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
      let open Workspace in
      let rrecs, brecs = recurrences in
      let banded_end bands =
        List.iter bands ~f:(fill_end rrecs ws);
        let spec_cols = List.map bands ~f:(fun b -> b.cols) in
        ws.emission <- brecs.spec_final_e spec_cols ws
      in
      banded_end (banded_middle (c.warmup + 1))

  end (* Bands *)

end (* ForwardMultipleGen *)

(* Instantiate the actual passes over the 2 implemented probability rings. *)
module ForwardS = ForwardSingleGen(MultiplicativeProbability)
module ForwardSLogSpace = ForwardSingleGen(LogProbabilities)

module ForwardM = ForwardMultipleGen(MultiplicativeProbability)
module ForwardMLogSpace = ForwardMultipleGen (LogProbabilities)

type t =
  { align_date      : string
  ; number_alleles  : int
  ; aset            : (module Alleles.Set)
  ; alleles         : string list
  ; merge_map       : (string * string) list
  ; emissions_a     : (BaseState.t * int) CAM.t array
  ; increment_a     : (Alleles.set * int) list array
  }

let construct input selectors =
  if not (Alleles.Input.imputed input) then
    invalid_argf "Allele input MUST be imputed!"
  else begin
    let open Mas_parser in
    Alleles.Input.construct input >>= fun (mp, merge_map) ->
      let nalt_elems =
        mp.alt_elems
        |> List.sort ~cmp:(fun (n1, _) (n2, _) -> Alleles.compare n1 n2)
        |> Alleles.Selection.apply_to_assoc selectors
      in
      let alleles = mp.reference :: List.map ~f:fst nalt_elems in
      let allele_index = Alleles.index alleles in
      let module Aset = Alleles.MakeSet (struct let index = allele_index end) in
      let aset = (module Aset : Alleles.Set) in
      let module Cm = CAM.Make(Aset) in
      let emissions_a, position_map =
        initialize_base_array_and_position_map aset mp.reference mp.ref_elems
      in
      let aaa = add_alternate_allele aset mp.reference ~position_map in
      List.iter ~f:(fun (allele, altseq) -> aaa allele altseq emissions_a) nalt_elems;
      let emissions_a =
        (* TODO: Move the CAM logic up into the construction algorithms *)
        Array.map emissions_a ~f:(fun l ->
          List.map l ~f:(fun (b, s) -> (s, b)) |> Cm.of_list)
      in
      let increment_a = Array.make (Array.length emissions_a - 1) Cm.empty in
      Array.iteri emissions_a ~f:(fun i s ->
        if i = 0 then () else
          Cm.map s ~f:(fun (_b, v) -> v)
          |> Cm.iter ~f:(fun s g ->
              let k = i + g in
              increment_a.(k) <- Cm.add s i increment_a.(k)));
      Ok { align_date = mp.align_date
         ; number_alleles = Alleles.length allele_index
         ; aset
         ; alleles
         ; merge_map
         ; emissions_a
         ; increment_a
         }
  end

let save_pphmm t =
  let fname = Filename.temp_file ~temp_dir:"." "pphmm" "" in
  let oc = open_out fname in
  Marshal.to_channel oc t [Marshal.Closures];
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

type proc =
  { output_ws_array : unit -> float array                     (* Allocate the right size *)
  ; doit            : bool -> string -> float array -> unit pass_result
  ; best_alleles    : unit -> (float * string) list
  ; best_positions  : unit -> (float * int) list
  ; per_allele_llhd : unit -> float array                     (* Pass'es output. *)
  }


(* TODO: expose this 5 if it becomes useful *)
let lg5 a i lst =
  largest 5 a i lst

(* Single, for one allele, forward pass *)
let setup_single_allele_forward_pass ?insert_p ?max_number_mismatches
  t read_length allele =
  let { emissions_a; aset; alleles; _ } = t in
  let module Aset = (val aset : Alleles.Set) in
  (* Recover allele's base sequence, aka, emission. *)
  if not (List.exists ~f:((=) allele) alleles) then
    invalid_argf "%s not found among the list of alleles!" allele
  else
    let module Cm = CAM.Make(Aset) in
    let allele_set = Aset.singleton allele in
    let allele_a =
      Array.fold_left emissions_a ~init:[] ~f:(fun acc cm ->
        match Cm.get_value allele_set cm with
        | None         -> acc           (* In gap. *)
        | Some (bs, _) -> bs :: acc)
      |> List.rev
      |> Array.of_list
    in
    let best_positions final =
      Array.fold_left final ~init:(0, [])
        ~f:(fun (p, acc) l ->
          (p + 1, lg5 l p acc))
      |> snd
    in
    let ref_length = Array.length allele_a in
    let ws = ForwardSLogSpace.Workspace.generate ~ref_length ~read_length in
    (* For comparison against all of the alleles we want to have the same
       transition probabilities, that depends upon the reference size.
       Therefore we'll use the reference size of all the alleles; size of
       emission_a.  Perhaps we want to expose a parameter to switch to the
       'local' reference size; size of allele_a. *)
    let transition_ref_length = Array.length emissions_a in
    let pass =
      ForwardSLogSpace.full ?insert_p ?max_number_mismatches
        ~transition_ref_length ~read_length ws allele_a
    in
    let doit rc rd rd_errors =
      let read = access rc rd rd_errors in
      pass read
    in
    let best_alleles () = [ws.emission, allele] in
    let best_positions () = best_positions ws.final in
    let per_allele_llhd () = [| ws.emission |] in
    let output_ws_array () = [| LogProbabilities.one |] in
    { doit ; best_alleles ; best_positions ; per_allele_llhd
    ; output_ws_array
    }

(* Return
  1. a function to process one read
  2. the workspace that it uses
  3. a function to combine results *)
let setup_single_pass ?band ?insert_p ?max_number_mismatches read_length t =
  let { number_alleles; emissions_a; increment_a; aset; alleles; _ } = t in
  let ref_length = Array.length emissions_a in
  let tm = Phmm.TransitionMatrix.init ~ref_length read_length in
  let module AS = (val aset : Alleles.Set) in
  let module F = ForwardMLogSpace(AS) in
  let r, br = F.recurrences ?insert_p tm read_length number_alleles in
  let ws = F.Workspace.generate number_alleles ref_length read_length in
  let last_read_index = read_length - 1 in
  let best_alleles emissions =
    Array.to_list emissions
    |> List.fold_left2 alleles
        ~init:[] ~f:(fun acc allele emission -> lg5 emission allele acc)
  in
  let best_positions final =
    Array.fold_left final ~init:(0, [])
      ~f:(fun (p, acc) fcam ->
        (p + 1, lg5 (F.cam_max fcam) p acc))
    |> snd
  in
  let reference i = emissions_a.(i) in
  let best_alleles () = best_alleles ws.emission in
  let best_positions () = best_positions ws.final in
  let per_allele_llhd () = ws.emission in
  let output_ws_array () = F.per_allele_emission_arr t.number_alleles in
  let normal () =
    let doit rc rd rd_errors =
      (* We do not have to clear the workspace, since a full pass will
         overwrite all elements of the workspace.
        F.Workspace.clear ws;*)
      let read = access rc rd rd_errors in
      match max_number_mismatches with
      | None ->
          F.Regular.full ws r ~reference ~read;
          Completed ()
      | Some number_mismatches ->
          let of_entry = F.cam_max_cell in
          let filter = F.Fc.match_filter ~number_mismatches of_entry in
          try
            F.Regular.full_f ~filter ws r ~reference ~read;
            Completed ()
          with PastThreshold msg ->
            Filtered msg
    in
    { doit ; best_alleles ; best_positions ; per_allele_llhd
    ; output_ws_array
    }
  in
  let banded c =
    let doit rc rd rd_errors =
      (* Clear the forward/final array since the banded pass algorithm relies
         on unfilled cells to indicate boundaries (places where we use
         heuristics). *)
      F.Workspace.clear ws;
      let read = access rc rd rd_errors in
      F.Regular.pass ws r ~reference ~read ~rows:c.warmup;
      F.Bands.pass c emissions_a increment_a ws (r, br) last_read_index read;
      Completed ()
    in
    { doit ; best_alleles ; best_positions ; per_allele_llhd
    ; output_ws_array
    }
  in
  match band with
  | None                                    -> normal ()
  | Some c when c.warmup >= last_read_index -> normal ()
  | Some c                                  -> banded c

(**** Aggregations of the passes. *)

type stats =
  { per_allele : (float * string) list
  ; positions  : (float * int) list
  }

let stats_to_strings ?(sep='\t') s =
  let l_to_s fmt l =
    String.concat ~sep:";" (List.map l ~f:(fun (l,a) -> sprintf fmt  a l))
  in
  let al_to_s l = l_to_s "%s:%0.2f" l in
  let pl_to_s l = l_to_s "%d:%0.2f" l in
  sprintf "%s%c%s" (al_to_s s.per_allele) sep (pl_to_s s.positions)

type 'a both_rc_states =
  { regular     : 'a pass_result
  ; complement  : 'a pass_result
  }

type mapped_stats = stats both_rc_states

(* The output is ordered so that the left most stat is the most likely.
   We add characters 'R', 'C', 'F' as little annotations. *)
let mapped_stats_to_string ?(sep='\t') ms =
  match ms.regular, ms.complement with
  | Completed r, Completed c ->
    if fst (List.hd_exn r.positions) > fst (List.hd_exn c.positions) then
      sprintf "R%c%s%cC%c%s"
        sep (stats_to_strings ~sep r)
        sep sep (stats_to_strings ~sep c)
    else
      sprintf "C%c%s%cR%c%s"
        sep (stats_to_strings ~sep c)
        sep sep (stats_to_strings ~sep r)
  | Completed r, Filtered m ->
      sprintf "R%c%s%cF%c%s" sep (stats_to_strings ~sep r) sep sep m
  | Filtered m, Completed c ->
      sprintf "C%c%s%cF%c%s" sep (stats_to_strings ~sep c) sep sep m
  | Filtered mr, Filtered mc ->
      sprintf "F%c%s%cF%c%s" sep mr sep sep mc

let pass_result_to_stat = function
  | Completed s -> fst (List.hd_exn s.positions)
  | Filtered _  -> neg_infinity

let best_stat ms =
  max (pass_result_to_stat ms.regular) (pass_result_to_stat ms.complement)

let both pass_res_to_stat proc read read_prob =
  let w rc =
    match proc.doit rc read read_prob with
    | Filtered m   -> Filtered m
    | Completed () -> Completed (pass_res_to_stat rc proc)
  in
  let regular = w false in
  let complement = w true in
  { regular; complement }

let mapper =
  both (fun _rc proc ->
          { per_allele = proc.best_alleles ()
          ; positions = proc.best_positions ()
          })

let compare_emissions e1 e2 =
  let r1 = Array.fold_left e1 ~init:neg_infinity ~f:max in
  let r2 = Array.fold_left e2 ~init:neg_infinity ~f:max in
  r1 >= r2

type reducer_stat =
  { comp : bool
  ; likelihood : float array
  }

let most_likely_between_rc = function
  | Filtered mr, Filtered mc -> Filtered (sprintf "Both! %s, %s" mr mc)
  | Filtered _r, Completed c -> Completed c
  | Completed r, Filtered _c -> Completed r
  | Completed r, Completed c ->
    if compare_emissions r.likelihood c.likelihood then
      Completed r
    else
      Completed c

let reducer proc ?(check_rc=true) read read_prob =
  if check_rc then begin
    let mp =
      both (fun comp proc ->
        { comp; likelihood = proc.per_allele_llhd () })
        proc read read_prob
    in
    most_likely_between_rc (mp.regular, mp.complement)
  end else begin
    match proc.doit false read read_prob with
    | Filtered m    ->
        Filtered m
    | Completed ()  ->
        Completed { comp = false ; likelihood = proc.per_allele_llhd () }
  end

let single_allele_forward_pass ?insert_p ?max_number_mismatches mode pt
  read_length allele =
  let proc =
    setup_single_allele_forward_pass ?insert_p ?max_number_mismatches pt
    read_length allele
  in
  match mode with
  | `Mapper   -> `Mapper (mapper proc)
  | `Reducer  -> `Reducer (proc.output_ws_array (), reducer proc)

let forward_pass ?insert_p ?max_number_mismatches ?band mode t read_length =
  let proc =
    setup_single_pass ?insert_p ?max_number_mismatches ?band read_length t in
  match mode with
  | `Mapper   -> `Mapper (mapper proc)
  | `Reducer  -> `Reducer (proc.output_ws_array (), reducer proc)

(* Combine results from multiple loci. *)
module MultiLoci = struct

  type best_state =
    { name  : string
    ; rspr  : reducer_stat pass_result
    ; maxl  : float
    }

  let max_arr = Array.fold_left ~init:neg_infinity ~f:max

  let pr_to_llhd = function
    | Filtered _  -> neg_infinity
    | Completed c -> max_arr c.likelihood

  let to_bs (name, rspr) =
    { name; rspr; maxl = pr_to_llhd rspr}

  let best_bs bs1 bs2 =
    (*printf "%s %f vs %s %f\n%!" bs1.name bs1.maxl bs2.name bs2.maxl; *)
    if bs1.maxl >= bs2.maxl then bs1 else bs2

  let best = function
    | []          -> invalid_argf "Can't select best from empty!"
    | h :: t ->
        let init = to_bs h in
        List.fold_left t ~init ~f:(fun bs p ->
          best_bs bs (to_bs p))

  let add_log_likelihoods ~into = function
    | Filtered _  -> ()                     (* Should we warn about ignoring a read? *)
    | Completed r ->
        let n = Array.length into in
        for i = 0 to n - 1 do
          into.(i) <- into.(i) +. r.likelihood.(i)
        done

  let merge current state =
    let b = best current in
    List.iter state ~f:(fun (n, into) ->
      if n = b.name then add_log_likelihoods ~into b.rspr);
    state

end (* MultiLoci *)

let forward_pass_m ?insert_p ?max_number_mismatches ?band mode tlst read_length =
  let make_passes = setup_single_pass ?insert_p ?max_number_mismatches ?band read_length in
  let passes = List.map_snd tlst ~f:make_passes in
  match mode with
  | `Mapper ->
      `Mapper (fun read read_prob ->
        List.map_snd passes ~f:(fun p -> mapper p read read_prob))
  | `Reporter ->
      `Reporter ( List.map_snd passes ~f:(fun proc -> proc.output_ws_array ())
                , fun read read_prob state ->
                    let current =
                      List.map_snd passes ~f:(fun p -> reducer p read read_prob)
                    in
                    MultiLoci.merge current state)
