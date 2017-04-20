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
    | []          -> invalid_arg "Rlele.encode: empty list"
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

  let expand_into_array ~f ~update ret rl =
    let rec fill_value i length v =
      for j = i to i + length - 1 do ret.(j) <- update ret.(i) v done
    in
    let rec loop i = function
      | []                    -> ()
      | { value; length} :: t -> fill_value i length (f value);
                                 loop (i + length) t
    in
    fill_value 0 rl.hd.length (f rl.hd.value);
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
let initialize_base_array_and_position_map ali reference ref_elems =
  let open Mas_parser in
  let ref_set () = Alleles.Set.singleton ali reference in
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
let add_alternate_allele ali reference ~position_map allele allele_instr arr =
  let base_and_offset b o ((bp,bo), _) = b = bp && o = bo in
  let add_to_base_state i b o =
    (*printf "adding base at %d %c %d\n" i (BaseState.to_char b) o; *)
    match List.find arr.(i) ~f:(base_and_offset b o) with
    | None              -> let s = Alleles.Set.singleton ali allele in
                           (*printf "single cell at %d for %s \n"  i allele; *)
                           arr.(i) <- ((b, o), s) :: arr.(i)
    | Some ((ab,ro), s) -> (*printf "At: %d %s to %s\n"  i allele (Alleles.Set.to_string ali s); *)
                           ignore (Alleles.Set.set ali s allele)
  in
  let is_reference_set (_, s) = Alleles.Set.is_set ali s reference in
  let add_to_reference_set offset start end_ =
    (*printf "adding reference set at %d %d %d\n" offset start end_; *)
    let rec loop i offset =
      if i = end_ then offset else begin
        match List.find arr.(i) ~f:is_reference_set with
        | None              -> (* Reference in gap -> so are we! *)
                               loop (i + 1) (-1)  (* Offset becomes -1 after 1st use! *)
        | Some ((rb,ro), s) ->
            if ro = offset then begin
              ignore (Alleles.Set.set ali s allele);
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

(* Probability Ring where we perform the forward pass calculation. *)
module type Ring = sig

  type t
  val to_string : t -> string
  val zero : t
  val one  : t

  val ( + ) : t -> t -> t
  val ( * ) : t -> t -> t

  (* Special constructs necessary for the probabilistic logic. *)
  (* Convert constant probabilities. *)
  val constant : float -> t

  (* Scale a probability be a third. *)
  val times_one_third : float -> t

  (* Complement probability. *)
  val complement_probability : float -> t

end (* Ring *)

type emissions = ((BaseState.t * int) * Alleles.Set.t) list
(* For every k there are 3 possible states. *)

type 'a cell =
  { match_  : 'a
  ; insert  : 'a
  ; delete  : 'a
  }

type 'a entry = (Alleles.Set.t * 'a cell) list
type 'a final_entry = (Alleles.Set.t * 'a) list

type workspace =
  { mutable forward             : float entry array array
  ; mutable final               : float final_entry array
  ; mutable per_allele_emission : float array
  }

let generate_workspace number_alleles bigK read_size =
  let just_zeros n = Array.make n 0. in
  { forward             = Array.init bigK ~f:(fun _ -> Array.make read_size [])
  ; final               = Array.make bigK []
  ; per_allele_emission = just_zeros number_alleles
  }

type 'a fwd_recurrences =
  { start     :  char -> float -> emissions -> 'a entry
  ; first_row : 'a entry array array -> Alleles.index -> char -> float ->
                  emissions -> i:int -> 'a entry
  ; middle    : 'a entry array array -> Alleles.index -> char -> float ->
                  emissions -> i:int -> k:int -> 'a entry
  (* Doesn't use the delete section. *)
  ; end_      : 'a entry array array -> int -> 'a final_entry

  ; emission  : 'a final_entry array -> Alleles.index -> 'a array

(* Combine emission results *)
  ; combine   : into:'a array -> 'a array -> unit
  }

(* Union, tail recursive. *)
let mutate_or_add value allele_set lst =
  let rec loop acc = function
    | (s, v) :: t when v = value -> acc @ (Alleles.Set.union s allele_set, v) :: t
    | h :: t                     -> loop (h :: acc) t
    | []                         -> (allele_set, value) :: acc
  in
  loop [] lst

(* let mutate_or_add value new_allele_set assoc =
  let added =
    List.fold assoc ~init:false ~f:(fun added (into, v) ->
      if added then
        added
      else if v = value then begin
        Alleles.Set.unite ~into new_allele_set;
        true
      end else
        false)
  in
  if added then
    assoc
  else
    (Alleles.Set.copy new_allele_set, value) :: assoc *)

let debug_set_assoc = ref false

let set_assoc ali to_find slst =
  if !debug_set_assoc then printf "set_assoc %d\n" (List.length slst);
  let to_s = Alleles.Set.to_string ali in
  let rec loop to_find acc = function
    | []          -> invalid_argf "Still missing! %s after looking in: %s"
                      (to_s to_find) (List.map slst ~f:(fun (s,_) -> to_s s) |> String.concat ~sep:"\n\t")
    | (s, v) :: t ->
        let inter = Alleles.Set.inter s to_find in
        let still_to_find = Alleles.Set.diff to_find s in
        if !debug_set_assoc then
          printf "For\t%s\n\
                  in\t%s\n\
                  found\t%s\n\
                  butnot\t%s\n"
            (to_s to_find)
            (to_s s)
            (to_s inter)
            (to_s still_to_find);
        if inter = to_find then                         (* Found everything *)
          (to_find, v) :: acc
        else if Alleles.Set.cardinal inter = 0 then       (* Found nothing. *)
          loop to_find acc t
        else                                            (* Found something. *)
          loop still_to_find ((inter, v) :: acc) t
  in
  loop to_find [] slst

module ForwardGen (R : Ring) = struct

  (* TODO. Avoid the `float_of_int (Phred_score.to_int c) /. -10.` round trip
      between converting to log10p and then back to log10, and just use char's
      instead for the quality calc. *)
  let to_match_prob base base_error =
    let compare_against c =
      if base = c then
        R.complement_probability base_error
      else
        R.times_one_third base_error
    in
    let open BaseState in
    function
    | A -> compare_against 'A'
    | C -> compare_against 'C'
    | G -> compare_against 'G'
    | T -> compare_against 'T'

  let per_allele_emission_arr len =
    Array.make len R.one

  let recurrences tm ~insert_prob read_size =

    let open R in                       (* Opening R shadows '+' and '*' below*)
    let t_s_m = constant (tm `StartOrEnd `Match) in
    let t_s_i = constant (tm `StartOrEnd `Insert) in
    let t_m_m = constant (tm `Match `Match) in
    let t_i_m = constant (tm `Insert `Match) in
    let t_d_m = constant (tm `Delete `Match) in

    let t_m_i = constant (tm `Match `Insert) in
    let t_i_i = constant (tm `Insert `Insert) in

    let t_m_d = constant (tm `Match `Delete) in
    let t_d_d = constant (tm `Delete `Delete) in

    let t_m_s = constant (tm `Match `StartOrEnd) in
    let t_i_s = constant (tm `Insert `StartOrEnd) in

    let start_i = t_s_i * insert_prob in
    { start   = begin fun base base_error emissions ->
                  List.fold_left emissions ~init:[]
                    ~f:(fun acc ((b, offset), allele_set) ->
                          let emp = to_match_prob base base_error b in
                          mutate_or_add (offset, emp)
                          allele_set acc)
                  |> List.map ~f:(fun (allele_set, (_offset, emissionp)) ->
                        allele_set,
                        { match_ = emissionp * t_s_m
                        ; insert = start_i
                        ; delete = zero
                        })
                end
    ; first_row = begin fun fm ali base base_error emissions ~i ->
                    let set_assoc = set_assoc ali in
                    List.fold_left emissions ~init:[]
                      ~f:(fun acc ((b, offset), allele_set) ->
                            let emp = to_match_prob base base_error b in
                            mutate_or_add (offset, emp) allele_set acc)
                    |> List.fold_left ~init:[] ~f:(fun acc (allele_set, (offset, emp)) ->
                        let inserts = set_assoc allele_set fm.(0).(i-1) in
                        (*printf "at k: %d i: %d offset:%d %d \n%!" 0 i offset j; *)
                          List.fold_left inserts ~init:acc ~f:(fun acc (als, c) ->
                            mutate_or_add
                              { match_ = emp * ( t_m_m * zero
                                              + t_i_m * zero
                                              + t_d_m * zero)
                              ; insert = insert_prob * (t_m_i * c.match_ + t_i_i * c.insert)
                              ; delete = t_m_d * zero + t_d_d * zero
                              } als acc))
                  end
    ; middle  = begin fun fm ali base base_error emissions ~i ~k ->
                  let set_assoc = set_assoc ali in
                  List.fold_left emissions ~init:[]
                    ~f:(fun acc ((b, offset), allele_set) ->
                          let emp = to_match_prob base base_error b in
                          mutate_or_add (offset, emp) allele_set acc)
                  |> List.fold_left ~init:[] ~f:(fun acc (allele_set, (offset, emp)) ->
                      let inserts = set_assoc allele_set fm.(k).(i-1) in
                      (*printf "at k: %d i: %d offset:%d %d \n%!" k i offset j; *)
                      let ks = Pervasives.(+) k offset in
                      let matches = fm.(ks).(i-1) in
                      let deletes = fm.(ks).(i) in
                      List.fold_left inserts ~init:acc
                        ~f:(fun init (insert_s, insert_c) ->
                            List.fold_left (set_assoc insert_s matches) ~init
                              ~f:(fun init (match_s, match_c) ->
                                  List.fold_left (set_assoc match_s deletes) ~init
                                    ~f:(fun acc (delete_s, delete_c) -> (* at delete_s intersects all 3!*)
                                          mutate_or_add
                                            { match_ = emp * ( t_m_m * match_c.match_
                                                             + t_i_m * match_c.insert
                                                             + t_d_m * match_c.delete)
                                            ; insert = insert_prob * ( t_m_i * insert_c.match_
                                                                     + t_i_i * insert_c.insert)
                                            ; delete = t_m_d * delete_c.match_
                                                     + t_d_d * delete_c.delete
                                            } delete_s acc))))
                end
    ; end_    = begin fun fm k ->
                  let fc = fm.(k).(read_size-2) in
                  List.map fc ~f:(fun (allele_set, c) ->
                    allele_set, (c.match_ * t_m_s + c.insert * t_i_s))
                end
    ; emission  = begin fun final ali ->
                    let ret = Alleles.Map.make ali zero in
                    Array.iter final ~f:(fun l ->
                      List.iter l ~f:(fun (allele_set, v) ->
                        Alleles.Map.update_from allele_set ~f:((+) v) ret));
                    Alleles.Map.to_array ret
                  end
    ; combine   = begin fun ~into em ->
                    Array.iteri em ~f:(fun i e -> into.(i) <- into.(i) * e)
                  end
    }

end (* ForwardGen *)

module MutliplicativeProbability = struct
  type t = float
  let zero  = 0.
  let one   = 1.
  let ( + ) = ( +. )
  let ( * ) = ( *. )

  let constant x = x

  let complement_probability p =
    1. -. p

  let times_one_third p =
    p /. 3.

  let to_string = sprintf "%f"

end (* MutliplicativeProbability *)

module Forward = ForwardGen(MutliplicativeProbability)

module LogProbabilities = struct

  let zero  = neg_infinity

  let one   = 0.  (* log10 1. *)

  let exp10 x = 10. ** x

  let ( * ) lx ly = lx +. ly

  let ( + ) lx ly =
         if lx = neg_infinity then ly
    else if ly = neg_infinity then lx
    else if lx > ly           then lx +. log10 (1. +. exp10 (ly -. lx))
    else (* lx < ly *)             ly +. log10 (1. +. exp10 (lx -. ly))

  type t = float

  let to_string = sprintf "%f"
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

module ForwardLogSpace = ForwardGen (LogProbabilities)

type t =
  { align_date    : string
  ; allele_index  : Alleles.index
  ; merge_map     : (string * string) list
  ; emissions_a   : emissions array
  }

let construct input selectors =
  if not (Alleles.Input.imputed input) then
    invalid_argf "Allele input MUST be imputed!"
  else
    let open Mas_parser in
    Alleles.Input.construct input >>= fun (mp, merge_map) ->
      let nalt_elems =
        mp.alt_elems
        |> List.sort ~cmp:(fun (n1, _) (n2, _) -> Alleles.compare n1 n2)
        |> Alleles.Selection.apply_to_assoc selectors
      in
      let alleles = mp.reference :: List.map ~f:fst nalt_elems in
      let allele_index = Alleles.index alleles in
      let ems, position_map = initialize_base_array_and_position_map
        allele_index mp.reference mp.ref_elems
      in
      let aaa = add_alternate_allele allele_index mp.reference ~position_map in
      List.iter ~f:(fun (allele, altseq) -> aaa allele altseq ems) nalt_elems;
      Ok { align_date = mp.align_date
         ; allele_index
         ; merge_map
         ; emissions_a = ems
         }

let debug_ref = ref false

let save_workspace ws =
  let fname = Filename.temp_file ~temp_dir:"." "forward_workspace" "" in
  let oc = open_out fname in
  Marshal.to_channel oc ws [];
  close_out oc;
  printf "Saved workspace to %s\n" fname

let load_workspace fname =
  let ic = open_in fname in
  let ws : workspace = Marshal.from_channel ic in
  close_in ic;
  ws

let generate_workspace_conf c read_size =
  let bigK = Array.length c.emissions_a in
  let number_alleles = Alleles.length c.allele_index in
  generate_workspace number_alleles bigK read_size

let forward_pass ?(logspace=true) ?ws { allele_index; emissions_a } read_size =
  let bigK = Array.length emissions_a in
  let number_alleles = Alleles.length allele_index in
  let tm = Phmm.TransitionMatrix.init ~ref_length:bigK read_size in
  let insert_prob = 0.25 in
  let recurrences =
    (if logspace then ForwardLogSpace.recurrences else Forward.recurrences)
      tm ~insert_prob read_size
  in
  let ws =
    match ws with
    | Some w -> w
    | None   -> generate_workspace number_alleles bigK read_size
  in
  fun ~into read read_prob ->
    (* special case the first row. *)
    ws.forward.(0).(0) <-
      recurrences.start (String.get_exn read 0) read_prob.(0) emissions_a.(0);
    for i = 1 to read_size - 1 do
      let base = String.get_exn read i in
      let base_prob = read_prob.(i) in
      ws.forward.(0).(i) <-
        recurrences.first_row ws.forward allele_index base base_prob ~i emissions_a.(0)
    done;
    (* All other rows. *)
    for k = 1 to bigK - 1 do
      let ek = emissions_a.(k) in
      ws.forward.(k).(0) <-
        recurrences.start (String.get_exn read 0) read_prob.(0) ek;
      for i = 1 to read_size - 1 do
        let base = String.get_exn read i in
        let base_prob = read_prob.(i) in
        ws.forward.(k).(i) <-
          recurrences.middle ws.forward allele_index base base_prob ~i ~k ek
      done
    done;
    for k = 0 to bigK - 1 do
      ws.final.(k) <- recurrences.end_ ws.forward k
    done;
    ws.per_allele_emission <- recurrences.emission ws.final allele_index;
    recurrences.combine ~into ws.per_allele_emission;
    if !debug_ref then save_workspace ws;
    into

let setup ?ws ~logspace t read_size =
  let perform_forward_pass = forward_pass ?ws ~logspace t read_size in
  let output_array =
    let number_alleles = Alleles.length t.allele_index in
    if logspace then
      ForwardLogSpace.per_allele_emission_arr number_alleles
    else
      Forward.per_allele_emission_arr number_alleles
  in
  perform_forward_pass, output_array
