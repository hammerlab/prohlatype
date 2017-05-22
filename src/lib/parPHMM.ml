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
    | []          -> invalid_arg "Rlel.encode: empty list"
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

let list_map_snd lst ~f =
  List.map lst ~f:(fun (k, v) -> k, f v)

let list_map2_snd l1 l2 ~f =
  List.map2 l1 l2 ~f:(fun (k1, v1) (k2, v2) ->
    assert (k1 = k2);
    (k1, f v1 v2))

type set = Alleles.set

(* CAM= Compressed Allele Map

   Since our PHMM is parameterized by alleles, we have to keep track many
   values on a per allele basis. This module aims to provide an abstraction
   for this map: allele -> 'a.

   It is different than the Alleles.Map module (and perhaps that one should be
   replaced or deprecated) because the implementation tries to be succinct, by
   compressing the values and avoiding O(n) (n = # of alleles) maps/folds.
   Specifically, the operations are performed on the unique values in the map.
*)
module type CAMs = sig

  type 'a t

  val allele_set_to_string : set -> string

  val empty : 'a t

  val to_string_full : ('a -> string) -> 'a t -> string

  val of_list : (set * 'a) list -> 'a t

  val singleton : set -> 'a -> 'a t

  val to_list : 'a t -> (set * 'a) list

  val domain : 'a t -> set

  val length : 'a t -> int

  val add : set -> 'a -> 'a t -> 'a t

  val join : 'a t -> 'a t -> 'a t

  val get : set -> 'a t -> 'a t option

  exception StillMissing of string

  val get_exn : set -> 'a t -> 'a t

  val iter : 'a t -> f:(set -> 'a -> unit) -> unit

  val iter_values : 'a t -> f:(int -> 'a -> unit) -> unit

  val fold : 'a t -> init:'b -> f:('b -> set -> 'a -> 'b) -> 'b

  (* Default to not bijective map *)
  val map : ?bijective:bool -> 'a t -> f:('a -> 'b) -> 'b t

  (* Not the perfect name for this function. *)
  val concat_map : 'a t -> f:(set -> 'a -> 'b t) -> 'b t

  val concat_map2 : 'a t -> by:'b t -> f:(set -> 'a -> 'b -> 'c t) -> 'c t

  val concat_map2_partial : 'a t -> by:'b t -> f:(set -> 'a -> 'b -> 'c t) ->
    missing:(set -> 'a -> 'c t) -> 'c t

  val map2 : 'a t -> 'b t -> f:('a -> 'b -> 'c) -> 'c t

  val map3 : 'a t -> 'b t -> 'c t -> f:('a -> 'b -> 'c -> 'd) -> 'd t

  val map4 : 'a t -> 'b t -> 'c t -> 'd t -> f:('a -> 'b -> 'c -> 'd -> 'e) -> 'e t

  val init_everything : 'a -> 'a t

  val map2_partial : 'a t -> by:'b t -> missing:(set -> 'a -> 'c t) ->
    f:('a -> 'b -> 'c) -> 'c t

  val map3_partial : 'a t ->
    by1:'b t -> missing1:(set -> 'a -> 'b t) ->
    by2:'c t -> missing2:(set -> 'a -> 'b -> 'c t) ->
    f:('a -> 'b -> 'c -> 'd) -> 'd t

  val partition_map : 'a t -> f:(set -> 'a -> [< `Fst of 'b | `Snd of 'c ]) ->
     'b t * 'c t

end (* CAMs *)

module MakeCam (AS : Alleles.Set) = struct

  type 'a t = (set * 'a) list

  let empty = []

  let allele_set_to_string s = AS.to_human_readable s

  let to_string t =
    String.concat ~sep:"\n\t"
      (List.map t ~f:(fun (s,_v) ->
        sprintf "%s" (allele_set_to_string s)))

  let to_string_full v_to_s t =
    String.concat ~sep:"\n\t"
      (List.map t ~f:(fun (s,v) ->
        sprintf "%s:%s" (allele_set_to_string s) (v_to_s v)))

  (* let mutate_or_add assoc new_allele_set value =
    let added =
      List.fold assoc ~init:false ~f:(fun added (into, v) ->
        if added then
          added
        else if v = value then begin
          AS.unite ~into new_allele_set;
          true
        end else
          false)
    in
    if added then
      assoc
    else
      (AS.copy new_allele_set, value) :: assoc *)

  (* Union, tail recursive. *)
  let mutate_or_add lst ((alleles, value) as p) =
    let rec loop acc = function
      | (s, v) :: t when v = value -> acc @ (AS.union s alleles, v) :: t
      | h :: t                     -> loop (h :: acc) t
      | []                         -> p :: acc
    in
    loop [] lst

  let add alleles v l = mutate_or_add l (alleles,v)

  let join l1 l2 = List.fold_left l1 ~init:l2 ~f:mutate_or_add

  let of_list l = List.fold_left l ~init:[] ~f:mutate_or_add

  let singleton s a = [s,a]

  let to_list l = l

  let domain = function
    | []             -> AS.init ()
    | (init, _) :: t -> List.fold_left t ~init ~f:(fun u (s, _) -> AS.union u s)

  let length = List.length

  exception StillMissing of string

  let still_missingf fmt =
    ksprintf (fun s -> raise (StillMissing s)) fmt

  let set_assoc_exn to_find t =
    let rec loop to_find acc = function
      | []          -> still_missingf "%s after looking in: %s"
                        (allele_set_to_string to_find) (to_string t)
      | (s, v) :: t ->
          let inter, still_to_find, same_intersect, no_intersect =
            AS.inter_diff to_find s in
          if same_intersect then begin                      (* Found everything *)
            (to_find, v) :: acc
          end else if no_intersect then begin                 (* Found nothing. *)
            loop to_find acc t
          end else begin                                    (* Found something. *)
            loop still_to_find ((inter, v) :: acc) t
          end
    in
    loop to_find [] t

  let set_assoc to_find t =
    try Some (set_assoc_exn to_find t)
    with (StillMissing _) -> None

  let get_exn = set_assoc_exn

  let get = set_assoc

  let iter l ~f = List.iter l ~f:(fun (a, s) -> f a s)

  let iter_values l ~f =
    List.iter l ~f:(fun (a, v) ->
      AS.iter_set_indices a ~f:(fun i -> f i v))

  let fold l ~init ~f = List.fold_left l ~init ~f:(fun b (a, s) -> f b a s)

  let set_assoc_with_mg to_find slst ~missing ~g ~init =
    let rec loop to_find acc = function
      | []          -> add to_find (missing to_find) acc
      | (s, v) :: t ->
          let inter, still_to_find, same_intersect, no_intersect =
            AS.inter_diff to_find s
          in
          if same_intersect then begin                      (* Found everything *)
            add to_find (g v) acc
          end else if no_intersect then begin                 (* Found nothing. *)
            loop to_find acc t
          end else begin                                    (* Found something. *)
            let nacc = add inter (g v) acc in
            loop still_to_find nacc t
          end
    in
    loop to_find init slst

  let map ?bijective l ~f =
    match bijective with
    | Some true         ->                                            (* O(n) *)
      list_map_snd ~f:(fun v -> f v) l
    | Some false | None ->                                          (* O(n^2) *)
      List.fold_left l ~init:[] ~f:(fun acc (s, v) -> add s (f v) acc)

  let set_assoc_k ?n ?missing to_find t ~k ~init =
    let rec loop to_find acc = function
      | []          -> begin match missing with
                       | None -> still_missingf "%s%s after looking in: %s"
                                  (Option.value ~default:"" n)
                                  (allele_set_to_string to_find) (to_string t)
                       | Some m -> m to_find acc
                       end
      | (s, v) :: t ->
          let inter, still_to_find, same_intersect, no_intersect =
            AS.inter_diff to_find s
          in
          if same_intersect then begin                      (* Found everything *)
            k to_find v acc
          end else if no_intersect then begin                 (* Found nothing. *)
            loop to_find acc t
          end else begin                                    (* Found something. *)
            let nacc = k inter v acc in
            loop still_to_find nacc t
          end
    in
    loop to_find init t

  let absorb_k t ~init ~f = List.fold_left t ~init ~f

  let absorb t ~init = absorb_k t ~init ~f:mutate_or_add

  let concat_map l ~f =
    List.fold_left l ~init:[] ~f:(fun init (s, a) -> absorb (f s a) ~init)

  (* The order of set arguments matters for performance. Better to fold over
     the longer list and lookup (set_assoc_k) into the shorter one. Remember
     that the lookup requires a Allele.inter_diff per item! Perhaps it makes
     sense to keep track of the length to avoid O(n) lookups and then
     automatically re-order functional arguments as necessary?

     Probably just need a better data structure. *)
  let concat_map2 l ~by ~f =
    (*printf "%d %d\n" (List.length l) (List.length by); *)
    List.fold_left l ~init:[] ~f:(fun init (s, a) ->
      set_assoc_k s by ~init ~k:(fun intersect b init ->
        absorb (f intersect a b) ~init))

  let concat_map2_partial l ~by ~f ~missing =
    List.fold_left l ~init:[] ~f:(fun init (s, a) ->
      set_assoc_k s by ~init
        ~k:(fun intersect b init -> absorb (f intersect a b) ~init)
        ~missing:(fun sm init -> absorb ~init (missing sm a)))

  let map2 l1 l2 ~f =
    List.fold_left l1 ~init:[] ~f:(fun init (s, a) ->
      set_assoc_k s l2 ~init ~k:(fun intersect b acc ->
        mutate_or_add acc (intersect, f a b)))

  let map3 l1 l2 l3 ~f =
    List.fold_left l1 ~init:[] ~f:(fun init (is1, a) ->
      set_assoc_k ~n:"1" is1 l2 ~init ~k:(fun is2 b init ->
        set_assoc_k ~n:"2" is2 l3 ~init ~k:(fun intersect c acc ->
          mutate_or_add acc (intersect, f a b c))))

  let map4 l1 l2 l3 l4 ~f =
    List.fold_left l1 ~init:[] ~f:(fun init (is1, a) ->
      set_assoc_k is1 l2 ~init ~k:(fun is2 b init ->
        set_assoc_k is2 l3 ~init ~k:(fun is3 c init ->
          set_assoc_k is3 l4 ~init ~k:(fun intersect d acc ->
            mutate_or_add acc (intersect, f a b c d)))))

  let map2_partial l ~by ~missing ~f =
    List.fold_left l ~init:[] ~f:(fun init (s, a) ->
      set_assoc_k s by ~init
        ~k:(fun intercept b acc -> mutate_or_add acc (intercept, f a b))
        ~missing:(fun sm init -> absorb ~init (missing sm a)))

  let map3_partial l ~by1 ~missing1 ~by2 ~missing2 ~f =
    List.fold_left l ~init:[] ~f:(fun init (is1, a) ->
      let k is2 b init =
        let k2 intercept c acc = mutate_or_add acc (intercept, f a b c) in
        set_assoc_k is2 by2 ~init ~k:k2
          ~missing:(fun sm init ->
            absorb_k (missing2 sm a b) ~init ~f:(fun init (s, b) -> k2 s b init))
      in
      set_assoc_k is1 by1 ~init ~k
        ~missing:(fun sm init ->
          absorb_k (missing1 sm a) ~init ~f:(fun init (s, b) -> k s b init)))

  let init_everything v =
    let nothing = AS.init () in
    [AS.complement nothing, v]

  let partition_map l ~f =
    let rec loop bs cs = function
      | []          -> bs, cs
      | (s, a) :: t ->
          match f s a with
          | `Fst b -> loop ((s, b) :: bs) cs t
          | `Snd c -> loop bs ((s, c) :: cs) t
    in
    loop [] [] l

end (* MakeCam *)

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

module MultiplicativeProbability = struct
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

end (* MultiplicativeProbability *)

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

(* For every k there are 3 possible states. *)
type 'a cell =
  { match_  : 'a
  ; insert  : 'a
  ; delete  : 'a
  }

let cell_to_string f c =
  sprintf "{match_: %s; insert: %s; delete: %s}"
    (f c.match_) (f c.insert) (f c.delete)

type 'a recurrences =
  { start   : 'a -> 'a cell
  ; top_row : 'a -> 'a cell -> 'a cell
  ; middle  : 'a -> insert_c:('a cell)
                 -> delete_c:('a cell)
                 -> match_c:('a cell)
                 -> 'a cell
  ; end_    : 'a cell -> 'a
  }

(** What are the values and equations that determine how probabilities are
    calculated in the forward pass. *)
module ForwardCalcs  (R : Ring) = struct

  (* TODO. Avoid the `float_of_int (Phred_score.to_int c) /. -10.` round trip
      between converting to log10p and then back to log10, and just use char's
      instead for the quality calc. *)
  let to_match_prob (base, base_error) =
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

  let g tm ~insert_p read_size =

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

    let start_i = insert_p * t_s_i in
    { start   = begin fun emission_p ->
                  { match_ = emission_p * t_s_m
                  ; insert = start_i
                  ; delete = zero
                  }
                end
    ; top_row = begin fun emission_p insert_c ->
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
                  (*let () =
                      printf "--------%s \n\tmatch_: %s\n\tinsert: %s\n\tdelete: %s\n\tafter : %s\n"
                        (R.to_string emission_p)
                        (cell_to_string R.to_string match_c)
                        (cell_to_string R.to_string insert_c)
                        (cell_to_string R.to_string delete_c)
                        (cell_to_string R.to_string r)
                  in*)
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

end (* ForwardCalcs *)

module MakeWorkspace (CAM: CAMs) = struct

  type 'a entry = 'a cell CAM.t
  type 'a final_entry = 'a CAM.t

  type 'a t =
    { mutable forward             : 'a entry array array
    ; mutable final               : 'a final_entry array
    ; mutable per_allele_emission : 'a array
    }

  let generate number_alleles bigK read_size =
    let just_zeros n = Array.make n 0. in
    { forward             = Array.init bigK ~f:(fun _ -> Array.make read_size CAM.empty)
    ; final               = Array.make bigK CAM.empty
    ; per_allele_emission = just_zeros number_alleles
    }

  let clear ws =
    let bigK = Array.length ws.final in
    let rs   = Array.length ws.forward.(0) in
    let numA = Array.length ws.per_allele_emission in
    ws.forward <- Array.init bigK ~f:(fun _ -> Array.make rs CAM.empty);
    ws.final   <- Array.make bigK CAM.empty;
    Array.fill ws.per_allele_emission ~pos:0 ~len:numA 0.

  let save ws =
    let fname = Filename.temp_file ~temp_dir:"." "forward_workspace" "" in
    let oc = open_out fname in
    Marshal.to_channel oc ws [];
    close_out oc;
    printf "Saved workspace to %s\n" fname

  let load fname =
    let ic = open_in fname in
    let ws : float t = Marshal.from_channel ic in
    close_in ic;
    ws

end (* MakeWorkspace *)

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
  { start_column  : int     (* How many columns of the forward space to fill
                              before using a banded pass. *)
  ; number        : int     (* How many bands to calculate. This logic seems
                              inelegant in the sense that if we reduce the
                              calculation from the full forward space to just
                              this number of bands, why not reduce it further
                              when some bands, inevitably, will have far less
                              probability mass. *)
  ; width         : int     (* How many rows of the band to calculate. *)
  }

let band_default =
  { start_column  = 10
  ; number        = 5
  ; width         = 3
  }

module ForwardGen (R : Ring)(Aset: Alleles.Set) = struct

  module CAM = MakeCam(Aset)

  (* Eh... not the best place for it. *)
  let cam_max = CAM.fold ~init:(neg_infinity) ~f:(fun m _s v -> max m v)

  module Workspace = MakeWorkspace(CAM)

  type base_emissions = (BaseState.t * int) CAM.t

  (* offset and emission probabilties *)
  type 'a oeps = (int * 'a) CAM.t
  type 'a emission_map = 'a oeps PosMap.t

  (* Pairing the observation makes it easier to abstract the regular vs
    reverse complement access pattern. Leaving this as a pair to avoid
    redundant pairing/unpairing.

    I'll use obsp (observation pair) as the variable name. *)
  type obs = char * float

  type 'a fwr =
    { f_start   :  obs -> base_emissions -> 'a Workspace.entry
    ; f_top_row : 'a Workspace.t
                      -> obs  -> base_emissions
                      -> i:int -> k:int
                      -> 'a Workspace.entry
    ; f_middle  : 'a Workspace.t
                      -> obs -> base_emissions
                      -> i:int -> k:int
                      -> 'a Workspace.entry

   (* This isn't the greatest design.... *)
    ; middle_emissions : obs -> base_emissions -> 'a oeps

    ; banded    : 'a Workspace.t -> 'a oeps
                    -> ?prev_col:('a cell) -> ?cur_col:('a Workspace.entry)
                    -> i:int -> k:int -> 'a Workspace.entry

   (* Doesn't use the delete section. *)
    ; f_end_    : 'a Workspace.t -> int -> 'a Workspace.final_entry

    ; emission  : ?spec_rows:(int list list) -> 'a Workspace.t -> int -> 'a array

   (* Combine emission results over indepdent observations
    ; combine   : into:'a array -> 'a array -> unit
   *)
    }

  let per_allele_emission_arr len =
    Array.make len R.one

  module Fc = ForwardCalcs(R)

  let recurrences tm ~insert_p read_size =
    let open Workspace in
    let r = Fc.g tm ~insert_p read_size in

   (* TODO: I could imagine some scenario's where it makes sense to cache,
       precompute or memoize this calculation. The # of base errors isn't
       that large (<100) and there are only 4 bases. So we could be performing
       the same lookup. *)
    let to_em_set obsp emissions =
      CAM.map emissions ~f:(fun (b, offset) ->
        offset, Fc.to_match_prob obsp b)
    in
    { f_start = begin fun obsp emissions ->
                  to_em_set obsp emissions
                  |> CAM.map ~bijective:true
                      ~f:(fun (_offset, emissionp) -> r.start emissionp)
                end
    ; f_top_row = begin fun ws obsp emissions ~i ~k ->
                    to_em_set obsp emissions
                    |> CAM.map2 (ws.forward.(k).(i-1))
                        ~f:(fun insert_c (_offset, emission_p) ->
                              r.top_row emission_p insert_c)
                 end
    ; f_middle = begin fun ws obsp emissions ~i ~k ->
                  let inserts = ws.forward.(k).(i-1) in
                  let ems = to_em_set obsp emissions in
                  CAM.concat_map2 inserts ~by:ems   (* ORDER matters for performance! *)
                      ~f:(fun inters insert_c (offset, emission_p) ->
                            let ks = Pervasives.(+) k offset in
                            let matches = ws.forward.(ks).(i-1) in
                            let deletes = ws.forward.(ks).(i) in
                            let insertsi = CAM.singleton inters insert_c in
                            (* inserti should come before other 2 for performance. *)
                            CAM.map3 insertsi deletes matches
                                ~f:(fun insert_c delete_c match_c ->
                                    r.middle emission_p ~insert_c ~delete_c ~match_c))
               end
    ; middle_emissions = to_em_set
    ; banded  = begin fun ws allele_ems ?prev_col ?cur_col ~i ~k ->
                  let with_insert inters (offset, emission_p) insert_c =
                    let calc insert_c delete_c match_c =
                      r.middle emission_p ~insert_c ~delete_c ~match_c
                    in
                    let ks = Pervasives.(+) k offset in
                    let matches = ws.forward.(ks).(i-1) in
                    let deletes = ws.forward.(ks).(i) in
                    let insertsi = CAM.singleton inters insert_c in
                    CAM.map3_partial insertsi
                      ~by1:deletes
                      ~missing1:(fun missing_deletes _insert_c ->
                          let default = CAM.singleton missing_deletes Fc.zero_cell in
                          Option.value_map ~default cur_col ~f:(fun as_ ->
                            Option.value (CAM.get missing_deletes as_) ~default))
                      ~by2:matches
                      ~missing2:(fun missing_matches _insert_c _delete_c ->
                        CAM.singleton missing_matches
                          (Option.value prev_col ~default:Fc.zero_cell))
                      ~f:calc
                  in
                  let inserts = ws.forward.(k).(i-1) in
                    CAM.concat_map2_partial allele_ems ~by:inserts
                      ~missing:(fun missing_inserts ep_pair ->
                          match prev_col with
                          | None -> invalid_argf "At %d %d looking for inserts still missing %s"
                                      k i (CAM.allele_set_to_string missing_inserts)
                          | Some v -> with_insert missing_inserts ep_pair v)
                      ~f:with_insert
                end
    ; f_end_  = begin fun ws k ->
                  CAM.map ~bijective:true ws.forward.(k).(read_size-1) ~f:r.end_
                end
    ; emission  = begin fun ?spec_rows ws len ->
                    let open R in
                    let ret = Array.make len zero in
                    let update_cam l =
                      CAM.iter_values l ~f:(fun i v -> ret.(i) <- ret.(i) + v)
                    in
                    let () =
                      match spec_rows with
                      | None   -> Array.iter ws.final ~f:update_cam
                      | Some l -> List.iter l ~f:(fun rows ->
                                    List.iter rows ~f:(fun k ->
                                      update_cam ws.final.(k)))
                    in
                    ret
                  end
(*    ; combine   = begin fun ~into em ->
                    let open R in
                    Array.iteri em ~f:(fun i e -> into.(i) <- into.(i) * e)
                  end *)
    }

    (* Regular (fill in everything) forward pass. *)
    module Regular = struct

      open Workspace

      let pass ws recurrences emissions_a rows columns a =
        (* special case the first row. *)
        ws.forward.(0).(0) <-
          recurrences.f_start (a 0) emissions_a.(0);
        for i = 1 to columns do
          ws.forward.(0).(i) <-
            recurrences.f_top_row ws (a i) ~i ~k:0 emissions_a.(0)
        done;
        (* All other rows. *)
        for k = 1 to rows do
          let ek = emissions_a.(k) in
          ws.forward.(k).(0) <- recurrences.f_start (a 0) ek;
          for i = 1 to columns do
            ws.forward.(k).(i) <-
              recurrences.f_middle ws (a i) ~i ~k ek
          done
        done

      let final ws recurrences rows ~number_alleles =
        for k = 0 to rows do
          ws.final.(k) <- recurrences.f_end_ ws k
        done;
        ws.per_allele_emission <- recurrences.emission ws number_alleles

      let both ws recurrences emissions_a ~rows ~columns ~number_alleles a =
        pass ws recurrences emissions_a rows columns a;
        final ws recurrences rows number_alleles

  end (* Regular *)

  (* Banded Pass logic **)
  module Bands = struct

    (* 1. Identify bands. *)
    let select_specific_band_indices ws c =
      let open Workspace in
      let lg = largest c.number in        (* compares by match_, insert, delete *)
      let ev = CAM.init_everything [] in
      Array.fold_left ws.forward ~init:(0, ev) ~f:(fun (k, acc) r ->
        let nacc =
          CAM.map2_partial acc ~by:r.(c.start_column)
            ~f:(fun lst c -> lg c k lst)
            ~missing:(fun s v -> CAM.singleton s v)
        in
        k + 1, nacc)
      |> snd

    let expand_allele_set_map l =
      CAM.to_list l
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
      CAM.concat_map inds ~f:(fun s (bv, ilst) ->
        let i = List.hd_exn ilst in
        CAM.get_exn s emissions.(i)
        |> CAM.map ~bijective:true ~f:(fun (_bs, o) ->
            if o = min_int then
              (bv, ilst)
            else
              (bv, (i + o) :: ilst)))

    let find_indices_below increments inds =
      let bigK = Array.length increments in
      CAM.concat_map inds ~f:(fun s (bv, ilst) ->
        let i = List.hd_exn ilst in
        if i = bigK then
          CAM.singleton s (bv, ilst)
        else
          CAM.get_exn s increments.(i)
          |> CAM.map ~bijective:true ~f:(fun o ->
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
      let lnds = CAM.map inds ~bijective:true ~f:(fun st -> st, [to_index st]) in
      let ai = find_indices_above_n c.width emissions_a lnds in
      let bi = find_indices_below_n c.width increment_a lnds in
      (* tl_exn -> drop the center band, so we don't duplicate it. *)
      CAM.map2 ai bi ~f:(fun (st, a) (_st2, b) -> st, a @ (List.tl_exn (List.rev b)))

    (* This step (see get_exn) partitions the current bands based on the
      last cell's. Since these value are already away from the optimal point
      in the band (that would be width above), I'm not certain that this
      is necessary. We're using this value as just an approximation in lieu
      of filling in the forward matrix. Specifically, we can have a less strict
      cell equality test, so that the allele sets are joined together.

      This is a general notion to consider; when are cells close enough
      (difference just due to numerical rounding) that it isn't worth the
      split. *)
    let lookup_previous_values ws col bands =
      let open Workspace in
      CAM.concat_map bands ~f:(fun s (bv, rows) ->
        let end_row = List.last rows |> Option.value_exn ~msg:"empty rows!" in
        CAM.get_exn s ws.forward.(end_row).(col)
        |> CAM.map ~bijective:true ~f:(fun lv ->
              (rows, bv, lv)))

    type 'a t =
      { rows        : int list
      ; best_value  : 'a cell
      ; last_value  : 'a cell         (* Doesn't have to occur at end_row *)
      ; alleles     : set
      } (*The order determines comparison. *)

    let to_string sc t =
      sprintf "rows: %s\tbv: [m: %s; i: %s; d: %s]\tlv: [m: %s; i: %s; d: %s]\n\t\ts: %s"
        (String.concat ~sep:";" (List.map t.rows ~f:(sprintf "%d")))
        (sc t.best_value.match_)
        (sc t.best_value.insert)
        (sc t.best_value.delete)
        (sc t.last_value.match_)
        (sc t.last_value.insert)
        (sc t.last_value.delete)
        (Aset.to_human_readable t.alleles)

    (* TODO: Should we shift these down 1 ala next_band ? *)
    let setup emissions_a increment_a c ws =
      select_specific_band_indices ws c
      |> expand_allele_set_map
      |> group_by_allele_value
      (* We have to keep the Allele.Set bands separate, not in an
        CAM.t to avoid overlaps. *)
      |> List.map ~f:(fun p ->
          CAM.of_list [p]
          |> to_bands emissions_a increment_a c ~to_index:(fun (bv, i) -> i)
          |> lookup_previous_values ws c.start_column
          |> CAM.to_list
          |> List.map ~f:(fun (alleles, (rows, (best_value, _), last_value)) ->
              { rows; best_value; last_value; alleles}))
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
        3. Keep track of the calculated rows for an allele. This allows us to
          adaptively figure out the rows of the next pass by moving width away
          from the best_row. See [find_next_row_from_fill_state].
      *)
    type fill_state =
      { best_row : int          (* where is the best, lowest match likelihood row *)
      ; best_c   : float cell
      ; worse    : int          (* Number of likelihoods < than best_c.match_ *)
      ; last_c   : float cell
      ; nrows    : int list     (* Where we're calculating. Since there might be
                                  gaps, the width needs to look inside this list
                                  for the next start/end_row *)
      }

    let init_fill_state row cell =
      { best_row  = row
      ; best_c    = cell
      ; worse     = 0
      ; last_c    = cell
      ; nrows     = [row]
      }

    let update_fill_state row fs cell =
      if cell.match_ > fs.best_c.match_ then
        { best_row = row
        ; best_c   = cell
        ; worse    = 0
        ; last_c   = cell
        ; nrows    = row :: fs.nrows
        }
      else
        { fs with worse  = fs.worse + 1
                ; last_c = cell
                ; nrows  = row :: fs.nrows
        }

    (* rows are in reverse, descending, order! *)
    let to_next_rows width best_row rows =
      let rec find_best acc = function
        | []     -> invalid_argf "Didn't find best row."
        | h :: t ->
            if h = best_row then
              (* TODO: These can silently take less than width. *)
              let s = List.take t width in
              let e = List.take acc width in
              (List.rev s) @ (h :: e)
            else
              find_best (h :: acc) t
      in
      find_best [] rows

    let find_next_row_from_fill_state c fs =
      to_next_rows c.width fs.best_row fs.nrows

    let next_band emissions_a increment_a c fs_map =
      CAM.concat_map fs_map ~f:(fun alleles fs ->
        (* Shift the band, by adjusting around best_row,  for next column *)
        CAM.get_exn alleles increment_a.(fs.best_row)
        (* Now fill in the width. *)
        |> to_bands emissions_a increment_a c ~to_index:(fun nr -> nr)
        |> CAM.map ~bijective:true ~f:(fun (_br,rows) -> (rows, fs.best_c, fs.last_c)))
      |> CAM.to_list
      |> List.map ~f:(fun (alleles, (rows, best_value, last_value)) ->
          { rows ; alleles ; best_value ; last_value })

    let fill_next emissions_a increment_a c recurrences em_map ws obsp i b col_values =
      let open Workspace in
      (*let () =
        let base, base_prob = obsp in
        printf "current bands for %c %f rows:%s at %d \n\t%s\n"
          base base_prob
          (String.concat ~sep:";" (List.map b.rows ~f:(sprintf "%d"))) i
            (to_string (sprintf "%f") b)
      in*)
      let cur_col = CAM.get b.alleles col_values in
      let update ?cur_col emp k alleles =
        let em_values, nem_map =
          try
            let emv = PosMap.find k emp in
            emv, emp
          with Not_found ->
            let es = recurrences.middle_emissions obsp emissions_a.(k) in
            let nemp = PosMap.add ~key:k ~data:es emp in
            es, nemp
        in
        let allele_emissions = CAM.get_exn alleles em_values in
        let entry =
          recurrences.banded ws allele_emissions
            (* Poor design: No harm in adding this as banded will only use this
              value in the missing case. So we're not going to track that we're
              at the right row. *)
            ~prev_col:b.last_value ?cur_col
            ~i ~k
        in
        ws.forward.(k).(i) <- CAM.join entry ws.forward.(k).(i);
        nem_map, entry
      in
      match b.rows with
      | []                -> invalid_argf "empty rows"
      | start_row :: trows -> begin
          let nem_map, first_entry = update ?cur_col em_map start_row b.alleles in
          let state =
            CAM.map first_entry ~bijective:true
              ~f:(init_fill_state start_row)
          in
          let update_fill_state prev nk cur =
            CAM.map2_partial prev ~by:cur
              ~f:(update_fill_state nk)
              ~missing:(fun s p -> CAM.singleton s p)
          in
          let rec loop em_map acc fill_state not_full_alleles cur_col = function
            | []        -> invalid_argf "empty row, was there only one row?"
            | k :: rows ->
                let nem_map, entry = update ~cur_col em_map k not_full_alleles in
                let new_fill_state = update_fill_state fill_state k entry in
                if rows <> [] then                  (* Still have remaining rows to fill. *)
                  loop nem_map acc new_fill_state not_full_alleles entry rows
                else begin                          (* Filled all that we're supposed to. *)
                  let full, not_full_state =
                    CAM.partition_map new_fill_state ~f:(fun _s fs ->
                      if fs.worse >= c.width then `Fst fs else `Snd fs)
                  in
                  let full_bands = next_band emissions_a increment_a c full in
                  (*let () =
                    printf "new bands for k:%d at %d: %d\n" k i (List.length full_bands);
                    List.iter full_bands ~f:(fun b ->
                      printf "\t%s\n" (to_string (sprintf "%f") b))
                  in*)
                  let nacc = full_bands @ acc in
                  if CAM.length not_full_state = 0 ||     (* Nothing left to fill -> Done *)
                    k = Array.length increment_a then                     (* Reached end! *)
                    nacc, nem_map
                  else begin
                    CAM.fold not_full_state ~init:(nacc, nem_map)
                      ~f:(fun init alleles state ->
                            (*printf "not_full_recurse %d %s in %s\n%!"
                              k (Alleles.Set.to_human_readable alleles)
                                (Alleles.Set.to_human_readable (CAM.domain increment_a.(k))); *)
                            CAM.get_exn alleles increment_a.(k)
                            |> CAM.fold ~init ~f:(fun (acc, em_map) alleles2 next_row ->
                                loop em_map acc (CAM.singleton alleles2 state) alleles2
                                entry [next_row]))
                  end
                end
          in
          loop nem_map [] state b.alleles first_entry trows
      end

    let fill_end recurrences ws b =
      let open Workspace in
      List.iter b.rows ~f:(fun k ->
        ws.final.(k) <- CAM.join (recurrences.f_end_ ws k) ws.final.(k))

    let pass c emissions_a increment_a number_alleles ws recurrences last_read_index a =
      (* order matters for passing along last_col *)
      let first_bands = setup emissions_a increment_a c ws |> List.sort ~cmp:compare in
      (*printf "first_bands %d \n" (List.length first_bands);
      List.iter first_bands ~f:(fun t -> printf "\t%s\n" (to_string (sprintf "%f") t)); *)
      let banded_middle bands start_column =
        let rec loop bands i =
          let new_bands_to_flatten, _last_em_map, _last_col_values =
            List.fold_left bands ~init:([], PosMap.empty, CAM.empty)
              ~f:(fun (acc, em_map, col_values) b ->
                    let nb, nem_map =
                      fill_next emissions_a increment_a c recurrences em_map ws (a i) i b col_values
                    in
                    let ncol_values =
                      List.map nb ~f:(fun t -> t.alleles, t.last_value)
                      |> CAM.of_list
                    in
                    nb :: acc, nem_map, ncol_values)
          in
          if i = last_read_index then
            bands (* We need these bands for end_ *)
          else begin
            let new_bands =
              List.flatten new_bands_to_flatten
              (* The default comparator will sort first by rows (first field),
                and within the int lists, the comparison is by the values,
                with smaller length lists taking precedence. *)
              |> List.sort ~cmp:compare
          in
          (*printf "bands at %d %d \n" i (List.length new_bands);
          List.iter new_bands ~f:(fun t -> printf "\t%s\n" (to_string (sprintf "%f") t)); *)
          loop new_bands (i + 1)
          end
        in
        loop bands start_column
      in
      let open Workspace in
      let banded_end bands =
        List.iter bands ~f:(fill_end recurrences ws);
        let spec_rows = List.map bands ~f:(fun b -> b.rows) in
        ws.per_allele_emission <- recurrences.emission ~spec_rows ws number_alleles
      in
      banded_end (banded_middle first_bands (c.start_column + 1))

  end (* Bands *)

end (* ForwardGen *)

module Forward = ForwardGen(MultiplicativeProbability)
module ForwardLogSpace = ForwardGen (LogProbabilities)

type t =
  { align_date      : string
  ; number_alleles  : int
  ; aset            : (module Alleles.Set)
  ; alleles         : string list
  ; merge_map       : (string * string) list
  ; emissions_a     : (Alleles.set * (BaseState.t * int)) list array
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
      let module CAM = MakeCam(Aset) in
      let emissions_a, position_map =
        initialize_base_array_and_position_map aset mp.reference mp.ref_elems
      in
      let aaa = add_alternate_allele aset mp.reference ~position_map in
      List.iter ~f:(fun (allele, altseq) -> aaa allele altseq emissions_a) nalt_elems;
      let emissions_a =
        (* TODO: Move the CAM logic up into the construction algorithms *)
        Array.map emissions_a ~f:(fun l ->
          List.map l ~f:(fun (b, s) -> (s, b)) |> CAM.of_list)
      in
      let increment_a = Array.make (Array.length emissions_a - 1) CAM.empty in
      Array.iteri emissions_a ~f:(fun i s ->
        if i = 0 then () else
          CAM.map s ~f:(fun (_b, v) -> v)
          |> CAM.iter ~f:(fun s g ->
              let k = i + g in
              increment_a.(k) <- CAM.add s i increment_a.(k)));
      Ok { align_date = mp.align_date
         ; number_alleles = Alleles.length allele_index
         ; aset
         ; alleles
         ; merge_map
         ; emissions_a
         ; increment_a
         }
  end

let debug_ref = ref false

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

let compare_emissions e1 e2 =
  let r1 = Array.fold_left e1 ~init:neg_infinity ~f:max in
  let r2 = Array.fold_left e2 ~init:neg_infinity ~f:max in
  r1 >= r2

let access rc read read_prob =
  let m = Array.length read_prob - 1 in
  if rc then
    fun i ->
      complement (String.get_exn read (m - i))
      , Array.get read_prob (m - i)
  else
    fun i -> String.get_exn read i , Array.get read_prob i

type mapped_stats =
  { regular     : (float * string) list
  ; rpositions  : (float * int) list
  ; complement  : (float * string) list
  ; cpositions  : (float * int) list
  }

let mapped_stats_to_string ?(sep='\t') ms =
  let l_to_s fmt l =
    String.concat ~sep:";" (List.map l ~f:(fun (l,a) -> sprintf fmt  a l))
  in
  let al_to_s l = l_to_s "%s:%0.2f" l in
  let pl_to_s l = l_to_s "%d:%0.2f" l in
  if fst (List.hd_exn ms.rpositions) > fst (List.hd_exn ms.cpositions) then
    sprintf "R %s%c%s%c%s%c%s"
      (al_to_s ms.regular)    sep
      (pl_to_s ms.rpositions) sep
      (al_to_s ms.complement) sep
      (pl_to_s ms.cpositions)
  else
    sprintf "C %s%c%s%c%s%c%s"
      (al_to_s ms.complement) sep
      (pl_to_s ms.cpositions) sep
      (al_to_s ms.regular)    sep
      (pl_to_s ms.rpositions)

let best_stat ms =
  max (fst (List.hd_exn ms.rpositions)) (fst (List.hd_exn ms.cpositions))

(*** Full Forward Pass *)

type proc =
  { output_ws_array : unit -> float array                     (* Allocate the right size *)
  ; doit            : bool -> string -> float array -> unit
  ; best_alleles    : unit -> (float * string) list
  ; best_positions  : unit -> (float * int) list
  ; per_allele_llhd : unit -> float array                     (* Pass'es output. *)
  ; save_workspace  : unit -> unit
  }

(* Return
  1. a function to process one read
  2. the workspace that it uses
  3. a function to combin results *)
let setup_single_pass ?band read_size t =
  let { number_alleles; emissions_a; increment_a; aset; alleles; _ } = t in
  if !debug_ref then save_pphmm t;
  let bigK = Array.length emissions_a in
  let tm = Phmm.TransitionMatrix.init ~ref_length:bigK read_size in
  let insert_p = 0.25 in
  let module AS = (val aset : Alleles.Set) in
  let module F = ForwardLogSpace(AS) in
  let r = F.recurrences tm ~insert_p read_size in
  let ws = F.Workspace.generate number_alleles bigK read_size in
  let last_row = bigK - 1 in
  let last_read_index = read_size - 1 in
  let lg5 a i lst = largest 5 a i lst in  (* TODO: expose this 5 if it becomes useful *)
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
  let best_alleles () = best_alleles ws.F.Workspace.per_allele_emission in
  let best_positions () = best_positions ws.F.Workspace.final in
  let per_allele_llhd () = ws.F.Workspace.per_allele_emission in
  let save_workspace () = F.Workspace.save ws in
  let output_ws_array () = F.per_allele_emission_arr t.number_alleles in
  let normal () =
    { doit =
        begin fun rc read read_prob ->
          let a = access rc read read_prob in
          F.Regular.both ws r emissions_a a
            ~rows:last_row ~columns:last_read_index ~number_alleles
        end
    ; best_alleles
    ; best_positions
    ; per_allele_llhd
    ; save_workspace
    ; output_ws_array
    }
  in
  let banded c =
    { doit =
        begin fun rc read read_prob ->
          let a = access rc read read_prob in
          (* clear the forward/final array since the banded pass algorithm relies on
                unfilled cells to indicate boundaries; places where we use heuristics.*)
          F.Workspace.clear ws;
          F.Regular.pass ws r emissions_a last_row c.start_column a;
          F.Bands.pass c emissions_a increment_a number_alleles
            ws r last_read_index a
        end
    ; best_alleles
    ; best_positions
    ; per_allele_llhd
    ; save_workspace
    ; output_ws_array
    }
  in
  match band with
  | None                                          -> normal ()
  | Some c when c.start_column >= last_read_index -> normal ()
  | Some c                                        -> banded c

let mapper pass read read_prob =
  pass.doit false read read_prob;                                 (* Regular. *)
  let regular     = pass.best_alleles () in
  let rpositions  = pass.best_positions () in
  pass.doit true read read_prob;                                (* Complement *)
  let complement  = pass.best_alleles () in
  let cpositions  = pass.best_positions () in
  { regular; rpositions; complement; cpositions}

let reducer pass ?(check_rc=true) read read_prob =
  if check_rc then begin
    pass.doit false read read_prob;                               (* Regular. *)
    let regular = Array.copy (pass.per_allele_llhd ()) in
    pass.doit true read read_prob;                              (* Complement *)
    let complement = pass.per_allele_llhd () in
    if compare_emissions regular complement then begin
      regular
    end else begin
      complement
    end
  end else begin
    pass.doit false read read_prob;
    pass.save_workspace ();
    pass.per_allele_llhd ()
  end

let forward_pass ?band mode t read_size =
  let pass = setup_single_pass ?band read_size t in
  match mode with
  | `Mapper   -> `Mapper (mapper pass)
  | `Reducer  -> `Reducer (pass.output_ws_array (), reducer pass)

module Many = struct

  type obs =
    { regular     : float array
    ; complement  : float array
    }

  let max_arr = Array.fold_left ~init:neg_infinity ~f:max

  type best_state =
    { name  : string
    ; maxl  : float
    ; llhd  : float array
    }

  let to_bs (name, {regular; complement}) =
    let rm = max_arr regular in
    let cm = max_arr complement in
    if rm >= cm then
      { name; maxl = rm ; llhd = regular}
    else
      { name; maxl = cm ; llhd = complement}

  let best_bs bs1 bs2 =
    (*printf "%s %f vs %s %f\n%!" bs1.name bs1.maxl bs2.name bs2.maxl; *)
    if bs1.maxl >= bs2.maxl then bs1 else bs2

  let best = function
    | []          -> invalid_argf "Can't select best from empty!"
    | h :: t ->
        let init = to_bs h in
        List.fold_left t ~init ~f:(fun bs p ->
          best_bs bs (to_bs p))

  let add_log_likelihoods ~into nl =
    let n = Array.length into in
    for i = 0 to n - 1 do
      into.(i) <- into.(i) +. nl.(i)
    done

  let merge current state =
    let b = best current in
    List.iter state ~f:(fun (n, into) ->
      if n = b.name then
        add_log_likelihoods ~into  b.llhd
      else ());
    state

end (* Many *)

let reporter pass read read_prob =
  pass.doit false read read_prob;                                 (* Regular. *)
  let regular = Array.copy (pass.per_allele_llhd ()) in
  pass.doit true read read_prob;                                (* Complement *)
  let complement = pass.per_allele_llhd () in
  { Many.regular; complement}


let forward_pass_m ?band mode tlst read_size =
  let passes =
    List.map tlst ~f:(fun (name,t) -> name, setup_single_pass ?band read_size t)
  in
  match mode with
  | `Mapper ->
      `Mapper (fun read read_prob ->
        List.map passes ~f:(fun (n, p) -> n, mapper p read read_prob))
  | `Reporter ->
      `Reporter ( List.map passes ~f:(fun (name, pass) -> name, pass.output_ws_array ())
                , fun read read_prob state ->
                  let current = List.map passes ~f:(fun (n, p) ->
                    n, reporter p read read_prob) in
                  Many.merge current state)
