(* CAM = Compressed Allele Map

   Since our PHMM is parameterized by alleles, we have to keep track many
   values on a per allele basis. This module aims to provide an abstraction
   for this map: allele -> 'a.

   It is different than the Alleles.Map module (and perhaps that one should be
   replaced or deprecated) because the implementation tries to be succinct, by
   compressing the values and avoiding O(n) (n = # of alleles) maps/folds.
   Specifically, the operations are performed on the unique values in the map.
*)

open Util

type set = Alleles.set

module type M = sig

  type 'a t

  val allele_set_to_string : set -> string

  val empty : 'a t

  val to_string_full : ('a -> string) -> 'a t -> string

  val of_list : ?eq:('a -> 'a -> bool) -> (set * 'a) list -> 'a t

  val singleton : set -> 'a -> 'a t

  val to_list : 'a t -> (set * 'a) list

  val domain : 'a t -> set

  val length : 'a t -> int

  val add : ?eq:('a -> 'a -> bool) -> set -> 'a -> 'a t -> 'a t

  val join : ?eq:('a -> 'a -> bool) -> 'a t -> 'a t -> 'a t

  val get : set -> 'a t -> 'a t option

  val get_value : set -> 'a t -> 'a option

  exception StillMissing of string

  val get_exn : set -> 'a t -> 'a t

  val iter : 'a t -> f:(set -> 'a -> unit) -> unit

  val iter_values : 'a t -> f:(int -> 'a -> unit) -> unit

  val fold : 'a t -> init:'b -> f:('b -> set -> 'a -> 'b) -> 'b

  (* Default to not bijective map *)
  val map : ?bijective:bool -> 'a t -> f:('a -> 'b) -> 'b t

  (* Not the perfect name for this function. *)
  val concat_map : ?eq:('b -> 'b -> bool)
                  -> 'a t
                  -> f:(set -> 'a -> 'b t)
                  -> 'b t

  val concat_map2 : ?eq:('c -> 'c -> bool)
                    -> 'a t
                    -> by:'b t
                    -> f:(set -> 'a -> 'b -> 'c t)
                    -> 'c t

  val concat_map2_partial : ?eq:('c -> 'c -> bool)
                            -> 'a t
                            -> by:'b t
                            -> f:(set -> 'a -> 'b -> 'c t)
                            -> missing:(set -> 'a -> 'c t)
                            -> 'c t

  val map2 : ?eq:('c -> 'c -> bool)
            -> 'a t
            -> 'b t
            -> f:('a -> 'b -> 'c)
            -> 'c t

  val map3 : ?eq:('d -> 'd -> bool)
            -> 'a t
            -> 'b t
            -> 'c t
            -> f:('a -> 'b -> 'c -> 'd)
            -> 'd t

  val init_everything : 'a -> 'a t

  val map2_partial : ?eq:('c -> 'c -> bool)
                    -> 'a t
                    -> by:'b t
                    -> missing:(set -> 'a -> 'c t)
                    -> f:('a -> 'b -> 'c)
                    -> 'c t

  val map3_partial : ?eq:('d -> 'd -> bool)
                    -> 'a t
                    -> by1:'b t
                    -> missing1:(set -> 'a -> 'b t)
                    -> by2:'c t
                    -> missing2:(set -> 'a -> 'b -> 'c t)
                    -> f:('a -> 'b -> 'c -> 'd)
                    -> 'd t

  val partition_map : 'a t -> f:(set -> 'a -> [< `Fst of 'b | `Snd of 'c ]) ->
     'b t * 'c t

end (* M *)

module Make (AS : Alleles.Set) : M = struct

  type 'a t = (set * 'a) list

  let empty = []

  let allele_set_to_string s = AS.to_human_readable s

  let to_string = function
    | []  -> "empty"
    | lst -> String.concat ~sep:"\n\t"
                (List.map lst ~f:(fun (s,_v) ->
                  sprintf "%s" (allele_set_to_string s)))

  let to_string_full v_to_s = function
    | []  -> "empty"
    | lst -> String.concat ~sep:"\n\t"
                (List.map lst ~f:(fun (s,v) ->
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
  let mutate_or_add ?eq lst ((alleles, value) as p) =
    let eq = Option.value eq ~default:(=) in
    let rec loop acc = function
      | (s, v) :: t when eq v value -> acc @ (AS.union s alleles, v) :: t
      | h :: t                      -> loop (h :: acc) t
      | []                          -> p :: acc
    in
    loop [] lst

  let add ?eq alleles v l = mutate_or_add ?eq l (alleles,v)

  let join ?eq l1 l2 = List.fold_left l1 ~init:l2 ~f:(mutate_or_add ?eq)

  let of_list ?eq l = List.fold_left l ~init:[] ~f:(mutate_or_add ?eq)

  let singleton s a = [s,a]

  let to_list l = l

  let domain = function
    | []             -> AS.init ()
    | (init, _) :: t -> List.fold_left t ~init ~f:(fun u (s, _) -> AS.union u s)

  let length = List.length

  exception StillMissing of string

  let still_missingf fmt =
    ksprintf (fun s -> raise (StillMissing s)) fmt

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

  let set_assoc_exn to_find t =
    set_assoc_k to_find t
      ~init:[]
      ~k:(fun to_find v acc -> (to_find, v) :: acc)
      ~missing:(fun to_find t -> still_missingf "%s after looking in: %s"
                        (allele_set_to_string to_find) (to_string t))

  let set_assoc to_find t =
    try Some (set_assoc_exn to_find t)
    with (StillMissing _) -> None

  let get_exn = set_assoc_exn

  let get = set_assoc

  let get_value s t =
    Option.bind (get s t) ~f:(function
      | [_,v] -> Some v
      | _     -> None)

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
      List.map_snd ~f:(fun v -> f v) l
    | Some false | None ->                                          (* O(n^2) *)
      List.fold_left l ~init:[] ~f:(fun acc (s, v) -> add s (f v) acc)

  let absorb_k t ~init ~f = List.fold_left t ~init ~f

  let absorb ?eq t ~init = absorb_k t ~init ~f:(mutate_or_add ?eq)

  let concat_map ?eq l ~f =
    List.fold_left l ~init:[] ~f:(fun init (s, a) -> absorb ?eq (f s a) ~init)

  (* Unroll the list so that there's less of it to search. *)
  let set_assoc_u ?n ?missing to_find t ~k ~init =
    let rec loop to_find acc macc = function
      | []          -> begin match missing with
                       | None -> still_missingf "%s%s after looking in: %s"
                                  (Option.value ~default:"" n)
                                  (allele_set_to_string to_find) (to_string t)
                       | Some m -> m to_find acc macc
                       end
      | (s, v) :: tl ->
          let inter, still_to_find, left_in_s,
              found_everything, no_intersect, left = AS.split3 to_find s
          in
          if found_everything then begin                    (* Found everything *)
            (* Where we end up using set_assoc_u, it turns out that the cost of
               reversing macc to preserver order doesn't improve performance. *)
            let nmacc = macc @ (if left then (left_in_s, v) :: tl else tl) in
            (k to_find v acc), nmacc
          end else if no_intersect then begin                 (* Found nothing. *)
            loop to_find acc ((s,v) :: macc) tl
          end else begin                                    (* Found something. *)
            let nacc = k inter v acc in
            let nmacc = (left_in_s, v) :: macc in
            loop still_to_find nacc nmacc tl
          end
    in
    loop to_find init [] t

(*  let concat_map2 ?eq l ~by ~f =
    List.fold_left l ~init:[] ~f:(fun init (s, a) ->
      set_assoc_k s by ~init ~k:(fun intersect b init ->
        absorb ?eq (f intersect a b) ~init))*)

  let concat_map2 ?eq l ~by ~f =
    List.fold_left l ~init:([], by) ~f:(fun (init, nby) (s, a) ->
      set_assoc_u s nby ~init
        ~k:(fun intersect b init    -> absorb ?eq (f intersect a b) ~init))
    |> fst

(*  let concat_map2_partial ?eq l ~by ~f ~missing =
    List.fold_left l ~init:[] ~f:(fun init (s, a) ->
      set_assoc_k s by ~init
        ~k:(fun intersect b init -> absorb ?eq (f intersect a b) ~init)
        ~missing:(fun sm init -> absorb ?eq ~init (missing sm a))) *)

  let concat_map2_partial ?eq l ~by ~f ~missing =
    List.fold_left l ~init:([], by) ~f:(fun (init, nby) (s, a) ->
      set_assoc_u s nby ~init
        ~k:(fun intersect b init -> absorb ?eq (f intersect a b) ~init)
        ~missing:(fun sm init nby -> absorb ?eq ~init (missing sm a), nby))
    |> fst

  (*let map2 ?eq l1 l2 ~f =
    List.fold_left l1 ~init:[] ~f:(fun init (s, a) ->
      set_assoc_k s l2 ~init ~k:(fun intersect b acc ->
        mutate_or_add ?eq acc (intersect, f a b))) *)

  let map2 ?eq l1 l2 ~f =
    List.fold_left l1 ~init:([],l2) ~f:(fun (init, nl2) (s, a) ->
      set_assoc_u s nl2 ~init ~k:(fun intersect b acc ->
        mutate_or_add ?eq acc (intersect, f a b)))
    |> fst

  let map3 ?eq l1 l2 l3 ~f =
    List.fold_left l1 ~init:[] ~f:(fun init (is1, a) ->
      set_assoc_k ~n:"1" is1 l2 ~init ~k:(fun is2 b init ->
        set_assoc_k ~n:"2" is2 l3 ~init ~k:(fun intersect c acc ->
          mutate_or_add ?eq acc (intersect, f a b c))))

(* TODO: For map3 it is TOO costly to construct nl2 & nl3. This asymmetry
   really demonstrates the difficulty of making this lookup efficient.
  let map3 ?eq l1 l2 l3 ~f =
    List.fold_left l1 ~init:([],l2,l3) ~f:(fun (acc, nl2, nl3) (is1, a) ->
      set_assoc_u ~n:"1" is1 nl2 ~init:(acc,nl3) ~k:(fun is2 b (init, nl3) ->
        set_assoc_u ~n:"2" is2 nl3 ~init ~k:(fun intersect c acc ->
          mutate_or_add ?eq acc (intersect, f a b c)))
      |> fun ((acc,ml3),ml2) -> (acc,ml2,ml3))
    |> fun (x, _, _) -> x
    *)

  let map2_partial ?eq l ~by ~missing ~f =
    List.fold_left l ~init:[] ~f:(fun init (s, a) ->
      set_assoc_k s by ~init
        ~k:(fun intercept b acc -> mutate_or_add ?eq acc (intercept, f a b))
        ~missing:(fun sm init -> absorb ?eq ~init (missing sm a)))

  let map3_partial ?eq l ~by1 ~missing1 ~by2 ~missing2 ~f =
    List.fold_left l ~init:[] ~f:(fun init (is1, a) ->
      let k is2 b init =
        let k2 intercept c acc = mutate_or_add ?eq acc (intercept, f a b c) in
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

end (* Make *)
