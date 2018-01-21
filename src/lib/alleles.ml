
open Util

type allele = string [@@deriving eq, ord]
let compare = compare_allele
let equal = equal_allele

type index =
  { size      : int
  ; to_index  : int StringMap.t
        [@equal StringMap.equal ~cmp:(fun (a:int) b -> a = b)]
        [@compare StringMap.compare ~cmp:Int.compare]
  ; to_allele : string array
  }
  [@@deriving eq, ord]

let empty_index =
  { size      = 0
  ; to_index  = StringMap.empty
  ; to_allele = [||]
  }

let length { size; _} = size

let index lst =
  let to_allele = Array.of_list (List.sort ~cmp:compare_allele lst) in
  let size, to_index =
    Array.fold_left to_allele ~init:(0, StringMap.empty )
      ~f:(fun (i, sm) a -> (i + 1, StringMap.add a i sm))
  in
  { size; to_index; to_allele }

let to_alleles { to_allele; _ } = Array.to_list to_allele

module CompressNames = struct

  let rejoin_colon = String.concat ~sep:":"

  module CompressIntSequencess = struct

    type t =
      | And of (t * t)
      | From of (int * int)
      | Just of int

    let rec to_last = function
      | And  (a, b) -> sprintf "%s,%s" (to_last a) (to_last b)
      | From (a, b) -> sprintf "%0.2d-%0.2d" a b
      | Just n      -> sprintf "%0.2d" n

    let rec last_val = function
      | And (_, v)  -> last_val v
      | From (_, l) -> l
      | Just n      -> n

    let rec add c v =
      match c with
      | Just n       -> if n + 1 = v then From (n,v) else And (c, Just v)
      | From (a, b)  -> if b + 1 = v then From (a, v) else And (c, Just v)
      | And (a, b)   -> And (a, add b v)

  end (* CompressIntSequencess *)

  (* TODO: Unify this with the Nomenclature suffix logic. *)
  let parse_last s =
    try Some (int_of_string s)
    with Failure _ -> None     (* for the ":5N" cases *)

  let to_comparable a =
    let l = Nomenclature.colon_split a in
    let n = List.length l in
    match List.split_n l (n - 1) with
    | [], _       -> `Special a         (* n = 0 *)
    | _,  []      -> eprintf "Odd allele length %d. Unexpected List.split_n \
                              behavior." n;
                     `Special a
    | d,  l :: [] ->
        begin
          match parse_last l with
          | None    -> `Special a
          | Some li -> `Template (d, li)
        end
    | _,  l :: _  -> eprintf "Asked to split at penultimate! %d. \
                              Unexpected List.split_n behavior." n;
                     `Special a


  let rec split_all =
    let comp = function
      | `Special s        -> `Special s
      | `Template (d, li) -> `Compress (d, [li])
    in
    function
    | []      -> []
    | h :: [] -> (comp h) :: []
    | h :: t  ->
        let rec loop cur acc = function
          | []      -> List.rev (cur :: acc)
          | h :: t  ->
              match h with
              | `Template (d, li) -> begin
                  match cur with
                  | `Compress (dc, ls) when d = dc
                                  -> loop (`Compress (d, li :: ls)) acc t
                  | `Special _
                  | `Compress _   -> loop (comp h) (cur :: acc) t
                  end
              | `Special _    -> begin
                  (* try to keep compressable targets current *)
                  match cur with
                  | `Special _  -> loop (comp h) (cur :: acc) t
                  | `Compress _ -> loop cur (comp h :: acc) t
                  end
        in
        loop (comp h) [] t

  let compress_int_list lst =
    let open CompressIntSequencess in
    match List.sort lst ~cmp:(fun (c1 : int) c2 -> Pervasives.compare c1 c2) with
    | []      -> ""  (* error instead? *)
    | h :: tl ->
        List.fold_left tl ~init:(Just h) ~f:add
        |> to_last

  let f lst =
    List.map lst ~f:to_comparable
    |> split_all
    |> fun clst ->
        let rec loop acc = function
          | []                      -> List.rev acc
          | `Special s :: tl        -> loop (s :: acc) tl
          | `Compress (t, l) :: tl  ->
              let ns = rejoin_colon (t @ [ compress_int_list l]) in
              loop (ns :: acc) tl
        in
        loop [] clst

end  (* CompressNames *)

module Set : sig

  type t

  val empty : t

  val equal : t -> t -> bool

  val compare : t -> t -> int

  val copy : t -> t

  val init : index -> t

  (** [set set allele] will make sure that [allele] is
      in [set], specifically [is_set set allele] will be [true].
      [set] is mutated. *)
  val set : t -> allele -> t

  (** [singleton allele] will create an edge set with just [allele]. *)
  val singleton : index -> allele -> t

  (** [clear set allele] will make sure that [allele] is not
      in [set], specifically [is_set set allele] will be
      [false]. [set] is mutated. *)
  val clear : t -> allele -> t

  (** [is_set set allele] is [allele] in [set]. *)
  val is_set : t -> allele -> bool

  val fold : f:('a -> allele -> 'a) -> init:'a -> t -> 'a

  val fold_set_indices : f:('a -> int -> 'a) -> init:'a -> t -> 'a

  val iter : f:(allele -> unit) -> t -> unit

  val iter_set_indices : f:(int -> unit) -> t -> unit

  (** [cardinal set] returns the number of alleles found in [set]. *)
  val cardinal : t -> int

  (** [union e1 e2] will return an edge set with all alleles found in
      [e1] and/or [e2]. *)
  val union : t -> t -> t

  val unite : into:t -> t -> unit

  (** [inter e1 e2] will return an edge set with alleles found in both
      [e1] and [e2]. *)
  val inter : t -> t -> t

  (** [diff e1 e2] will return an edge set with alleles found in [e1] but not
      in [e2]. *)
  val diff : t -> t -> t

  val inter_diff : t -> t -> t * t * bool * bool

  (** [complement set] returns a set of all the alleles not in [set].*)
  val complement : t -> t

  (** Construct a string of all the alleles found in the edge set. *)
  val to_string : ?compress:bool -> t -> string

  (** Construct a string with all of the alleles {b not} found in the edge set.

      @param prefix will determine if the complement string is prefixed
        (ie. "Complement of ") *)
  val complement_string : ?compress:bool -> ?prefix:string -> t -> string

  (** [to_human_readable] uses heuristics to pick a shorter string
      representation of the edge set.

      @param compress Use allele run compression
        (ex. ["A*01:01"; "A*01:02"; .. "A*01:05"] -> "A*01:01-05"), defaults to true.

      @param max_length Trim the string, defaults to first 500 characters.
      @param complement Whether to consider printing the complement string
        (defaults to `Yes):
        [`No]       -> Never print the complement string.
        [`Yes]      -> If the complement set has fewer alleles then print it with
                        "C. of" as a prefix
        [`Prefix s] -> If the complement set has fewer alleles then print it with
                        [s] as a prefix
      *)
  val to_human_readable : ?compress:bool -> ?max_length:int ->
    ?complement:[ `No | `Prefix of string | `Yes] -> t -> string

end (* Set *) = struct

  type t =
    { i : index
    ; b : Bitvector.t
    }
  [@@deriving eq, ord]

  let empty =
    { i = empty_index
    ; b = Bitvector.empty
    }

  let copy s =
    { s with b = Bitvector.copy s.b }

  let init index =
    { b = Bitvector.create ~size:index.size false
    ; i = index
    }

  let set s allele =
    Bitvector.set s.b (StringMap.find allele s.i.to_index);
    s

  let unite ~into from =
    Bitvector.union_into ~into:into.b from.b

  let singleton index allele =
    let s = init index in
    set s allele

  let clear s allele =
    Bitvector.reset s.b (StringMap.find allele s.i.to_index);
    s

  let is_set s allele =
    Bitvector.get s.b (StringMap.find allele s.i.to_index)

  let fold_set_indices ~f ~init s =
    let r = ref init in
    Bitvector.iter_true s.b (fun i -> r := f !r i);
    !r

  let fold ~f ~init s =
    fold_set_indices s ~init ~f:(fun a i -> f a s.i.to_allele.(i))

  let iter ~f s =
    fold ~init:() ~f:(fun () a -> f a) s

  let iter_set_indices ~f s =
    fold_set_indices ~init:() ~f:(fun () i -> f i) s

  let cardinal s = Bitvector.cardinal s.b

  let union s t =
    { s with b = Bitvector.union s.b t.b }

  let inter s t =
    { s with b = Bitvector.inter s.b t.b }

  let diff s t =
    { s with b = Bitvector.diff s.b t.b }

  let inter_diff s t =
    let bi, bd, same, nots = Bitvector.inter_diff s.b t.b in
    { s with b = bi}
    , { s with b = bd }
    , same
    , nots

  let complement s =
    { s with b = Bitvector.negate s.b }

  let to_string ?(compress=false) s =
    fold ~f:(fun a s -> s :: a) ~init:[] s
    |> List.rev
    |> (fun l -> if compress then CompressNames.f l else l)
    |> String.concat ~sep:";"

  let complement_string ?(compress=false) ?prefix s =
    Array.fold_left s.i.to_allele ~init:(0, [])
        ~f:(fun (i, acc) a -> if Bitvector.get s.b i
                              then (i + 1, acc)
                              else (i + 1, a :: acc))
    |> snd
    |> List.rev
    |> (fun l -> if compress then CompressNames.f l else l)
    |> String.concat ~sep:";"
    |> function
        | ""  -> "Nothing"
        | s   -> match prefix with
                 | None    -> s
                 | Some cp -> cp ^ s

  let to_human_readable ?(compress=true) ?(max_length=1000) ?(complement=`Yes) s =
    let make_shorter =
      let p = cardinal s in
      if p = s.i.size then
        "Everything"
      else if p = 0 then
        "Nothing"
      else if p <= s.i.size / 2 then
        to_string ~compress s
      else
        match complement with
        | `Yes           -> complement_string ~compress ~prefix:"C. of " s
        | `Prefix prefix -> complement_string ~compress ~prefix s
        | `No            -> to_string ~compress s
    in
    String.take make_shorter ~index:max_length

end (* Set *)

module Map : sig

  type 'a t

  val to_array : 'a t -> 'a array

  (** [make default_value]. *)
  val make : index -> 'a -> 'a t

  (** [init f]. *)
  val init : index -> (allele -> 'a) -> 'a t

  (** [get map allele]. *)
  val get : 'a t -> allele -> 'a

  (** [update2 source dest u] update the value of [dest] with the result of
      [u] and the values of the same allele of [source] and [dest]. *)
  val update2 : source:'a t -> dest:'b t -> ('a -> 'b -> 'b) -> unit

  (** [fold f init map] fold over all alleles found in the [map]. *)
  val fold : f:('a -> 'b -> allele -> 'a) -> init:'a -> 'b t -> 'a

  (** [iter f map] iter over all allele assignments in the [map]. *)
  val iter : f:('b -> allele -> unit) -> 'b t -> unit

  (** [map f cm] return a new map based on applying [f] to each element
       of [cm]. *)
  val map : f:('a -> allele -> 'b) -> 'a t -> 'b t

  (** [fold_wa f init map] fold over all values found in the [map],
      without regard for the allele key (wa = without allele). *)
  val fold_wa : f:('a -> 'b -> 'a) -> init:'a -> 'b t -> 'a

  (** [iter_wa f map] iter over all allele assignments in the [map],
      without regard for the allele key (wa = without allele). *)
  val iter_wa : f:('b -> unit) -> 'b t -> unit

  (** [map_wa f m] maps the values of the map [m] without regard to the allele
      key (wa = without allele) .*)
  val map_wa : f:('a -> 'b) -> 'a t -> 'b t

  (** [map2_wa f m1 m2] map from two maps. *)
  val map2_wa : f:('a -> 'b -> 'c) -> 'a t -> 'b t -> 'c t

  (** [values_assoc m] compress (invert) the values found in [m] into
      an association list.
  val values_assoc : 'a map -> ('a * Set.set) list
*)

  val update_from : Set.t -> f:('a -> 'a) -> 'a t -> unit

  (** [update_from_and_fold set f init map] updates and folds over the [set]
      elements of [map].*)
  val update_from_and_fold : Set.t -> f:('a -> 'b -> 'b * 'a) -> init:'a -> 'b t -> 'a

  (** [choose m] return a single element from the map. *)
  val choose : 'a t -> 'a

end (* Map *) = struct

  type 'a t =
    { i : index
    ; a : 'a array
    }

  let to_array m = m.a

  let make index e =
    { a = Array.make index.size e
    ; i = index
    }

  let init index f =
    { a = Array.init index.size ~f:(fun i -> f index.to_allele.(i))
    ; i = index
    }

  let get m a =
    m.a.(StringMap.find a m.i.to_index)

  let update2 ~source ~dest f =
    for i = 0 to Array.length source.a - 1 do
      dest.a.(i) <- f source.a.(i) dest.a.(i)
    done

  let map2_wa ~f m1 m2 =
    { m1 with a = Array.init m1.i.size ~f:(fun i -> f m1.a.(i) m2.a.(i))}

  let fold ~f ~init m =
    let s = ref init in
    for i = 0 to m.i.size - 1 do
      s := f !s m.a.(i) m.i.to_allele.(i)
    done;
    !s

  let iter ~f m =
    fold ~init:() ~f:(fun () v a -> f v a) m

  let map ~f m =
    { m with a = Array.mapi m.a ~f:(fun i v -> f v m.i.to_allele.(i)) }

  let fold_wa ~f ~init m  = Array.fold_left m.a ~f ~init

  let iter_wa ~f m = Array.iter m.a ~f

  let map_wa ~f m =
    { m with a = Array.map m.a ~f }

  let update_from s ~f m =
    Set.fold_set_indices s ~init:() ~f:(fun () i -> m.a.(i) <- f m.a.(i))

  let update_from_and_fold s ~f ~init m =
    Set.fold_set_indices s ~init ~f:(fun acc i ->
      let nv, nacc = f acc m.a.(i) in
      m.a.(i) <- nv;                    (* Use in 1st position ala fold_left. *)
      nacc)

  let choose m =
    m.a.(0)

end (* Map *)

(* Different logic of how we choose which alleles to include in the analysis. *)
module Selectors = struct

  (* TODO: The selector steps need to come before imputation! *)
  type t =
    | Regex of string
    | Specific of string
    | Number of int
    | DoNotIgnoreSuffixed
    [@@ deriving eq, ord, show ]
    (* When ppx_deriving >4.1 hits:
    [@@ deriving eq, ord, show { with_path = false }] *)

  let filenameable s =
    String.to_character_list s
    |> List.map ~f:(function
      | '*'  -> "s"
      | ':'  -> "c"        (* This is a weird MAC-Finder behavior.*)
      | '/'  -> "fs"
      | '\\' -> "bs"       (* Will we ever have windows users? *)
      | c    -> sprintf "%c" c)
    |> String.concat ~sep:""

  (* A bit more concise than show and easier for filenames.*)
  let to_string = function
    | Regex r             -> sprintf "R%s" (filenameable r)
    | Specific s          -> sprintf "S%s" (filenameable s)
    | Number n            -> sprintf "N%d" n
    | DoNotIgnoreSuffixed -> "DNIS"

  let string_of_list =
    string_of_list ~sep:"_" ~f:to_string

  let is_do_not_ignored_suffixed = function
    | DoNotIgnoreSuffixed -> true
    | Regex _
    | Specific _
    | Number _            -> false

  let matching_regex r alts =
    let p = Re_posix.compile_pat r in
    List.filter alts ~f:(fun a -> Re.execp p a.MSA.Parser.allele)

  let specific_allele s alts =
    List.filter alts ~f:(fun a -> a.MSA.Parser.allele = s)

  let remove_suffixed alts =
    List.fold_left alts ~init:[] ~f:(fun acc a ->
      match Nomenclature.trim_suffix a.MSA.Parser.allele with
      | Ok (_a, None)         -> a :: acc
      | Ok (_a, Some _suffix) -> acc                  (* Ignore suffixed! *)
      | Error m               -> failwith m)

  let apply_to_alt_elems ?(sort_to_nomenclature_order=true) lst =
    fun alts ->
      let suffixed, filters =
        List.partition lst ~f:is_do_not_ignored_suffixed
      in
      let nalts =
        match suffixed with
        | DoNotIgnoreSuffixed :: _ -> alts
        | _                        -> remove_suffixed alts
      in
      match filters with
      | []    -> nalts
      | apply ->
        List.fold_left apply ~init:[] ~f:(fun acc t ->
            match t with
            | Regex r    -> (matching_regex r nalts) @ acc
            | Specific s -> (specific_allele s nalts) @ acc
            | Number n   -> List.take nalts n @ acc
            | DoNotIgnoreSuffixed -> assert false        (* filtered! *))
      |> fun l ->
          if sort_to_nomenclature_order then
            MSA.Parser.sort_alts_by_nomenclature l
          else
            l

  let apply_to_mp lst mp =
    let open MSA.Parser in
    { mp with alt_elems = apply_to_alt_elems lst mp.alt_elems }

end (* Selectors. *)

(* Where do we get the allele information? *)
module Input = struct

  type t =
    | AlignmentFile of
        { path      : string              (* path to file (ex. ../alignments/A_nuc.txt) *)
        ; selectors : Selectors.t list
        ; distance  : Distances.logic option
        (* How to measure distances, if specified than we impute. *)
        }
    | MergeFromPrefix of
        { prefix_path : string                   (* path to prefix (ex ../alignments/A) *)
        ; selectors   : Selectors.t list
        ; drop_sv     : bool                                    (* drop splice variants *)
        ; distance    : Distances.logic     (* How to measure distances, always impute! *)
        }

  let alignment ?(selectors=[]) ?distance path =
    AlignmentFile { path; selectors; distance }

  let merge ?(selectors=[]) ?(drop_sv=true)
    ?(distance=Distances.WeightedPerSegment) prefix_path =
    MergeFromPrefix {prefix_path; selectors; distance; drop_sv}

  let is_imputed = function
    | AlignmentFile { distance; _ } ->
        begin match distance with
        | None    -> false
        | Some _  -> true
        end
    | MergeFromPrefix _ -> true

  let to_short_fname_prefix = function
    | AlignmentFile { path; _ } ->
        (Filename.chop_extension (Filename.basename path))
    | MergeFromPrefix { prefix_path; } ->
        (Filename.basename prefix_path)

  let to_string = function
    | AlignmentFile { path; selectors; distance } ->
        sprintf "%s_%s_%s"
          (Filename.chop_extension (Filename.basename path))
          (Selectors.string_of_list selectors)
          (match distance with
           | None -> "not-imputed"
           | Some dl -> Distances.show_logic dl)
    | MergeFromPrefix { prefix_path; selectors; drop_sv; distance } ->
        sprintf "Merged_%s_%s_%s_%s"
          (Filename.basename prefix_path)
          (Selectors.string_of_list selectors)
          (if drop_sv then "dropped-slice-variants" else "kept-slice-variants")
          (Distances.show_logic distance)

  module MergeConstruction = struct

    (* IMGT alignment format, (nor the XML) provide a robust way to detect splice
      variants and their respective start/stop positions with respect to the
      global alignmnet. In the past, I've hard-coded transformations that perform
      the appropriate adjustment. This solution isn't feasible in the long term,
      and is abandonded to the comments above. This association list
      (locus -> alleles) provides a quick way to just remove those alleles from
      the graph.

      To try creating the graph without them set ~drop_known_splice_variants to
      false *)
    let splice_variants =
      let open Nomenclature in
      [ A, [ "A*01:11N" ; "A*29:01:01:02N" ; "A*26:01:01:03N" ; "A*03:01:01:02N" ]
      ; B, [ "B*15:01:01:02N"; "B*44:02:01:02S"; "B*07:44N" ]
      ; C, [ "C*04:09N" ]
      ; E, []
      ; F, []
      ; G, []
      ; H, []
      ; J, []
      ; K, []
      ; L, []
      ; T, []
      ; V, []
      ; W, []
      ]

    let splice_variant_filter p =
      match List.Assoc.get p splice_variants with
      | None      -> eprintf "Warning: splice variant filter not implemented for %s\n"
                      (Nomenclature.show_locus p);
                     fun l -> l
      | Some set  -> List.filter ~f:(fun a -> not (List.mem a.MSA.Parser.allele ~set))

    let select_and_impute selectors dl mp =
      Alter_MSA.Impute.do_it dl (Selectors.apply_to_mp selectors mp)

    (* As of IMGT 3.30 HLA-P doesn't have a P_nuc.txt (all of the data is form
     * gDNA).  * Asking for a merge shouldn't fail because the file doesn't
     * exist.
     * Furthermore HLA-E has an inconsistency between the Exon sequences. We'll
     * use just gen.
     * *)
    let handle_special path selectors dl =
      let open MSA.Parser in
      select_and_impute selectors dl (from_file path)

    let do_it ?(drop_known_splice_variants=true) prefix selectors dl =
      let open MSA.Parser in
      if Filename.basename prefix = "P" then begin
        eprintf "HLA-P doesn't have a nuc, ignoring request.\n";
        handle_special (prefix ^ "_gen.txt") selectors dl
      end else if Filename.basename prefix = "E" then begin
        eprintf "HLA-E has an inconsistency between gen and nuc, ignoring request.\n";
        handle_special (prefix ^ "_gen.txt") selectors dl
      end else
        let gen_mp = from_file (prefix ^ "_gen.txt") in
        let nuc_mp = from_file (prefix ^ "_nuc.txt") in
        if gen_mp.reference <> nuc_mp.reference then
          error "References don't match %s vs %s" gen_mp.reference nuc_mp.reference
        else if gen_mp.align_date <> nuc_mp.align_date then
          error "Align dates don't match %s vs %s" gen_mp.align_date nuc_mp.align_date
        else
          (* Check if it is even worthwhile to merge, do the gen and nuc mp's
           * have the same alleles? If they do, just impute the gen. *)
          let an { allele; _ } = allele in
          let gen_ans = string_set_of_list (List.map gen_mp.alt_elems ~f:an) in
          let nuc_ans = string_set_of_list (List.map nuc_mp.alt_elems ~f:an) in
          if StringSet.equal gen_ans nuc_ans then
            Selectors.apply_to_mp selectors gen_mp
            |> Alter_MSA.Impute.do_it dl
          else
            let gen_mp, nuc_mp =
              if drop_known_splice_variants then begin
                let f = splice_variant_filter gen_mp.locus in
                { gen_mp with alt_elems = f gen_mp.alt_elems},
                { nuc_mp with alt_elems = f nuc_mp.alt_elems}
              end else
                gen_mp, nuc_mp
            in
            let gen_mp = Selectors.apply_to_mp selectors gen_mp in
            let nuc_mp = Selectors.apply_to_mp selectors nuc_mp in
            Alter_MSA.Merge.do_it ~gen_mp ~nuc_mp dl

  end (* MergeConstruction *)

  let construct = function
    | AlignmentFile { path; selectors; distance } ->
        let mp = MSA.Parser.from_file path in
        let smp = Selectors.apply_to_mp selectors mp in
        begin match distance with
        | None    -> Ok smp
        | Some dl -> Alter_MSA.Impute.do_it dl smp
        end
    | MergeFromPrefix
        { prefix_path; selectors; drop_sv; distance } ->
          MergeConstruction.do_it prefix_path selectors distance

end (* Input *)
