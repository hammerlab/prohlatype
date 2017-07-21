
open Util

type allele = string [@@deriving eq, ord]
let compare = compare_allele
let equal = equal_allele

type index =
  { size      : int
  ; to_index  : int StringMap.t
  ; to_allele : string array
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

  let split_colon = String.split ~on:(`Character ':')
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
    let l = split_colon a in
    let n = List.length l in
    match List.split_n l (n - 1) with
    | [], _       -> `Special a         (* n = 0 *)
    | _,  []      -> invalid_argf "odd length %d" n
    | d,  l :: [] ->
        begin
          match parse_last l with
          | None    -> `Special a
          | Some li -> `Template (d, li)
        end
    | _,  l :: _  -> invalid_argf "Split at the end! %d" n

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

type set = Fixed_width.t

let set_compare = Pervasives.compare  (*?*)

let set_equals = (=) (*?*)

let set_empty = Fixed_width.empty

module type Index = sig
  val index : index
end

module type Set = sig

  val copy : set -> set

  val init : unit -> set

  (** [set set allele] will make sure that [allele] is
      in [set], specifically [is_set set allele] will be [true].
      [set] is mutated. *)
  val set : set -> allele -> set

  (** [singleton allele] will create an edge set with just [allele]. *)
  val singleton : allele -> set

  (** [clear set allele] will make sure that [allele] is not
      in [set], specifically [is_set set allele] will be
      [false]. [set] is mutated. *)
  val clear : set -> allele -> set

  (** [is_set set allele] is [allele] in [set]. *)
  val is_set : set -> allele -> bool

  val fold : f:('a -> allele -> 'a) -> init:'a -> set -> 'a

  val fold_set_indices : f:('a -> int -> 'a) -> init:'a -> set -> 'a

  val iter : f:(allele -> unit) -> set -> unit

  val iter_set_indices : f:(int -> unit) -> set -> unit

  (** [cardinal set] returns the number of alleles found in [set]. *)
  val cardinal : set -> int

  (** [union e1 e2] will return an edge set with all alleles found in
      [e1] and/or [e2]. *)
  val union : set -> set -> set

  val unite : into:set -> set -> unit

  (** [inter e1 e2] will return an edge set with alleles found in both
      [e1] and [e2]. *)
  val inter : set -> set -> set

  (** [diff e1 e2] will return an edge set with alleles found in [e1] but not
      in [e2]. *)
  val diff : set -> set -> set

  val inter_diff : set -> set -> set * set * bool * bool

  val split3 : set -> set -> set * set * set * bool * bool * bool

  (** [complement set] returns a set of all the alleles not in [set].*)
  val complement : set -> set

  (** Construct a string of all the alleles found in the edge set. *)
  val to_string : ?compress:bool -> set -> string

  (** Construct a string with all of the alleles {b not} found in the edge set.

      @param prefix will determine if the complement string is prefixed
        (ie. "Complement of ") *)
  val complement_string : ?compress:bool -> ?prefix:string -> set -> string

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
    ?complement:[ `No | `Prefix of string | `Yes] -> set -> string

end (* Set *)

module MakeSet (I: Index) : Set = struct

  module Fw = Fixed_width.Make(struct let size = I.index.size end)

  let copy = Fw.copy

  let init () = Fw.create false

  let set s allele =
    Fw.set s (StringMap.find allele I.index.to_index);
    s

  let unite ~into from = Fw.union_into ~into from

  let singleton allele =
    let s = init () in
    set s allele

  let clear s allele =
    Fw.reset s (StringMap.find allele I.index.to_index);
    s

  let is_set s allele =
    Fw.get s (StringMap.find allele I.index.to_index)

  let fold_set_indices ~f ~init s =
    let r = ref init in
    Fw.iter_true s (fun i -> r := f !r i);
    !r

  let fold ~f ~init s =
    fold_set_indices s ~init ~f:(fun a i -> f a I.index.to_allele.(i))

  let iter ~f s =
    fold ~init:() ~f:(fun () a -> f a) s

  let iter_set_indices ~f s =
    fold_set_indices ~init:() ~f:(fun () i -> f i) s

  let cardinal = Fw.cardinal

  let union = Fw.union

  let inter = Fw.inter

  let diff = Fw.diff

  let inter_diff = Fw.inter_diff

  let split3 = Fw.split3

  let complement = Fw.negate

  let to_string ?(compress=false) s =
    fold ~f:(fun a s -> s :: a) ~init:[] s
    |> List.rev
    |> (fun l -> if compress then CompressNames.f l else l)
    |> String.concat ~sep:";"

  let complement_string ?(compress=false) ?prefix s =
    Array.fold_left I.index.to_allele ~init:(0, [])
        ~f:(fun (i, acc) a -> if Fw.get s i
                              then (i + 1, acc)
                              else (i + 1, a :: acc))
    |> snd
    |> List.rev
    |> (fun l -> if compress then CompressNames.f l else l)
    |> String.concat ~sep:";"
    |> function
        | ""  -> invalid_argf "Complement of everything!"
        | s   -> match prefix with
                 | None    -> s
                 | Some cp -> cp ^ s

  let to_human_readable ?(compress=true) ?(max_length=1000) ?(complement=`Yes) s =
    let make_shorter =
      let p = cardinal s in
      if p = I.index.size then
        "Everything"
      else if p = 0 then
        "Nothing"
      else if p <= I.index.size / 2 then
        to_string ~compress s
      else
        match complement with
        | `Yes           -> complement_string ~compress ~prefix:"C. of " s
        | `Prefix prefix -> complement_string ~compress ~prefix s
        | `No            -> to_string ~compress s
    in
    String.take make_shorter ~index:max_length

end (* MakeSet *)

type 'a map = 'a array

module type Map = sig

  val to_array : 'a map -> 'a array

  (** [make default_value]. *)
  val make : 'a -> 'a map

  (** [init f]. *)
  val init : (allele -> 'a) -> 'a map

  (** [get map allele]. *)
  val get : 'a map -> allele -> 'a

  (** [update2 source dest u] update the value of [dest] with the result of
      [u] and the values of the same allele of [source] and [dest]. *)
  val update2 : source:'a map -> dest:'b map -> ('a -> 'b -> 'b) -> unit

  (** [fold f init map] fold over all alleles found in the [map]. *)
  val fold : f:('a -> 'b -> allele -> 'a) -> init:'a -> 'b map -> 'a

  (** [iter f map] iter over all allele assignments in the [map]. *)
  val iter : f:('b -> allele -> unit) -> 'b map -> unit

  (** [map f cm] return a new map based on applying [f] to each element
       of [cm]. *)
  val map : f:('a -> allele -> 'b) -> 'a map -> 'b map

  (** [fold_wa f init map] fold over all values found in the [map],
      without regard for the allele key (wa = without allele). *)
  val fold_wa : f:('a -> 'b -> 'a) -> init:'a -> 'b map -> 'a

  (** [iter_wa f map] iter over all allele assignments in the [map],
      without regard for the allele key (wa = without allele). *)
  val iter_wa : f:('b -> unit) -> 'b map -> unit

  (** [map_wa f m] maps the values of the map [m] without regard to the allele
      key (wa = without allele) .*)
  val map_wa : f:('a -> 'b) -> 'a map -> 'b map

  (** [map2_wa f m1 m2] map from two maps. *)
  val map2_wa : f:('a -> 'b -> 'c) -> 'a map -> 'b map -> 'c map

  (** [values_assoc m] compress (invert) the values found in [m] into
      an association list.
  val values_assoc : 'a map -> ('a * Set.set) list
*)

  val update_from : set -> f:('a -> 'a) -> 'a map -> unit

  (** [update_from_and_fold set f init map] updates and folds over the [set]
      elements of [map].*)
  val update_from_and_fold : set -> f:('a -> 'b -> 'b * 'a) -> init:'a -> 'b map -> 'a

  (** [choose m] return a single element from the map. *)
  val choose : 'a map -> 'a

end (* Map *)

module MakeMap (I: Index) : Map = struct

  module Fw = Fixed_width.Make(struct let size = I.index.size end)

  let to_array a = a

  let make e =
    Array.make I.index.size e

  let init f =
    Array.init I.index.size ~f:(fun i -> f I.index.to_allele.(i))

  let get m a =
    m.(StringMap.find a I.index.to_index)

  let update2 ~source ~dest f =
    for i = 0 to Array.length source - 1 do
      dest.(i) <- f source.(i) dest.(i)
    done

  let map2_wa ~f m1 m2 =
    Array.init (Array.length m1) ~f:(fun i -> f m1.(i) m2.(i))

  let fold ~f ~init amap =
    let s = ref init in
    for i = 0 to I.index.size - 1 do
      s := f !s amap.(i) I.index.to_allele.(i)
    done;
    !s

  let iter ~f amap =
    fold ~init:() ~f:(fun () m a -> f m a) amap

  let map ~f amap =
    Array.mapi amap ~f:(fun i c -> f amap.(i) (I.index.to_allele.(i)))

  let fold_wa = Array.fold_left

  let iter_wa = Array.iter

  let map_wa ~f amap =
    Array.map amap ~f

(*
  let values_assoc amap =
    Array.fold_left amap ~init:(0, []) ~f:(fun (i, asc) v ->
      let a = I.index.to_allele.(i) in
      try
        let s, rest = remove_and_assoc v asc in
        (* TODO: make these modules recursive to allow maps to see inside sets *)
        (i + 1, (v, Set.set s a) :: rest)
      with Not_found ->
        (i + 1, (v, Set.singleton a) ::asc))
    |> snd
    *)

  let update_from s ~f m =
    Fw.iter_true s (fun i -> m.(i) <- f m.(i))

  let update_from_and_fold s ~f ~init m =
    let r = ref init in
    Fw.iter_true s (fun i ->
      let nm, nacc = f !r m.(i) in (* Use in 1st position ala fold_left. *)
      r := nacc;
      m.(i) <- nm);
    !r

  let choose m =
    m.(0)

end (* MakeMap *)

(* Different logic of how we choose which alleles to include in the analysis. *)
module Selectors = struct

  (* Defined in this order to so that to apply multiple selections, we can sort
     the list and apply Regex's first to generate possibilities and then use
     reducing selectors like Without and Number afterwards. *)
  (* TODO: The selector steps need to come before imputation! *)
  type t =
    | Regex of string
    | Specific of string
    | Without of string
    | Number of int
    | DoNotIgnoreSuffixed
    [@@ deriving eq, ord, show ]
    (* When ppx_deriving >4.1 hits:
    [@@ deriving eq, ord, show { with_path = false }] *)

  (* A bit more concise than show and easier for filenames.*)
  let to_string = function
    | Regex r             -> sprintf "R%s" (Digest.string r |> Digest.to_hex)
    | Specific s          -> sprintf "S%s" s
    | Without e           -> sprintf "W%s" e
    | Number n            -> sprintf "N%d" n
    | DoNotIgnoreSuffixed -> "DNIS"

  let list_to_string l =
    String.concat ~sep:"_" (List.map l ~f:(to_string))

  let sort_by_nomenclature lst =
    List.map lst ~f:(fun (a, s) -> Nomenclature.parse_to_resolution_exn a, a, s)
    |> List.sort ~cmp:(fun (n1, _,_) (n2,_,_) -> Nomenclature.compare n1 n2)
    |> List.map ~f:(fun (_n, a, s) -> (a, s))

  let apply_to_assoc ?(sort_to_nomenclature_order=true) lst =
    let sorted = List.sort ~cmp:compare lst in
    fun assoc ->
      let assoc =
        if List.mem DoNotIgnoreSuffixed ~set:sorted then
          assoc
        else
          List.fold_left assoc ~init:[] ~f:(fun acc (allele, v) ->
            match Nomenclature.trim_suffix allele with
            | Ok (_a, None)         -> (allele, v) :: acc
            | Ok (_a, Some _suffix) -> acc                  (* Ignore suffixed! *)
            | Error m               -> failwith m)
          (*|> List.rev *)
      in
      List.fold_left sorted ~init:assoc ~f:(fun acc -> function
        | Regex r             -> let p = Re_posix.compile_pat r in
                                 List.filter acc ~f:(fun (allele, _) -> Re.execp p allele)
        | Specific s          -> List.filter acc ~f:(fun (allele, _) -> allele = s)
        | Without e           -> List.filter acc ~f:(fun (allele, _) -> allele <> e)
        | Number n            -> List.take acc n
        | DoNotIgnoreSuffixed -> acc        (* No-op at this point *))
      |> fun l ->
          if sort_to_nomenclature_order then sort_by_nomenclature l else l

  let apply_to_mp lst mp =
    let open MSA.Parser in
    { mp with alt_elems = apply_to_assoc lst mp.alt_elems }

end (* Selectors. *)

(* Where do we get the allele information? *)
module Input = struct

  type t =
    | AlignmentFile of
        { path      : string    (* path to file (ex. ../alignments/A_nuc.txt) *)
        ; selectors : Selectors.t list
        ; impute    : bool                                         (* impute? *)
        }
    | MergeFromPrefix of
        { prefix_path : string         (* path to prefix (ex ../alignments/A) *)
        ; selectors   : Selectors.t list
        ; drop_sv     : bool                          (* drop splice variants *)
        ; distance    : Distances.logic          (* how to measure distances. *)
        ; impute      : bool                                       (* impute? *)
        }

  let alignment ?(selectors=[]) ~impute path =
    AlignmentFile { path; impute ; selectors }

  let merge ?(selectors=[]) ?(drop_sv=true) ~distance ~impute prefix_path =
    MergeFromPrefix {prefix_path; distance; impute; selectors; drop_sv}

  let is_imputed = function
    | AlignmentFile { impute; _ }
    | MergeFromPrefix { impute; _} -> impute

  let to_string = function
    | AlignmentFile { path; selectors; impute } ->
        sprintf "AF_%s_%s_%b"
          (Filename.chop_extension (Filename.basename path))
          (Selectors.list_to_string selectors)
          impute
    | MergeFromPrefix { prefix_path; selectors; drop_sv; distance; impute} ->
        sprintf "MGD_%s_%s_%b_%s_%b"
          (Filename.basename prefix_path)
          (Selectors.list_to_string selectors)
          drop_sv
          (Distances.show_logic distance)
          impute

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
      [ "A", [ "A*01:11N" ; "A*29:01:01:02N" ; "A*26:01:01:03N" ; "A*03:01:01:02N" ]
      ; "B", [ "B*15:01:01:02N"; "B*44:02:01:02S"; "B*07:44N" ]
      ; "C", [ "C*04:09N" ]
      ]

    let splice_variant_filter p =
      match List.Assoc.get p splice_variants with
      | None      -> eprintf "Splice variant filter not implemented for %s\n" p;
                     fun l -> l
      | Some set  -> List.filter ~f:(fun (a, _) -> not (List.mem a ~set))

    let do_it ?verbose ?(drop_known_splice_variants=true) prefix selectors dl =
      let open MSA.Parser in
      let gen_mp = from_file (prefix ^ "_gen.txt") in
      let nuc_mp = from_file (prefix ^ "_nuc.txt") in
      if gen_mp.reference <> nuc_mp.reference then
        error "References don't match %s vs %s" gen_mp.reference nuc_mp.reference
      else if gen_mp.align_date <> nuc_mp.align_date then
        error "Align dates don't match %s vs %s" gen_mp.align_date nuc_mp.align_date
      else
        let gen_mp, nuc_mp =
          if drop_known_splice_variants then begin
            let f = splice_variant_filter (Filename.basename prefix) in
            { gen_mp with alt_elems = f gen_mp.alt_elems},
            { nuc_mp with alt_elems = f nuc_mp.alt_elems}
          end else
            gen_mp, nuc_mp
        in
        let gen_mp = Selectors.apply_to_mp selectors gen_mp in
        let nuc_mp = Selectors.apply_to_mp selectors nuc_mp in
        Alter_MSA.Merge.do_it_spec ?verbose ~gen_mp ~nuc_mp dl

  end (* MergeConstruction *)

  let construct = function
    | AlignmentFile { path; selectors; impute} ->
        let mp = MSA.Parser.from_file path in
        let smp = Selectors.apply_to_mp selectors mp in
        if impute then
          Alter_MSA.Impute.do_it smp
        else
          Ok (smp, []) (* empty merge_map *)
    | MergeFromPrefix
        { prefix_path; selectors; drop_sv; distance; impute } ->
          if impute then
            failwith "Imputation and merging not implemented!"
          else
            MergeConstruction.do_it prefix_path selectors distance

end (* Input *)
