
open Util
module BitSet = Batteries.BitSet
module Enum = Batteries.Enum

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

let to_alleles { to_allele; _ } = Array.to_list to_allele

module Set = struct

  type t = BitSet.t

  let copy = BitSet.copy

  let init {size; _} =
    BitSet.create size

  let set { to_index; _ } s allele =
    BitSet.set s (StringMap.find allele to_index);
    s

  let unite ~into from = BitSet.unite into from

  let singleton index allele =
    let s = init index in
    set index s allele

  let clear { to_index; _ } s allele =
    BitSet.unset s (StringMap.find allele to_index);
    s

  let is_set { to_index; _ } s allele =
    BitSet.mem s (StringMap.find allele to_index)

  let fold { to_allele; } ~f ~init s =
    Enum.fold (fun a i -> f a to_allele.(i)) init (BitSet.enum s)

  let iter index ~f s =
    fold index ~init:() ~f:(fun () a -> f a) s

  let cardinal = BitSet.count

  let empty = BitSet.empty
  let compare = BitSet.compare
  let equals = BitSet.equal

  let union = BitSet.union

  let inter = BitSet.inter

  let diff = BitSet.diff

  let complement {size; _} t =
    let c = BitSet.copy t in
    for i = 0 to size - 1 do BitSet.toggle c i done;
    c

  let to_string ?(compress=false) index s =
    fold index ~f:(fun a s -> s :: a) ~init:[] s
    |> List.rev
    |> (fun l -> if compress then CompressNames.f l else l)
    |> String.concat ~sep:";"

  let complement_string ?(compress=false) ?prefix { to_allele; _} s =
    Array.fold_left to_allele ~init:(0, [])
        ~f:(fun (i, acc) a -> if BitSet.mem s i
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

  let to_human_readable ?(compress=true) ?(max_length=500) ?(complement=`Yes) t s =
    let make_shorter =
      if BitSet.count s = t.size then
        "Everything"
      else if BitSet.count s = 0 then
        "Nothing"
      else if BitSet.count s <= t.size / 2 then
        to_string ~compress t s
      else
        match complement with
        | `Yes           -> complement_string ~compress ~prefix:"C. of " t s
        | `Prefix prefix -> complement_string ~compress ~prefix t s
        | `No            -> to_string ~compress t s
    in
    String.take make_shorter ~index:max_length

end (* Set *)

module Map = struct

  type 'a t = 'a array

  let to_array a = a

  let make { size; _} e =
    Array.make size e

  let init { size; to_allele; _} f =
    Array.init size ~f:(fun i -> f to_allele.(i))

  let get { to_index; _} m a =
    m.(StringMap.find a to_index)

  let cardinal = Array.length

  let update2 ~source ~dest f =
    for i = 0 to Array.length source - 1 do
      dest.(i) <- f source.(i) dest.(i)
    done

  let map2_wa ~f m1 m2 =
    Array.init (Array.length m1) ~f:(fun i -> f m1.(i) m2.(i))

  let fold { to_allele; size; _} ~f ~init amap =
    let s = ref init in
    for i = 0 to size - 1 do
      s := f !s amap.(i) to_allele.(i)
    done;
    !s

  let iter i ~f amap =
    fold i ~init:() ~f:(fun () m a -> f m a) amap

  let map { to_allele; _} ~f amap =
    Array.mapi amap ~f:(fun i c -> f amap.(i) (to_allele.(i)))

  let fold_wa = Array.fold_left

  let iter_wa = Array.iter

  let map_wa ~f amap =
    Array.map amap ~f

  let values_assoc index amap =
    Array.fold_left amap ~init:(0, []) ~f:(fun (i, asc) v ->
      let a = index.to_allele.(i) in
      try
        let s, rest = remove_and_assoc v asc in
        (* TODO: make these modules recursive to allow maps to see inside sets *)
        (i + 1, (v, Set.set index s a) :: rest)
      with Not_found ->
        (i + 1, (v, Set.singleton index a) ::asc))
    |> snd

  let update_from s ~f m =
    Enum.iter (fun i -> m.(i) <- f m.(i))
      (BitSet.enum s)

  let update_from_and_fold s ~f ~init m =
    Enum.fold (fun acc i ->
      let nm, nacc = f acc m.(i) in (* Use in 1st position ala fold_left. *)
      m.(i) <- nm;
      nacc) init (BitSet.enum s)

  let choose m =
    m.(0)

end (* Map *)

(* Different logic of how we choose which alleles to include in the analysis. *)
module Selection = struct

  (* Defined in this order to so that to apply multiple selections, we can sort
     the list and apply Regex's first to generate possibilities and then use
     reducing selectors like Without and Number afterwards. *)
  type t =
    | Regex of string
    | Specific of string
    | Without of string
    | Number of int
    [@@ deriving eq, ord, show ]
    (* When ppx_deriving >4.1 hits:
    [@@ deriving eq, ord, show { with_path = false }] *)

  (* A bit more concise than show and easier for filenames.*)
  let to_string = function
    | Regex r     -> sprintf "R%s" (Digest.string r |> Digest.to_hex)
    | Specific s  -> sprintf "S%s" s
    | Without e   -> sprintf "W%s" e
    | Number n    -> sprintf "N%d" n

  let list_to_string l =
    String.concat ~sep:"_" (List.map l ~f:(to_string))

  let apply_to_assoc lst =
    let sorted = List.sort ~cmp:compare lst in
    fun assoc ->
      List.fold_left sorted ~init:assoc ~f:(fun acc -> function
        | Regex r    -> let p = Re_posix.compile_pat r in
                        List.filter acc ~f:(fun (allele, _) -> Re.execp p allele)
        | Specific s -> List.filter acc ~f:(fun (allele, _) -> allele = s)
        | Without e  -> List.filter acc ~f:(fun (allele, _) -> allele <> e)
        | Number n   -> List.take acc n)

end (* Selection. *)

(* Where do we get the allele information? *)
module Input = struct

  type t =
    | AlignmentFile
          of (string            (* path to file (ex. ../alignments/A_nuc.txt) *)
             * bool)                                               (* impute? *)
    | MergeFromPrefix
          of (string                   (* path to prefix (ex ../alignments/A) *)
             * Distances.logic                   (* how to measure distances. *)
             * bool)                                               (* impute? *)

  let imputed = function
    | AlignmentFile (_, i)      -> i
    | MergeFromPrefix (_, _, i) -> i

  let to_string = function
    | AlignmentFile (path, i)   -> sprintf "AF_%s_%b"
                                    (Filename.chop_extension
                                      (Filename.basename path))
                                      i
    | MergeFromPrefix (p, d, i) -> sprintf "MGD_%s_%s_%b"
                                    (Filename.basename p)
                                    (Distances.show_logic d)
                                    i

  let construct = function
    | AlignmentFile (file, impute) ->
        let mp = Mas_parser.from_file file in
        if impute then
          Merge_mas.naive_impute mp
        else
          Ok (mp, []) (* empty merge_map *)
    | MergeFromPrefix (prefix, distance_logic, impute) ->
        if impute then
          failwith "Imputation and merging not implemented!"
        else
          Merge_mas.do_it prefix distance_logic

end (* Input *)
