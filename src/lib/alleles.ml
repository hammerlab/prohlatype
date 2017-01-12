
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

  let init {size; _} =
    BitSet.create (size - 1)

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

module MSet = struct

  (* The order of the indices should be sorted! *)
  type t = int list list

  let of_set s =
    BitSet.enum s
    |> Enum.fold (fun acc i -> i :: acc) []
    |> List.rev
    |> fun a -> a :: []

  let complement {size; _} lst =
    List.map lst ~f:BitSet.of_list
    |> List.reduce  ~f:BitSet.union
    |> begin function
        | None   -> []
        | Some b ->
            let d = size - 1 in
            let rec loop acc i =
              if i > d then
                [ List.rev acc ]
              else if BitSet.mem b i then
                loop acc i
              else
                loop (i :: acc) (i + 1)
            in
            loop [] 0
        end

  let rec cardinal = function
    | []     -> 0
    | l :: t -> List.length l + cardinal t

  let union_separate l1 l2 =
    l1 @ l2

  let to_set { size; _} l =
    let s = BitSet.create (size - 1) in
    List.iter l ~f:(List.iter ~f:(BitSet.set s));
    s

  let to_human_readable ?compress ?max_length ?complement i t =
    let s = to_set i t in
    Set.to_human_readable ?compress ?max_length ?complement i s

end (* MSet *)

module Map = struct

  type 'a t = 'a array

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

  let update_from_and_fold s ~f ~init m =
    Enum.fold (fun acc i ->
      let nm, nacc = f acc m.(i) in (* Use in 1st position ala fold_left. *)
      m.(i) <- nm;
      nacc) init (BitSet.enum s)

  let update_from_and_fold_mset s ~f ~init m =
    List.fold_left s ~init ~f:(fun acc a ->
      List.fold_left a ~init:acc ~f:(fun acc i ->
        let nm, nacc = f acc m.(i) in (* Use in 1st position ala fold_left. *)
        m.(i) <- nm;
        nacc))

end (* Map *)
