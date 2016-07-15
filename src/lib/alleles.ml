
open Util

type allele = string [@@deriving eq, ord]
let compare = compare_allele
let equal = equal_allele
let to_string a = a

module SMap = Map.Make (struct
    type t = string
    let compare = compare_allele
  end)

type index =
  { size      : int
  ; to_index  : int SMap.t
  ; to_allele : string array
  }

let index lst =
  let to_allele = Array.of_list (List.sort ~cmp:compare_allele lst) in
  let size, to_index =
    Array.fold_left to_allele ~init:(0, SMap.empty )
      ~f:(fun (i, sm) a -> (i + 1, SMap.add a i sm))
  in
  { size; to_index; to_allele }

module CompressNames = struct

  let split_semi = String.split ~on:(`Character ':')

  type last_state =
    | And of (last_state * last_state)
    | From of (int * int)
    | Just of int

  let rec to_last = function
    | And  (a, b) -> sprintf "%s,%s" (to_last a) (to_last b)
    | From (a, b) -> sprintf "%d-%d" a b
    | Just n      -> sprintf "%d" n

  let rec last_val = function
    | And (_, v)  -> last_val v
    | From (_, l) -> l
    | Just n      -> n

  let rec add c v =
    match c with
    | Just n       -> if n + 1 = v then From (n,v) else And (c, Just v)
    | From (a, b)  -> if b + 1 = v then From (a, v) else And (c, Just v)
    | And (a, b)   -> And (a, add b v)

  let parse_last s =
    try Some (int_of_string s)
    with Failure _ -> None     (* for the ":5N" cases *)

  let rec f = function
    | []      -> []
    | h :: [] -> [h]
    | h :: t  ->
        let to_string ts lint =
          String.concat ~sep:":" (ts @ [to_last lint])
        in
        let rec loop ts lint acc = function
          | []      -> List.rev ((to_string ts lint) :: acc)
          | h :: t  ->
              let ss = split_semi h in
              let n  = List.length ss in
              let reset () =
                  let current = List.rev ((to_string ts lint) :: acc) in
                  let rest = h :: t in
                  current @ f rest
              in
              match List.split_n ss (n - 1) with
              | _,   []        -> reset ()
              | nts, last :: _ ->
                  match parse_last last with
                  | None          -> reset ()
                  | Some last_int ->
                      if nts = ts then
                        loop ts (add lint last_int) acc t
                      else
                        loop nts (Just last_int) (to_string ts lint :: acc) t
        in
        let ss = split_semi h in
        let n  = List.length ss in
        match List.split_n ss (n - 1) with
        | _,  []        -> h :: (f t)
        | ts, last :: _ ->
            match parse_last last with
            | None    -> h :: (f t)
            | Some l  -> loop ts (Just l) [] t

end  (* CompressNames *)

module Set = struct

  type t = BitSet.t

  let init {size; _} =
    BitSet.create (size - 1)

  let set { to_index; _ } s allele =
    BitSet.set s (SMap.find allele to_index)

  let singleton index allele =
    let s = init index in
    set index s allele;
    s

  let clear { to_index; _ } s allele =
    BitSet.unset s (SMap.find allele to_index)

  let is_set { to_index; _ } s allele =
    BitSet.is_set s (SMap.find allele to_index)

  let fold { to_allele; } ~f ~init s =
    Enum.fold (fun i a -> f a to_allele.(i)) init (BitSet.enum s)

  let iter index ~f s =
    fold index ~init:() ~f:(fun () a -> f a) s

  let empty = BitSet.empty
  let compare = BitSet.compare
  let equals = BitSet.equals

  let union = BitSet.union

  let inter = BitSet.inter

  let complement {size; _} t =
    let c = BitSet.copy t in
    for i = 0 to size - 1 do BitSet.toggle c i done;
    c

  let is_empty t =
    BitSet.count t = 0

  let any t =
    BitSet.count t > 0

  let all {size; _} t =
    BitSet.count t = size

  let to_string ?(compress=false) index s =
    fold index ~f:(fun a s -> s :: a) ~init:[] s
    |> List.rev
    |> (fun l -> if compress then CompressNames.f l else l)
    |> String.concat ~sep:" "

  let complement_string ?(compress=false) ?complement_prefix { to_allele; _} s =
    Array.fold_left to_allele ~init:(0, [])
        ~f:(fun (i, acc) a -> if BitSet.is_set s i
                              then (i + 1, acc)
                              else (i + 1, a :: acc))
    |> snd
    |> List.rev
    |> (fun l -> if compress then CompressNames.f l else l)
    |> String.concat ~sep:" "
    |> function
        | ""  -> invalid_argf "Complement of everything!"
        | s   -> match complement_prefix with
                 | None    -> s
                 | Some cp -> cp ^ s

  let to_human_readable ?compress ?(max_length=500) ?complement_prefix t s =
    let make_shorter =
      if BitSet.count s = t.size then
        "Everything"
      else if BitSet.count s > t.size / 2 then
        let complement_prefix =
          Option.value complement_prefix ~default:"C. of"
        in
        complement_string ?compress ~complement_prefix t s
      else
        to_string ?compress t s
    in
    String.take make_shorter ~index:max_length

end (* Set *)

module Map = struct

  type 'a t = 'a array

  let make { size; _} e =
    Array.make size e

  let init { size; to_allele; _} f =
    Array.init size ~f:(fun i -> f to_allele.(i))

  let get { to_index; _} m a =
    m.(SMap.find a to_index)

  (* let map { to_allele; _} s f =
    Array.mapi to_allele ~f:(fun i a -> f (BitSet.is_set s i) a)*)

  let update_spec { to_index; _} m a f =
    let j = SMap.find a to_index in
    m.(j) <- f m.(j)

  let update_all s m f =
    for j = 0 to Array.length m - 1 do
      m.(j) <- f (BitSet.is_set s j) m.(j)
    done

  let update_from s m f =
    Enum.iter (fun i -> m.(i) <- f m.(i)) (BitSet.enum s)

  let update2 m1 m2 f =
    for i = 0 to Array.length m1 - 1 do
      m2.(i) <- f m1.(i) m2.(i)
    done

  let fold { to_allele; size; _} ~f ~init amap =
    let s = ref init in
    for i = 0 to size - 1 do
      s := f !s amap.(i) to_allele.(i)
    done;
    !s

end (* Map *)
