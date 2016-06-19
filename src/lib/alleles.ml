
open Util

module SMap = Map.Make (struct
    type t = string
    let compare = String.compare
  end)

type index =
  { size      : int
  ; to_index  : int SMap.t
  ; to_allele : string array
  }

let index lst =
  let to_allele = Array.of_list (List.sort ~cmp:compare lst) in
  let size, to_index =
    Array.fold_left to_allele ~init:(0, SMap.empty )
      ~f:(fun (i, sm) a -> (i + 1, SMap.add a i sm))
  in
  { size; to_index; to_allele }

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

  let empty = BitSet.empty
  let compare = BitSet.compare
  let equals = BitSet.equals

  let union = BitSet.union

  let inter = BitSet.inter

  let to_string index s =
    fold index ~f:(fun a s -> s :: a) ~init:[] s
    |> List.rev
    |> String.concat ~sep:" "

  let complement_string ?(annotate=true) { to_allele; _} s =
    Array.fold_left to_allele ~init:(0, [])
        ~f:(fun (i, acc) a -> if BitSet.is_set s i
                              then (i + 1, acc)
                              else (i + 1, a :: acc))
    |> snd
    |> List.rev
    |> String.concat ~sep:" "
    |> (fun s -> if annotate then s ^ "Complement of " else s)

  let to_human_readable t s =
    if BitSet.count s > t.size / 2 then
      complement_string t s
        else
      to_string t s

end

module Map = struct

  type 'a t = 'a array

  let make { size; _} e =
    Array.make size e

  (*let init { size; _} s p =
    Array.init size ~f:(fun i -> p (BitSet.is_set s i))

  let map { to_allele; _} s f =
    Array.mapi to_allele ~f:(fun i a -> f (BitSet.is_set s i) a)*)

  let update_from s m f =
    Enum.iter (fun i -> m.(i) <- f m.(i)) (BitSet.enum s)

  let fold { to_allele; size; _} ~f ~init amap =
    let s = ref init in
    for i = 0 to size - 1 do
      s := f !s amap.(i) to_allele.(i)
    done;
    !s

end
