
open Util

module SMap = Map.Make (struct type t = string let compare = String.compare end)
module IMap = Map.Make (struct type t = int let compare = Int.compare end)

type ds =
  { size      : int
  ; to_index  : int SMap.t
  ; to_allele : string IMap.t
  }

let construct lst =
  List.sort ~cmp:compare lst
  |> List.fold_left ~init:(0, SMap.empty, IMap.empty)
      ~f:(fun (i, sm, im) a -> (i + 1, SMap.add a i sm, IMap.add i a im))
  |> fun (size, to_index, to_allele) ->
      { size; to_index; to_allele }

let empty {size; _} =
  BitSet.create (size - 1)

let allele_list_to_bitset ({ to_index; _ } as ds) allele_lst =
  let bt = empty ds in
  List.iter allele_lst
    ~f:(fun allele -> BitSet.set bt (SMap.find allele to_index))

let set_allele { to_index; _ } bt allele =
  BitSet.set bt (SMap.find allele to_index)

let singleton ds allele =
  let bt = empty ds in
  set_allele ds bt allele;
  bt

let is_set { to_index; _ } bt allele =
  BitSet.is_set bt (SMap.find allele to_index)

let union = BitSet.union 

let inter = BitSet.inter

let to_string { to_allele; _} bt =
  IMap.fold (fun i s a -> if BitSet.is_set bt i then s :: a else a)
    to_allele []
  |> String.concat ~sep:"; "


