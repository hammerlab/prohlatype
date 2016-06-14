
open Util

module SMap = Map.Make (struct type t = string let compare = String.compare end)
module IMap = Map.Make (struct type t = int let compare = Int.compare end)

type ds =
  { size      : int
  ; to_index  : int SMap.t
  ; to_allele : string array
  }

let construct lst =
  let to_allele = Array.of_list (List.sort ~cmp:compare lst) in
  let size, to_index =
    Array.fold_left to_allele ~init:(0, SMap.empty )
      ~f:(fun (i, sm) a -> (i + 1, SMap.add a i sm))
  in
  { size; to_index; to_allele }

let empty {size; _} =
  BitSet.create (size - 1)

let allele_list_to_bitset ds allele_lst =
  let bt = empty ds in
  List.iter allele_lst
    ~f:(fun allele -> BitSet.set bt (SMap.find allele ds.to_index))

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
  BitSet.enum bt
  |> Enum.fold (fun i a -> to_allele.(i) :: a) []
  |> List.rev
  |> String.concat ~sep:" "

let complement_string { to_allele; _} bt =
  Array.fold_left to_allele ~init:(0, [])
      ~f:(fun (i, acc) a -> if BitSet.is_set bt i
                            then (i + 1, acc)
                            else (i + 1, a :: acc))
  |> snd
  |> List.rev
  |> String.concat ~sep:" "
  |> ( ^ ) "Complement of "

let to_human_readable t bt =
  let first_try =
    if BitSet.count bt > t.size / 2 then
      complement_string t bt
        else
      to_string t bt
  in
  if String.length first_try < 100 then
    first_try
  else
    first_try
      (*
    Printf.sprintf "%d" (Hashtbl.hash first_try)
         *)

    (*
  IMap.fold (fun i s a -> if BitSet.is_set bt i then s :: a else a)
    to_allele []
  |> String.concat ~sep:"; " *)


