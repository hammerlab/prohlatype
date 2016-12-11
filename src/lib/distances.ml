(** Measure distances between different alleles from (at the moment) parsed
    multiple alignment files. *)

open Util
open MoreLabels

module Amap = Map.Make (struct
  type t = string
  let compare = compare
end)

module TrieDistances = struct

  let init_trie elems =
    let open Nomenclature in
    list_fold_ok elems ~init:Trie.empty ~f:(fun trie s ->
      parse s >>= fun (_gene, (allele_resolution, suffix_opt)) ->
        Ok (Trie.add allele_resolution suffix_opt trie))

  let f targets candidates =
    let open Nomenclature in
    let just_candidates = Amap.bindings candidates |> List.map ~f:fst in
    init_trie just_candidates >>= fun trie ->
      Amap.bindings targets |> list_fold_ok ~init:Amap.empty ~f:(fun m (ta, _) ->
        parse ta >>= fun (gene, (allele_resolution, suffix_opt)) ->
          let closest_allele_res = Trie.nearest allele_resolution trie in
          let closest_allele_str = resolution_and_suffix_opt_to_string ~gene closest_allele_res in
          Ok (Amap.add ~key:ta ~data:[closest_allele_str] m))

end

let compute targets candidates = function
  | `Trie -> TrieDistances.f targets candidates
