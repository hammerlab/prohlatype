(** Measure distances between different alleles from (at the moment) parsed
    multiple alignment files. *)

open Util

module TrieDistances = struct

  let init_trie elems =
    let open Nomenclature in
    list_fold_ok elems ~init:Trie.empty ~f:(fun trie s ->
      parse s >>= fun (_gene, (allele_resolution, suffix_opt)) ->
        Ok (Trie.add allele_resolution suffix_opt trie))

  let f ~targets ~candidates =
    let open Nomenclature in
    let just_candidates = StringMap.bindings candidates |> List.map ~f:fst in
    init_trie just_candidates >>= fun trie ->
      StringMap.bindings targets |> list_fold_ok ~init:StringMap.empty ~f:(fun m (ta, _) ->
        parse ta >>= fun (gene, (allele_resolution, suffix_opt)) ->
          let closest_allele_res = Trie.nearest allele_resolution trie in
          let closest_allele_str = resolution_and_suffix_opt_to_string ~gene closest_allele_res in
          Ok (StringMap.add ~key:ta ~data:[(closest_allele_str, 1.)] m))

end

module AverageExon = struct

  let against_mask ~init ~f =
    List.fold_left ~init ~f:(fun a o ->
      match o with | None -> a | Some (mismatches, ref_len) ->
        f a ~mismatches ~ref_len)

  let apply_mask ~init ~f =
    let open Option in
    List.fold_left2 ~init:(Some init) ~f:(fun aopt b c ->
      aopt >>= fun a ->
        match b,c  with
        | None,   None   -> Some a
        | None,   _
        | _,      None   -> None
        | Some _, Some (mismatches, ref_len) ->
            Some (f a ~mismatches ~ref_len))

  let dist ~normalize tlen =
    0.,
    (fun a ~mismatches ~ref_len ->
      let w = ref_len /. tlen in
      if normalize then
        a +. mismatches *. w
      else
        a +. mismatches)

  let f ref_allele reference ~targets ~candidates =
    let open Mas_parser in
    StringMap.mapi targets ~f:(fun al1 allele1 ->
      let dist_to_ref = allele_distances ~reference ~allele:allele1 in
      let ref_mask =
        List.map dist_to_ref ~f:(fun s ->
          match s.allele_relationships with
          | Full _ -> Some (float s.mismatches, float s.reference_seq_length)
          | _      -> None)
      in
      let tlen = against_mask ref_mask ~init:0.
          ~f:(fun a ~mismatches ~ref_len -> a +. ref_len)
      in
      let dist_init, dist_f = dist ~normalize:true tlen in
      let ref_diff = against_mask ~init:dist_init ~f:dist_f ref_mask in
      let all_distances =
        StringMap.fold candidates ~init:[] ~f:(fun ~key:al2 ~data:allele2 acc ->
          let dlst =
            Mas_parser.allele_distances_between
              ~reference ~allele1 ~allele2
            |> List.map ~f:(fun s ->
                match s.allele_relationships with
                | (Full _), (Full _) ->
                    Some (float s.mismatches, float s.reference_seq_length)
                | _                  ->
                    None)
          in
          match apply_mask ~init:dist_init ~f:dist_f ref_mask dlst with
          | None      -> acc
          | Some dist -> (al2, dist) :: acc)
      in
      let with_reference = (ref_allele, ref_diff) :: all_distances in
      List.sort with_reference ~cmp:(fun (_,d1) (_,(d2:float)) -> compare d1 d2))

end

type logic =
  | Trie
  | AverageExon
  [@@deriving show]

let compute ref_allele reference ~targets ~candidates = function
  | Trie        -> TrieDistances.f ~targets ~candidates
  | AverageExon -> Ok (AverageExon.f ref_allele reference ~targets ~candidates)
