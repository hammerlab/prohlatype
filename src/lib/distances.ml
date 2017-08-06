(** Measure distances between different alleles from (at the moment) parsed
    Multiple Sequence Alignment files. *)

open Util

module Trie_distances = struct

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

end (* Trie_distances *)

module Weighted_per_segment = struct

  let against_mask ~init ~f =
    List.fold_left ~init ~f:(fun a -> function
      | None -> a
      | Some (mismatches, ref_len) -> f a ~mismatches ~ref_len)

  let apply_mask ~init ~f =
    let open Option in
    List.fold_left2 ~init:(Some init) ~f:(fun aopt b c ->
      aopt >>= fun a ->
        match b, c with
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

  let one ~reference ~reference_sequence ~candidates ~allele =
    let open MSA.Segments in
    let dist_to_ref = distances ~reference:reference_sequence ~allele in
    let ref_mask =
      List.map dist_to_ref ~f:(fun s ->
        match s.relationship with
        | Full _ -> Some (float s.mismatches, float s.seq_length)
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
          distances_between ~reference:reference_sequence ~allele1:allele ~allele2
          |> List.map ~f:(fun s ->
              match s.relationship with
              | (Full _), (Full _) ->
                  Some (float s.mismatches, float s.seq_length)
              | _                  ->
                  None)
        in
        match apply_mask ~init:dist_init ~f:dist_f ref_mask dlst with
        | None      -> acc
        | Some dist -> (al2, dist) :: acc)
    in
    let with_reference = (reference, ref_diff) :: all_distances in
    List.sort with_reference ~cmp:(fun (_,d1) (_,(d2:float)) -> compare d1 d2)

  let f ~reference ~reference_sequence ~targets ~candidates =
    let c = one ~reference ~reference_sequence ~candidates in
    StringMap.mapi targets ~f:(fun _al1 allele -> c allele)

end (* Weighted_per_segment *)

module Reference = struct

  let one ~reference ~reference_sequence ~candidates ~allele =
    let is_ref, isn't =
      StringMap.bindings candidates
      |> List.partition ~f:(fun (al, _seq) -> al = reference)
    in
    List.map is_ref ~f:(fun _ -> reference, 0.0)
    @ List.map isn't ~f:(fun (a, _) -> a, infinity)

  let f ~reference ~reference_sequence ~targets ~candidates =
    StringMap.map targets ~f:(fun _s ->
      one ~reference ~reference_sequence ~candidates
        ~allele:("Allele sequence doesn't matter", []))

end (* Reference *)

type logic =
  | Reference
  | Trie
  | WeightedPerSegment
  [@@deriving show]

type alignment_sequence = string MSA.alignment_sequence

let one ~reference ~reference_sequence ~allele ~candidates = function
  | Reference          ->
      let _aname, aseq = allele in
      Ok (Reference.one ~reference ~reference_sequence ~candidates
            ~allele:aseq)
  | Trie               ->
      let aname, aseq = allele in
      let targets = StringMap.singleton aname aseq in
      Trie_distances.f ~targets ~candidates >>= fun m ->
        Ok (StringMap.find aname m)
  | WeightedPerSegment ->
      let _aname, aseq = allele in
      Ok (Weighted_per_segment.one ~reference ~reference_sequence ~candidates
            ~allele:aseq)

let compute ~reference ~reference_sequence ~targets ~candidates = function
  | Reference ->
      Ok (Reference.f ~reference ~reference_sequence ~targets ~candidates)
  | Trie               ->
      Trie_distances.f ~targets ~candidates
  | WeightedPerSegment ->
      Ok (Weighted_per_segment.f ~reference ~reference_sequence ~targets ~candidates)
