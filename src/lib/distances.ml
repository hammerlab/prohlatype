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
        parse ta >>= fun (locus, (allele_resolution, suffix_opt)) ->
          let closest_allele_res = Trie.nearest allele_resolution trie in
          let closest_allele_str = to_string ~locus closest_allele_res in
          Ok (StringMap.add ~key:ta ~data:[(closest_allele_str, 1.)] m))

end (* Trie_distances *)

module Weighted_per_segment = struct

  let debug = ref false

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

  let one ~reference ~reference_sequence ~candidates ~allele ~allele_name =
    let open MSA.Segments in
    distances ~reference:reference_sequence ~allele >>= fun dist_to_ref ->
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
      StringMap.bindings candidates
      |> list_fold_ok ~init:[] ~f:(fun acc (al2, allele2) ->
          if !debug then
            printf "Calculating weighted differences for %s vs %s\n" allele_name al2;
          distances_between ~reference:reference_sequence ~allele1:allele ~allele2
          >>| List.map ~f:(fun s ->
                match s.relationship with
                | (Full _), (Full _) ->
                    Some (float s.mismatches, float s.seq_length)
                | _                  ->
                    None)
          >>= fun dlst ->
            match apply_mask ~init:dist_init ~f:dist_f ref_mask dlst with
            | None      -> Ok acc
            | Some dist -> Ok ((al2, dist) :: acc))
      >>= fun all_distances ->
        let with_reference = (reference, ref_diff) :: all_distances in
        Ok (List.sort with_reference
                ~cmp:(fun (_,d1) (_,(d2:float)) -> compare d1 d2))

  let f ~reference ~reference_sequence ~targets ~candidates =
    let c = one ~reference ~reference_sequence ~candidates in
    StringMap.bindings targets
    |> list_fold_ok ~init:[] ~f:(fun acc (allele_name, allele) ->
      c ~allele_name ~allele >>= fun d -> Ok ((allele_name, d) :: acc))
    >>| string_map_of_assoc

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
      let allele_name, aseq = allele in
      Weighted_per_segment.one ~reference ~reference_sequence ~candidates
            ~allele_name ~allele:aseq

type arg =
  { reference : string
  ; reference_sequence : alignment_sequence
  ; targets : alignment_sequence StringMap.t
  ; candidates : alignment_sequence StringMap.t
  }

let compute { reference; reference_sequence; targets; candidates } = function
  | Reference ->
      Ok (Reference.f ~reference ~reference_sequence ~targets ~candidates)
  | Trie               ->
      Trie_distances.f ~targets ~candidates
  | WeightedPerSegment ->
      Weighted_per_segment.f ~reference ~reference_sequence ~targets ~candidates
