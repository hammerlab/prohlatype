(** Compute distances between MHC alleles. *)
open Prohlatype

let app_name = "allele_distances"

let to_distance_targets_and_candidates ?alignment ?merge () =
  let open MSA.Parser in
  match merge, alignment with
  | (Some prefix), _  ->
      let gen = from_file (prefix ^ "_gen.txt") in
      let nuc = from_file (prefix ^ "_nuc.txt") in
      let gen_alleles = List.map ~f:(fun a -> a.allele) gen.alt_elems in
      Ok (Alter_MSA.Merge.to_distance_arguments nuc gen_alleles)
  | None, (Some af)   ->
      let mp = from_file af in
      let targets =
        mp.alt_elems 
        |> List.map ~f:(fun alt -> (alt.allele, alt.seq))
        |> string_map_of_assoc 
      in
      Ok { Distances.reference = mp.reference
         ; Distances.reference_sequence = mp.ref_elems
         ; targets
         ; candidates = targets           (* compute distances to all alleles *)
         }
  | None, None        ->
      Error "Either a file or merge argument must be specified"

let allele_distances alignment merge distance_logic =
  to_distance_targets_and_candidates ?alignment ?merge () >>=
    begin fun args ->
      Distances.compute args distance_logic >>= fun dmap ->
        printf "allele, closests alleles \n";
        Ok (StringMap.iter dmap ~f:(fun ~key ~data ->
            let allst = List.map data ~f:(fun (s,d) -> sprintf "%s,%f" s d) in
            printf "%s,%s\n" key (String.concat ~sep:"," allst)))
  end
  |> function
      | Error e -> eprintf "%s" e; 1
      | Ok ()   -> 0

let () =
  let open Cmdliner in
  let open Cmdline_options in
  let allele_distances =
    let version = "0.0.0" in
    let doc = "Compute distances between HLA alleles." in
    let bug =
      sprintf "Browse and report new issues at <https://github.com/hammerlab/%s"
        repo
    in
    let man =
      [ `S "AUTHORS"
      ; `P "Leonid Rozenberg <leonidr@gmail.com>"
      ; `Noblank
      ; `S "BUGS"
      ; `P bug
      ]
    in
    Term.(const allele_distances
            (* Allele set construction args *)
            $ alignment_arg $ merge_arg
            (* How to measure distances. *)
            $ defaulting_distance_flag
        , info app_name ~version ~doc ~man)
  in
  match Term.eval allele_distances with
  | `Ok n            -> exit n
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
