(** Compute distances between MHC alleles. *)
open Util

let app_name = "allele_distances"

let allele_distances alignment_file_opt merge_opt distance_logic =
  Common_options.to_distance_targets_and_candidates alignment_file_opt merge_opt >>=
    begin fun (reference, reference_sequence, targets, candidates) ->
      Distances.compute ~reference ~reference_sequence ~targets ~candidates
        distance_logic >>= fun dmap ->
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
  let open Common_options in
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
            $ file_arg $ merge_arg
            (* How to measure distances. *)
            $ defaulting_distance_flag
        , info app_name ~version ~doc ~man)
  in
  match Term.eval allele_distances with
  | `Ok n            -> exit n
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
