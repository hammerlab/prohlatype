(** Compute distances between MHC alleles. *)
open Prohlatype
open Cmdline_options
open Cmdliner

let app_name = "allele_distances"

let to_distance_targets_and_candidates file_or_prefix =
  let open MSA.Parser in
  match file_or_prefix with
  | `Prefix p  ->
      let gen = from_file (p ^ "_gen.txt") in
      let nuc = from_file (p ^ "_nuc.txt") in
      let gen_alleles = List.map ~f:(fun a -> a.allele) gen.alt_elems in
      Alter_MSA.Merge.to_distance_arguments nuc gen_alleles
  | `File path ->
      let mp = from_file path in
      let targets =
        mp.alt_elems
        |> List.map ~f:(fun alt -> (alt.allele, alt.seq))
        |> string_map_of_assoc
      in
      { Distances.reference = mp.reference
      ; Distances.reference_sequence = mp.ref_elems
      ; targets
      ; candidates = targets           (* compute distances to all alleles *)
      }

let distance_calculation_error = 2
let distance_calculation_error_exit_info =
  Term.exit_info ~doc:"distance calculation error"
    distance_calculation_error

let allele_distances file_or_prefix distance_logic =
  let args = to_distance_targets_and_candidates file_or_prefix in
  let dmap_res = Distances.compute args distance_logic in
  match dmap_res with
  | Error m -> errored distance_calculation_error "%s" m
  | Ok dmap ->
    printf "allele, closests alleles \n";
    StringMap.iter dmap ~f:(fun ~key ~data ->
      let allst = List.map data ~f:(fun (s,d) -> sprintf "%s,%f" s d) in
      printf "%s,%s\n" key (String.concat ~sep:"," allst));
    Term.exit_status_success

let () =
  let allele_distances =
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
            $ (file_or_prefix_arg ~action:"use")
            (* How to measure distances. *)
            $ defaulting_distance_flag
        , info app_name ~version ~doc ~man
              ~exits:(distance_calculation_error_exit_info :: default_exits))
  in
  Term.(exit_status (eval allele_distances))
