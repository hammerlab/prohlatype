
open Util
open Common_options

let app_name = "align2fasta"

let fail_on_parse a =
  match Nomenclature.parse a with
  | Ok (g, r) -> g, r
  | Error e   -> failwith e

let against_mp ?width mp merge_assoc out =
  let open MSA in
  let open MSA.Parser in
  let r = reference_sequence mp in
  let reference = mp.ref_elems in
  let gene, ref_res = fail_on_parse mp.reference in
  let check_for_merge =
    match merge_assoc with
    | []  -> fun _ -> " "
    | ma  -> let smap = string_map_of_assoc ma in
             begin fun a ->
               try
                 let info = StringMap.find a smap in
                 sprintf " %s " (Alter_MSA.info_to_string info)
               with Not_found -> " "
             end
  in
  let oc = open_out out in
  try
    List.map mp.alt_elems ~f:(fun (a, allele) ->
        let g, r = fail_on_parse a in
          if g <> gene then failwithf "Different genes: %s vs %s" gene g;
          r, a, allele_sequence ~reference ~allele ())
    |> fun l -> ((ref_res, mp.reference, r) :: l)
    |> List.sort ~cmp:(fun (r1,_,_) (r2,_,_) -> Nomenclature.compare_by_resolution r1 r2)
    |> List.iter ~f:(fun (_, a, s) ->
      fprintf oc ">%s%slength: %d\n" a (check_for_merge a) (String.length s);
      print_line ?width oc s);
    close_out oc
  with e ->
    close_out oc;
    raise e

let convert
  (* output destination. *)
  ofile
  (* output option *)
  width
  (* input *)
  alignment_file merge_file
  (* optional distance to trigger imputation, merging *)
  distance
  (* selectors *)
  regex_list specific_list without_list number_alleles
  do_not_ignore_suffixed_alleles
  =
  let selectors =
    aggregate_selectors ~regex_list ~specific_list ~without_list
      ?number_alleles ~do_not_ignore_suffixed_alleles
  in
  to_allele_input ?alignment_file ?merge_file ?distance ~selectors >>= fun i ->
    Alleles.Input.construct i >>= fun (mp, merge_assoc) ->
      let ofiledefault = Alleles.Input.to_short_fname_prefix i in
      let out = sprintf "%s.fasta" (Option.value ofile ~default:ofiledefault) in
      Ok (against_mp ?width mp merge_assoc out)

let () =
  let open Cmdliner in
  let output_fname_arg =
    let docv = "FILE" in
    let doc  = "Output file name, defaults to \"(input file)_.fasta\"."
    in
    Arg.(value & opt (some string) None & info ~doc ~docv ["o"; "output"])
  in
  let width_arg =
    let docv = "POSITIVE INTEGER" in
    let doc  = "Sequence, per-line, output width." in
    Arg.(value
          & opt (some positive_int) None
          & info ~doc ~docv ["output-width"])
  in
  let convert =
    let version = "0.0.0" in
    let doc = "Transform MHC IMGT alignments to fasta." in
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
    Term.(const convert
          $ output_fname_arg
          $ width_arg
          (* Allele information source *)
          $ alignment_arg $ merge_arg $ optional_distance_flag
          (* Allele selectors *)
          $ regex_arg $ allele_arg $ without_arg $ num_alt_arg
          $ do_not_ignore_suffixed_alleles_flag

        , info app_name ~version ~doc ~man)
  in
  match Term.eval convert with
  | `Ok (Ok ())      -> exit 0
  | `Ok (Error e)    -> failwith e
  | `Error _         -> failwith "error"
  | `Version | `Help -> exit 0
