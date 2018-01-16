
open Prohlatype
open Cmdline_options
open Cmdliner

let app_name = "align2fasta"

let fail_on_parse a =
  match Nomenclature.parse a with
  | Ok (l, r) -> l, r
  | Error e   -> failwith e

let alters_to_string = function
  | []  -> " "
  | lst -> sprintf " %s "
              (string_of_list lst ~sep:"," ~f:MSA.Alteration.to_string)

let with_oc out f =
  let oc = open_out out in
  try
    f oc;
    close_out oc
  with e ->
    close_out oc;
    raise e

let single_mp ?width oc mp =
  let open MSA in
  let open MSA.Parser in
  let r = reference_sequence mp in
  let reference = mp.ref_elems in
  let locus, ref_res = fail_on_parse mp.reference in
  List.map mp.alt_elems ~f:(fun a ->
      let l, r = fail_on_parse a.allele in
      if l <> locus then
        failwithf "Different loci: %s vs %s"
          (Nomenclature.show_locus locus)
          (Nomenclature.show_locus l);
      let aseq = allele_sequence ~reference ~allele:a.seq () in
      r, a.allele, aseq, a.alters)
  |> fun l -> ((ref_res, mp.reference, r, []) :: l)
  |> List.sort ~cmp:(fun (r1,_,_,_) (r2,_,_,_) -> Nomenclature.compare_by_resolution r1 r2)
  |> List.iter ~f:(fun (_, a, s, alters) ->
      fprintf oc ">%s%slength: %d\n"
        a (alters_to_string alters) (String.length s);
      print_line ?width oc s)

let against_mps ?width out mplst =
  with_oc out (fun oc -> List.iter ~f:(single_mp ?width oc) mplst)

let invalid_alignment_file_error = 4
let invalid_alignment_file_error_exit_info =
  Term.exit_info ~doc:"alignment file error" invalid_alignment_file_error

let convert
  (* output destination. *)
  ofile
  (* output option *)
  width
  (* optional distance to trigger imputation, merging *)
  distance
  (* selectors *)
  regex_list
  specific_list
  number_alleles
  do_not_ignore_suffixed_alleles
  (* input *)
  only_class1
  gen_nuc_mgd
  file_prefix_or_directory
  =
  let selectors =
    aggregate_selectors ~regex_list ~specific_list ~number_alleles
      ~do_not_ignore_suffixed_alleles
  in
  let ailst =
    file_prefix_or_directory_to_allele_input_list
      ?distance ~selectors
      file_prefix_or_directory
      only_class1 gen_nuc_mgd
  in
  match ailst with
  | []      -> errored Term.exit_status_cli_error "No input sent"
  | h :: t  ->
      let prefix =
        match ofile with
        | Some v  -> v
        | None    ->
            let ofiledefault = Alleles.Input.to_short_fname_prefix h in
            (* Warn if t isn't empty, ie there is more than one input. *)
            if t <> [] then begin
              eprintf "There is more than one input gene input and no output \
                       file. Writing to file %s.fasta"
                ofiledefault;
              ofiledefault
            end else
              ofiledefault
      in
      let out = sprintf "%s.fasta" prefix in
      let mps =
        list_fold_ok (h :: t) ~init:[] ~f:(fun acc i ->
            Alleles.Input.construct i >>= fun mp -> Ok (mp :: acc))
      in
      match mps with
      | Error m  -> errored invalid_alignment_file_error "%s" m
      | Ok mplst -> against_mps ?width out (List.rev mplst);
                    Term.exit_status_success

let () =
  let open Cmdliner in
  let output_fname_arg =
    let docv = "FILE" in
    let doc  = "Output file name, defaults to \"(input file basename).fasta\"."
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
    let doc = "Transform IMGT/HLA's alignments to fasta." in
    let description =
      [ `P  "$(tname) program transforms IMGT/HLA's per locus alignment \
            files to FASTA files. The alignments are the text files found in \
            the alignments folder and start with
            \"HLA-* Genomic Sequence Alignments\" header."

      ; `P "The IMGT/HLA database already distributes FASTA files for the \
            MHC genes in their database, so why create a tool to recreate \
            them? One minor reason is that, rarely, the alignment information \
            may be inconsistent with the FASTA files. Because downstream tools \
            in this project rely upon the alignment information as the \
            ground-truth, this tool may be helpful."

      ; `P (sprintf
           "The raison d'être for this tool is to invoke powerful processing \
            algorithms in this project to impute the missing parts of alleles \
            in a given gene. At the time of the creation of this project we \
            have incomplete sequence information for many of the (sometimes \
            thousands) alleles. To alleviate the sparseness one may specify \
            distance arguments (%s, %s, or %s) that will invoke imputation \
            logic for the missing segments. One may also ask for merged loci \
            which will combine the alleles whose data is derived from gDNA \
            with those derived from cDNA. Finally, one may tailor the list of \
            gene's and/or alleles to include in the output."
              trie_argument
              weighted_per_segment_argument
              reference_distance_argument)

      ; `P (allele_selector_paragraph "one can compare sequences between")
      ]
    in
    let examples =
      [ `P "To create a fasta for HLA-A from genetic DNA sequenced alleles \
            where the missing data (typically at the ends) has been imputed \
            by borrowing from the allele closest in the nomenclature \
            trie and write to a specific file:"
      ; `Pre (sprintf
          "%s path-to-IMGTHLA/alignments/A_gen.txt --%s -o hla_a_trie_imputed"
           app_name trie_argument)

      ; `P "To create a fasta for HLA-C merging gDNA (C_gen.txt) and \
            cDNA (C_nuc.txt) derivied sequences \
            and using the reference sequence \
            to fill in missing segments:"
      ; `Pre (sprintf
          "%s --%s -o output path-to-IMGTHLA/alignments/C"
          app_name reference_distance_argument)

      ; `P "To cast the widest possible net; all class I, merge the gDNA and \
            cDNA and impute missing segments."
      ; `Pre (sprintf
          "%s --%s -o output path-to-IMGTHLA/alignments"
          app_name merged_argument)
      ]
    in
    let man =
      let open Manpage in
      [ `S s_description
      ; `Blocks description
      ; `S s_examples
      ; `Blocks examples
      ; `S s_bugs
      ; `P bugs
      ; `S s_authors
      ; `P author
      ]
    in
    Term.(const convert
          $ output_fname_arg
          $ width_arg
          (* Allele information source *)
          $ optional_distance_flag
          (* Allele selectors *)
          $ regex_arg
          $ specific_arg
          $ number_alleles_arg
          $ do_not_ignore_suffixed_alleles_flag
          (* Input *)
          $ only_class1_directory_flag
          $ gen_nuc_merged_flag
          $ (file_prefix_or_directory_arg "process")
        , info app_name
            ~version
            ~doc
            ~man
            ~exits:(invalid_alignment_file_error_exit_info :: default_exits))
  in
  Term.(exit_status (eval convert))
