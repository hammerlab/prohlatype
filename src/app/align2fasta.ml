
open Util
open Common_options

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

let convert
  (* output destination. *)
  ofile
  (* output option *)
  width
  (* input *)
  class1_gen_dir
  class1_nuc_dir
  class1_mgd_dir
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
  let opt_to_lst = Option.value_map ~default:[] ~f:(fun s -> [s]) in
  let alignment_files, merge_files =
    class_selectors class1_gen_dir class1_nuc_dir class1_mgd_dir
      (opt_to_lst alignment_file)
      (opt_to_lst merge_file)
  in
  to_allele_inputs ~alignment_files ~merge_files ?distance ~selectors >>= function
    |  []    -> Error "No input sent"
    | h :: t ->
      let ofiledefault = Alleles.Input.to_short_fname_prefix h in
      let out = sprintf "%s.fasta" (Option.value ofile ~default:ofiledefault) in
      list_fold_ok (h :: t) ~init:[] ~f:(fun acc i ->
          Alleles.Input.construct i >>= fun mp -> Ok (mp :: acc))
        >>= fun mplst -> Ok (against_mps ?width out (List.rev mplst))

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
    let doc = "Transform IMGT/HLA's alignments to fasta." in
    let bugs =
      sprintf "Browse and report new issues at <https://github.com/hammerlab/%s>"
              repo
    in
    let description =
      let open Common_options in
      [ `P (sprintf
           "%s is a program that transforms IMGT/HLA's per locus alignment \
            files to FASTA files. The alignments are the text files found in \
            the alignments folder and start with
            \"HLA-A Genomic Sequence Alignments\"." app_name)

      ; `P "The IMGT-HLA database already distributes FASTA files for the \
            MHC genes in their database, so why create a tool to recreate \
            them? One minor reason is that, rarely, the alignment information \
            may be inconsistent with the FASTA files. Because downstream tools \
            in this project rely upon the alignment information as the \
            ground-truth, this tool may be helpful."

      ; `P (sprintf
           "The raison d'être for this tool is to invoke powerful processing \
            algorithms in this project to impute the missing parts of alleles \
            in a given gene. At the time of the creation of this project we \
            have incomplete sequence information for all of the (sometimes \
            thousands) alleles. To alleviate the sparseness one may specify \
            distance arguments (%s, %s, or %s) that will invoke imputation \
            logic for the missing segments. One may also ask for merged loci \
            which will combine the alleles whose data is derived from gDNA \
            with those derived from cDNA. Finally, one may tailor the list of \
            gene's and/or alleles to include in the output."
              trie_argument
              weighted_per_segment_argument
              reference_distance_argument)
      ]
    in
    let examples =
      [ `P "To create a fasta for HLA-A from genetic DNA sequenced alleles \
            where the missing data (typically at the ends) has been imputed \
            by looking borrowing from the allele closest in the nomenclature \
            trie:"
      ; `Pre (sprintf
          "\t%s --alignment path-to-IMGTHLA/alignments/A_gen.txt --trie -o hla_a_trie_imputed.fasta"
           app_name)
      ; `P "To create a fasta for all class I genes (A,B,C) that $(b,merges) \
            their gDNA (ex. C_gen.txt) and the cDNA (ex. C_nuc.txt) derived \
            alleles using a weighted by segment (exon's) distance similarity \
            to fill in missing intron/UTR segments:"
      ; `Pre (sprintf
          "\t%s--class1-mgd path-to-IMGTHLA/alignments --weighted-segment -o hla_class_I_ws_imputed_mgd.fasta "
          app_name)
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
      ; `P "Leonid Rozenberg <leonidr@gmail.com>"
      ]
    in
    Term.(const convert
          $ output_fname_arg
          $ width_arg
          (* Allele information source *)
          $ class1gen_arg
          $ class1nuc_arg
          $ class1mgd_arg
          $ alignment_arg
          $ merge_arg
          $ optional_distance_flag
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
