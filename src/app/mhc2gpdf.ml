
open Prohlatype

let to_filename_and_graph_args
  (* Allele information source *)
    ?alignment_file ?merge_file ?distance
  (* Allele selectors *)
    ~regex_list
    ~specific_list
    ~without_list
    ?number_alleles
    ~do_not_ignore_suffixed_alleles
  (* Graph modifiers. *)
    ~join_same_sequence =
    let selectors =
      Cmdline_options.aggregate_selectors ~regex_list ~specific_list ~without_list
        ?number_alleles ~do_not_ignore_suffixed_alleles
    in
    Cmdline_options.to_allele_input ?alignment_file ?merge_file ?distance ~selectors
      >>= fun input ->
            let arg = { Ref_graph.join_same_sequence } in
            let graph_arg = Cache.graph_args ~arg ~input in
            let option_based_fname = Cache.graph_args_to_string graph_arg in
            Ok (option_based_fname, graph_arg)


let construct
  (* input *)
  alignment_file merge_file distance
  (* output *)
  ofile
  (* allele selection. *)
  regex_list specific_list without_list number_alleles
  do_not_ignore_suffixed_alleles
  (* Construction args *)
  not_join_same_seq
  (* output configuration. *)
  notshort no_pdf no_open skip_disk_cache max_edge_char_length
  not_human_edges not_compress_edges not_compress_start
  not_insert_newlines =
  let open Ref_graph in
  let open Cache in
  let join_same_sequence = (not not_join_same_seq) in
  let fname_cargs_result =
    to_filename_and_graph_args
      ?alignment_file ?merge_file ?distance
      ~specific_list ~regex_list ~without_list ?number_alleles
      ~do_not_ignore_suffixed_alleles
      ~join_same_sequence
  in
  match fname_cargs_result with
  | Error msg ->
      eprintf "%s" msg;
      1
  | Ok (option_based_fname, cargs) ->
      let ofile = Option.value ofile ~default:option_based_fname in
      let short = not notshort in
      let pdf   = not no_pdf in
      let open_ = not no_open in
      let max_length = max_edge_char_length in
      let human_edges = not not_human_edges in
      let compress_edges = not not_compress_edges in
      let compress_start = not not_compress_start in
      let insert_newlines = not not_insert_newlines in
      Cache.graph ~skip_disk_cache cargs
      |> Ref_graph.output ~compress_edges ~compress_start ~insert_newlines
                            ~human_edges ?max_length ~short ~pdf ~open_ ofile

let app_name = "mhc2gpdf"

let () =
  let open Cmdliner in
  let open Cmdline_options in
  let output_fname_arg =
    let docv = "FILE" in
    let doc  = "Output file name, defaults to a file \"(input file)_(number of \
                non reference alleles).[pdf|dot]\" ."
    in
    Arg.(value & opt (some string) None & info ~doc ~docv ["o"; "output"])
  in
  let not_short_flag =
    let doc = "Whether to turn off short labels on sequence nodes.\n\
               By default sequences longer than 10 are compressed to the first \
               4, ellipsis and then last 3 (ex ACGT...TGA). Providing this \
               flag will create nodes with the full sequence."
    in
    Arg.(value & flag & info ~doc ["not-short"])
  in
  let no_pdf_flag =
    let doc = "Whether to create a pdf output with the graphviz program dot.\n\
               By default, after a dot file has been written, this program \
               calls dot on it to render it to pdf. Provide this flag to skip \
               that step."
    in
    Arg.(value & flag & info ~doc ["no-pdf"])
  in
  let no_open_flag =
    let doc = "Whether to open the pdf file.\n\
               By default if a pdf file is created this program will try to \
               open it with the open system command. Specify this flag to turn \
               off that behavior."
    in
    Arg.(value & flag & info ~doc ["no-open"])
  in
  let not_human_edges_flag =
    let doc = "Do not use heuristics to shorten edge names." in
    Arg.(value & flag & info ~doc ["not-human-edges"])
  in
  let max_edge_char_length_flag =
    let docv = "POSITIVE INTEGER" in
    let doc = "Specify a max number of characters to label the edges. (ex 100)" in
    Arg.(value & opt (some positive_int) None
               & info ~doc ~docv ["max-edge-char-length"])
  in
  let do_not_compress_edges_flag =
    let doc = "Do not run encode the alleles along the edges." in
    Arg.(value & flag & info ~doc ["do-not-compress-edges"])
  in
  let do_not_compress_start_flag =
    let doc = "Do not run encode the alleles inside a merged start box." in
    Arg.(value & flag & info ~doc ["do-not-compress-start"])
  in
  let do_not_insert_newlines_flag =
    let doc = "Do not insert newlines into a list of alleles." in
    Arg.(value & flag & info ~doc ["do-not-insert-newline-into-alleles"])
  in
  let construct =
    let version = "0.0.0" in
    let doc = "Transform MHC IMGT alignments to pdf graphs." in
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
    Term.(const construct
            (* input files *)
            $ alignment_arg $ merge_arg $ optional_distance_flag
            (* output file *)
            $ output_fname_arg
            (* allele selection. *)
            $ regex_arg
            $ allele_arg
            $ without_arg
            $ num_alt_arg
            $ do_not_ignore_suffixed_alleles_flag
            (* construction args. *)
            $ do_not_join_same_sequence_paths_flag
            (* output configuration *)
            $ not_short_flag $ no_pdf_flag $ no_open_flag
            $ no_cache_flag
            $ max_edge_char_length_flag
            $ not_human_edges_flag
            $ do_not_compress_edges_flag
            $ do_not_compress_start_flag
            $ do_not_insert_newlines_flag
        , info app_name ~version ~doc ~man)
  in
  match Term.eval construct with
  | `Ok n    -> exit n
  | `Error _ -> failwith "error"
  | `Version | `Help -> exit 0
