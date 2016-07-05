
open Util
open Common_options

let construct ofile alignment_file num_alt_to_add allele_list notshort no_pdf
      no_open skip_disk_cache max_edge_char_length not_human_edges_flag =
  let open To_graph in
  let open Cache in
  let option_based_fname, cargs =
    to_filename_and_graph_args alignment_file num_alt_to_add allele_list
  in
  let ofile = Option.value ofile ~default:option_based_fname in
  let short = not notshort in
  let pdf   = not no_pdf in
  let open_ = not no_open in
  let max_length = max_edge_char_length in
  let human_edges = not not_human_edges_flag in
  Cache.graph ~skip_disk_cache cargs (*{ alignment_file; which } *)
  |> Ref_graph.output ~human_edges ?max_length ~short ~pdf ~open_ ofile

let app_name = "mhc2gpdf"

let () =
  let open Cmdliner in
  let output_fname_arg =
    let docv = "file" in
    let doc  = "Output file name, defaults to a file (input file)_(number of \
                non reference alleles).[pdf|dot]."
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
    let docv = "positive integer" in
    let doc = "Specify a max number of characters to label the edges. (ex 100)" in
    Arg.(value & opt (some positive_int) None & info ~doc ~docv ["max-edge-char-length"])
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
            $ output_fname_arg $ file_arg
            $ num_alt_arg $ allele_arg
            $ not_short_flag $ no_pdf_flag $ no_open_flag
            $ no_cache_flag
            $ max_edge_char_length_flag
            $ not_human_edges_flag
        , info app_name ~version ~doc ~man)
  in
  match Term.eval construct with
  | `Ok n -> exit n
  | `Error _ -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
