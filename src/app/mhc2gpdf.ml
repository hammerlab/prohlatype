
open Util

let construct ofile alignment_file num_alt_to_add allele_list notshort no_pdf no_open =
  let open To_graph in
  let base = Filename.basename alignment_file |> Filename.chop_extension in
  let option_based_fname, which =
    match allele_list with
    | []  ->
        begin
          match num_alt_to_add with
          | None   -> sprintf "%s_all" base, None
          | Some n -> sprintf "%s_%d" base n, (Some (NumberOfAlts n))
        end
    | lst -> sprintf "%s_%d" base (List.length lst), (Some (SpecificAlleles lst))
  in
  let ofile = Option.value ofile ~default:option_based_fname in
  let short = not notshort in
  let pdf   = not no_pdf in
  let open_ = not no_open in
  construct_from_file { alignment_file; which }
  |> Ref_graph.output ~short ~pdf ~open_ ofile

let app_name = "mhc2gpdf"
let repo = "supreme-tribble"

let () =
  let open Cmdliner in
  let file_arg =
    let docv = "file" in
    let def  = "./A_gen.txt" in
    let doc  =
      sprintf "File to lookup IMGT allele alignments. Defaults to %s." def
    in
    Arg.(value & opt file "A_gen.txt" & info ~doc ~docv ["f"; "file"])
  in
  let output_fname_arg =
    let docv = "file" in
    let doc  = "Output file name, defaults to a file (input file)_(number of \
                non reference alleles).[pdf|dot]."
    in
    Arg.(value & opt (some string) None & info ~doc ~docv ["o"; "output"])
  in
  let num_alt_arg =
    let docv = "positive integer" in
    let doc  = "Number of alternate alleles to add to the graph. \
                If not specified, all of the alternate alleles in a file are \
                added. Take note that rendering (converting to pdf with dot) \
                a graph with a lot (> 50) of paths could be time consuming."
    in
    Arg.(value & opt (some int) None & info ~doc ~docv ["n"; "num-alt"])
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
  let allele_arg =
    let docv = "allele name" in
    let doc = "Specify specfic alternate alleles to add to the graph.\
               By default the behavior is to add [num_alt] alleles (in \
               alphabetical order) to the graph. Use this option to construct \
               a graph with specific alleles (ex. A*01:02). Specify multiple \
               alleles with multiple arguments. This flag will override \
               [num_alt]."
    in
    Arg.(value & opt_all string [] & info ~doc ~docv ["a"; "allele"])
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
        , info app_name ~version ~doc ~man)
  in
  match Term.eval construct with
  | `Ok n -> exit n
  | `Error _ -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
