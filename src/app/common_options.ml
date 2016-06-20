(* Common arguments shared between programs. *)

open Util
open Cmdliner

let repo = "prohlatype"

let file_arg =
  let docv = "file" in
  let def  = "./A_gen.txt" in
  let doc  =
    sprintf "File to lookup IMGT allele alignments. Defaults to %s." def
  in
  Arg.(value & opt file "A_gen.txt" & info ~doc ~docv ["f"; "file"])

let num_alt_arg =
  let docv = "positive integer" in
  let doc  = "Number of alternate alleles to add to the graph. \
              If not specified, all of the alternate alleles in a file are \
              added."
  in
  Arg.(value & opt (some int) None & info ~doc ~docv ["n"; "num-alt"])

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

let no_cache_flag =
  let doc =
    sprintf "Do not use a disk cache (in %s sub directory of the current \
             directory) to search for previously (and then save) constructed \
             graphs."
      Cache.dir
  in
  Arg.(value & flag & info ~doc ["no-cache"])

let to_filename_and_graph_args alignment_file num_alt_to_add allele_list =
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
    | lst -> sprintf "%s_spec_%d" base (Hashtbl.hash lst), (Some (SpecificAlleles lst))
  in
  option_based_fname, { Cache.alignment_file; Cache.which }

