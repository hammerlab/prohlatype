(* Common arguments shared between programs. *)

open Util
open Cmdliner

let repo = "prohlatype"

let positive_int_parser =
  fun s ->
    try
      let d = Scanf.sscanf s "%d" (fun x -> x) in
      if d <= 0 then
        `Error (s ^ " is not a positive integer")
      else
        `Ok d
    with Scanf.Scan_failure msg ->
      `Error msg

let positive_int =
  positive_int_parser, (fun frmt -> Format.fprintf frmt "%d")

let file_arg =
  let docv = "FILE" in
  let doc  = sprintf "File to lookup IMGT allele alignments." in
  Arg.(value & opt (some file) None & info ~doc ~docv ["f"; "file"])

let merge_arg =
  let parser path =
    let s = Filename.basename path in
    let n = path ^ "_nuc.txt" in
    let g = path ^ "_gen.txt" in
    if not (List.mem ~set:Merge_mas.supported_genes s) then
      `Error ("gene not supported: " ^ s)
    else if not (Sys.file_exists n) then
      `Error ("Nuclear alignment file doesn't exist: " ^ n)
    else if not (Sys.file_exists g) then
      `Error ("Genetic alignment file doesn't exist: " ^ n)
    else
      `Ok s   (* Return s, and do appending later, the prefix is more useful. *)
  in
  let convrtr = parser, (fun frmt -> Format.fprintf frmt "%s") in
  let docv = sprintf "[%s]" (String.concat ~sep:"|" Merge_mas.supported_genes) in
  let doc  =
    sprintf "Construct a merged (gDNA and cDNA) graph of the specified \
             prefix path. Currently only supports %s genes. The argument must \
             be a path to files with $(docv)_nuc.txt and $(docv)_gen.txt. \
             Overrides file arguments."
      (String.concat ~sep:", " Merge_mas.supported_genes)
  in
  Arg.(value & opt (some convrtr) None & info ~doc ~docv ["m"; "merge"])

let num_alt_arg =
  let docv = "POSITIVE INTEGER" in
  let doc  = "Number of alternate alleles to add to the graph. \
              If not specified, all of the alternate alleles in the alignment \
              file are added."
  in
  Arg.(value & opt (some positive_int) None & info ~doc ~docv ["n"; "num-alt"])

let allele_arg =
  let docv = "STRING" in
  let doc  = "Specify specfic alternate alleles to add to the graph. \
              By default the behavior is to add all (or [num-alt] if \
              specified) alleles (in alphabetical order) to the graph. Use \
              this option to construct a graph with specific alleles (ex. \
              A*01:02). One can repeat this argument to specify multiple \
              alleles. This argument will override the [num-alt] argument \
              but may be used in tandem with the allele-regex arguments."
  in
  Arg.(value & opt_all string [] & info ~doc ~docv ["a"; "allele"])

let allele_regex_arg =
  let docv = "REGEX" in
  let doc  = "Specify specfic alternate alleles to add to the graph via a \
              regex. This option is similar to allele, but lets the user \
              specify a class (specifically a pattern) of alleles to include \
              in the graph, without having to add them individually via the \
              allele argument. This argument will override [num-alt] but may be \
              used in tandem with allele flag. The regex format is POSIX, and \
              the '*' used in HLA allele names must be properly escaped from \
              the command line: ex. -ar \"A\\*02:03\""
  in
  Arg.(value & opt_all string [] & info ~doc ~docv ["ar"; "allele-regex"])

let remove_reference_flag =
  let doc  = "Remove the reference allele from the graph. The reference \
              allele is the one that is listed first in the alignments file. \
              Graphs are currently constructed based upon their \"diff\" to \
              the reference as represented in the alignments file. Therefore, \
              the original reference sequence must be a part of the graph \
              during construction. Specifying this flag will remove it from \
              the graph."
  in
  Arg.(value & flag & info ~doc ["no-reference"])

let no_cache_flag =
  let doc =
    sprintf "Do not use a disk cache (in %s sub directory of the current \
             directory) to search for previously (and then save) constructed \
             graphs."
      Cache.dir
  in
  Arg.(value & flag & info ~doc ["no-cache"])

let do_not_join_same_sequence_paths_flag =
  let doc = "Do not join same sequence paths; remove overlapping nodes at the same position, \
             in the string graph."
  in
  Arg.(value & flag & info ~doc ["do-not-join-same-sequence-paths"])

let to_filename_and_graph_args ?alignment_file ?merge_file num_alt_to_add
  ~allele_list ~allele_regex_list ~join_same_sequence ~remove_reference =
  let open Ref_graph in
  let to_ofname_and_cache_args base alignment_file =
    let option_based_fname, which =
      match allele_list, allele_regex_list with
      | [], [] ->
          begin
            match num_alt_to_add with
            | None      -> sprintf "%s_%b_%b_all" base join_same_sequence remove_reference
                          , None
            | Some n    -> sprintf "%s_%b_%b_%d" base join_same_sequence remove_reference n
                          , (Some (NumberOfAlts n))
          end
      | specific, regex -> let sr_dig =
                              String.concat (specific @ regex)
                              |> Digest.string
                              |> Digest.to_hex
                          in
                          sprintf "%s_spec_%s_%b" base sr_dig remove_reference
                          , (Some (Alleles { specific; regex }))
    in
    option_based_fname
    , { Cache.alignment_file
      ; Cache.which
      ; Cache.join_same_sequence
      ; Cache.remove_reference
      }
  in
  match alignment_file, merge_file with
  | _, (Some prefix)          ->
      let base = Filename.basename prefix ^ "_mgd" in
      Ok (to_ofname_and_cache_args base (MergeFromPrefix prefix))
  | (Some alignment_file), _  ->
      let base = Filename.basename alignment_file |> Filename.chop_extension in
      Ok (to_ofname_and_cache_args base (AlignmentFile alignment_file))
  | None, None                ->
      Error "Either a file or merge argument must be specified"

let verbose_flag =
  let doc = "Print progress messages to stdout." in
  Arg.(value & flag & info ~doc ["v"; "verbose"])

let kmer_size_arg =
  let default = 10 in
  let docv = "POSITIVE INTEGER" in
  let doc =
    sprintf "Number of consecutive nucleotides to use consider in K-mer
              index construction. Defaults to %d." default
  in
  Arg.(value & opt positive_int default & info ~doc ~docv ["k"; "kmer-size"])

let fastq_file_arg =
  let docv = "FILE" in
  let doc = "Fastq formatted DNA reads file." in
  Arg.(non_empty & pos_all file [] & info ~doc ~docv [])

let num_reads_arg =
  let docv = "POSITIVE INTEGER" in
  let doc = "Number of reads to take from the front of the FASTA file" in
  Arg.(value & opt (some positive_int) None & info ~doc ~docv ["reads"])

