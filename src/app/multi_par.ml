(* Typing via a Parametric PHMM. *)
open Util

let app_name = "multi_par"

let (//) = Filename.concat

let to_read_size_dependent
  (* Allele information source *)
    ~alignment_files ~merge_files ~distance ~impute
    ~skip_disk_cache =
    let als =
      List.map alignment_files ~f:(Common_options.input_alignments ~impute)
    in
    let mls =
      List.map merge_files ~f:(Common_options.input_merges ~distance ~impute)
    in
    match als @ mls with
    | []     -> Error "Neither a merge nor alignment file specified!"
    | inputs ->
      let selectors = [] in           (* Selectors NOT supported, on purpose! *)
      Ok (fun read_size ->
        List.map inputs ~f:(fun input ->
          let name = Alleles.Input.to_string input in
          let par_phmm_arg = Cache.par_phmm_args ~input ~selectors ~read_size in
          let pt = Cache.par_phmm ~skip_disk_cache par_phmm_arg in
          name, pt))

module Pdml = ParPHMM_drivers.Mulitple_loci

let to_comp conf read_size rp mode =
  let ptlst =
    time (sprintf "Setting up ParPHMM transitions with %d read_size" read_size)
      (fun () -> rp read_size)
  in
  time "Allocating forward pass workspaces"
    (fun () -> Pdml.init conf read_size ptlst mode)

let to_apply conf mode =
  fun acc fqi ->
    match acc with
    | `Setup rp ->
        let read_size = String.length fqi.Biocaml_unix.Fastq.sequence in
        let c = to_comp conf read_size rp mode in
        `Set (Pdml.apply c fqi; c)
    | `Set c ->
        `Set (Pdml.apply c fqi; c)

let to_paired conf mode =
  fun acc fq1 fq2 ->
    match acc with
    | `Setup rp ->
        let read_size = String.length fq1.Biocaml_unix.Fastq.sequence in
        let c = to_comp conf read_size rp mode in
        `Set (Pdml.paired c fq1 fq2; c)
    | `Set c ->
        `Set (Pdml.paired c fq1 fq2; c)

let across_fastq conf mode
  (* Fastq specific arguments. *)
  ?number_of_reads ~specific_reads file init =
  let f = to_apply conf mode in
  try
    Fastq.fold ?number_of_reads ~specific_reads ~init file ~f
    |> function
        | `Setup _  -> eprintf "Didn't find any reads."
        | `Set c    -> Pdml.output c stdout
  with ParPHMM_drivers.Read_error_parsing e ->
    eprintf "%s" e

let across_paired ~finish_singles conf mode
    (* Fastq specific arguments. *)
    ?number_of_reads ~specific_reads file1 file2 init =
  let f = to_paired conf mode in
  try
    begin
      if finish_singles then
        let ff = to_apply conf mode in
        let fs = ff in
        Fastq.fold_paired_both ?number_of_reads ~specific_reads file1 file2
          ~init ~f ~ff ~fs
      else
        Fastq.fold_paired ?number_of_reads ~specific_reads file1 file2 ~init ~f
    end
    |> function
        | `BothFinished o
        | `FinishedSingle o
        | `OneReadPairedFinished (_, o)
        | `StoppedByFilter o
        | `DesiredReads o ->
            match o with
            | `Setup _ -> eprintf "Didn't find any reads."
            | `Set c   -> Pdml.output c stdout
  with ParPHMM_drivers.Read_error_parsing e ->
    eprintf "%s" e

let type_
  (* Allele information source *)
    class1_gen_dir
    alignment_files merge_files distance
  (* Process *)
    skip_disk_cache
  (* What to do? *)
    fastq_file_lst number_of_reads specific_reads
    do_not_finish_singles
  (* options *)
    insert_p
    do_not_past_threshold_filter
    max_number_mismatches
    read_size_override
    not_band
    warmup
    number
    width
    likelihood_first
  (* how are we typing *)
    map_depth
    not_incremental_pairs
    mode
    forward_accuracy_opt
    =
  Option.value_map forward_accuracy_opt ~default:()
    ~f:(fun fa -> ParPHMM.dx := fa);
  let impute   = true in
  let band     =
    if not_band then
      None
    else
      Some { ParPHMM.warmup; number; width }
  in
  let alignment_files, merge_files =
    match class1_gen_dir with
    | None   -> alignment_files, merge_files
    | Some d -> [ d // "A_gen.txt"
                ; d // "B_gen.txt"
                ; d // "C_gen.txt"
                ], []
  in
  let past_threshold_filter = not do_not_past_threshold_filter in
  let incremental_pairs = not not_incremental_pairs in
  let finish_singles = not do_not_finish_singles in
  let conf = Pdml.conf ~insert_p ?band ?max_number_mismatches
                    ~past_threshold_filter ~incremental_pairs ()
  in
  let need_read_size_r =
    to_read_size_dependent
      ~alignment_files ~merge_files ~distance ~impute
      ~skip_disk_cache
  in
  let mode =
    match mode with
    | `Reducer -> `Reducer likelihood_first
    | `Mapper  -> `Mapper map_depth
  in
  match need_read_size_r with
  | Error e           -> eprintf "%s" e
  | Ok need_read_size ->
    let init =
      match read_size_override with
      | None   -> `Setup need_read_size
      | Some r -> `Set (to_comp conf r need_read_size mode)
    in
    begin match fastq_file_lst with
    | []              -> invalid_argf "Cmdliner lied!"
    | [fastq]         -> across_fastq conf mode
                            ?number_of_reads ~specific_reads fastq init
    | [read1; read2]  -> across_paired ~finish_singles conf mode
                            ?number_of_reads ~specific_reads read1 read2 init
    | lst             -> invalid_argf "More than 2, %d fastq files specified!"
                           (List.length lst)
    end

let () =
  let open Cmdliner in
  let open Common_options in
  let files_arg =
    let docv = "FILE" in
    let doc  = "File to lookup IMGT allele alignments. The alleles found in this \
                file will initially define the set of alleles to be used for \
                named locus. Supply multiple file (or merge) arguments to \
                compute the likelihood of a read against each locus and add to \
                the most likely (highest likelihood)." in
    Arg.(value & opt_all file [] & info ~doc ~docv ["f"; "file"])
  in
  let merges_arg =
    let parser_ path =
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
        `Ok path  (* Return path, and do appending later, the prefix is more useful. *)
    in
    let convrtr = parser_, (fun frmt -> Format.fprintf frmt "%s") in
    let docv = sprintf "[%s]" (String.concat ~sep:"|" Merge_mas.supported_genes) in
    let doc  =
      sprintf "Construct a merged (gDNA and cDNA) graph of the specified \
              prefix path. Currently only supports %s genes. The argument must \
              be a path to files with $(docv)_nuc.txt and $(docv)_gen.txt. \
              Combines with the file arguments to determine the set of loci to \
              type at the same time.. The set of alleles is defined by the \
              ones in the nuc file."
        (String.concat ~sep:", " Merge_mas.supported_genes)
    in
    Arg.(value & opt_all convrtr [] & info ~doc ~docv ["m"; "merge"])
  in
  let max_number_mismatches_arg =
    let docv = "POSITIVE INTEGER" in
    let doc = "Setup a filter on the reads to cancel evaluation once we've \
               seen this many mismatches." in
    Arg.(value & opt (some positive_int) None & info ~doc ~docv
          ["max-mismatches"])
  in
  let mode_flag =
    let open Arg in
    let modes =
      [ ( `Reducer
        , info ~doc:"Aggregate read data into per allele likelihoods, \
                    default. Each read's likelihood is added only to the most \
                    likely loci's."
          [ "reducer" ])
      ; ( `Mapper
        , info ~doc:
            (sprintf "Map each read into the most likely alleles and \
                     positions, across loci. Specify the %S argument to change
                     the number of elements that are reported."
              map_depth_argument)
          [ "map" ])
      ]
    in
    value & vflag `Reducer  modes
  in
  let class1gen_arg =
    let docv  = "DIRECTORY" in
    let doc   = "Short-cut argument that expands the given dir to look for \
                  A_gen.txt, B_gen.txt and C_gen.txt. This overrides any \
                 alignment files or merge arguments." in
    Arg.(value & opt (some dir) None & info ~doc ~docv ["class1-gen"])
  in
  let do_not_use_incremental_pairs_flag =
    let doc   = "By default when performing paired typing, we use the first \
                 (as specified by the read found in the first file on the \
                 command line) pair, to determine which loci is the most \
                 likely. This is faster but potentially less accurate. \
                 To use both pairs pass this flag. "
    in
    Arg.(value & flag & info ~doc ["do-not-incremental-pair"])
  in
  let type_ =
    let version = "0.0.0" in
    let doc = "Use a Parametric Profile Hidden Markov Model of HLA allele to \
               type fastq samples." in
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
    Term.(const type_
            (* Allele information source *)
            $ class1gen_arg
            $ files_arg $ merges_arg $ distance_flag
            (* What to do ? *)
            $ no_cache_flag
            (* What are we typing *)
            $ fastq_file_arg $ num_reads_arg $ specific_read_args
            $ do_not_finish_singles_flag
            (* options. *)
            $ insert_probability_arg
            $ do_not_past_threshold_filter_flag
            $ max_number_mismatches_arg
            $ read_size_override_arg
            $ not_band_flag
            $ band_warmup_arg
            $ number_bands_arg
            $ band_width_arg
            $ likelihood_first_flag
            (* How are we typing *)
            $ map_depth_arg
            $ do_not_use_incremental_pairs_flag
            $ mode_flag
            $ forward_pass_accuracy_arg
            (* $ map_allele_arg
            $ filter_flag $ multi_pos_flag $ stat_flag $ likelihood_error_arg
              $ upto_kmer_hood_arg
              $ allto_kmer_hood_arg
            (* Output *)
            $ print_top_flag
              $ do_not_normalize_flag
              $ bucket_arg
              $ error_output_flag
              $ reduce_resolution_arg *)
        , info app_name ~version ~doc ~man)
  in
  match Term.eval type_ with
  | `Ok ()           -> exit 0
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
