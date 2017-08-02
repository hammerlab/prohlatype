(* Typing via a Parametric PHMM. *)
open Util

let app_name = "par_type"

let to_read_size_dependent
  (* Allele information source *)
    ?alignment_file ?merge_file ~distance
  (* Allele selectors *)
    ~regex_list
    ~specific_list
    ~without_list
    ?number_alleles
    ~do_not_ignore_suffixed_alleles
  (* Cache logic. *)
    ~skip_disk_cache
    =
    let open Common_options in
    let selectors =
      aggregate_selectors ~regex_list ~specific_list
        ~without_list ?number_alleles ~do_not_ignore_suffixed_alleles
    in
    to_allele_input ?alignment_file ?merge_file ~distance ~selectors
      >>= fun input ->
          Ok (fun read_size ->
                let par_phmm_args = Cache.par_phmm_args ~input ~read_size in
                Cache.par_phmm ~skip_disk_cache par_phmm_args)

module BoilerPlate (S : ParPHMM_drivers.S) = struct

  module D = ParPHMM_drivers.Make_single_loci(S)

  let init conf rp read_size opt =
    let pt =
      time (sprintf "Setting up ParPHMM transitions with %d read_size" read_size)
        (fun () -> rp read_size)
    in
    time "Allocating forward pass workspace"
      (fun () -> D.init conf read_size pt opt)

  let single conf opt =
    fun acc fqi ->
      match acc with
      | `Setup rp ->
          let read_size = String.length fqi.Biocaml_unix.Fastq.sequence in
          let t, s = init conf rp read_size opt in
          D.merge t s (D.single t fqi);
          `Set (t, s)
      | `Set (t,s) ->
          D.merge t s (D.single t fqi);
          acc

  let paired conf opt =
    fun acc fq1 fq2 ->
      match acc with
      | `Setup rp ->
          let read_size = String.length fq1.Biocaml_unix.Fastq.sequence in
          let t, s = init conf rp read_size opt in
          D.merge t s (D.paired t fq1 fq2);
          `Set (t, s)
      | `Set (t, s) ->
          D.merge t s (D.paired t fq1 fq2);
          acc

  let across_fastq conf opt
      (* Fastq specific arguments. *)
      ?number_of_reads ~specific_reads file init =
    let f = single conf opt in
    try
      Fastq.fold ?number_of_reads ~specific_reads file ~init ~f
      |> function
          | `Setup _    -> eprintf "Didn't find any reads."
          | `Set (t, s) -> D.output t s stdout
    with ParPHMM_drivers.Fastq_items.Read_error_parsing e ->
      eprintf "%s" e

  let across_paired ~finish_singles conf opt
      (* Fastq specific arguments. *)
      ?number_of_reads ~specific_reads file1 file2 init =
    let f = paired conf opt in
    try
      begin
        if finish_singles then
          let fs = single conf opt in
          let ff = fs in
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
              | `Setup _    -> eprintf "Didn't find any reads."
              | `Set (t, s) -> D.output t s stdout
    with ParPHMM_drivers.Fastq_items.Read_error_parsing e ->
      eprintf "%s" e

end (* BoilerPlate *)

module Bv = BoilerPlate(ParPHMM_drivers.Viterbi)

let viterbi read_size_override need_read_size conf fastq_file_list
  number_of_reads specific_reads finish_singles =
  let init =
    match read_size_override with
    | None   -> `Setup need_read_size
    | Some r -> `Set (Bv.init conf need_read_size r ())
  in
  match fastq_file_list with
  | []              -> invalid_argf "Cmdliner lied!"
  | [fastq]         -> Bv.across_fastq conf ()
                            ?number_of_reads ~specific_reads fastq init
  | [read1; read2]  -> Bv.across_paired ~finish_singles conf ()
                            ?number_of_reads ~specific_reads read1 read2 init
  | lst             -> invalid_argf "More than 2, %d fastq files specified!"
                          (List.length lst)

module Bf = BoilerPlate(ParPHMM_drivers.Forward)

let forward opt read_size_override need_read_size conf fastq_file_list
  number_of_reads specific_reads finish_singles =
  let init =
    match read_size_override with
    | None   -> `Setup need_read_size
    | Some r -> `Set (Bf.init conf need_read_size r opt)
  in
  match fastq_file_list with
  | []              -> invalid_argf "Cmdliner lied!"
  | [fastq]         -> Bf.across_fastq conf opt
                          ?number_of_reads ~specific_reads fastq init
  | [read1; read2]  -> Bf.across_paired ~finish_singles conf opt
                          ?number_of_reads ~specific_reads read1 read2 init
  | lst             -> invalid_argf "More than 2, %d fastq files specified!"
                        (List.length lst)

let type_
  (* Allele information source *)
    alignment_file merge_file distance
  (* Allele selectors *)
    regex_list specific_list without_list number_alleles
    do_not_ignore_suffixed_alleles
    allele
  (* Process *)
    skip_disk_cache
  (* What to do? *)
    fastq_file_list number_of_reads specific_reads
    do_not_finish_singles
  (* options *)
    insert_p
    do_not_past_threshold_filter
    max_number_mismatches
    read_size_override
    not_check_rc
    not_band
    warmup
    number
    width
    likelihood_first
    zygosity_report_size
  (* how are we typing *)
    map_depth
    mode
    forward_accuracy_opt
    =
  Option.value_map forward_accuracy_opt ~default:()
    ~f:(fun fa -> ParPHMM.dx := fa);
  to_read_size_dependent
    ?alignment_file ?merge_file ~distance
    ~regex_list ~specific_list ~without_list ?number_alleles
    ~do_not_ignore_suffixed_alleles
    ~skip_disk_cache
  |> function
      | Error e           -> eprintf "%s" e       (* Construction args don't make sense.*)
      | Ok need_read_size ->
        (* Rename some arguments to clarify logic. *)
        let band     =
          if not_band then
              None
            else
              Some { ParPHMM.warmup; number; width }
          in
        let check_rc = not not_check_rc in
        let past_threshold_filter = not do_not_past_threshold_filter in
        let finish_singles = not do_not_finish_singles in
        let conf =
          ParPHMM_drivers.single_conf ?allele ~insert_p ?band
            ?max_number_mismatches ~past_threshold_filter ~check_rc ()
        in
        match mode with
        | `Viterbi ->
            viterbi read_size_override need_read_size conf fastq_file_list
              number_of_reads specific_reads finish_singles
        | `Forward ->
            let opt = { ParPHMM_drivers.likelihood_first
                      ; zygosity_report_size
                      ; report_size = map_depth
                      } in
            forward opt read_size_override need_read_size conf fastq_file_list
              number_of_reads specific_reads finish_singles

let () =
  let open Cmdliner in
  let open Common_options in
  let not_check_rc_flag =
    let doc = "Do not check the reverse complement." in
    Arg.(value & flag & info ~doc ["not-check-rc"])
  in
  let specific_allele_argument = "allele" in
  let spec_allele_arg =
    let docv = "ALLELE" in
    let doc  = "Use a faster mode where we measure the likelihood for just \
                the passed allele. The allele must be found in the alignment \
                or merge file." in
    Arg.(value & opt (some string) None
               & info ~doc ~docv [specific_allele_argument])
  in
  let mode_flag =
    let open Arg in
    let modes =
      [ ( `Viterbi, info ~doc:
            (sprintf "Map each read into the viterbi decoded path against a \
                      specified allele. Specify the allele with %S \
                      argument, defaults to using the reference of the loci."
              specific_allele_argument)
          [ "viterbi" ])
      ; ( `Forward, info ~doc:
          (sprintf "Default mode: aggregate read data into per allele \
            likelihoods. This mode reports the aggregated per allele \
            likelihoods as well as a subset of all possible zygosity
            likelihoods (control the number with %S). Furthermore, for each \
            read we'll report the most likely alleles, their likely emission \
            and the position with the largest emission. Specify the %S \
            argument to change the number of elements that are reported."
              zygosity_report_size_argument
              map_depth_argument)
          [ "forward" ])
      ]
    in
    value & vflag `Forward modes
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
            $ alignment_arg $ merge_arg $ defaulting_distance_flag
            (* Allele selectors *)
            $ regex_arg $ allele_arg $ without_arg $ num_alt_arg
            $ do_not_ignore_suffixed_alleles_flag
            $ spec_allele_arg
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
            $ not_check_rc_flag
            $ not_band_flag
            $ band_warmup_arg
            $ number_bands_arg
            $ band_width_arg
            $ likelihood_first_flag
            $ zygosity_report_size_arg
            (* How are we typing *)
            $ map_depth_arg
            $ mode_flag
            $ forward_pass_accuracy_arg
            (* $ map_allele_arg
            $ filter_flag $ multi_pos_flag $ stat_flag $ likelihood_error_arg
              $ do_not_check_rc_flag
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
