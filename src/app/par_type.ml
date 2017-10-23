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

module Pd = ParPHMM_drivers
module Bf = Pd.Sequential(Pd.Forward)

let forward read_length_override need_read_length conf fastq_file_list
  number_of_reads specific_reads finish_singles =
  let init =
    match read_length_override with
    | None   -> `Setup need_read_length
    | Some r -> `Set (Bf.init need_read_length conf r)
  in
  match fastq_file_list with
  | []              -> invalid_argf "Cmdliner lied!"
  | [fastq]         -> Bf.across_fastq conf
                          ?number_of_reads ~specific_reads fastq init
  | [read1; read2]  -> Bf.across_paired ~finish_singles conf
                          ?number_of_reads ~specific_reads read1 read2 init
  | lst             -> invalid_argf "More than 2, %d fastq files specified!"
                        (List.length lst)

module Bv = Pd.Sequential(Pd.Viterbi)

let viterbi read_length_override need_read_length conf fastq_file_list
  number_of_reads specific_reads finish_singles =
  let init =
    match read_length_override with
    | None   -> `Setup need_read_length
    | Some r -> `Set (Bv.init need_read_length conf r)
  in
  match fastq_file_list with
  | []              -> invalid_argf "Cmdliner lied!"
  | [fastq]         -> Bv.across_fastq conf
                            ?number_of_reads ~specific_reads fastq init
  | [read1; read2]  -> Bv.across_paired ~finish_singles conf
                            ?number_of_reads ~specific_reads read1 read2 init
  | lst             -> invalid_argf "More than 2, %d fastq files specified!"
                          (List.length lst)


module Pf = Pd.Parallel(Pd.Forward)

let p_forward read_length_override need_read_length conf fastq_file_list
  number_of_reads specific_reads finish_singles nprocs =
  match read_length_override with
  | None   -> invalid_argf "Must specify read size for pallel!"
  | Some r -> let state = Pf.init need_read_length conf r in
              begin match fastq_file_list with
              | []              -> invalid_argf "Cmdliner lied!"
              | [fastq]         -> Pf.across_fastq conf
                                      ?number_of_reads ~specific_reads ~nprocs
                                      fastq state
              | [read1; read2]  -> Pf.across_paired conf
                                      ?number_of_reads ~specific_reads ~nprocs
                                      read1 read2 state
              | lst             -> invalid_argf "More than 2, %d fastq files specified!"
                                    (List.length lst)
              end


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
    read_length_override
    not_check_rc
  (* band logic *)
    band_warmup_arg
    band_number_arg
    band_radius_arg
  (* outputting logic. *)
    allele_depth
    likelihood_report_size
    zygosity_report_size
    zygosity_non_zero_value
    per_reads_report_size
    output_format
  (* how are we typing *)
    split
    mode
    not_prealigned
    forward_accuracy_opt
    number_processes_opt
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
      | Ok need_read_length ->
        (* Rename some arguments to clarify logic. *)
        let band =
          let open Option in
          band_warmup_arg >>= fun warmup ->
            band_number_arg >>= fun number ->
              band_radius_arg >>= fun radius ->
                Some { ParPHMM.warmup; number; radius }
        in
        let check_rc = not not_check_rc in
        let prealigned_transition_model = not not_prealigned in
        let past_threshold_filter = not do_not_past_threshold_filter in
        let finish_singles = not do_not_finish_singles in
        let output_opt =
          { ParPHMM_drivers.allele_depth
          ; output_format
          ; depth =
            { ParPHMM_drivers.Output.num_likelihoods = likelihood_report_size
            ; num_zygosities         = Common_options.to_num_zygosities
                                        ~zygosity_non_zero_value
                                        ~zygosity_report_size
            ; num_per_read           = per_reads_report_size
            }
          }
        in
        let conf =
          ParPHMM_drivers.single_conf ?allele ~insert_p ?band ?split
            ?max_number_mismatches ~prealigned_transition_model
            ~past_threshold_filter ~check_rc ~output_opt ()
        in
        match mode with
        | `Viterbi ->
            viterbi read_length_override need_read_length conf fastq_file_list
              number_of_reads specific_reads finish_singles
        | `Forward ->
            begin match number_processes_opt with
            | None  ->
                forward read_length_override need_read_length conf fastq_file_list
                  number_of_reads specific_reads finish_singles
            | Some n ->
                p_forward read_length_override need_read_length conf fastq_file_list
                  number_of_reads specific_reads finish_singles n
            end

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
      sprintf "Browse and report new issues at https://github.com/hammerlab/%s"
        repo
    in
    let man =
      [ `S "BUGS"
      ; `P bug
      ; `S "AUTHORS"
      ; `P "Leonid Rozenberg <leonidr@gmail.com>"
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
            $ band_warmup_arg
            $ number_bands_arg
            $ band_radius_arg
            $ allele_depth_arg
            $ likelihood_report_size_arg
            $ zygosity_report_size_arg
            $ zygosity_non_zero_value_arg
            $ per_reads_report_size_arg
            $ output_format_flag
            (* How are we typing *)
            $ split_arg
            $ mode_flag
            $ not_prealigned_flag
            $ forward_pass_accuracy_arg
            $ number_processes_arg
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
