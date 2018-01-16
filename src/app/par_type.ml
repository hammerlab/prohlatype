(* Typing via a Parametric PHMM. *)
open Prohlatype
open Cmdline_options
open Cmdliner

let app_name = "par_type"

let to_read_size_dependent
  (* Allele information source *)
    file_or_prefix
    ~distance
  (* Allele selectors *)
    ~regex_list
    ~specific_list
    ~number_alleles
    ~do_not_ignore_suffixed_alleles
  (* Cache logic. *)
    ~skip_disk_cache
    =
    let selectors =
      aggregate_selectors ~regex_list ~specific_list ~number_alleles
        ~do_not_ignore_suffixed_alleles
    in
    let input = file_or_prefix_to_allele_input ~distance ~selectors file_or_prefix in
    begin fun read_size ->
      let par_phmm_args = Cache.par_phmm_args ~input ~read_size in
      Cache.par_phmm ~skip_disk_cache par_phmm_args
    end

module Pd = ParPHMM_drivers

(*Wrap Execution *)
module We (W : Pd.Worker) = struct

  module E = Pd.Sequential(W)

  let f ~log_oc ~data_oc read_length_override need_read_length conf
    number_of_reads specific_reads finish_singles fastq_file_list =
    let init =
      match read_length_override with
      | None   -> `Setup need_read_length
      | Some r -> `Set (E.init log_oc need_read_length conf r)
    in
    at_most_two_fastqs
      ~single:(E.across_fastq ~log_oc ~data_oc conf
                  ?number_of_reads ~specific_reads init)
      ~paired:(E.across_paired ~log_oc ~data_oc ~finish_singles conf
                  ?number_of_reads ~specific_reads init)
      fastq_file_list

end (* We *)

module Wef = We(Pd.Forward)
module Wev = We(Pd.Viterbi)

module Pf = Pd.Parallel(Pd.Forward)

let p_forward ~log_oc ~data_oc read_length_override need_read_length conf
  number_of_reads specific_reads finish_singles nprocs fastq_file_list =
  match read_length_override with
  | None   -> errored Term.exit_status_cli_error
                  "Must specify read size for parallel mode."
  | Some r -> let state = Pf.init log_oc need_read_length conf r in
              at_most_two_fastqs
                ~single:(Pf.across_fastq ~log_oc ~data_oc conf
                          ?number_of_reads ~specific_reads ~nprocs state)
                ~paired:(Pf.across_paired ~log_oc ~data_oc conf
                          ?number_of_reads ~specific_reads ~nprocs state)
                fastq_file_list

let type_
  (* Allele information source *)
    file_or_prefix
    distance
  (* Allele selectors *)
    regex_list specific_list number_alleles
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
    do_not_check_rc
  (* outputting logic. *)
    allele_depth
    likelihood_report_size
    zygosity_report_size
    zygosity_non_zero_value
    per_reads_report_size
    output_format
    output
  (* how are we typing *)
    split
    viterbi
    not_prealigned
    forward_accuracy_opt
    number_processes_opt
    =
  Option.value_map forward_accuracy_opt ~default:()
    ~f:(fun fa -> ParPHMM.dx := fa);
  let log_oc, data_oc = setup_oc output output_format in
  let need_read_length =
    to_read_size_dependent
      file_or_prefix
      ~distance
      ~regex_list ~specific_list ~number_alleles
      ~do_not_ignore_suffixed_alleles
      ~skip_disk_cache
  in
  (* Rename some arguments to clarify logic. *)
  let check_rc = not do_not_check_rc in
  let prealigned_transition_model = not not_prealigned in
  let past_threshold_filter = not do_not_past_threshold_filter in
  let finish_singles = not do_not_finish_singles in
  let output_opt =
    { ParPHMM_drivers.allele_depth
    ; output_format
    ; depth =
      { ParPHMM_drivers.Output.num_likelihoods = likelihood_report_size
      ; num_zygosities         = to_num_zygosities
                                  ~zygosity_non_zero_value
                                  ~zygosity_report_size
      ; num_per_read           = per_reads_report_size
      }
    }
  in
  let conf =
    ParPHMM_drivers.single_conf ?allele ~insert_p ?split
      ?max_number_mismatches ~prealigned_transition_model
      ~past_threshold_filter ~check_rc ~output_opt ()
  in
  if viterbi then
    Wev.f ~log_oc ~data_oc read_length_override need_read_length conf
      number_of_reads specific_reads finish_singles fastq_file_list
  else (* forward *)
    match number_processes_opt with
    | None  ->
        Wef.f ~log_oc ~data_oc  read_length_override need_read_length
          conf number_of_reads specific_reads finish_singles
          fastq_file_list
    | Some n ->
        p_forward ~log_oc ~data_oc  read_length_override need_read_length
          conf number_of_reads specific_reads finish_singles n
          fastq_file_list

let () =
  let do_not_check_rc_flag =
    let doc = "Do not check the reverse complement alignment of a read. \
               This option is mostly relevant when investigating a specific \
               reads behavior." in
    Arg.(value & flag & info ~doc ["do-not-check-reverse-complement"])
  in
  let specific_allele_argument = "allele" in
  let spec_allele_arg =
    let docv = "ALLELE" in
    let doc  = "Use a faster mode where we measure the likelihood for just \
                the passed allele. The allele must be found in the alignment \
                file or merged set. " in
    Arg.(value & opt (some string) None
               & info ~doc ~docv [specific_allele_argument])
  in
  let viterbi_flag_str = "viterbi" in
  let viterbi_flag =
    let doc = sprintf
      "Map each read into the viterbi decoded path against a specified allele. \
       Specify the allele with --%s argument, defaults to using the reference \
       of the loci."
      specific_allele_argument
    in
    Arg.(value & flag & info ~doc [viterbi_flag_str])
  in
  let type_ =
    let doc = "Use a Parametric Profile Hidden Markov Model of HLA alleles to \
               type FASTQ samples." in
    let description =
      [ `P "$(tname) uses sequence information found in FASTQ files to build a \
            posterior distribution over a given HLA gene's alleles, thereby \
            providing typing information."

      ; `P "We use a Profile Hidden Markov Model's (PHMM) forward pass's \
            emission as a likelihood function in Bayes theorem to compute the \
            posterior. Since computing such a forward for $(b,each) allele \
            would be computationally infeasible, we parameterize the \
            internal transitions such that we are able to compute the \
            emission for $(b,all) alleles with one pass."

      ; `P (sprintf
          "This tool is purposefully limited (as opposed to $(b,multi_par)) \
          to working with one gene. As has been discussed in the preprint, \
          $(i,Prohlatype: A Probabilistic Framework for HLA Typing), this \
          is insufficient for general typing. \
          Consequently, this tool exists mostly for debugging and \
          analytical purposes. \
          It has options to perform a forward pass for a single allele, \
          without the parameterization (--%s), and it has logic to perform \
          a Viterbi pass that allows one to see the most likely alignment \
          of a read to an allele (--%s)."
          specific_allele_argument
          viterbi_flag_str)

      ; `P (sprintf
            "In the default mode, $(tname) will aggregate read data found \
            in the passed FASTQ files into diploid zygosity and \
            per allele likelihoods. The tool will report only a subset of all \
            possible zygosity likelihoods (control the number with --%s). \
            Furthermore, for each read $(tname) will report the most likely \
            alleles, their likely emission and the position with the largest \
            emission. \
            The --%s argument can change the number of elements that \
            are reported for each read."
              zygosity_report_size_argument
              map_depth_argument)

      ; `P (allele_selector_paragraph "one could PHMMs")
      ]
    in
    let examples =
      [ `P "Type samples based on merged HLA-A:"
      ; `Pre (sprintf "%s path-to-IMGTHLA/alignments/A samples.fastq" app_name)

      ;  `P "Type samples based on cDNA HLA-A and paired reads:"
      ; `Pre (sprintf
          "%s path-to-IMGTHLA/alignments/A_nuc samples_1.fastq samples_2.fastq"
          app_name)

      ; `P "Compute the likelihood of one read versus one allele:"
      ; `Pre (sprintf
          "%s --%s A*30:106 --sr name-of-read \
           path-to-IMGTHLA/alignments/A_gen.txt samples.fastq"
          app_name
          specific_allele_argument)

      ; `P "See the viterbi path that read took:"
      ; `Pre (sprintf
          "%s --%s A*30:106 --sr name-of-read --%s \
           path-to-IMGTHLA/alignments/A_gen.txt samples.fastq"
          app_name
          specific_allele_argument
          viterbi_flag_str)
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
    Term.(const type_
            (* Allele information source *)
            $ (file_or_prefix_arg ~action:"use as reference")
            $ defaulting_distance_flag
            (* Allele selectors *)
            $ regex_arg $ specific_arg $ number_alleles_arg
            $ do_not_ignore_suffixed_alleles_flag
            $ spec_allele_arg
            (* What to do ? *)
            $ (no_cache_flag "PHMM")
            (* What are we typing *)
            $ fastq_file_arg
            $ num_reads_arg
            $ specific_read_args
            $ do_not_finish_singles_flag
            (* options. *)
            $ insert_probability_arg
            $ do_not_past_threshold_filter_flag
            $ max_number_mismatches_arg
            $ read_size_override_arg
            $ do_not_check_rc_flag
            (*$ band_warmup_arg
            $ number_bands_arg
            $ band_radius_arg *)
            $ allele_depth_arg
            $ likelihood_report_size_arg
            $ zygosity_report_size_arg
            $ zygosity_non_zero_value_arg
            $ per_reads_report_size_arg
            $ output_format_flag
            $ output_arg
            (* How are we typing *)
            $ split_arg
            $ viterbi_flag
            $ not_prealigned_flag
            $ forward_pass_accuracy_arg
            $ number_of_processors_arg
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
        , info app_name
             ~version
             ~doc
             ~man
             ~exits:default_exits)
  in
  Term.(exit_status (eval type_))
