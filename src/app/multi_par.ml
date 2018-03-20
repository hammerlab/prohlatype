(* Typing via a Parametric PHMM. *)
open Prohlatype
open Cmdline_options
open Cmdliner

let app_name = "multi_par"

let to_read_size_dependent
  (* Allele information source *)
    file_prefix_or_directory
    only_class1
    gen_nuc_mgd
    ~distance
  (* Cache logic. *)
    ~skip_disk_cache
    =
    let inputs =
      (* We do NOT accept selectors. *)
      file_prefix_or_directory_to_allele_input_list ~distance ~selectors:[]
        file_prefix_or_directory only_class1 gen_nuc_mgd
    in
    begin fun read_size ->
      List.map inputs ~f:(fun input ->
        let par_phmm_args = Cache.par_phmm_args ~input ~read_size in
        match Cache.par_phmm ~skip_disk_cache par_phmm_args with
        | Ok phmm -> phmm
        | Error e -> raise (Phmm_construction_error e))
    end

module Pd = ParPHMM_drivers

module Sequential = Pd.Sequential(Pd.Multiple_loci)

module Parallel = Pd.Parallel(Pd.Multiple_loci)

let type_
  (* Input *)
    file_prefix_or_directory
    only_class1
    gen_nuc_mgd
  (* Process *)
    distance
    skip_disk_cache
  (* What to do? *)
    fastq_file_lst
    number_of_reads
    specific_reads
    do_not_finish_singles
  (* options *)
    number_processes_opt
    insert_p
    do_not_past_threshold_filter
    max_number_mismatches
    read_size_override
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
    not_prealigned
    forward_accuracy_opt
    =
  let commandline = String.concat ~sep:" " (Array.to_list Sys.argv) in
  Option.value_map forward_accuracy_opt ~default:()
    ~f:(fun fa -> Probability.dx := fa);
  let log_oc, data_oc = setup_oc output output_format in
  let past_threshold_filter = not do_not_past_threshold_filter in
  let prealigned_transition_model = not not_prealigned in
  let finish_singles = not do_not_finish_singles in
  let output_opt =
    { Pd.allele_depth
    ; output_format
    ; depth =
        { Pd.Output.num_likelihoods = likelihood_report_size
        ; num_zygosities            = to_num_zygosities
                                        ~zygosity_non_zero_value
                                        ~zygosity_report_size
        ; num_per_read              = per_reads_report_size
        }
    }
  in
  let conf =
    Pd.multiple_conf ~insert_p ?max_number_mismatches ?split
      ~prealigned_transition_model ~past_threshold_filter
      ~output_opt
      commandline
  in
  try
    let need_read_size =
      to_read_size_dependent ~distance ~skip_disk_cache
        file_prefix_or_directory only_class1 gen_nuc_mgd
    in
    match number_processes_opt with
    | None  ->
        let init =
          match read_size_override with
          | None   -> `Setup need_read_size
          | Some r -> `Set (Sequential.init log_oc need_read_size conf r)
        in
        at_most_two_fastqs fastq_file_lst
          ~single:(Sequential.across_fastq ~log_oc ~data_oc conf
                    ?number_of_reads ~specific_reads init)
          ~paired:(Sequential.across_paired ~log_oc ~data_oc ~finish_singles conf
                    ?number_of_reads ~specific_reads init)
    | Some nprocs ->
        begin match read_size_override with
        | None -> errored Term.exit_status_cli_error
                    "Must specify read size override in parallel mode."
        | Some r ->
            let state = Parallel.init log_oc need_read_size conf r in
            at_most_two_fastqs fastq_file_lst
              ~single:(Parallel.across_fastq ~log_oc ~data_oc conf
                        ?number_of_reads ~specific_reads ~nprocs state)
              ~paired:(Parallel.across_paired ~log_oc ~data_oc conf
                        ?number_of_reads ~specific_reads ~nprocs state)
        end
  with (Phmm_construction_error e) ->
    errored phmm_construction_error "%s" e

let () =
  let type_ =
    let doc = "Use a Parametric Profile Hidden Markov Model of HLA allele to \
               type fastq samples." in
    let description =
      [ `P "$(tname) uses sequence information found in FASTQ files to build a \
            posterior distribution over HLA gene's alleles, thereby \
            providing typing information."

      ; `P "We use a Profile Hidden Markov Model's (PHMM) forward pass's \
            emission as a likelihood function in Bayes theorem to compute the \
            posterior. Since computing such a forward for $(b,each) allele \
            would be computationally infeasible, we parameterize the \
            internal transitions such that we are able to compute the \
            emission for $(b,all) alleles with one pass."

      ; `P "This tool (as opposed to $(b,par_type)) \
            designed to work with multiple gene, as discussed in the \
            preprint: $(i,Prohlatype: A Probabilistic Framework for HLA \
            Typing), and is intended for robust typing."

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
              per_reads_report_size_argument)

      ]
    in
    let examples =
      [ `P "Type samples based on merged HLA-A:"
      ; `Pre (sprintf "%s path-to-IMGTHLA/alignments/A sample.fastq" app_name)

      ;  `P "Type samples based on cDNA HLA-A and paired reads:"
      ; `Pre (sprintf
          "%s path-to-IMGTHLA/alignments/A_nuc sample_1.fastq sample_2.fastq"
          app_name)

      ; `P "The most common usage case:"
      ; `Pre (sprintf
          "%s path-to-IMGTHLA/alignments/ sample_1.fastq sample_2.fastq -o sample_output"
          app_name)

      ; `P "Consider all allele information via merging, for just traditional \
            class I genes and split the passes:"
      ; `Pre (sprintf
          "%s --%s --%s --%s 4 -o sample_output
          path-to-IMGTHLA/alignments/ sample_1.fastq sample_2.fastq"
          merged_argument
          only_class1_flag
          split_argument
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
      ; `P author
      ]
    in
    Term.(const type_
            (* Allele information source *)
            $ (file_prefix_or_directory_arg "use as reference")
            $ only_class1_directory_flag
            $ gen_nuc_merged_flag
            (* What to do ? *)
            $ defaulting_distance_flag
            $ (no_cache_flag "PHMM")
            (* What are we typing *)
            $ fastq_file_arg
            $ num_reads_arg
            $ specific_read_args
            $ do_not_finish_singles_flag
            (* options. *)
            $ number_of_processors_arg
            $ insert_probability_arg
            $ do_not_past_threshold_filter_flag
            $ max_number_mismatches_arg
            $ read_size_override_arg
            (* output logic. *)
            $ allele_depth_info_arg
            $ likelihood_report_size_arg
            $ zygosity_report_size_arg
            $ zygosity_non_zero_value_arg
            $ per_reads_report_size_arg
            $ output_format_flag
            $ output_arg
            (* How are we typing *)
            $ split_arg
            $ not_prealigned_flag
            $ forward_pass_accuracy_arg
         , info app_name
              ~version
              ~doc
              ~man
              ~exits:(phmm_construction_exit_info :: default_exits))
  in
  Term.(exit_status (eval type_))
