(* Typing via a Parametric PHMM. *)
open Util

let app_name = "multi_par"

let to_read_size_dependent
  (* Allele information source *)
  ~alignment_files ~merge_files ~distance
  (* Cache management *)
  ~skip_disk_cache =
  Common_options.to_allele_inputs ~alignment_files ~merge_files ~distance
    ~selectors:[] (* Selectors NOT supported, on purpose! *)
    >>= fun inputs ->
      Ok (fun read_size ->
          List.map inputs ~f:(fun input ->
            let par_phmm_arg = Cache.par_phmm_args ~input ~read_size in
            Cache.par_phmm ~skip_disk_cache par_phmm_arg))

module Pd = ParPHMM_drivers

module Sequential = Pd.Sequential(Pd.Multiple_loci)

module Parallel = Pd.Parallel(Pd.Multiple_loci)

let type_
  (* Allele information source *)
    class1
    full_class1
    gen_nuc_mgd
    alignment_files
    merge_files
    distance
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
    output
  (* how are we typing *)
    split
    not_prealigned
    forward_accuracy_opt
    number_processes_opt
    =
  Option.value_map forward_accuracy_opt ~default:()
    ~f:(fun fa -> ParPHMM.dx := fa);
  let log_oc, data_oc = Common_options.setup_oc output output_format in
  let band =
    let open Option in
      band_warmup_arg >>= fun warmup ->
        band_number_arg >>= fun number ->
          band_radius_arg >>= fun radius ->
            Some { ParPHMM.warmup; number; radius }
  in
  let alignment_files, merge_files =
    Common_options.class_selectors class1 full_class1 gen_nuc_mgd
      alignment_files merge_files
  in
  let past_threshold_filter = not do_not_past_threshold_filter in
  let prealigned_transition_model = not not_prealigned in
  let finish_singles = not do_not_finish_singles in
  let output_opt =
    { Pd.allele_depth
    ; output_format
    ; depth =
        { Pd.Output.num_likelihoods = likelihood_report_size
        ; num_zygosities            = Common_options.to_num_zygosities
                                        ~zygosity_non_zero_value
                                        ~zygosity_report_size
        ; num_per_read              = per_reads_report_size
        }
    }
  in
  let conf =
    Pd.multiple_conf ~insert_p ?band ?max_number_mismatches ?split
      ~prealigned_transition_model ~past_threshold_filter
      ~output_opt
      ()
  in
  let need_read_size_r =
    to_read_size_dependent
      ~alignment_files ~merge_files ~distance ~skip_disk_cache
  in
  match need_read_size_r with
  | Error e           -> eprintf "%s" e
  | Ok need_read_size ->
      begin
        match number_processes_opt with
        | None  ->
            let init =
              match read_size_override with
              | None   -> `Setup need_read_size
              | Some r -> `Set (Sequential.init log_oc need_read_size conf r)
            in
            begin match fastq_file_lst with
            | []              -> invalid_argf "Cmdliner lied!"
            | [fastq]         -> Sequential.across_fastq ~log_oc ~data_oc conf
                                    ?number_of_reads ~specific_reads fastq init
            | [read1; read2]  -> Sequential.across_paired ~log_oc ~data_oc ~finish_singles conf
                                    ?number_of_reads ~specific_reads read1 read2 init
            | lst             -> invalid_argf "More than 2, %d fastq files specified!"
                                  (List.length lst)
            end
        | Some nprocs ->
            let r = Option.value_exn read_size_override
                      ~msg:"Must specify read size override in parallel mode"
            in
            let state = Parallel.init log_oc need_read_size conf r in
            begin match fastq_file_lst with
            | []              -> invalid_argf "Cmdliner lied!"
            | [fastq]         -> Parallel.across_fastq ~log_oc ~data_oc conf
                                    ?number_of_reads ~specific_reads ~nprocs
                                    fastq state
            | [read1; read2]  -> Parallel.across_paired ~log_oc ~data_oc conf
                                    ?number_of_reads ~specific_reads ~nprocs
                                    read1 read2 state
            | lst             -> invalid_argf "More than 2, %d fastq files specified!"
                                  (List.length lst)
            end
      end

let () =
  let open Cmdliner in
  let open Common_options in
  let alignments_arg =
    let docv = "FILE" in
    let doc  = "File to lookup IMGT allele alignments. The alleles found in this \
                file will initially define the set of alleles to be used for \
                named locus. Supply multiple file (or merge) arguments to \
                compute the likelihood of a read against each locus and add to \
                the most likely (highest likelihood)." in
    Arg.(value & opt_all file [] & info ~doc ~docv ["alignments"])
  in
  let max_number_mismatches_arg =
    let docv = "POSITIVE INTEGER" in
    let doc = "Setup a filter on the reads to cancel evaluation once we've \
               seen this many mismatches." in
    Arg.(value & opt (some positive_int) None & info ~doc ~docv
          ["max-mismatches"])
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
            $ class1_directory_arg
            $ full_class1_directory_arg
            $ gen_nuc_merged_flag
            $ alignments_arg $ merges_arg $ defaulting_distance_flag
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
            $ band_warmup_arg
            $ number_bands_arg
            $ band_radius_arg
            (* output logic. *)
            $ allele_depth_arg
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
            $ number_processes_arg
         , info app_name ~version ~doc ~man)
  in
  match Term.eval type_ with
  | `Ok ()           -> exit 0
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
