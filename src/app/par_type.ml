(* Typing via a Parametric PHMM. *)
open Util

let app_name = "par_type"

let time s f =
  let n = Sys.time () in
  let r = f () in
  Printf.printf "%s total running time in seconds: %f\n%!" s (Sys.time () -. n);
  r

let to_read_size_dependent
  (* Allele information source *)
  ?alignment_file ?merge_file ~distance ~impute
  (* Allele selectors *)
    ~regex_list
    ~specific_list
    ~without_list
    ?number_alleles
    ~do_not_ignore_suffixed_alleles
    ~skip_disk_cache
    =
    Common_options.to_input ?alignment_file ?merge_file ~distance ~impute ()
      >>= fun input ->
        let selectors =
          Common_options.aggregate_selectors ~regex_list ~specific_list
            ~without_list ?number_alleles ~do_not_ignore_suffixed_alleles
        in
        Ok (fun read_size ->
            let par_phmm_args = Cache.par_phmm_args ~input ~selectors ~read_size in
            Cache.par_phmm ~skip_disk_cache par_phmm_args)

exception PPE of string

let add_log_likelihoods n =
  let open ParPHMM in
  fun ~into ->
    function
    | Filtered _   -> ()
    | Completed rs ->
      for i = 0 to n - 1 do into.(i) <- into.(i) +. rs.likelihood.(i) done

let to_reduce_update_f perform_forward_pass add_ll update_me fqi =
  time (sprintf "updating on %s" fqi.Biocaml_unix.Fastq.name) (fun () ->
    let open Core_kernel.Std in
    match Fastq.phred_log_probs fqi.Biocaml_unix.Fastq.qualities with
    | Result.Error e       -> raise (PPE (Error.to_string_hum e))
    | Result.Ok read_probs ->
        let nl = perform_forward_pass fqi.Biocaml_unix.Fastq.sequence read_probs in
        add_ll ~into:update_me nl;
        update_me)

let to_map_update_f f acc fqi =
  time (sprintf "updating on %s" fqi.Biocaml_unix.Fastq.name) (fun () ->
    let open Core_kernel.Std in
    match Fastq.phred_log_probs fqi.Biocaml_unix.Fastq.qualities with
    | Result.Error e       -> raise (PPE (Error.to_string_hum e))
    | Result.Ok read_probs ->
        let t = f fqi.Biocaml_unix.Fastq.sequence read_probs in
        (fqi.Biocaml_unix.Fastq.name, t) :: acc)

type 'a g =
  { f         : 'a -> Biocaml_unix.Fastq.item -> 'a
  ; mutable s : 'a
  ; fin       : 'a -> unit
  }

let proc_g = function
  | `Mapper g -> begin fun fqi ->
                  g.s <- g.f g.s fqi;
                  `Mapper g
                 end
  | `Reducer g -> begin fun fqi ->
                  g.s <- g.f g.s fqi;
                  `Reducer g
                 end

let output_mapper lst =
  printf "Reads: %d\n" (List.length lst);
  List.sort lst ~cmp:(fun (_n1, ms1) (_n2, ms2) ->
    ParPHMM.(compare (best_stat ms1) (best_stat ms2)))
  |> List.rev
  |> List.iter ~f:(fun (n, s) ->
    printf "%s\t%s\n" n (ParPHMM.mapped_stats_to_string ~sep:'\t' s))

let to_set ?insert_p ?max_number_mismatches ~past_threshold_filter
  mode specific_allele ~check_rc ?band rp read_size =
  let pt =
    time (sprintf "Setting up ParPHMM transitions with %d read_size" read_size)
      (fun () -> rp read_size)
  in
  let g =
    match specific_allele with
    | None ->
      let add_ll = add_log_likelihoods pt.ParPHMM.number_alleles in
      time (sprintf "Allocating forward pass workspaces")
        (fun () ->
          let r =
            ParPHMM.forward_pass mode ?insert_p ?max_number_mismatches
              ~past_threshold_filter ?band pt read_size
          in
          match r with
          | `Reducer (u, f) ->
              `Reducer
                { f   = to_reduce_update_f (f ~check_rc) add_ll
                ; s   = u
                ; fin = ParPHMM.output (`Pt pt)
                }
          | `Mapper m ->
              `Mapper
                { f   = to_map_update_f m
                ; s   = []
                ; fin = output_mapper
                })
    | Some allele ->
        let add_ll = add_log_likelihoods 1 in
        time (sprintf "Allocating forward pass workspaces")
          (fun () ->
            let r =
              ParPHMM.single_allele_forward_pass ?insert_p ?max_number_mismatches
              mode pt read_size allele
            in
            match r with
            | `Reducer (u, f) ->
                `Reducer
                  { f   = to_reduce_update_f (f ~check_rc) add_ll
                  ; s   = u
                  ; fin = ParPHMM.output (`Allele allele)
                  }
            | `Mapper m ->
                `Mapper
                  { f   = to_map_update_f m
                  ; s   = []
                  ; fin = output_mapper
                  })
  in
  `Set g

let fin = function
  | `Setup _ -> eprintf "Didn't fine any reads."
  | `Set (`Mapper g)  -> g.fin g.s
  | `Set (`Reducer g) -> g.fin g.s

let across_fastq ?insert_p ?max_number_mismatches ~past_threshold_filter
    ?number_of_reads ~specific_reads ~check_rc
    specific_allele ?band mode file init =
  try
    Fastq.fold ?number_of_reads ~specific_reads file ~init
      ~f:(fun acc fqi ->
            match acc with
            | `Setup rp ->
                let read_size = String.length fqi.Biocaml_unix.Fastq.sequence in
                let `Set g = to_set ?insert_p ?max_number_mismatches
                                ~past_threshold_filter
                                mode specific_allele ~check_rc ?band rp read_size in
                `Set (proc_g g fqi)
            | `Set g ->
                `Set (proc_g g fqi))
    |> fin
  with PPE e ->
    eprintf "%s" e

let type_
  (* Allele information source *)
    alignment_file merge_file distance not_impute
  (* Allele selectors *)
    regex_list specific_list without_list number_alleles
    do_not_ignore_suffixed_alleles
    specific_allele
  (* Process *)
    skip_disk_cache
  (* What to do? *)
    fastq_file_lst number_of_reads specific_reads
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
  (* how are we typing *)
    map
    forward_accuracy_opt
    =
  Option.value_map forward_accuracy_opt ~default:()
    ~f:(fun fa -> ParPHMM.dx := fa);
  let impute   = not not_impute in
  let check_rc = not not_check_rc in
  let band     =
    if not_band then
      None
    else
      Some { ParPHMM.warmup; number; width }
  in
  let need_read_size_r =
    to_read_size_dependent
      ?alignment_file ?merge_file ~distance ~impute
      ~regex_list ~specific_list ~without_list ?number_alleles
      ~do_not_ignore_suffixed_alleles
      ~skip_disk_cache
  in
  let past_threshold_filter = not do_not_past_threshold_filter in
  match need_read_size_r with
  | Error e           -> eprintf "%s" e
  | Ok need_read_size ->
    let mode = if map then `Mapper else `Reducer in
    let init =
      match read_size_override with
      | None   -> `Setup need_read_size
      | Some r -> to_set ~insert_p ?max_number_mismatches
                    ~past_threshold_filter
                    mode specific_allele ~check_rc ?band need_read_size r
    in
    begin match fastq_file_lst with
    | []              -> invalid_argf "Cmdliner lied!"
    | [read1; read2]  -> invalid_argf "implement pairs!"
    | [fastq]         -> across_fastq ~insert_p ?max_number_mismatches
                          ~past_threshold_filter
                          ?number_of_reads ~specific_reads ~check_rc ?band
                          specific_allele mode fastq init
    | lst             -> invalid_argf "More than 2, %d fastq files specified!" (List.length lst)
    end

let () =
  let open Cmdliner in
  let open Common_options in
  let not_check_rc_flag =
    let doc = "Do not check the reverse complement." in
    Arg.(value & flag & info ~doc ["not-check-rc"])
  in
  let do_not_impute_flag =
    let doc  = "Do NOT fill in the missing segments of alleles with an \
                iterative algorithm that picks the closest allele with full \
                length. The default behavior is to impute the sequences, \
                so that (1) the per-allele likelihoods are comparable \
                (the default behavior of unknown sequences biases reads \
                to map correctly in that region) and (2) this limits the \
                actual branching of transitions inside the PPHMM."
    in
    Arg.(value & flag & info ~doc ["do-not-impute"])
  in
  let spec_allele_arg =
    let docv = "ALLELE" in
    let doc  = "Use a faster mode where we measure the likelihood for just \
                the passed allele. The allele must be found in the alignment \
                or merge file." in
    Arg.(value & opt (some string) None
               & info ~doc ~docv ["allele"])
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
            $ file_arg $ merge_arg $ distance_flag $ do_not_impute_flag
            (* Allele selectors *)
            $ regex_arg $ allele_arg $ without_arg $ num_alt_arg
            $ do_not_ignore_suffixed_alleles_flag
            $ spec_allele_arg
            (* What to do ? *)
            $ no_cache_flag
            (* What are we typing *)
            $ fastq_file_arg $ num_reads_arg $ specific_read_args
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
            (* How are we typing *)
            $ map_flag
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
