(* Report viterbi path of PHMM. *)
open Util

let app_name = "viterbi"

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

let to_read_proc perform_viterbi_pass acc fqi =
  time (sprintf "updating on %s" fqi.Biocaml_unix.Fastq.name)
    (fun () ->
      let open Core_kernel.Std in
      match Fastq.phred_log_probs fqi.Biocaml_unix.Fastq.qualities with
      | Result.Error e       -> raise (PPE (Error.to_string_hum e))
      | Result.Ok read_probs ->
          let res = perform_viterbi_pass fqi.Biocaml_unix.Fastq.sequence read_probs in
          (fqi.Biocaml_unix.Fastq.name, res)  :: acc)

(* This is unnecessarily parameterized at the moment. *)
type 'a comp =
  { f         : 'a -> Biocaml_unix.Fastq.item -> 'a
  ; mutable s : 'a
  ; fin       : 'a -> unit
  }

let apply c fqi =
  c.s <- c.f c.s fqi;
  c

let write_paths_to_stdout =
  let open ParPHMM in
  List.iter ~f:(fun (read_name, (rc, emission, path_list)) ->
    let { reference; read } = path_list_to_strings path_list in
    printf "%s:\t %b: %f %d\n%s\n"
      read_name rc emission (List.length path_list)
        (manual_comp_display ~width:120 ~labels:("ref ", "read ")
            reference read))

let to_comp ?insert_p allele rp read_size =
  let pt =
    time (sprintf "Setting up ParPHMM transitions with %d read_size" read_size)
      (fun () -> rp read_size)
  in
  time (sprintf "Allocating viterbi pass workspaces")
    (fun () ->
      let f = ParPHMM.setup_single_allele_viterbi_pass ?insert_p pt read_size
                ~allele in
      { f   = to_read_proc f
      ; s   = []
      ; fin = write_paths_to_stdout
      })

let fin = function
  | `Setup _  -> eprintf "Didn't find any reads."
  | `Set c    -> c.fin c.s

let across_fastq ?insert_p ?number_of_reads ~specific_reads allele file init =
  try
    Fastq.fold ?number_of_reads ~specific_reads file ~init
      ~f:(fun acc fqi ->
            match acc with
            | `Setup rp ->
                let read_size = String.length fqi.Biocaml_unix.Fastq.sequence in
                let c = to_comp ?insert_p allele rp read_size in
                `Set (apply c fqi)
            | `Set c ->
                `Set (apply c fqi))
    |> fin
  with PPE e ->
    eprintf "%s" e

let viterbi
  (* Allele information source *)
    alignment_file merge_file distance not_impute
  (* Allele selectors *)
    regex_list specific_list without_list number_alleles
    do_not_ignore_suffixed_alleles
  (* Which allele ? *)
    specific_allele
  (* Process *)
    skip_disk_cache
  (* What to do? *)
    fastq_file_lst number_of_reads specific_reads
  (* options *)
    insert_p
    read_size_override
  (* how are we typing *)
    forward_accuracy_opt
    =
  Option.value_map forward_accuracy_opt ~default:()
    ~f:(fun fa -> ParPHMM.dx := fa);
  let impute   = not not_impute in
  let need_read_size_r =
    to_read_size_dependent
      ?alignment_file ?merge_file ~distance ~impute
      ~regex_list ~specific_list ~without_list ?number_alleles
      ~do_not_ignore_suffixed_alleles
      ~skip_disk_cache
  in
  match need_read_size_r with
  | Error e           -> eprintf "%s" e
  | Ok need_read_size ->
    let init =
      match read_size_override with
      | None   -> `Setup need_read_size
      | Some r -> `Set (to_comp ~insert_p specific_allele need_read_size r)
    in
    begin match fastq_file_lst with
    | []              -> invalid_argf "Cmdliner lied!"
    | [read1; read2]  -> invalid_argf "implement pairs!"
    | [fastq]         -> across_fastq ~insert_p ?number_of_reads ~specific_reads 
                          specific_allele fastq init
    | lst             -> invalid_argf "More than 2, %d fastq files specified!" (List.length lst)
    end

let () =
  let open Cmdliner in
  let open Common_options in
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
    let doc  = "What allele to report viterbi path" in
    Arg.(required & pos 0 (some string) None 
                  & info ~doc ~docv [])
  in
  (* Different than the one from Common_options. *)
  let fastq_file_arg =
    let docv = "FASTQ FILE" in
    let doc = "Fastq formatted DNA reads file, only one file per sample. \
               List paired end reads as 2 sequential files." in
    Arg.(non_empty & pos_right 0 file [] & info ~doc ~docv [])
  in
  let viterbi =
    let version = "0.0.0" in
    let doc = "Use a Parametric Profile Hidden Markov Model of HLA allele to \
               print out the Viterbi decoding of reads to HLA loci." in
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
    Term.(const viterbi
            (* Allele information source *)
            $ file_arg $ merge_arg $ distance_flag $ do_not_impute_flag
            (* Allele selectors *)
            $ regex_arg $ allele_arg $ without_arg $ num_alt_arg
            $ do_not_ignore_suffixed_alleles_flag
            (* Which allele? *)
            $ spec_allele_arg
            (* What to do ? *)
            $ no_cache_flag
            (* What are we typing *)
            $ fastq_file_arg $ num_reads_arg $ specific_read_args
            (* options. *)
            $ insert_probability_arg
            $ read_size_override_arg
            (* How are we typing *)
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
  match Term.eval viterbi with
  | `Ok ()           -> exit 0
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
