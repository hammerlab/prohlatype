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
    ~skip_disk_cache
    =
    Common_options.to_input ?alignment_file ?merge_file ~distance ~impute ()
      >>= fun input ->
        let selectors =
          regex_list @ specific_list @ without_list @
            (match number_alleles with | None -> [] | Some s -> [s])
        in
        Ok (fun read_size ->
            let par_phmm_args = Cache.par_phmm_args ~input ~selectors ~read_size in
            Cache.par_phmm ~skip_disk_cache par_phmm_args)

exception PPE of string

let to_reduce_update_f ~check_rc ~logspace perform_forward_pass update_me fqi =
  time (sprintf "updating on %s" fqi.Biocaml_unix.Fastq.name) (fun () ->
    let open Core_kernel.Std in
    let p = if logspace then Fastq.phred_log_probs else Fastq.phred_probabilities in
    match p fqi.Biocaml_unix.Fastq.qualities with
    | Result.Error e       -> raise (PPE (Error.to_string_hum e))
    | Result.Ok read_probs ->
      perform_forward_pass ~into:update_me fqi.Biocaml_unix.Fastq.sequence read_probs)

let to_map_update_f ~logspace f acc fqi =
  time (sprintf "updating on %s" fqi.Biocaml_unix.Fastq.name) (fun () ->
    let open Core_kernel.Std in
    let p = if logspace then Fastq.phred_log_probs else Fastq.phred_probabilities in
    match p fqi.Biocaml_unix.Fastq.qualities with
    | Result.Error e       -> raise (PPE (Error.to_string_hum e))
    | Result.Ok read_probs ->
        let t = f fqi.Biocaml_unix.Fastq.sequence read_probs in
        (fqi.Biocaml_unix.Fastq.name, t) :: acc)

type 'a g =
  { f         : 'a -> Biocaml_unix.Fastq.item -> 'a
  ; fin       : 'a -> unit
  ; mutable s : 'a
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


let to_set ~map ~check_rc ?band ~logspace rp read_size =
  let pt =
    time (sprintf "Setting up ParPHMM transitions with %d read_size" read_size)
      (fun () -> rp read_size)
  in
  let alleles = Alleles.to_alleles pt.ParPHMM.allele_index in
  let g =
    time (sprintf "Allocating forward pass workspaces")
      (fun () ->
        match ParPHMM.setup ~map ?band ~logspace pt read_size with
        | `Reducer (f, u) ->
            `Reducer
              { f   = to_reduce_update_f ~check_rc ~logspace (f ~check_rc)
              ; s   = u
              ; fin = fun final_likelihoods ->
                  List.mapi alleles ~f:(fun i a -> (final_likelihoods.(i), a))
                  |> List.sort ~cmp:(fun (l1,_) (l2,_) -> compare l2 l1)
                  |> List.iter ~f:(fun (l,a) ->
                      let v = (*if logspace then 10. ** l else*) l in
                      printf "%10s\t%0.20f\n" a v)
              }
        | `Mapper m ->
            `Mapper
              { f   = to_map_update_f ~logspace m
              ; s   = []
              ; fin = begin fun lst ->
                        printf "got %d\n" (List.length lst);
                        List.sort lst ~cmp:(fun (_n1, ms1) (_n2, ms2) ->
                          ParPHMM.(compare (best_stat ms1) (best_stat ms2)))
                        |> List.rev
                        |> List.iter ~f:(fun (n, s) ->
                          printf "%s\t%s\n" n (ParPHMM.mapped_stats_to_string ~sep:'\t' s))
                      end
              })
  in
  `Set g

let fin = function
  | `Setup _ -> eprintf "Didn't fine any reads."
  | `Set (`Mapper g)  -> g.fin g.s
  | `Set (`Reducer g) -> g.fin g.s



let across_fastq ~map ?number_of_reads ~specific_reads ~check_rc ?band ~logspace file init =
  try
    Fastq.fold ?number_of_reads ~specific_reads file ~init
      ~f:(fun acc fqi ->
            match acc with
            | `Setup rp ->
                let read_size = String.length fqi.Biocaml_unix.Fastq.sequence in
                let `Set g = to_set ~map ~check_rc ?band ~logspace rp read_size in
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
  (* Process *)
    skip_disk_cache
  (* What to do? *)
    fastq_file_lst number_of_reads specific_reads
  (* options *)
    read_size_override
    not_logspace
    not_check_rc
    not_band
    start_column
    number
    width
  (* how are we typing *)
    map
    =
  let logspace = not not_logspace in
  let impute   = not not_impute in
  let check_rc = not not_check_rc in
  let band     =
    if not_band then
      None
    else
      Some { ParPHMM.Bands.start_column; number; width }
  in
  let need_read_size_r =
    to_read_size_dependent
      ?alignment_file ?merge_file ~distance ~impute
      ~regex_list ~specific_list ~without_list ?number_alleles
      ~skip_disk_cache
  in
  match need_read_size_r with
  | Error e           -> eprintf "%s" e
  | Ok need_read_size ->
    let init =
      match read_size_override with
      | None   -> `Setup need_read_size
      | Some r -> to_set ~map ~check_rc ?band ~logspace need_read_size r
    in
    begin match fastq_file_lst with
    | []              -> invalid_argf "Cmdliner lied!"
    | [read1; read2]  -> invalid_argf "implement pairs!"
    | [fastq]         -> across_fastq ~map ?number_of_reads ~specific_reads ~check_rc
                            ?band ~logspace fastq init
    | lst             -> invalid_argf "More than 2, %d fastq files specified!" (List.length lst)
    end

let () =
  let open Cmdliner in
  let open Common_options in
  let read_size_override_arg =
    let docv = "POSITIVE INTEGER" in
    let doc = "Override the number of bases to calculate the likelihood over, \
               instead of using the number of bases in the FASTQ." in
    Arg.(value & opt (some positive_int) None & info ~doc ~docv ["read-size"])
  in
  let not_logspace_flag =
    let doc = "Do not perform the calculation in \"log space\"." in
    Arg.(value & flag & info ~doc ["not-logspace"])
  in
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
  (* Band arguments *)
  let not_band_flag =
    let doc = "Calculate the full forward pass matrix." in
    Arg.(value & flag & info ~doc ["do-not-band"])
  in
  let band_start_column_arg =
    let default = ParPHMM.Bands.(default.start_column) in
    let docv = "POSITIVE INTEGER" in
    let doc =
      sprintf "At which column in the forward pass to compute bands instead \
               of the full pass. Defaults to: %d." default
    in
    Arg.(value
          & opt positive_int default
          & info ~doc ~docv ["band-start-column"])
  in
  let number_bands_arg =
    let default = ParPHMM.Bands.(default.number) in
    let docv = "POSITIVE INTEGER" in
    let doc  =
      sprintf "Number of bands to calculate. Defaults to %d" default
    in
    Arg.(value
          & opt positive_int default
          & info ~doc ~docv ["number_bands"])
  in
  let band_width_arg =
    let default = ParPHMM.Bands.(default.width) in
    let docv = "greater than 1" in
    let doc  =
      sprintf "Width of bands to calculate. Defaults to %d. Must be greater than 1." default
    in
    Arg.(value
          & opt greater_than_one default
          & info ~doc ~docv ["band-width"])
  in
  let map_flag =
    let doc = "Map (do not reduce) the typing of the individual reads." in
    let docv = "This switch turns the typing logic into more of a diagnostic \
                mode. The end user can see how individual reads are typed."
    in
    Arg.(value & flag & info ~doc ~docv ["map"])
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
            (* What to do ? *)
            $ no_cache_flag
            (* What are we typing *)
            $ fastq_file_arg $ num_reads_arg $ specific_read_args
            (* options. *)
            $ read_size_override_arg
            $ not_logspace_flag
            $ not_check_rc_flag
            $ not_band_flag
            $ band_start_column_arg
            $ number_bands_arg
            $ band_width_arg
            (* How are we typing *)
            $ map_flag
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
