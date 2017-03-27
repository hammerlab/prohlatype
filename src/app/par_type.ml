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
  ?alignment_file ?merge_file ~distance
  (* Allele selectors *)
    ~regex_list
    ~specific_list
    ~without_list
    ?number_alleles
    ~skip_disk_cache
    =
    Common_options.to_input ?alignment_file ?merge_file ~distance () >>= fun input ->
      let selectors =
        regex_list @ specific_list @ without_list @
          (match number_alleles with | None -> [] | Some s -> [s])
      in
      Ok (fun read_size ->
          let par_phmm_args = Cache.par_phmm_args ~input ~selectors ~read_size in
          Cache.par_phmm ~skip_disk_cache par_phmm_args)

exception PPE of string

let to_update_f final_run_len fp update_me fqi =
  time (sprintf "updating on %s" fqi.Biocaml_unix.Fastq.name) (fun () ->
    let open Core_kernel.Std in
    match Fastq.phred_probabilities fqi.Biocaml_unix.Fastq.qualities with
    | Result.Error e       -> raise (PPE (Error.to_string_hum e))
    | Result.Ok read_probs ->
      let fe = fp fqi.Biocaml_unix.Fastq.sequence read_probs in
      let likelihoods = ParPHMM.to_allele_arr fe final_run_len in
      Array.iteri likelihoods ~f:(fun i l -> update_me.(i) <- update_me.(i) *. l);
      update_me)

let to_set rp read_size =
  let pt =
    time (sprintf "Setting up ParPHMM transitions with %d read_size" read_size)
      (fun () -> rp read_size)
  in
  let alleles = pt.ParPHMM.alleles in
  let n = List.length alleles in
  let update_me = Array.make n 1. in
  let fp = time (sprintf "Allocating forward pass workspace")
    (fun () -> ParPHMM.forward_pass pt.ParPHMM.conf read_size)
  in
  let update = to_update_f pt.ParPHMM.conf.ParPHMM.final_run_len fp in
  `Set (alleles, update, update_me)

let across_fastq ?number_of_reads ~specific_reads file init =
  try
    Fastq.fold ?number_of_reads ~specific_reads file ~init
      ~f:(fun acc fqi ->
            match acc with
            | `Setup rp ->
                let read_size = String.length fqi.Biocaml_unix.Fastq.sequence in
                let `Set (alleles, update, update_me) =
                  to_set rp read_size
                in
                `Set (alleles, update, update update_me fqi)
            | `Set (alelles, update, ret) ->
                `Set (alelles, update, update ret fqi))
    |> function
        | `Setup _ -> eprintf "Didn't fine any reads."
        | `Set (alleles, _, final_likelihoods) ->
            List.mapi alleles ~f:(fun i a -> (final_likelihoods.(i), a))
            |> List.sort ~cmp:(fun (l1,_) (l2,_) -> compare l2 l1)
            |> List.iter ~f:(fun (l,a) -> printf "%10s\t%0.6f\n" a l)
  with PPE e ->
    eprintf "%s" e

let type_
  (* Allele information source *)
  alignment_file merge_file distance
  (* Allele selectors *)
    regex_list specific_list without_list number_alleles
  (* Process *)
    skip_disk_cache
  (* What to do? *)
    fastq_file_lst number_of_reads specific_reads
  (* options *)
    read_size_override
    =
  let need_read_size_r =
    to_read_size_dependent ?alignment_file ?merge_file ~distance
      ~regex_list ~specific_list ~without_list ?number_alleles
      ~skip_disk_cache
  in
  match need_read_size_r with
  | Error e           -> eprintf "%s" e
  | Ok need_read_size ->
    let init =
      match read_size_override with
      | None   -> `Setup need_read_size
      | Some r -> to_set need_read_size r
    in
    begin match fastq_file_lst with
    | []              -> invalid_argf "Cmdliner lied!"
    | [read1; read2]  -> invalid_argf "implement pairs!"
    | [fastq]         -> across_fastq ?number_of_reads ~specific_reads fastq init
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
            $ file_arg $ merge_arg $ distance_flag
            (* Allele selectors *)
            $ regex_arg $ allele_arg $ without_arg $ num_alt_arg
            (* What to do ? *)
            $ no_cache_flag
            (* What are we typing *)
            $ fastq_file_arg $ num_reads_arg $ specific_read_args
            (* options. *)
            $ read_size_override_arg
            (* How are we typing
            $ map_arg $ map_allele_arg
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
