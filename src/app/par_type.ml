(* Typing via a Parametric PHMM. *)
open Util

let app_name = "par_type"
let time s f =
  let n = Sys.time () in
  let r = f () in
  Printf.printf "%-10s:%f\n" s (Sys.time () -. n);
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

let across_fastq ?number_of_reads ~specific_reads file read_size_dep =
  let module E = struct exception E of string end in
  try
    Fastq.fold ?number_of_reads ~specific_reads file ~init:(`Setup read_size_dep)
      ~f:(fun acc fqi ->
            match acc with
            | `Setup rp ->
                let read_size = String.length fqi.Biocaml_unix.Fastq.sequence in
                let pt = rp read_size in (* Invalid_args on error. *)
                let n = List.length pt.ParPHMM.alleles in
                let update_me = Array.make n 1. in
                let update update_me fqi =
                  match Fastq.phred_probabilities fqi.Biocaml_unix.Fastq.qualities with
                  | Core_kernel.Std.Result.Error e -> raise (E.E (Core_kernel.Std.Error.to_string_hum e))
                  | Core_kernel.Std.Result.Ok read_probs ->
                    let fp =
                      ParPHMM.forward_pass pt.ParPHMM.conf
                        fqi.Biocaml_unix.Fastq.sequence read_probs
                    in
                    let likelihoods = ParPHMM.to_allele_arr fp pt.ParPHMM.conf.ParPHMM.final_run_len in
                    Array.iteri likelihoods ~f:(fun i l -> update_me.(i) <- update_me.(i) *. l);
                    update_me
                in
                `Set (pt.ParPHMM.alleles, update, update update_me fqi)
            | `Set (alelles, update, ret) ->
                `Set (alelles, update, update ret fqi))
    |> function
        | `Setup _ -> eprintf "Didn't fine any reads."
        | `Set (alleles, _, final_likelihoods) ->
            List.mapi alleles ~f:(fun i a -> (final_likelihoods.(i), a))
            |> List.sort ~cmp:(fun (l1,_) (l2,_) -> compare l2 l1)
            |> List.iter ~f:(fun (l,a) -> printf "%10s\t%0.6f\n" a l)
  with E.E e ->
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
    = 
  let to_read_size =
    to_read_size_dependent ?alignment_file ?merge_file ~distance
      ~regex_list ~specific_list ~without_list ?number_alleles
      ~skip_disk_cache
  in
  match to_read_size with
  | Error e         -> eprintf "%s" e
  | Ok to_read_size ->
    begin match fastq_file_lst with
    | []              -> invalid_argf "Cmdliner lied!"
    | [read1; read2]  -> invalid_argf "implement pairs!"
    | [fastq]         -> across_fastq ?number_of_reads ~specific_reads fastq to_read_size
    | lst             -> invalid_argf "More than 2, %d fastq files specified!" (List.length lst)
    end
  
let () =
  let open Cmdliner in
  let open Common_options in
  let type_ =
    let version = "0.0.0" in
    let doc = "Use a Parametric Profile Hidden Markov Model of HLA allele to type fastq samples." in
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
