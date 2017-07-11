(* Typing via a Parametric PHMM. *)
open Util

let app_name = "multi_par"

let (//) = Filename.concat

let time s f =
  let n = Sys.time () in
  let r = f () in
  Printf.printf "%s total running time in seconds: %f\n%!" s (Sys.time () -. n);
  r

let to_read_size_dependent
  (* Allele information source *)
    ~alignment_files ~merge_files ~distance ~impute
    ~skip_disk_cache =
    let als =
      List.map alignment_files ~f:(Common_options.input_alignments ~impute)
    in
    let mls =
      List.map merge_files ~f:(Common_options.input_merges ~distance ~impute)
    in
    match als @ mls with
    | []     -> Error "Neither a merge nor alignment file specified!"
    | inputs ->
      let selectors = [] in           (* Selectors NOT supported, on purpose! *)
      Ok (fun read_size ->
        List.map inputs ~f:(fun input ->
          let name = Alleles.Input.to_string input in
          let par_phmm_arg = Cache.par_phmm_args ~input ~selectors ~read_size in
          let pt = Cache.par_phmm ~skip_disk_cache par_phmm_arg in
          name, pt))

exception PPE of string

let ppe s = raise (PPE s)
let ppef fmt =
  ksprintf ppe fmt

(* TODO: write out type annotation for f, otherwise it is too abstract *)
let to_update_f read_size ~f acc fqi =
  time (sprintf "updating on %s" fqi.Biocaml_unix.Fastq.name) (fun () ->
    let open Core_kernel.Std in
    let ls = String.length fqi.Biocaml_unix.Fastq.sequence in
    let lq = String.length fqi.Biocaml_unix.Fastq.qualities in
    if ls <> read_size then
      ppef "Sequence length %d not configured length:%d, skipping." ls read_size
    else if lq <> read_size then
      ppef "Sequence length %d not configured length:%d, skipping." ls read_size
    else
      match Fastq.phred_log_probs fqi.Biocaml_unix.Fastq.qualities with
      | Result.Error e       -> ppe (Error.to_string_hum e)
      | Result.Ok read_probs -> f acc ~name:fqi.Biocaml_unix.Fastq.name
                                      ~seq:fqi.Biocaml_unix.Fastq.sequence
                                      ~read_probs)

type 'a g =
  { f         : 'a -> Biocaml_unix.Fastq.item -> 'a
  ; fin       : 'a -> unit
  ; mutable s : 'a
  }

let apply = function
  | `Mapper g -> begin fun fqi ->
                  try
                    g.s <- g.f g.s fqi;
                    `Mapper g
                  with PPE e ->
                    eprintf "for %s: %s" fqi.Biocaml_unix.Fastq.name e;
                    `Mapper g
                 end
  | `Reporter g -> begin fun fqi ->
                    try
                      g.s <- g.f g.s fqi;
                      `Reporter g
                  with PPE e ->
                    eprintf "for %s: %s" fqi.Biocaml_unix.Fastq.name e;
                    `Reporter g
                 end

let sort_mapped_output lst =
  let open ParPHMM in
  let highest_llhd l =
    List.map l ~f:(fun (_name, s) -> best_stat s)
    |> List.reduce ~f:max
  in
  List.sort lst ~cmp:(fun (_rn1, nlst1) (_rn2, nlst2) ->
    let sl1 = highest_llhd nlst1 in
    let sl2 = highest_llhd nlst2 in
    compare sl2 sl1 (* higher is better *))

let output_mapped lst =
  List.iter lst ~f:(fun (read, lst) ->
    printf "%s\n\t%s\n"
      read
      (String.concat ~sep:"\n\t"
        (List.map lst ~f:(fun (name, ms) ->
          sprintf "%s\t%s" name
            (ParPHMM.mapped_stats_to_string ~sep:'\t' ms)))))

let to_set ?insert_p ~past_threshold_filter ?max_number_mismatches ?band mode rp read_size =
  let ptlst =
    time (sprintf "Setting up ParPHMM transitions with %d read_size" read_size)
      (fun () -> rp read_size)
  in
  let g =
    time (sprintf "Allocating forward pass workspaces")
      (fun () ->
        let r =
          ParPHMM.forward_pass_m ?insert_p ?max_number_mismatches ?band
            ~past_threshold_filter mode ptlst read_size
        in
        match r with
        | `Mapper m ->
            `Mapper
              { f   = to_update_f read_size
                          ~f:(fun l ~name ~seq ~read_probs ->
                                (name, m seq read_probs) :: l)
              ; s   = []
              ; fin = begin fun lst ->
                        printf "Reads: %d\n" (List.length lst);
                        output_mapped (sort_mapped_output lst)
                      end
              }
        | `Reporter (lst, r) ->
            `Reporter
              { f   = to_update_f read_size
                          ~f:(fun l ~name ~seq ~read_probs ->
                                r seq read_probs l)
              ; s   = lst
              ; fin = begin fun lst ->
                        List.iter2 ptlst lst ~f:(fun (name, pt) (_name, llhd) ->
                          printf "%s\n" name;
                          ParPHMM.output (`Pt pt) llhd)
                      end
              })
  in
  `Set g

let fin = function
  | `Setup _ -> eprintf "Didn't find any reads."
  | `Set (`Mapper g)  -> g.fin g.s
  | `Set (`Reporter g) -> g.fin g.s

let across_fastq ?insert_p ?max_number_mismatches
  ~past_threshold_filter
  ?number_of_reads ~specific_reads ?band mode file init =
  try
    Fastq.fold ?number_of_reads ~specific_reads ~init file
      ~f:(fun acc fqi ->
            match acc with
            | `Setup rp ->
                let read_size = String.length fqi.Biocaml_unix.Fastq.sequence in
                let `Set g =
                  to_set ?insert_p ~past_threshold_filter
                    ?max_number_mismatches ?band mode rp read_size
                in
                `Set (apply g fqi)
            | `Set g ->
                `Set (apply g fqi))
    |> fin
  with PPE e ->
    eprintf "%s" e

let type_
  (* Allele information source *)
    class1_gen_dir
    alignment_files merge_files distance
  (* Process *)
    skip_disk_cache
  (* What to do? *)
    fastq_file_lst number_of_reads specific_reads
  (* options *)
    insert_p
    do_not_past_threshold_filter
    max_number_mismatches
    read_size_override
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
  let impute   = true in
  let band     =
    if not_band then
      None
    else
      Some { ParPHMM.warmup; number; width }
  in
  let alignment_files, merge_files =
    match class1_gen_dir with
    | None   -> alignment_files, merge_files
    | Some d -> [ d // "A_gen.txt"
                ; d // "B_gen.txt"
                ; d // "C_gen.txt"
                ], []
  in
  let past_threshold_filter = not do_not_past_threshold_filter in
  let need_read_size_r =
    to_read_size_dependent
      ~alignment_files ~merge_files ~distance ~impute
      ~skip_disk_cache
  in
  let mode = match map with | Some n -> `Mapper n | None -> `Reporter in
  match need_read_size_r with
  | Error e           -> eprintf "%s" e
  | Ok need_read_size ->
    let init =
      match read_size_override with
      | None   -> `Setup need_read_size
      | Some r -> to_set ~insert_p ?max_number_mismatches
                    ~past_threshold_filter ?band mode need_read_size r
    in
    begin match fastq_file_lst with
    | []              -> invalid_argf "Cmdliner lied!"
    | [read1; read2]  -> invalid_argf "implement pairs!"
    | [fastq]         -> across_fastq ~insert_p ?max_number_mismatches
                            ~past_threshold_filter
                            ?number_of_reads ~specific_reads
                            ?band mode fastq init
    | lst             -> invalid_argf "More than 2, %d fastq files specified!" (List.length lst)
    end

let () =
  let open Cmdliner in
  let open Common_options in
  let files_arg =
    let docv = "FILE" in
    let doc  = "File to lookup IMGT allele alignments. The alleles found in this \
                file will initially define the set of alleles to be used for \
                named locus. Supply multiple file (or merge) arguments to \
                compute the likelihood of a read against each locus and add to \
                the most likely (highest likelihood)." in
    Arg.(value & opt_all file [] & info ~doc ~docv ["f"; "file"])
  in
  let merges_arg =
    let parser_ path =
      let s = Filename.basename path in
      let n = path ^ "_nuc.txt" in
      let g = path ^ "_gen.txt" in
      if not (List.mem ~set:Merge_mas.supported_genes s) then
        `Error ("gene not supported: " ^ s)
      else if not (Sys.file_exists n) then
        `Error ("Nuclear alignment file doesn't exist: " ^ n)
      else if not (Sys.file_exists g) then
        `Error ("Genetic alignment file doesn't exist: " ^ n)
      else
        `Ok path  (* Return path, and do appending later, the prefix is more useful. *)
    in
    let convrtr = parser_, (fun frmt -> Format.fprintf frmt "%s") in
    let docv = sprintf "[%s]" (String.concat ~sep:"|" Merge_mas.supported_genes) in
    let doc  =
      sprintf "Construct a merged (gDNA and cDNA) graph of the specified \
              prefix path. Currently only supports %s genes. The argument must \
              be a path to files with $(docv)_nuc.txt and $(docv)_gen.txt. \
              Combines with the file arguments to determine the set of loci to \
              type at the same time.. The set of alleles is defined by the \
              ones in the nuc file."
        (String.concat ~sep:", " Merge_mas.supported_genes)
    in
    Arg.(value & opt_all convrtr [] & info ~doc ~docv ["m"; "merge"])
  in
  let max_number_mismatches_arg =
    let docv = "POSITIVE INTEGER" in
    let doc = "Setup a filter on the reads to cancel evaluation once we've \
               seen this many mismatches." in
    Arg.(value & opt (some positive_int) None & info ~doc ~docv
          ["max-mismatches"])
  in
  let class1gen_arg =
    let docv  = "DIRECTORY" in
    let doc   = "Short-cut argument that expands the given dir to look for \
                  A_gen.txt, B_gen.txt and C_gen.txt. This overrides any \
                 alignment files or merge arguments." in
    Arg.(value & opt (some dir) None & info ~doc ~docv ["class1-gen"])
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
            $ class1gen_arg
            $ files_arg $ merges_arg $ distance_flag
            (* What to do ? *)
            $ no_cache_flag
            (* What are we typing *)
            $ fastq_file_arg $ num_reads_arg $ specific_read_args
            (* options. *)
            $ insert_probability_arg
            $ do_not_past_threshold_filter_flag
            $ max_number_mismatches_arg
            $ read_size_override_arg
            $ not_band_flag
            $ band_warmup_arg
            $ number_bands_arg
            $ band_width_arg
            (* How are we typing *)
            $ map_flag
            $ forward_pass_accuracy_arg
            (* $ map_allele_arg
            $ filter_flag $ multi_pos_flag $ stat_flag $ likelihood_error_arg
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
