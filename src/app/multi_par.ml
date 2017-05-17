(* Typing via a Parametric PHMM. *)
open Util

let app_name = "multi_par"

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
      | Result.Ok read_probs -> f acc fqi.Biocaml_unix.Fastq.name fqi.Biocaml_unix.Fastq.sequence read_probs)

type 'a g =
  { f         : 'a -> Biocaml_unix.Fastq.item -> 'a
  ; fin       : 'a -> unit
  ; mutable s : 'a
  }

let proc_g = function
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
      
let to_set ?band mode rp read_size =
  let ptlst =
    time (sprintf "Setting up ParPHMM transitions with %d read_size" read_size)
      (fun () -> rp read_size)
  in
  let g =
    time (sprintf "Allocating forward pass workspaces")
      (fun () ->
        match ParPHMM.forward_pass_m ?band mode ptlst read_size with
        | `Mapper m ->
            `Mapper
              { f   = to_update_f read_size ~f:(fun l n s r -> (n, m s r) :: l)
              ; s   = []
              ; fin = begin fun lst ->
                        printf "Reads: %d\n" (List.length lst);
                        output_mapped (sort_mapped_output lst)
                      end
              }
        | `Reporter (lst, ff) ->
            `Reporter
              { f   = to_update_f read_size ~f:(fun l _n s r -> ff s r l)
              ; s   = lst
              ; fin = begin fun lst ->
                        List.iter2 ptlst lst ~f:(fun (name, pt) (_name, llhd) ->
                          printf "%s\n" name;
                          List.mapi pt.ParPHMM.alleles ~f:(fun i a -> llhd.(i), a)
                          |> List.sort ~cmp:(fun (l1,_) (l2,_) -> compare l2 l1)
                          |> List.iter ~f:(fun (l,a) -> printf "%10s\t%0.20f\n" a l))
                      end
              })
  in
  `Set g

let fin = function
  | `Setup _ -> eprintf "Didn't fine any reads."
  | `Set (`Mapper g)  -> g.fin g.s
  | `Set (`Reporter g) -> g.fin g.s

let across_fastq ?number_of_reads ~specific_reads ?band mode file init =
  try
    Fastq.fold ?number_of_reads ~specific_reads ~init file
      ~f:(fun acc fqi ->
            match acc with
            | `Setup rp ->
                let read_size = String.length fqi.Biocaml_unix.Fastq.sequence in
                let `Set g = to_set ?band mode rp read_size in
                `Set (proc_g g fqi)
            | `Set g ->
                `Set (proc_g g fqi))
    |> fin
  with PPE e ->
    eprintf "%s" e

let type_
  (* Allele information source *)
    alignment_files merge_files distance
  (* Process *)
    skip_disk_cache
  (* What to do? *)
    fastq_file_lst number_of_reads specific_reads
  (* options *)
    read_size_override
    not_band
    start_column
    number
    width
  (* how are we typing *)
    map
    =
  let impute   = true in
  let band     =
    if not_band then
      None
    else
      Some { ParPHMM.start_column; number; width }
  in
  let need_read_size_r =
    to_read_size_dependent
      ~alignment_files ~merge_files ~distance ~impute
      ~skip_disk_cache
  in
  let mode = if map then `Mapper else `Reporter in
  match need_read_size_r with
  | Error e           -> eprintf "%s" e
  | Ok need_read_size ->
    let init =
      match read_size_override with
      | None   -> `Setup need_read_size
      | Some r -> to_set ?band mode need_read_size r
    in
    begin match fastq_file_lst with
    | []              -> invalid_argf "Cmdliner lied!"
    | [read1; read2]  -> invalid_argf "implement pairs!"
    | [fastq]         -> across_fastq ?number_of_reads ~specific_reads 
                            ?band mode fastq init
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
  (* Band arguments *)
  let not_band_flag =
    let doc = "Calculate the full forward pass matrix." in
    Arg.(value & flag & info ~doc ["do-not-band"])
  in
  let band_start_column_arg =
    let default = ParPHMM.(band_default.start_column) in
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
    let default = ParPHMM.(band_default.number) in
    let docv = "POSITIVE INTEGER" in
    let doc  =
      sprintf "Number of bands to calculate. Defaults to %d" default
    in
    Arg.(value
          & opt positive_int default
          & info ~doc ~docv ["number-bands"])
  in
  let band_width_arg =
    let default = ParPHMM.(band_default.width) in
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
            $ files_arg $ merges_arg $ distance_flag
            (* What to do ? *)
            $ no_cache_flag
            (* What are we typing *)
            $ fastq_file_arg $ num_reads_arg $ specific_read_args
            (* options. *)
            $ read_size_override_arg
            $ not_band_flag
            $ band_start_column_arg
            $ number_bands_arg
            $ band_width_arg
            (* How are we typing *)
            $ map_flag
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
