
open Util

let app_name = "type"

let sort_values_by_likelihood_assoc assoc =
  (* higher values first! *)
  List.sort assoc ~cmp:(fun (v1, a1) (v2, a2) ->
    let r = compare v2 v1 in
    if r = 0 then compare a2 a1 else r)

let sort_values_by_mismatches_assoc assoc =
  (* lower values first! *)
  List.sort assoc ~cmp:(fun (v1, a1) (v2, a2) ->
    let r = compare v1 v2 in
    if r = 0 then compare a2 a1 else r)

let output_values_assoc contents to_string aindex assoc =
  printf "%s\tset size\tset contents\n" contents;
  List.iter assoc ~f:(fun (w, a) ->
    printf "%s\t\t%d\t\t%s\n"
      (to_string w)
      (Alleles.Set.cardinal a)
      (insert_chars ['\t'; '\t'; '\n']
        (Alleles.Set.to_human_readable aindex ~max_length:120 ~complement:`No a)))

let default_error_fname = "typing_errors.log"

let report_errors ~all_options ~error_output elst =
  let oc, close =
    match error_output with
    | `InputPrefixed  -> let fname = sprintf "%s_typing_errors.log" all_options in
                         open_out fname, true
    | `DefaultFile    -> open_out default_error_fname, true
    | `Stderr         -> stderr, false
    | `Stdout         -> stdout, false
  in
  List.iter elst ~f:(fun (sae, fqi) ->
    fprintf oc "%s\n%s\n%s\n"
      fqi.Biocaml_unix.Fastq.name
      fqi.Biocaml_unix.Fastq.sequence
      (Path_inference.show_sequence_alignment_error sae));
  if close then close_out oc

let report_mismatches ~print_top g amap =
  let open Ref_graph in
  begin match print_top with
    | None ->
        Alleles.Map.values_assoc g.aindex amap
        |> sort_values_by_mismatches_assoc
    | Some n ->
        Alleles.Map.fold g.aindex amap ~init:[]
          ~f:(fun a v al -> (v, Alleles.Set.singleton g.aindex al) :: a)
        |> sort_values_by_mismatches_assoc
        |> fun l -> List.take l n
  end
  |> output_values_assoc "mismatches" (sprintf "%3d") g.aindex

let report_mismatches_list ~print_top g amap =
  let open Ref_graph in
  let mismatch_sum_across_lists =
    List.fold_left ~init:0 ~f:(fun init l ->
      List.fold_left l ~init ~f:(fun s (_, m) -> s + m))
  in
  Alleles.Map.fold g.aindex amap ~init:[] ~f:(fun acc l allele ->
    (mismatch_sum_across_lists l, l, allele) :: acc)
  |> List.sort ~cmp:(fun (s1, _, _) (s2, _, _) -> compare s1 s2)
  |> begin fun l ->
      match print_top with
      | None -> l
      | Some n -> List.take l n
     end
  |> List.iter ~f:(fun (_sum, l, allele) ->
        printf "%s\n" allele;
        List.iter l ~f:(fun l ->
          List.map l ~f:(fun (p, m) -> sprintf "(%3d,%2d)" p m)
          |> String.concat ~sep:"; "
          |> insert_chars ~every:100 ~token:';' ['\t'; '\n']  (* in reverse *)
          |> printf "\t%s\n"))

let report_likelihood ?reduce_resolution ?bucket ~print_top contents g amap =
  let open Ref_graph in
  match reduce_resolution with
  | None ->
    begin match print_top with
      | None ->
          begin
            match bucket with
            | None   -> Alleles.Map.values_assoc g.aindex amap
            | Some n ->
                let s,e =
                  Alleles.Map.fold_wa amap ~init:(infinity, neg_infinity)
                    ~f:(fun (s, e) v -> (min s v, max e v))
                in
                let nf = float n in
                let d = (e -. s) /. nf in
                let arr = Array.init n (fun _ -> Alleles.Set.init g.aindex) in
                Alleles.Map.iter g.aindex amap ~f:(fun v allele ->
                  let index = if v = e then n - 1 else truncate ((v -. s) /. d) in
                  arr.(index) <- Alleles.Set.set g.aindex arr.(index) allele);
                Array.mapi ~f:(fun i set -> s +. d *. (float i), set) arr
                |> Array.to_list
          end
          |> sort_values_by_likelihood_assoc
      | Some n ->
          Alleles.Map.fold g.aindex amap ~init:[]
            ~f:(fun a v al -> (v, Alleles.Set.singleton g.aindex al) :: a)
          |> sort_values_by_likelihood_assoc
          |> fun l -> List.take l n
    end
    |> output_values_assoc contents (sprintf "%0.3f") g.aindex
  | Some n ->
      let values =
        (* TODO: move this logic to Alleles.. *)
        Alleles.Map.fold g.Ref_graph.aindex ~init:[]
          ~f:(fun l v al -> (al, v) :: l) amap
      in
      let sort values =
        List.sort values ~cmp:(fun (a1, v1) (a2, v2) ->
          let r = compare v2 v1 in
          if r = 0 then compare a2 a1 else r)
      in
      let output l = List.iter l ~f:(fun (a, l) -> printf "%0.3f\t%s\n" l a) in
      match n with
      | 4 -> sort values |> output
      | 3 -> sort values |> Compress.compress ~level:`Three |> sort |> output
      | 2 -> sort values |> Compress.compress ~level:`Two |> sort |> output
      | 1 -> sort values |> Compress.compress ~level:`One |> sort |> output
      | _ -> invalid_argf "You promised me parsing! %d" n

let basename_without_extension f =
  let bn = Filename.basename f in
  try Filename.chop_extension bn
  with Invalid_argument _ -> bn

let type_
  (* Graph construction args. *)
  alignment_file
  merge_file
  distance
    num_alt_to_add allele_list allele_regex_list not_join_same_seq remove_reference
  (* Index *)
    k
  (* Process *)
    skip_disk_cache
  (* What are we typing. *)
    fastq_file_lst number_of_reads
  (* How are we typing *)
    filter multi_pos as_stat likelihood_error dont_check_rc max_distance
  (* Output *)
    print_top do_not_normalize bucket error_output reduce_resolution =
  let join_same_sequence = not not_join_same_seq in
  Common_options.to_filename_and_graph_args
    ?alignment_file ?merge_file ~distance num_alt_to_add
    ~allele_list ~allele_regex_list
    ~join_same_sequence ~remove_reference
    >>= begin fun (option_based_fname, cargs) ->
        let open Cache in
        let open Ref_graph in
        let g, idx = Cache.(graph_and_two_index ~skip_disk_cache { k ; graph_args = cargs}) in
        let new_as =
          match as_stat with
          | `List_mismatches_of_reads       -> `List_mismatches_of_reads
          | `Number_of_mismatches_of_reads  -> `Number_of_mismatches_of_reads
          | `Likelihood_of_reads            -> `Likelihood_of_reads likelihood_error
          | `LogLikelihood_of_reads         -> `LogLikelihood_of_reads likelihood_error
          | `Phred_likelihood_of_reads      -> `Phred_likelihood_of_reads
        in
        let fastq_fold_args =
          { Path_inference.number_of_reads = number_of_reads
          ; check_rc                       = Some (not dont_check_rc)
          ; max_distance                   = max_distance
          }
        in
        let ffs =
          sprintf "%s%s%s"
            (Option.value_map ~default:"" ~f:string_of_int number_of_reads)
            (if dont_check_rc then "" else "-check_rc-")
            (Option.value_map ~default:"" ~f:string_of_int max_distance);
        in
        let input_files, res =
          match fastq_file_lst with
          | []              -> invalid_argf "Cmdliner lied!"
          | [fastq_file]    ->
              basename_without_extension fastq_file,
              Path_inference.type_ ?filter ~as_:new_as fastq_fold_args g idx ~fastq_file
          | [read1; read2] ->
              sprintf "%s-%s"
                (basename_without_extension read1)
                (basename_without_extension read2),
              Path_inference.type_paired ?filter ~as_:new_as fastq_fold_args g idx read1 read2
          | lst             -> invalid_argf "More than 2, %d fastq files specified!" (List.length lst)
        in
        let all_options = sprintf "%s_%s_%s" option_based_fname ffs input_files in
        let re = report_errors ~all_options ~error_output in
        let () =
          match res with
          | `Number_of_mismatches_of_reads (error_list, amap) ->
              re error_list;
              report_mismatches ~print_top g amap
          | `List_mismatches_of_reads (error_list, amap)      ->
              re error_list;
              report_mismatches_list ~print_top g amap
          | `Likelihood_of_reads (error_list, amap)           ->
              re error_list;
              report_likelihood ?reduce_resolution ?bucket ~print_top
                "likelihood" g amap
          | `LogLikelihood_of_reads (error_list, amap)        ->
              re error_list;
              report_likelihood ?reduce_resolution ?bucket ~print_top
                "loglikelihood" g amap
          | `Phred_likelihood_of_reads (error_list, amap)     ->
              re error_list;
              report_likelihood ?reduce_resolution ?bucket ~print_top
                "phredlihood" g amap
        in
        Ok ()
  end
  |> function
      | Error msg -> eprintf "%s" msg; 1
      | Ok ()     -> 0

let () =
  let open Cmdliner in
  let open Common_options in
  let print_top_flag =
    let doc = "Print only the specified number (positive integer) of alleles" in
    Arg.(value & opt (some int) None & info ~doc ["print-top"])
  in
  let multi_pos_flag =
    let d = "How to aggregate multiple position matches: " in
    Arg.(value & vflag `Best
      [ `TakeFirst, info ~doc:(d ^ "take the first, as found in Index.") ["pos-take-first"]
      ; `Average,   info ~doc:(d ^ "average over all positions") ["pos-average"]
      ; `Best,      info ~doc:(d ^ "the best over all positions (default).") ["pos-best"]
      ])
  in
  let stat_flag =
    let d = "What statistics to compute over each sequences: " in
    Arg.(value & vflag `LogLikelihood_of_reads
      [ `LogLikelihood_of_reads,   info ~doc:(d ^ "log likelihood") ["log-likelihood"]
      ; `Likelihood_of_reads,      info ~doc:(d ^ "likelihood") ["likelihood"]
      ; `Number_of_mismatches_of_reads,      info ~doc:(d ^ "mismatches, that are then added then added together") ["mismatches"]
      ; `List_mismatches_of_reads,  info ~doc:(d ^ "list of mismatches") ["mis-list"]
      ; `Phred_likelihood_of_reads, info ~doc:(d ^ "sum of log likelihoods based of phred qualities") ["phred-llhd"]
      ])
  in
  let filter_flag =
    let doc  = "Filter, do not include in the likelihood calculation, sequences \
                where the highest number of mismatches is greater than the passed argument." in
    Arg.(value & opt (some int) None & info ~doc ["filter-matches"])
  in
  let do_not_normalize_flag =
    let doc  = "Do not normalize the per allele likelihoods to report accurate probabilities." in
    Arg.(value & flag & info ~doc ["do-not-normalize"])
  in
  let bucket_arg =
    let doc  = "When printing the allele likelihoods, aggregate the PDF in the \
                specified number of buckets." in
    Arg.(value & opt (some int) None & info ~doc ["bucket"])
  in
  let likelihood_error_arg =
    let default = 0.025 in
    let doc  =
      sprintf "Specify the error value used in likelihood calculations, defaults to %f"
        default
    in
    Arg.(value & opt float default & info ~doc ["likelihood-error"])
  in
  let error_output_flag =
    let doc dest =
      sprintf "Output errors such as sequences that don't match to %s. \
               By default output is written to %s." dest default_error_fname
    in
    Arg.(value & vflag `InputPrefixed
      [ `Stdout,        info ~doc:(doc "standard output") ["error-stdout"]
      ; `Stderr,        info ~doc:(doc "standard error") ["error-stderr"]
      ; `DefaultFile,   info ~doc:(doc "default filename") ["error-default"]
      ; `InputPrefixed, info ~doc:(doc "input prefixed") ["error-input-prefixed"]
      ])
  in
  let reduce_resolution_arg =
    let doc  = "Reduce the resolution of the PDF, to a lower number of \
                \"digits\". The general HLA allele nomenclature \
                has 4 levels of specificity depending on the number of colons \
                in the name. For example A*01:01:01:01 has 4 and A*01:95 \
                has 2. This argument specifies the number (1,2,or 3) of \
                digits to reduce results to. For example, specifying 2 will \
                choose the best (depending on metric) allele out of \
                A*02:01:01:01, A*02:01:01:03, ...  A*02:01:02, ... \
                A*02:01:100 for A*02:01. The resulting set will have at most \
                the specified number of digit groups, but may have less."
    in
    let one_to_three_parser s =
      match positive_int_parser s with
      | `Ok x when x = 1 || x = 2 || x = 3 -> `Ok x
      | `Ok x                              -> `Error (sprintf  "not 1 to 3: %d" x)
      | `Error e                           -> `Error e
    in
    let one_to_three = one_to_three_parser , (fun frmt -> Format.fprintf frmt "%d") in
    Arg.(value & opt (some one_to_three) None & info ~doc ["reduce-resolution"])
  in
  let do_not_check_rc_flag =
    Arg.(value & flag & info ["do-not-check-rc"])
  in
  let max_distance_arg =
    let docv = "NON-NEGATIVE INTEGER" in
    let doc  = "Distance away from base kmer to check for positions." in
    Arg.(value & opt (some non_negative_int) None & info ~doc ["kmer-distance"])
  in
  let type_ =
    let version = "0.0.0" in
    let doc = "Use HLA string graphs to type fastq samples." in
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
            (* Graph construction args *)
            $ file_arg $ merge_arg $ distance_flag
            $ num_alt_arg $ allele_arg $ allele_regex_arg
              $ do_not_join_same_sequence_paths_flag
              $ remove_reference_flag
            (* Index *)
            $ kmer_size_arg
            (* Process *)
            $ no_cache_flag
            (* What are we typing *)
            $ fastq_file_arg $ num_reads_arg
            (* How are we typing *)
            $ filter_flag $ multi_pos_flag $ stat_flag $ likelihood_error_arg
              $ do_not_check_rc_flag
              $ max_distance_arg
            (* Output *)
            $ print_top_flag
              $ do_not_normalize_flag
              $ bucket_arg
              $ error_output_flag
              $ reduce_resolution_arg
        , info app_name ~version ~doc ~man)
  in
  match Term.eval type_ with
  | `Ok n            -> exit n
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
