
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
        (Alleles.Set.to_human_readable aindex ~max_length:60 ~complement:`No a)))

let default_error_fname = "typing_errors.log"

let report_errors ~error_output elst =
  let oc, close =
    match error_output with
    | `DefaultFile -> open_out default_error_fname, true
    | `Stderr      -> stderr, false
    | `Stdout      -> stdout, false
  in
  List.iter elst ~f:(fun (sae, fqi) ->
    fprintf oc "%s\n%s\n" fqi.Biocaml_unix.Fastq.sequence
      (Path_inference.sequence_alignment_error_to_string sae));
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

let report_likelihood ?bucket ~print_top contents g amap =
  let open Ref_graph in
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
                Alleles.Set.set g.aindex arr.(index) allele);
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

let type_
  (* Graph construction args. *)
  alignment_file
    num_alt_to_add allele_list allele_regex_list not_join_same_seq remove_reference
  (* Index *)
    k
  (* Process *)
    skip_disk_cache
  (* What are we typing. *)
    fastq_file number_of_reads
  (* How are we typing *)
    filter multi_pos as_stat likelihood_error
  (* Output *)
    print_top do_not_normalize bucket error_output =
  let open Cache in
  let open Ref_graph in
  let _option_based_fname, graph_args =
    Common_options.to_filename_and_graph_args
      ~alignment_file num_alt_to_add
      ~allele_list ~allele_regex_list
        ~join_same_sequence:(not not_join_same_seq)
        ~remove_reference
  in
  let g, idx = Cache.graph_and_two_index ~skip_disk_cache { k ; graph_args } in
  let new_as =
    match as_stat with
    | `MismatchesList   -> `MismatchesList
    | `Mismatches       -> `Mismatches
    | `Likelihood       -> `Likelihood likelihood_error
    | `LogLikelihood    -> `LogLikelihood likelihood_error
    | `PhredLikelihood  -> `PhredLikelihood
  in
  match Path_inference.type_ ?filter ~as_:new_as g idx ?number_of_reads ~fastq_file with
  | `Mismatches (error_list, amap)      -> report_errors ~error_output error_list;
                                           report_mismatches ~print_top g amap
  | `MismatchesList (error_list, amap)  -> report_errors ~error_output error_list;
                                           report_mismatches_list ~print_top g amap
  | `Likelihood (error_list, amap)      -> report_errors ~error_output error_list;
                                           report_likelihood ?bucket ~print_top "likelihood" g amap
  | `LogLikelihood (error_list, amap)   -> report_errors ~error_output error_list;
                                           report_likelihood ?bucket ~print_top "loglikelihood" g amap
  | `PhredLikelihood (error_list, amap) -> report_errors ~error_output error_list;
                                           report_likelihood ?bucket ~print_top "phredlihood" g amap

       (*    let amap =
              if do_not_normalize then amap else
                let sum = Alleles.Map.fold_wa ~f:(+.) ~init:0. amap in
                Alleles.Map.map_wa ~f:(fun v -> v /. sum) amap
            in
            report_likelihood g amap do_not_bucket print_top
        | `LogLikelihood  ->
            let emap =
              if do_not_normalize then
                amap
              else
                let sum = Alleles.Map.fold_wa ~init:0. ~f:(fun s v -> v +. s) amap in
                Alleles.Map.map_wa ~f:(fun v -> v /. sum) amap
            in
            report_likelihood g emap do_not_bucket print_top
      end
  | `PhredLikelihood ->
      begin
        let open Core_kernel.Std in
        let init, f = Path_inference.multiple_phred ~verbose ~multi_pos ?early_stop g idx in
        let amap =
          Fastq.fold ?number_of_reads fastq_file ~init ~f:(fun amap i ->
            let q = i.Biocaml_unix.Fastq.qualities in
            match Fastq.phred_probabilities q with
            | Error e -> if verbose then printf "error\t%s: %s" q (Error.to_string_hum e); amap
            | Ok pro  ->
                let seq = i.Biocaml_unix.Fastq.sequence in
                match f amap (seq, pro) with
                | Error e -> if verbose then printf "error\t%s: %s\n" seq e; amap
                | Ok a    -> if verbose then printf "matched\t%s \n" seq; a)
        in
        let amap =
          if do_not_normalize then amap else
            let sum = Alleles.Map.fold_wa ~f:(+.) ~init:0. amap in
            Alleles.Map.map_wa ~f:(fun v -> v /. sum) amap
        in
        report_likelihood g amap do_not_bucket print_top
      end

*)
let () =
  let open Cmdliner in
  let open Common_options in
  let print_top_flag =
    let docv = "Print only most likely" in
    let doc = "Print only the specified number (positive integer) of alleles" in
    Arg.(value & opt (some int) None & info ~doc ~docv ["print-top"])
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
    Arg.(value & vflag `LogLikelihood
      [ `LogLikelihood,   info ~doc:(d ^ "log likelihood") ["log-likelihood"]
      ; `Likelihood,      info ~doc:(d ^ "likelihood") ["likelihood"]
      ; `Mismatches,      info ~doc:(d ^ "mismatches, that are then added then added together") ["mismatches"]
      ; `MismatchesList,  info ~doc:(d ^ "list of mismatches") ["mis-list"]
      ; `PhredLikelihood, info ~doc:(d ^ "sum of log likelihoods based of phred qualities") ["phred-llhd"]
      ])
  in
  let filter_flag =
    let docv = "Filter out sequences" in
    let doc  = "Filter, do not include in the likelihood calculation, sequences \
                where the highest number of mismatches is greater than the passed argument." in
    Arg.(value & opt (some int) None & info ~doc ~docv ["filter-matches"])
  in
  let do_not_normalize_flag =
    let docv = "Do not normalize the likelihoods" in
    let doc  = "Do not normalize the per allele likelihoods to report accurate probabilities." in
    Arg.(value & flag & info ~doc ~docv ["do-not-normalize"])
  in
  let bucket_arg =
    let docv = "Bucket the likelihoods." in
    let doc  = "When printing the allele likelihoods, aggregate the PDF in the \
                specified number of buckets." in
    Arg.(value & opt (some int) None & info ~doc ~docv ["bucket"])
  in
  let likelihood_error_arg =
    let default = 0.025 in
    let docv = "Override the likelihood error" in
    let doc  =
      sprintf "Specify the error value used in likelihood calculations, defaults to %f"
        default
    in
    Arg.(value & opt float default & info ~doc ~docv ["likelihood-error"])
  in
  let error_output_flag =
    let doc dest =
      sprintf "Output errors such as sequences that don't match to %s. \
               By default output is written to %s." dest default_error_fname
    in
    Arg.(value & vflag `DefaultFile
      [ `Stdout,      info ~doc:(doc "standard output") ["error-stdout"]
      ; `Stderr,      info ~doc:(doc "standard error") ["error-stderr"]
      ; `DefaultFile, info ~doc:(doc "default filename") ["error-default"]
      ])
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
            $ file_arg $ num_alt_arg $ allele_arg $ allele_regex_arg
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
            (* Output *)
            $ print_top_flag
              $ do_not_normalize_flag
              $ bucket_arg
              $ error_output_flag
        , info app_name ~version ~doc ~man)
  in
  match Term.eval type_ with
  | `Ok ()           -> exit 0
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
