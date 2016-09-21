(* How good is the alignment? *)

open Util
(*open Common *)

let repo = "prohlatype"
let app_name = "explore_alignment"

let mismatches_from ?(verbose=false) (gall, idx) seq =
  Index.lookup idx seq >>= fun lst ->
    if verbose then printf "%d positions %!" (List.length lst);
    List.fold_left lst ~init:(Ok [])
      ~f:(fun ac_err p -> ac_err >>= fun acc ->
            Alignment.compute_mismatches gall seq p >>= fun al ->
              Ok ((p, al) :: acc))
    >>= fun lst -> Ok (seq, lst)

let mismatches_from_lst (gall, idx) seq =
  Index.lookup idx seq >>= fun lst ->
      List.fold_left lst ~init:(Ok [])
        ~f:(fun ac_err p -> ac_err >>= fun acc ->
              Alignment.compute_mismatches_lst gall seq p >>= fun al ->
                Ok ((p, al) :: acc))
      >>= fun lst -> Ok (seq, lst)

let report_msm (gall, _idx) (seq, reslst) =
  List.iter reslst ~f:(fun (p, res) ->
    printf "At position: %s\n" (Index.show_position p);
    Alleles.Map.values_assoc gall.Ref_graph.aindex res
    |> List.sort ~cmp:compare
    |> List.iter ~f:(fun (n, als) ->
        printf "\t%d\t%s\n" n (Alleles.Set.to_human_readable gall.Ref_graph.aindex als)))

let sum_mismatches = List.fold_left ~init:0 ~f:(fun s (_, m) -> s + m)

let report_msm_lst (gall, _idx) (seq, reslst) =
  List.iter reslst ~f:(fun (p, res) ->
    printf "At position: %s\n" (Index.show_position p);
    Alleles.Map.values_assoc gall.Ref_graph.aindex res
    |> List.sort ~cmp:(fun (m1,_) (m2,_) -> compare (sum_mismatches m1) (sum_mismatches m2))
    |> List.iter ~f:(fun (nlst, als) ->
        let ma = Alleles.Set.min_elt gall.Ref_graph.aindex als in
        let as_seq = Ref_graph.sequence ~start:(`AtPos (p.Index.alignment + p.Index.offset)) ~stop:(`Length 100) gall ma |> unwrap_ok in
        printf "\t%d %s:\n%s\n\t%s\n\t\t%s\n"
          (sum_mismatches nlst)
          ma
          (manual_comp_display seq as_seq)
          (String.concat ~sep:";" (List.map nlst ~f:(fun (p,d) -> sprintf "(%d,%d)" p d)))
          (Alleles.Set.to_human_readable gall.Ref_graph.aindex als)))

let best_match_allele_map = Alleles.Map.fold_wa ~init:max_int ~f:min

let best_matches ars =
  List.filter_map ars ~f:(function
    | Error e   -> printf "failed to match %s\n" e; None
    | Ok (s,pl) ->
        Some (List.map pl ~f:(fun (p,a) ->
          let bm = best_match_allele_map a in
          (bm, s, p, a))))
  |> List.concat
  |> List.sort ~cmp:(fun (bm1, _, _, _) (bm2, _, _, _) ->
      compare bm1 bm2)

let mismatch_histogram verbose k file fastq_file number_of_reads width =
  let widthf = float width in
  printf "file: %s\n" file;
  let gip = Cache.(graph_and_two_index { k = k; g = graph_arg ~file () }) in
  let rs = Fastq_reader.reads_from_fastq ?number_of_reads fastq_file in
  printf "number of reads: %d\n" (List.length rs);
  let msms = List.mapi rs ~f:(fun i seq ->
              if verbose then printf ", %d: " i;
              mismatches_from ~verbose gip seq) in
  let best = best_matches msms in
  let histo_me = List.map best ~f:(fun (i, _, _, _) -> float i) |> Array.of_list in
  let hist = Oml.Statistics.Descriptive.histogram (`Width widthf) histo_me in
  print_newline ();
  Array.iter hist ~f:(fun (b,v) -> printf "%f \t %d\n" b v)

let () =
  let open Cmdliner in
  let open Common_options in
  let width_arg =
    let docv = "Histogram width" in
    let doc = "Specify the width of the reported mismatch histogram. Defaults to 10." in
    Arg.(value & opt (some positive_int) (Some 10)
        & info ~doc ~docv ["w"; "width"])
  in
  let mismatch_histogram_ =
    let version = "0.0.0" in
    let doc = "Report alignment of reads in a fastaq against HLA string graphs." in
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
    Term.(const mismatch_histogram
            $ verbose_flag
            $ kmer_size_arg
            $ file_arg
            $ fastq_file_arg
            $ num_reads_arg
            $ width_arg
        , info app_name ~version ~doc ~man)
  in
  match Term.eval mismatch_histogram_ with
  | `Ok ()           -> exit 0
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
