
open Util

let app_name = "type"

(*
let freader = Future.Reader.open_file file in
let fastq_rdr = Fastq.read freader in
Future.Pipe.iter fastq_rdr ~f:(fun oe ->
  let fastq_item = Or_error.ok_exn oe in
  let seq = fastq_item.Fastq.sequence in
  match Alignment.align ~mub:3 g kmt seq with
  | Error msg -> eprintf "error %s for seq: %s\n" msg seq
  | Ok als    -> List.iter als ~f:(Alignment.alignments_to_weights amap));
  *)

let likelihood ?(alph_size=4) ?(er=0.01) ~len mismatches =
  let lmp = log (er /. (float (alph_size - 1))) in
  let lcp = log (1. -. er) in
  let c = len - mismatches in
  exp ((float c) *. lcp +. (float mismatches) *. lmp)

let sort_values_by_likelihood_assoc =
  (* higher values first! *)
  List.sort ~cmp:(fun (v1, _) (v2, _) -> compare v2 v1)

let sort_values_by_mismatches_assoc =
  (* lower values first! *)
  List.sort ~cmp:(fun (v1, _) (v2, _) -> compare v1 v2)

let output_values_assoc aindex =
  List.iter ~f:(fun (w, a) ->
    printf "%0.8f\t%s\n" w
      (insert_chars ['\t'; '\t'; '\n']
        (Alleles.Set.to_human_readable aindex ~max_length:1000 ~complement:`No a)))

let type_ verbose alignment_file num_alt_to_add allele_list k skip_disk_cache
  fastq_file not_join_same_seq number_of_reads print_top multi_pos as_mismatches =
  let open Cache in
  let open Ref_graph in
  let option_based_fname, g =
    Common_options.to_filename_and_graph_args alignment_file num_alt_to_add
      allele_list (not not_join_same_seq)
  in
  let g, idx = Cache.graph_and_two_index ~skip_disk_cache { k ; g } in
  let as_likelihood = not as_mismatches in
  let init, f = Path_inference.multiple_fold ~verbose ~multi_pos ~as_likelihood g idx in
  let amap =
    (* This is backwards .. *)
    Fastq_reader.fold ?number_of_reads fastq_file ~init ~f:(fun amap seq ->
      if verbose then print_endline "--------------------------------";
      match f amap seq with
      | Error e -> if verbose then printf "error\t%s: %s\n" seq e; amap
      | Ok a    -> if verbose then printf "matched\t%s \n" seq; a)
  in
  if as_likelihood then begin
    let sum = Alleles.Map.fold g.aindex ~f:(fun s v _ -> s +. v) ~init:0. amap in
    let amap = Alleles.Map.map g.aindex ~f:(fun v _allele -> v /. sum) amap in
    match print_top with
    | None ->
        (* Round the values so that it is easier to display. *)
        Alleles.Map.map g.aindex ~f:(fun x _allele -> (ceil (x *. 10.)) /. 10.) amap
        |> Alleles.Map.values_assoc g.aindex
        |> sort_values_by_likelihood_assoc
        |> output_values_assoc g.aindex
    | Some n ->
        Alleles.Map.fold g.aindex amap ~init:[]
          ~f:(fun a v al -> (v, Alleles.Set.singleton g.aindex al) :: a)
        |> sort_values_by_likelihood_assoc
        |> fun l -> List.take l n
        |> output_values_assoc g.aindex
  end else begin
    match print_top with
    | None ->
        Alleles.Map.values_assoc g.aindex amap
        |> sort_values_by_mismatches_assoc
        |> output_values_assoc g.aindex
    | Some n ->
        Alleles.Map.fold g.aindex amap ~init:[]
          ~f:(fun a v al -> (v, Alleles.Set.singleton g.aindex al) :: a)
        |> sort_values_by_mismatches_assoc
        |> fun l -> List.take l n
        |> output_values_assoc g.aindex
  end

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
  let as_mismatches_flag =
    let docv = "Print mismatches" in
    let doc = "Print mismatches for read as opposed to converting to likelihood." in
    Arg.(value & flag & info ~doc ~docv ["as-mismatches"])
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
            $ verbose_flag
            $ file_arg $ num_alt_arg $ allele_arg $ kmer_size_arg $ no_cache_flag
            $ fastq_file_arg
            $ do_not_join_same_sequence_paths_flag
            $ num_reads_arg
            $ print_top_flag
            $ multi_pos_flag
            $ as_mismatches_flag
        , info app_name ~version ~doc ~man)
  in
  match Term.eval type_ with
  | `Ok ()           -> exit 0
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
