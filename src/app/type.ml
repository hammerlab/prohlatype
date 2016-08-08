
open Util
open Common_options

let app_name = "type"

let shitty_fastq_sequence_reader file ~f ~init =
  let ic = open_in file in
  let r = ref init in
  try
    let rec loop i =
      let line = input_line ic in
      if i mod 4 = 1 then r := f !r line;
      loop (i + 1)
    in
    loop 0
  with End_of_file ->
    close_in ic;
    !r

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

let sort_values_assoc =
  (* higher values first! *)
  List.sort ~cmp:(fun (v1, _) (v2, _) -> compare v2 v1)

let output_values_assoc aindex =
  List.iter ~f:(fun (w, a) ->
    printf "%0.8f\t%s\n" w
      (insert_chars ['\t'; '\t'; '\n']
        (Alleles.Set.to_human_readable aindex ~max_length:1000 ~complement:`No a)))

let type_ verbose alignment_file num_alt_to_add allele_list k skip_disk_cache
  fastq_file not_join_same_seq print_top =
  let open Cache in
  let open Ref_graph in
  let option_based_fname, g =
    to_filename_and_graph_args alignment_file num_alt_to_add allele_list
      (not not_join_same_seq)
  in
  let g, idx = Cache.graph_and_two_index ~skip_disk_cache { k ; g } in
  if verbose then
    printf " Got graph and index!\n%!";
  let amap = Alignment.init_alignment_map g.aindex in
  if verbose then
    printf " Aligning!\n%!";
  shitty_fastq_sequence_reader fastq_file ~init:() ~f:(fun () seq ->
    if verbose then
      printf "aligning: %s, %!" seq;
    match String.index_of_character seq 'N' with
    | Some _ -> if verbose then printf "skipping N!%!\n"
    | None   ->
      match Index.lookup idx seq with
      | Error m -> if verbose then printf "error looking up %s.\n" m
      | Ok []   -> if verbose then printf "empty positions. \n"
      | Ok lst  ->  (* TODO, more than one! *)
          let n = List.length lst in
          if verbose then printf " found %d alignments, averaging... \n" n;
          let weight = 1. /. (float n) in
          List.iter lst ~f:(fun p ->
            match Alignment.compute_mismatches g seq p with
            | Error m ->
                if verbose then
                  printf "error during mismatch "
            | Ok md ->
                let len = String.length seq in
                Alleles.Map.update2 md amap (fun m c -> c +. weight *. (likelihood ~len m))));

  let sum = Alleles.Map.fold g.aindex ~f:(fun s v _ -> s +. v) ~init:0. amap in
  let amap = Alleles.Map.map g.aindex ~f:(fun v _allele -> v /. sum) amap in
  match print_top with
  | None ->
      (* Round the values so that it is easier to display. *)
      Alleles.Map.map g.aindex ~f:(fun x _allele -> (ceil (x *. 10.)) /. 10.) amap
      |> Alleles.Map.values_assoc g.aindex
      |> sort_values_assoc
      |> output_values_assoc g.aindex
  | Some n ->
      Alleles.Map.fold g.aindex amap ~init:[]
        ~f:(fun a v al -> (v, Alleles.Set.singleton g.aindex al) :: a)
      |> sort_values_assoc
      |> fun l -> List.take l n
      |> output_values_assoc g.aindex

let () =
  let open Cmdliner in
  let kmer_size_arg =
    let default = 5 in
    let docv = "Kmer size" in
    let doc =
      sprintf "Number of consecutive nucleotides to use consider in K-mer
               index construction. Defaults to %d." default
    in
    Arg.(value & opt int default & info ~doc ~docv ["k"; "kmer-size"])
  in
  let fastq_file_arg =
    let docv = "Fastq samples" in
    let doc = "Fastq formatted DNA reads file." in
    Arg.(required & pos 0 (some file) None & info ~doc ~docv [])
  in
  let verbose_flag =
    let docv = "Be verbose" in
    let doc = "Print progress messages to stdout." in
    Arg.(value & flag & info ~doc ~docv ["v"; "verbose"])
  in
  let print_top_flag =
    let docv = "Print only most likely" in
    let doc = "Print only the specified number (positive integer) of alleles" in
    Arg.(value & opt (some int) None & info ~doc ~docv ["print-top"])
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
            $ print_top_flag
        , info app_name ~version ~doc ~man)
  in
  match Term.eval type_ with
  | `Ok ()           -> exit 0
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0

