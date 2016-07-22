
open Util
open Common_options

let app_name = "type"

let shitty_fastq_sequence_reader file f =
  let ic = open_in file in
  try
    let rec loop i =
      let line = input_line ic in
      if i mod 4 = 1 then f line;
      loop (i + 1)
    in
    loop 0
  with End_of_file ->
    close_in ic

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


let type_ alignment_file num_alt_to_add allele_list k skip_disk_cache fastq_file not_normalize =
  let open Cache in
  let open Ref_graph in
  let option_based_fname, g =
    to_filename_and_graph_args alignment_file num_alt_to_add allele_list
      (not not_normalize)
  in
  let g, idx = Cache.graph_and_two_index ~skip_disk_cache { k ; g } in
  printf " Got graph and index!\n%!";
  let amap = Alignment.init_alignment_map g.aindex in
  printf " Aligning!\n%!";
  shitty_fastq_sequence_reader fastq_file (fun seq ->
    Printf.printf "aligning: %s\n%!" seq;
    match String.index_of_character seq 'N' with
    | Some _ -> printf "skipping N!\n"
    | None   ->
      match Index.lookup idx seq with
      | Error m     -> printf "Error looking up %s\n" m
      | Ok []       -> printf "Empty for %s\n" seq
      | Ok (p :: _) ->  (* TODO, more than one! *)
          let md = Alignment.compute_mismatches g seq p in
          let len = String.length seq in
          Alleles.Map.update2 md amap (fun m c ->
            c +. likelihood ~len m));

  Alignment.most_likely g.aindex amap
  |> List.iter ~f:(fun (w, a) -> printf "%f \t %s\n" w a)

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
            $ file_arg $ num_alt_arg $ allele_arg $ kmer_size_arg $ no_cache_flag
            $ fastq_file_arg
            $ do_not_normalize_flag
        , info app_name ~version ~doc ~man)
  in
  match Term.eval type_ with
  | `Ok ()           -> exit 0
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0

