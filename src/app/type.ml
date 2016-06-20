
(* open Core.Std
open Biocaml_unix.Std
open Future_unix.Std *)
open Nonstd
open Util

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

let () =
  let cargs =
    { To_graph.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
    ; To_graph.which = None
    }
  in
  let aset, g = To_graph.(construct_from_file cargs) in
  printf " Got graph!\n%!";
  let kt = Graph_index.kmer_list ~k:5 g in
  printf " Got index!\n%!";
  let kmt = kt.Graph_index.full in
  let upenn_opti_res_dir = "~/Documents/projects/hlatyping/upenn/opti/merged" in
  let file = Filename.concat upenn_opti_res_dir "120013_TGACCA/120013_TGACCA_1fin.fastq" in
  let amap = Graph_alignment.init_alingment_map aset in
  printf " Aligning!\n%!";
  shitty_fastq_sequence_reader file (fun seq ->
    Printf.printf "aligning: %s\n%!" seq;
    match String.index_of_character seq 'N' with
    | Some _ -> printf "skipping N!\n"
    | None   ->
      match Graph_alignment.align ~mub:1 g kmt seq with
      | Error msg -> eprintf "error %s for seq: %s\n" msg seq
      | Ok als    -> List.iter als ~f:(Graph_alignment.alignments_to_weights amap));
  (*
  let freader = Future.Reader.open_file file in
  let fastq_rdr = Fastq.read freader in
  Future.Pipe.iter fastq_rdr ~f:(fun oe ->
    let fastq_item = Or_error.ok_exn oe in
    let seq = fastq_item.Fastq.sequence in
    match Graph_alignment.align ~mub:3 g kmt seq with
    | Error msg -> eprintf "error %s for seq: %s\n" msg seq
    | Ok als    -> List.iter als ~f:(Graph_alignment.alignments_to_weights amap));
    *)
  Graph_alignment.most_likely aset amap
  |> List.iter ~f:(fun (w, a) -> printf "%f \t %s\n" w a)

