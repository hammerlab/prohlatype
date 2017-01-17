(* Compare the Kmer table generate via Index vs a manually created one over
   the allele fasta file. *)


open Util
open Common

let kmer_table_from_fasta ~k file =
  let init = Kmer_table.init k (fun _ -> 0) in
  Fasta.fold file
  ~init ~f:(fun kt s ->
    let seq = s.Biocaml_unix.Fasta.sequence in
    Index.fold_over_kmers_in_string seq ~k ~init:kt
      ~f:(fun kt ss ->
            match ss.Index.length with
            | `Part _ -> kt
            | `Whole  ->
                Kmer_table.update succ kt
                  (Kmer_to_int.encode ~pos:ss.Index.index ~len:k seq);kt))

let comp ~k file =
  let fasta_kt = kmer_table_from_fasta ~k (to_fasta_file file) in
  let input = Ref_graph.AlignmentFile (to_alignment_file file) in
  let gm = Cache.(graph (graph_args ~input ())) in
  let graph_kt = Index.kmer_counts ~biological:true ~k gm in
  fasta_kt, graph_kt

let diff kt1 kt2 =
  let k = Kmer_table.k kt1 in
  if Kmer_table.k kt2 <> k then
    invalid_argf "Different k!"
  else
    let last_index = (1 lsl (k * 2)) - 1 in
    let rec loop i =
      if i = last_index then None else
        let v1 = Kmer_table.lookup kt1 i in
        let v2 = Kmer_table.lookup kt2 i in
        if v1 = 0 && v2 = 0 then
          loop (i + 1)
        else if v1 > 0 && v2 > 0 then
          loop (i + 1)
        else
          Some (i, v1, v2)
    in
    loop 0

let k = 10

let () =
  if !Sys.interactive then
    ()
  else
    let n = Array.length Sys.argv in
    let file = if n <= 1 then "A_nuc" else Sys.argv.(1) in
    let from_fasta, from_graph = comp ~k file in
    match diff from_fasta from_graph with
    | None             -> printf "same\n"
    | Some (i, p1, p2) -> printf "for %d %s, fasta %d, graph %d\n"
                            i (Kmer_to_int.decode ~k i) p1 p2;
                          exit (-1)
