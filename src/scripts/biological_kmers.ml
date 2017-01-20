(* Compare the Kmer table generate via Index vs a manually created one over
   the allele fasta file. *)


open Util
open Common

let kmer_table_from_fasta ~k file =
  let kt = Kmer_table.init k (fun _ -> 0) in
  let fasta_seqs =
    Fasta.fold file
    ~init:[] ~f:(fun acc s ->
      let seq = s.Biocaml_unix.Fasta.sequence in
      Index.fold_over_kmers_in_string seq ~k ~init:()
        ~f:(fun () ss ->
              match ss.Index.length with
              | `Part _ -> ()
              | `Whole  ->
                  Kmer_to_int.encode_N_tolerant ~pos:ss.Index.index ~len:k seq
                  |> List.iter ~f:(Kmer_table.update succ kt));
      (s.Biocaml_unix.Fasta.description, seq) :: acc)
  in
  kt, fasta_seqs

let create_kmer_counts ~k file =
  let fasta_kt, fasta_seqs = kmer_table_from_fasta ~k (to_fasta_file file) in
  let input = Ref_graph.AlignmentFile (to_alignment_file file) in
  let gm = Cache.(graph (graph_args ~input ())) in
  let graph_kt = Index.kmer_counts ~biological:true ~k gm in
  fasta_kt, fasta_seqs, graph_kt

let is_substring s ~sub =
  match String.index_of_string s ~sub with
  | None -> false
  | Some _ -> true

let header_has_allele_name ~header allele =
  is_substring header ~sub:allele

let known_allele_in_header ~header known_alleles =
  List.find known_alleles ~f:(header_has_allele_name ~header)

let diff kt1 kt2 fasta_seqs known_alleles =
  let k = Kmer_table.k kt1 in
  if Kmer_table.k kt2 <> k then
    invalid_argf "Different k!"
  else
    let last_index = (1 lsl (k * 2)) - 1 in
    let rec loop i =
      if i = last_index then
        None
      else
        let v1 = Kmer_table.lookup kt1 i in
        let v2 = Kmer_table.lookup kt2 i in
        if v1 = 0 && v2 = 0 then
          loop (i + 1)
        else if v1 > 0 && v2 > 0 then
          loop (i + 1)
        else
          check_known_diffs i (Kmer_to_int.decode ~k i) v1 v2 fasta_seqs
    and check_known_diffs i s v1 v2 = function
      | [] -> Some (sprintf "for %d Kmer: %s (fasta occ: %d, graph occ: %d) \
                            didn't find it in Fasta sequences!" i s v1 v2)
      | (header, seq) :: t ->
          if is_substring seq ~sub:s then
            match known_allele_in_header ~header known_alleles with
            | None ->
                Some (sprintf "for %d Kmer: %s (fasta occ: %d, graph occ: %d), \
                               is in %s's sequence but not known! Check round \
                               trip?"
                        i s v1 v2 header)
            | Some allele ->
                printf "for %d Kmer: %s (fasta occ: %d, graph occ: %d), is known diff for %s.\n"
                    i s v1 v2 allele;
                loop (i + 1)
          else
            check_known_diffs i s v1 v2 t
    in
    loop 0

let k = 10

let () =
  if !Sys.interactive then
    ()
  else
    let n = Array.length Sys.argv in
    let file, known_diffs =
      if n <= 1 then
        "A_nuc", []
      else
        Sys.argv.(1), Array.to_list (Array.sub Sys.argv 1 (n - 1))
    in
    let from_fasta, fasta_seqs, from_graph = create_kmer_counts ~k file in
    match diff from_fasta from_graph fasta_seqs known_diffs with
    | None     -> print_endline "same"
    | Some msg -> print_endline msg; exit (-1)
