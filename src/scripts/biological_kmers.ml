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
  let arg = Ref_graph.default_construction_arg in
  let input = Alleles.Input.AlignmentFile (to_alignment_file file) in
  let gm = Cache.(graph (graph_args ~input ~arg)) in
  let graph_kt = Index.kmer_counts ~biological:true ~k gm in
  let idx = Index.create ~k gm in
  fasta_kt, fasta_seqs, graph_kt, idx

let is_substring s ~sub =
  match String.index_of_string s ~sub with
  | None -> false
  | Some _ -> true

let header_has_allele_name ~header allele =
  is_substring header ~sub:allele

let known_allele_in_header ~header known_alleles =
  List.find known_alleles ~f:(header_has_allele_name ~header)

let diff ~fasta_table ~graph_table fasta_seqs known_alleles idx =
  let graph_wrong, fasta_wrong =
    List.partition_map known_alleles ~f:(fun s ->
      if String.get_exn s ~index:0 = 'f' then
        `Snd (String.drop s ~index:1)
      else
        `Fst s)
  in
  let k = Kmer_table.k fasta_table in
  if Kmer_table.k graph_table <> k then
    invalid_argf "Different k!"
  else
    let last_index = (1 lsl (k * 2)) - 1 in
    let rec loop i =
      if i = last_index then
        None
      else
        let v1 = Kmer_table.lookup fasta_table i in
        let v2 = Kmer_table.lookup graph_table i in
        if v1 = 0 && v2 = 0 then
          loop (i + 1)
        else if v1 > 0 && v2 > 0 then
          loop (i + 1)
        else if v1 = 0 then
          check_missing_from_fasta i (Kmer_to_int.decode ~k i) v1 v2
        else
          check_known_diffs i (Kmer_to_int.decode ~k i) v1 v2 fasta_seqs
    and check_known_diffs i s v1 v2 = function
      | [] ->
          let pos_msg =
            match Index.lookup idx s with
            | Error ts    ->
                sprintf "(Couldn't even find it in index: %s)"
                  (show_too_short ts)
            | Ok pos_lst  ->
                sprintf "(at: %s)"
                  (String.concat ~sep:","
                    (List.map pos_lst ~f:Index.show_position))
          in
          Some (sprintf "for %d Kmer: %s (fasta occ: %d, graph occ: %d, %s) \
                         didn't find it in Fasta sequences!" i s v1 v2 pos_msg)
      | (header, seq) :: t ->
          if is_substring seq ~sub:s then
            match known_allele_in_header ~header graph_wrong with
            | None ->
                Some (sprintf "for %d Kmer: %s (fasta occ: %d, graph occ: %d), \
                               is in %s's sequence but not known! Check round \
                               trip?"
                        i s v1 v2 header)
            | Some allele ->
                begin
                  printf "for %d Kmer: %s (fasta occ: %d, graph occ: %d), is \
                          known diff for %s.\n" i s v1 v2 allele;
                  loop (i + 1)
                end
          else
            check_known_diffs i s v1 v2 t
    and check_missing_from_fasta i s v1 v2 =
      (* Only thing we can verify is that one of the fasta_wrong doesn't have
         the kmer in the sequence! *)
      let is_there_at_least_one_seq_with_missing_kmer =
        List.map fasta_wrong ~f:(fun allele ->
            List.filter_map fasta_seqs ~f:(fun (header, seq) ->
              if header_has_allele_name ~header allele then
                None
              else if is_substring seq ~sub:s then
                None
              else
                Some allele))
      in
      begin match is_there_at_least_one_seq_with_missing_kmer with
          | []  -> Some (sprintf "for %d Kmer: %s (fasta occ: %d, graph occ: %d) \
                                  didn't find in known fasta errors."
                                    i s v1 v2)
          | ls  -> None
      end
    in
    loop 0

let k = 10

let () =
  if !Sys.interactive then
    ()
  else
    let n = Array.length Sys.argv in
    let file, known_diff_alleles =
      if n <= 1 then
        "A_nuc", []
      else
        Sys.argv.(1), Array.to_list (Array.sub Sys.argv 1 (n - 1))
    in
    let fasta_table, fasta_seqs, graph_table, idx = create_kmer_counts ~k file in
    match diff ~fasta_table ~graph_table fasta_seqs known_diff_alleles idx with
    | None     -> print_endline "same"
    | Some msg -> print_endline msg; exit (-1)
