
(* Verbose version of the kmer counting algorithm used for debugging. *)

let index_string s index =
  let (b, a) = String.split_at s ~index in
  b ^ "." ^ a

let verbose_counts ~k g =
  let init = Kmer_table.make k 0 in
  let f tbl (al, sequence) { Graph_index.index; length = `Whole } =
    let () = Printf.printf "At %d, %s whole at %d\n" al (index_string sequence index) index in
    let i = Kmer_to_int.encode sequence ~pos:index ~len:k in
    Kmer_table.update ((+) 1) tbl i;
    tbl
  in
  let extend (al, sequence) { Graph_index.index; length = `Part len } cur_opt =
    let () = Printf.printf "At %d, %s extend at %d len %d" al (index_string sequence index) index len in
    match cur_opt with
    | None     ->
        let nj = Kmer_to_int.encode sequence ~pos:index ~len in
        let () = Printf.printf " cur: None \t nj %d \n" nj in
        nj
    | Some ext ->
        let nj = Kmer_to_int.encode sequence ~pos:index ~len ~ext in
        let () = Printf.printf " cur: Some %d \t nj %d\n" ext nj in
        nj
  in
  let close tbl (al, sequence) { Graph_index.index; length = `Part len} ext =
    let () = Printf.printf "At %d, %s close at %d len %d ext %d\n" al (index_string sequence index) index len ext in
    let i = Kmer_to_int.encode sequence ~pos:index ~len ~ext in
    let () = Printf.printf "our i: %d\n" i in
    Kmer_table.update ((+) 1) tbl i;
    let () = Printf.printf "updated\n" in
    tbl
  in
  let fs = Graph_index.fold_over_kmers_in_graph ~k g ~f ~close ~extend ~init in
  fs.Graph_index.full


