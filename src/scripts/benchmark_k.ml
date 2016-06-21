
let dir = "../foreign/IMGTHLA/alignments/"
let files =
  [ "A_nuc.txt"
  ; "B_nuc.txt"
  ; "C_nuc.txt"
  ; "DMA_nuc.txt"
  ; "DMB_nuc.txt"
  ; "DOA_nuc.txt"
  ; "DOB_nuc.txt"
  ; "DPA1_nuc.txt"
  ; "DPB1_nuc.txt"
  ; "DPB2_nuc.txt"
  ; "DQA1_nuc.txt"
  ; "DQB1_nuc.txt"
  ; "DRA_nuc.txt"
  ; "DRB_nuc.txt"
  ]

let () =
  let to_cargs file =
    { To_graph.alignment_file = Filename.concat dir file
    ; To_graph.which = None }
  in
  let kmers_to_test = [ 5;6;7;8;9;10;11;12] in
  let doit k file =
    let carg = to_cargs file in
    let (_, g) = To_graph.construct_from_file carg in
    let kt = Ref_graph.kmer_counts g ~k in
    let dst = Kmer_table.distr kt in
    let ds_str =
      Array.to_list dst
      |> List.map string_of_int
      |> String.concat ","
    in
    Printf.printf "%s,%d,%s\n%!" file k ds_str
  in
  List.iter (fun file -> List.iter (fun k -> doit k file) kmers_to_test) files


