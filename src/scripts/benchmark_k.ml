(* Figure out the distribution of matches in K-mer tables
   constructed on our graphs. *)

open Common

let () =
  let kmers_to_test = [ 5;6;7;8;9;10;11;12] in
  let doit k file =
    let input = Alleles.Input.alignment (to_alignment_file file) ~impute:false in
    let arg = Ref_graph.default_construction_arg in
    let carg = Cache.graph_args ~input ~arg in
    let g = Cache.graph carg in
    let kt = Index.kmer_counts g ~biological:true ~k in
    let dst = Kmer_table.distr kt in
    let ds_str =
      Array.to_list dst
      |> List.map string_of_int
      |> String.concat ", "
    in
    Printf.printf "%s, %d, %s\n%!" file k ds_str
  in
  Printf.printf "file, k, freq_1, freq_2, freq_3, ... \n";
  List.iter (fun file -> List.iter (fun k -> doit k file) kmers_to_test)
    common_alignment_files
