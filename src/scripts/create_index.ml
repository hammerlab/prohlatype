
let cargs n =
  { Cache.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
  ; Cache.which = Some (Ref_graph.NumberOfAlts n)
  ; Cache.normalize = true
  } ;;

let all_args =
  { Cache.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
  ; Cache.which = None
  ; Cache.normalize = true
  } ;;

let k = 10 ;;

let gall, idxall = Cache.graph_and_two_index { Cache.k = k; Cache.g = all_args };;

let g_and_idx ?(k=10) n = (n, Cache.graph_and_two_index { Cache.k = k; Cache.g = cargs n })

let lst = [3; 5; 10; 15; 20; 50; 100; 150; 200; 250]
let idxs = List.map lst ~f:g_and_idx
