
let cargs n =
  { Cache.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
  ; Cache.which = Some (To_graph.NumberOfAlts n)
  } ;;

let all_args = { Cache.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
            ; Cache.which = None
            } ;;

let k = 10 ;;

(*let g10, aindex, idx = Cache.graph_and_two_index { Cache.k = k; Cache.g = cargs };;*)

let lst = [3; 5; 10; 15; 20; 50; 100; 150; 200; 250]

let idxs =
  List.map lst ~f:(fun n ->
    (n, Cache.graph_and_two_index { Cache.k = 10; Cache.g = cargs n }))
