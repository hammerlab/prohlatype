(* Load a graph and index. *)
open Common

let k = 10 ;;

let file = to_alignment_file "A_nuc"
let gall, idxall =
  Cache.(graph_and_two_index { k = k ; g = graph_arg ~file () }) ;;

let g_and_idx ?(k=10) ?file ?n () =
  let file = Option.value file ~default:(to_alignment_file "A_nuc") in
  let g2 =
    Cache.(graph_and_two_index { k = k; g = graph_arg ?n ~file () })
  in
  n, g2

let lst = [3; 5; 10; 15; 20; 50; 100; 150; 200; 250]
let idxs = List.map lst ~f:(fun n -> g_and_idx ~n:(Ref_graph.NumberOfAlts n) ())

