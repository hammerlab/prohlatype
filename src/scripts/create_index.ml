(* Load a graph and index. *)
open Common

let k = 10 ;;

let default_input = Ref_graph.AlignmentFile (to_alignment_file "A_nuc")
let gall, idxall =
  Cache.(graph_and_two_index { k = k ; graph_args = graph_args ~input:default_input () }) ;;

let g_and_idx ?(k=10) ?input ?n () =
  let input = Option.value input ~default:(Ref_graph.AlignmentFile (to_alignment_file "A_nuc")) in
  let g2 =
    Cache.(graph_and_two_index { k = k; graph_args = graph_args ?n ~input () })
  in
  n, g2

(*
let lst = [3; 5; 10; 15; 20; 50; 100; 150; 200; 250]
let idxs = List.map lst ~f:(fun n -> g_and_idx ~n:(Ref_graph.NumberOfAlts n) ())
*)
