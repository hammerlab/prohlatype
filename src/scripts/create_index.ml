(* Load a graph and index. *)
open Common

let k = 10 ;;

let default_input =
  Alleles.Input.AlignmentFile (to_alignment_file "A_gen", true) ;;

let gall, idxall =
  Cache.(graph_and_two_index { k = k
        ; graph_args = graph_args ~input:default_input
                          ~arg:Ref_graph.default_construction_arg  }) ;;

let g_and_idx ?(k=10) ?input ~arg () =
  let default = Alleles.Input.AlignmentFile (to_alignment_file "A_nuc", false) in
  let input = Option.value input ~default in
  Cache.(graph_and_two_index { k = k; graph_args = graph_args ~input ~arg })

  (*
let _, (gm, im) =
  g_and_idx ~input:(Alleles.Input.MergeFromPrefix
                      (to_merge_prefix "C"
                      , Distances.AverageExon
                      , false)) ()
                      *)
