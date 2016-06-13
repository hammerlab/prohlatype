
let cargs = { To_graph_ge.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
            ; To_graph_ge.which = None
            }
(* Old
let g2 = To_graph.(construct_from_file { cargs with which = Some (NumberOfAlts 2)})
let g3 = To_graph.(construct_from_file { cargs with which = Some (NumberOfAlts 3)})
let g20 = To_graph.(construct_from_file { cargs with which = Some (NumberOfAlts 20)})
let g_all = To_graph.(construct_from_file cargs)
   *)

let g2 = To_graph_ge.(construct_from_file { cargs with which = Some (NumberOfAlts 2)})
let g3 = To_graph_ge.(construct_from_file { cargs with which = Some (NumberOfAlts 3)})
let g20 = To_graph_ge.(construct_from_file { cargs with which = Some (NumberOfAlts 20)})
let g_all = To_graph_ge.(construct_from_file cargs)
