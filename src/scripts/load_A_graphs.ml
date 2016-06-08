
let cargs = { To_graph.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
            ; To_graph.which = None
            }
let g2 = To_graph.(construct_from_file { cargs with which = Some (NumberOfAlts 2)})
let g3 = To_graph.(construct_from_file { cargs with which = Some (NumberOfAlts 3)})
let g20 = To_graph.(construct_from_file { cargs with which = Some (NumberOfAlts 20)})
let g_all = To_graph.(construct_from_file cargs)
