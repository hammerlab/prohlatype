(* Load A graphs *)

open Common

let file = to_alignment_file "A_nuc"
let cff = Ref_graph.construct_from_file ~join_same_sequence:true ~remove_reference:false 
let g2 = Ref_graph.(cff ~which:(NumberOfAlts 2) file)
let g3 = Ref_graph.(cff ~which:(NumberOfAlts 3) file)
let g20 = Ref_graph.(cff ~which:(NumberOfAlts 20) file)
let g_all = Ref_graph.(cff file)
