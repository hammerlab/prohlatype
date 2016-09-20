(* Load A graphs *)

open Common

let file = to_alignment_file "A_nuc"
let g2 = Ref_graph.(construct_from_file ~which:(NumberOfAlts 2) ~join_same_sequence:true file)
let g3 = Ref_graph.(construct_from_file ~which:(NumberOfAlts 3) ~join_same_sequence:true file)
let g20 = Ref_graph.(construct_from_file ~which:(NumberOfAlts 20) ~join_same_sequence:true file)
let g_all = Ref_graph.(construct_from_file ~join_same_sequence:true file)
