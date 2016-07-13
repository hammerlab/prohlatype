

let file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
let g2 = To_graph.(construct_from_file ~which:(NumberOfAlts 2) file)
let g3 = To_graph.(construct_from_file ~which:(NumberOfAlts 3) file)
let g20 = To_graph.(construct_from_file ~which:(NumberOfAlts 20) file)
let g_all = To_graph.(construct_from_file file)


