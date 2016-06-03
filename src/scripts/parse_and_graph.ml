
let f = "A_nuc.txt" ;;
let d = "../foreign/IMGTHLA/alignments" ;;
let mp = Mas_parser.from_file (Filename.concat d f) ;;
let g = To_graph.construct ~num_alt_to_add:20 [] mp ;;
let g2 = To_graph.construct ~num_alt_to_add:2 [] mp ;;
let g3 = To_graph.construct ~num_alt_to_add:3 [] mp ;;


