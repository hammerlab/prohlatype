
let f = "A_nuc.txt" ;;
let d = "../foreign/IMGTHLA/alignments" ;;
let mp = Mas_parser.from_file (Filename.concat d f) ;;
let just_gaps = (function Mas_parser.Gap _ -> true | _ -> false) ;;
let just_bns = (function Mas_parser.Boundary _ -> true | _ -> false) ;;
let oa = "A*26:40";;
List.assoc oa mp.Mas_parser.alt_elems |> List.map ~f:Mas_parser.al_el_to_string |> List.iter ~f:(Printf.printf "%s\n") ;;
