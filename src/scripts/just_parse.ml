(* Old, deprecated, test of the parsing logic. *)
let f = "A_nuc.txt" ;;
let d = "../foreign/IMGTHLA/alignments" ;;
let mp = MSA.Parser.from_file (Filename.concat d f) ;;
let just_gaps = (function MSA.Gap _ -> true | _ -> false) ;;
let just_bns = (function MSA.Boundary _ -> true | _ -> false) ;;
let print_parsed oa =
  List.Assoc.get oa mp.MSA.Parser.alt_elems
  |> Option.map ~f:(fun l ->
    List.map l ~f:MSA.al_el_to_string
    |> List.iter ~f:(Printf.printf "%s\n")) ;;
