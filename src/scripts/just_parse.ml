(* Old, deprecated, test of the parsing logic. *)
let f = "A_nuc.txt" ;;
let d = "../foreign/IMGTHLA/alignments" ;;
let mp = MSA.Parser.from_file (Filename.concat d f) ;;
let just_gaps = (function MSA.Gap _ -> true | _ -> false) ;;
let just_bns = (function MSA.Boundary _ -> true | _ -> false) ;;

let print_alignment_sequence l =
    List.map l ~f:MSA.al_el_to_string
    |> List.iter ~f:(Printf.printf "%s\n") ;;

let mp_to_seq mp =
  let open MSA in
  let no_sequences = List.filter ~f:(fun a -> not (is_sequence a)) in
  let als =
    List.map mp.Parser.alt_elems ~f:(fun { Parser.allele; seq; _} ->
      ( allele
      , allele_sequences ~reference:mp.Parser.ref_elems ~allele:seq))
  in
  let rp =
    ( mp.Parser.reference
    , allele_sequences ~reference:mp.Parser.ref_elems
        ~allele:(no_sequences mp.Parser.ref_elems))
  in
  rp :: als

let mismatches s1 s2 =
  let n1 = String.length s1 in
  let n2 = String.length s2 in
  let m = min n1 n2 in
  let x = max n1 n2 - m in
  let rec loop index s =
    if index = m then s else
      if String.get s1 ~index = String.get s2 ~index then
        loop (index + 1) s
      else
        loop (index + 1) (s + 1)
  in
  loop 0 x

let search ?m ?n a1 a2 =
  let m = Option.value m ~default:(Array.length a1) in
  let n = Option.value n ~default:(Array.length a2) in
  let rec loop i j mm p =
    if i = m then
      mm, p
    else if j = n then
      loop (i + 1) 0 mm p
    else
      let al1, s1 = a1.(i) in
      let al2, s2 = a2.(j) in
      (*printf "comparing %s vs %s\n%!" al1 al2;*)
      let msms = mismatches s1 s2 in
      if msms < mm then
        loop i (j + 1) msms (al1, al2)
      else
        loop i (j + 1) mm p
  in
  loop 0 0 max_int ("","")
