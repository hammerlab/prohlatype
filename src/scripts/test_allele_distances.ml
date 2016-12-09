
open Util
open Common

let () =
  let open Mas_parser in
  let file = if Array.length Sys.argv < 2 then "A_gen" else Sys.argv.(0) in
  let mp = from_file (Common.to_alignment_file file) in
  List.iter mp.alt_elems ~f:(fun (al1, allele) ->
    let ds = allele_distances ~reference:mp.ref_elems ~allele
            |> List.map ~f:(function
                | MissingFromAllele d -> sprintf "%d" d
                | Partial p -> sprintf "%d" p.mismatches
                | Full f -> sprintf "%d" f.mismatches)
            |> String.concat ~sep:";"
    in
    printf "ref vs %s: %s\n" al1 ds;
    List.iter mp.alt_elems ~f:(fun (al2, allele2) ->
      let dvs =
        allele_distances_between ~reference:mp.ref_elems
          ~allele1:allele ~allele2
        |> List.map ~f:(fun (_,d) -> sprintf "%d" d)
        |> String.concat ~sep:";"
      in
      printf " %s vs %s: %s\n" al1 al2 dvs))
