
open Util
open Common

let () =
  if !Sys.interactive then () else
    let open Mas_parser in
    let n = Array.length Sys.argv in
    let file = if n < 2 then "A_gen" else Sys.argv.(1) in
    let mp = from_file (Common.to_alignment_file file) in
    List.iter mp.alt_elems ~f:(fun (al1, allele) ->
      let ds = allele_distances ~reference:mp.ref_elems ~allele
              |> List.map ~f:(fun d -> sprintf "%d" d.mismatches)
              |> String.concat ~sep:";"
      in
      printf "ref vs %s: %s\n" al1 ds;
      List.iter mp.alt_elems ~f:(fun (al2, allele2) ->
        let dvs =
          allele_distances_between ~reference:mp.ref_elems
            ~allele1:allele ~allele2
          |> List.map ~f:(fun d -> sprintf "%d" d.mismatches)
          |> String.concat ~sep:";"
        in
        printf " %s vs %s: %s\n" al1 al2 dvs))
