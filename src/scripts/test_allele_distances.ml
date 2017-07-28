
open Util
open Common

let () =
  if !Sys.interactive then () else
    let open MSA in
    let open Parser in
    let n = Array.length Sys.argv in
    let file = if n < 2 then "A_gen" else Sys.argv.(1) in
    let mp = from_file (Common.to_alignment_file file) in
    List.iter mp.alt_elems ~f:(fun (al1, allele) ->
      printf "ref vs %s:%!" al1;
      let ds = 
        Segments.distances ~reference:mp.ref_elems ~allele
        |> List.map ~f:(fun d -> sprintf "%d" d.Segments.mismatches)
        |> String.concat ~sep:";"
      in
      printf "%s\n" ds;
      List.iter mp.alt_elems ~f:(fun (al2, allele2) ->
        printf " %s vs %s:" al1 al2;
        let dvs =
          Segments.distances_between ~reference:mp.ref_elems
            ~allele1:allele ~allele2
          |> List.map ~f:(fun d -> sprintf "%d" d.Segments.mismatches)
          |> String.concat ~sep:";"
        in
        printf " %s\n" dvs))
