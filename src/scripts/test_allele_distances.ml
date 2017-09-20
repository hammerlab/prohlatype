
open Util
open Common

let () =
  if !Sys.interactive then () else
    let open MSA in
    let open Parser in
    let n = Array.length Sys.argv in
    let file = if n < 2 then "A_gen" else Sys.argv.(1) in
    let mp = from_file (Common.to_alignment_file file) in
    List.iter mp.alt_elems ~f:(fun alt1 ->
      printf "ref vs %s:%!" alt1.allele;
      let ds = 
        Segments.distances ~reference:mp.ref_elems ~allele:alt1.seq
        |> string_of_list ~sep:";" ~f:(fun d -> sprintf "%d" d.Segments.mismatches)
      in
      printf "%s\n" ds;
      List.iter mp.alt_elems ~f:(fun alt2 ->
        printf " %s vs %s:" alt1.allele alt2.allele;
        let dvs =
          Segments.distances_between ~reference:mp.ref_elems
            ~allele1:alt1.seq ~allele2:alt2.seq
          |> string_of_list ~sep:";" ~f:(fun d -> sprintf "%d" d.Segments.mismatches)
        in
        printf " %s\n" dvs))
