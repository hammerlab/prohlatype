(* Extract reads from a result that were particularly well aligned.
 *
 * These reads were used for the benchmarking test.
 * *)

let j =
  Post_analysis.of_json_file
    "FILL ME/205673_1.json"

let j_bl = Post_analysis.reads_by_loci j

let j_bl_A = assoc_exn Nomenclature.A j_bl
let j_bl_B = assoc_exn Nomenclature.B j_bl
let j_bl_C = assoc_exn Nomenclature.C j_bl

let select_reads l =
  let open Post_analysis in
  let open ParPHMM_drivers in
  let open Multiple_loci in
  List.filter_map l ~f:(fun (_allele, readname, pa_pr) ->
    match pa_pr with
    | Pr (FirstOrientedSecond { first_o = false; _}) -> Some readname
    | _ -> None)

