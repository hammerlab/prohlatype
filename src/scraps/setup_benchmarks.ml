(* Extract reads from a result that were particularly well aligned.
 *
 * These reads were used for the benchmarking test.
 * *)

let j =
  Post_analysis.of_json_file 
    "FILL ME/205673_1.json"

let j_bl = Post_analysis.reads_by_loci j

let j_bl_A = List.Assoc.get Nomenclature.A j_bl |> Option.value_exn ~msg:""
let j_bl_B = List.Assoc.get Nomenclature.B j_bl |> Option.value_exn ~msg:""
let j_bl_C = List.Assoc.get Nomenclature.C j_bl |> Option.value_exn ~msg:""

let select_reads l =
  let open Post_analysis in
  let open ParPHMM_drivers in
  let open Multiple_loci in
  List.filter_map l ~f:(fun (_allele, readname, pa_pr) ->
    match pa_pr with
    | Pr (FirstOrientedSecond { first_o = false; _}) -> Some readname
    | _ -> None)

