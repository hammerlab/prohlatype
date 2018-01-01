
open Util

let j =
  Post_analysis.of_json_file
    "FILL ME/res/2017_12_12_full_class1/004.json"

let j_bl =
  Post_analysis.reads_by_loci j
let j_bl_A =
  List.Assoc.get Nomenclature.A j_bl |> Option.value_exn ~msg:""
let j_bl_A1, j_bl_A2 =
  List.partition ~f:(fun (a, _, _) -> a = "A*26:08") j_bl_A

let specific_reads =
  List.map j_bl_A ~f:(fun (_, rn, _) -> rn)

let file1 = "FILL ME/fastqs/PGV_004_1.fastq"
let file2 = "FILL ME/fastqs/PGV_004_2.fastq"

let to_s = Biocaml_unix.Fastq.item_to_string

let () =
  let oc1 = open_out "pgv004_1.fastq" in
  let oc2 = open_out "pgv004_2.fastq" in
  Fastq.fold_paired_both
    ~specific_reads
    ~init:()
    ~f:(fun () r1 r2 -> fprintf oc1 "%s" (to_s r1); fprintf oc2 "%s" (to_s r2))
    ~ff:(fun () r -> fprintf oc1 "%s" (to_s r))
    ~fs:(fun () r -> fprintf oc2 "%s" (to_s r))
    file1
    file2
  |> ignore
