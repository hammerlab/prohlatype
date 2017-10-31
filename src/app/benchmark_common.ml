open Util

let (//) = Filename.concat
let imgthla_dir =
  try Sys.getenv "IMGTHLA_DIR"
  with _ -> "../foreign/IMGTHLA"

let to_alignment_file f = imgthla_dir // "alignments" // (f ^ ".txt")
let to_merge_prefix p = imgthla_dir // "alignments" // p
let to_fasta_file f = imgthla_dir // "fasta" // (f ^ ".fasta")

let sample_path = "Users/leonidrozenberg/Documents/projects/hlatyping/upenn/pro/fastqs/"
let sample = sample_path // "120013_1.fastq"

let setup_samples ?(seed=10) ?(n=100) () =
  let r = Fastq.all sample in
  let ar = Array.of_list r in
  Random.init seed;
  let par = Oml.Util.Array.permute ~copy:true ar in
  Array.sub par ~pos:0 ~len:n
(*
  let without_names = Array.map ~f:(fun fqi -> { fqi with Biocaml_unix.Fastq.name = "" }) in
  let as_sexp_array = Sexplib.Conv.sexp_of_array Biocaml_unix.Fastq.sexp_of_item without_names in
  Sexplib.Sexp.to_string as_sexp_array *)

let load_samples s =
  Sexplib.Sexp.of_string s
  |> Sexplib.Conv.array_of_sexp Biocaml_unix.Fastq.item_of_sexp

(*let sa1 = load_samples s1  *)

(* Constant parameters. *)

let prealigned_transition_model = true
let loci_to_test = [ "A"; "B"; "C"]
let read_length = 100
let number_of_alleles_to_test = 50
let number_of_samples = 30

let distance = Distances.WeightedPerSegment

let to_parPHMM_t gen s =
  let input =
    if gen then
      Alleles.Input.alignment ~distance
        (to_alignment_file (s ^ "_gen"))
    else
      Alleles.Input.merge ~distance
        (to_merge_prefix "A")
  in
  ParPHMM.construct input
  |> unwrap_ok

let gen_parPHMM_ts gen =
  List.map loci_to_test ~f:(to_parPHMM_t gen)


