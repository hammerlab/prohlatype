open Prohlatype
open Common

let setup_samples ?(seed=10) ?(n=100) ?fastq () =
  let ar = Array.of_list (Test_sample.all_reads ?fastq ()) in
  Random.init seed;
  let par = Oml.Util.Array.permute ~copy:true ar in
  Array.sub par ~pos:0 ~len:n

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
  |> unwrap

let gen_parPHMM_ts gen =
  List.map loci_to_test ~f:(to_parPHMM_t gen)

let reads ?fastq () =
  [ Test_sample.a_reads ?fastq ()
  ; Test_sample.b_reads ?fastq ()
  ; Test_sample.c_reads ?fastq ()
  ]
