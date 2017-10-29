(* Benchmark and compare performance of running a single allele PHMM against
   the parameterized PHMM. *)

open Util
open Benchmark_common

let time_aggregate ?(seed=11) parPHMM_t reads =
  let pta = ParPHMM.setup_single_pass ~prealigned_transition_model read_length parPHMM_t in
  Array.map reads ~f:(fun r ->
    let read = r.Biocaml_unix.Fastq.sequence in
    let read_errors = unwrap_ok (Fastq.phred_log_probs r.Biocaml_unix.Fastq.qualities) in
    let d =
      Digest.bytes (sprintf "%s%s%s"
        (Nomenclature.show_locus parPHMM_t.ParPHMM.locus)
        read
        r.Biocaml_unix.Fastq.qualities)
      |> Digest.to_hex
    in
    let name = sprintf "Testing %s %s %s"
      (Nomenclature.show_locus parPHMM_t.ParPHMM.locus)
      (Util.short_seq read) d
    in
    let f () = ignore (pta.ParPHMM.single ~read ~read_errors false) in
    name, f)
  |> Array.to_list

let () =
  let gen = true in
  let reads = setup_samples ~n:number_of_samples () in
  let tests =
    gen_parPHMM_ts gen
    |> List.map ~f:(fun pt -> time_aggregate pt reads)
    |> List.concat
  in
  let open Core_bench.Std in
  Core.Command.run (Bench.make_command
    (List.map tests ~f:(fun (name, f) -> Bench.Test.create ~name f)))
