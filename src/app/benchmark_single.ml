
open Util
open Benchmark_common

let time_individual ?(seed=12) parPHMM_t reads =
  let alleles =
    Random.init seed;
    Oml.Util.Array.permute ~copy:true parPHMM_t.ParPHMM.alleles
    |> Array.map ~f:fst
    |> Array.sub ~pos:0 ~len:number_of_alleles_to_test
    |> Array.to_list
  in
  let allele_str =
    string_of_list ~sep:"" ~f:(fun x -> x) alleles
  in
  Printf.printf "locus: %s, number alleles: %d, allele_str: %s\n"
    (Nomenclature.show_locus parPHMM_t.ParPHMM.locus)
    (List.length alleles)
    (Digest.bytes allele_str |> Digest.to_hex);
  List.map alleles ~f:(fun allele ->
    ParPHMM.setup_single_allele_forward_pass ~prealigned_transition_model
      read_length allele parPHMM_t)
  |> List.map ~f:(fun pta ->
    List.map reads ~f:(fun r ->
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
      name, f))
  |> List.concat

let () =
  let gen = true in
  let tests =
    gen_parPHMM_ts gen
    |> List.map2 ~f:(fun reads_by_loci pt -> time_individual pt reads_by_loci) reads
    |> List.concat
  in
  let open Core_bench.Std in
  Core.Command.run (Bench.make_command
    (List.map tests ~f:(fun (name, f) -> Bench.Test.create ~name f)))
