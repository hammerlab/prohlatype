(* Test that the Parametric PHMM has the correct ordering for some subset
 * of alleles as if they were calculated by a regular pass. *)

open Prohlatype
open Common

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

let reads ~n ?fastq () =
  [ Test_sample.a_reads ~n ?fastq ()
  ; Test_sample.b_reads ~n ?fastq ()
  ; Test_sample.c_reads ~n ?fastq ()
  ]

let setup_one_read parPHMM_t rd =
  let open ParPHMM in
  let module Pt = Partition_map in
  let pta = setup_single_pass ~prealigned_transition_model read_length parPHMM_t in
  let read = rd.Biocaml_unix.Fastq.sequence in
  let read_errors = unwrap (Fastq.phred_log_probs rd.Biocaml_unix.Fastq.qualities) in
  Pass_result.unwrap (pta.single ~read ~read_errors false);
  let aggregate_result = pta.per_allele_llhd () in
  let read_length = String.length read in
  let individual_assoc_arr =
    Array.map parPHMM_t.alleles ~f:(fun (allele, _alts) ->
      let sta = setup_single_allele_forward_pass
                  ~prealigned_transition_model read_length allele parPHMM_t
      in
      Pass_result.unwrap (sta.single ~read ~read_errors false);
      allele, Pt.get (sta.per_allele_llhd ()) 0)
  in
  let aggregate_assoc_arr =
    Array.map2 parPHMM_t.alleles (Pt.to_array aggregate_result)
      ~f:(fun (allele, _alts) lp -> allele, lp)
  in
  aggregate_assoc_arr, individual_assoc_arr

let describe individual_assoc aggregate_assoc =
  let open ParPHMM in
  let mgd = merge_assoc individual_assoc aggregate_assoc
      ~f:(fun allele m -> match m with
            | `First _  -> invalid_argf "Only had individual result for %s\n" allele
            | `Second _ -> invalid_argf "Only had aggregate result for %s\n" allele
            | `Both (il, al) -> Some (il, al))
  in
  printf "allele\tindiv lklhd\tagg lklhd\n";
  List.iter mgd ~f:(fun (allele, (il, al)) ->
      printf "%s\t%0.9f\t%0.9f\n" allele il al)

let describe_sorted individual_assoc aggregate_assoc =
  let open ParPHMM in
  let cmp (al1, l1) (al2, l2) =
    let lc = compare l1 l2 in
    if lc = 0 then
      compare al1 al2
    else
      lc
  in
  let is_st = List.sort ~cmp individual_assoc in
  let as_st = List.sort ~cmp aggregate_assoc in
  printf "allele\tindiv lklhd\tallele\tagg lklhd\n";
  List.iter2 is_st as_st ~f:(fun (al1, il) (al2, al) ->
      printf "%s\t%0.9f\t%s\t%0.9f\n" al1 il al2 al)

exception TestFailure of string

let test_one_read parPHMM_t loci rn (aggregate_assoc_arr, individual_assoc_arr) =
  printf "Testing relative PHMM calculations for: %s in %s \n" rn loci;
  let open ParPHMM in
  let open Oml in
  let module Sd = Statistics.Descriptive in
  let lp_asf x = Lp.as_float ~precision:5 x in
  let aggregate_assoc_arr = Array.map ~f:(fun (a, l) -> a, lp_asf l) aggregate_assoc_arr in
  let individual_assoc_arr = Array.map ~f:(fun (a, l) -> a, lp_asf l) individual_assoc_arr in
  let sa = Array.map ~f:snd in
  (* First check that their spearman correlation is ~ 1.0 *)
  let s =
    Sd.spearman (sa individual_assoc_arr) (sa aggregate_assoc_arr)
  in
  if not (Util.equal_floats ~d:1e-2 1.0 s) then begin
    eprintf "Loci %s Read %s likelihood Spearman not 1.0 but %f\n" loci rn s;
    describe_sorted (Array.to_list individual_assoc_arr)
      (Array.to_list aggregate_assoc_arr);
    raise (TestFailure "")
  end else begin
    let _n =
    Array.fold_left ~init:0 individual_assoc_arr
      ~f:(fun i (allele, il) ->
          let (_allele, al) = aggregate_assoc_arr.(i) in
          if not (Util.equal_floats ~d:1e-2 il al) then
            raise (TestFailure
              (sprintf "Likelihoods for %s not close enought, il %0.9f vs al %0.9f\n"
              allele il al))
          else
            i + 1)
    in
    ()
  end

let () =
  let n = 1 in
  try
    List.iter2 loci_to_test (reads ~n ())
      ~f:(fun loci reads ->
            let phmm = to_parPHMM_t true loci in
            List.iter reads ~f:(fun rd ->
                test_one_read phmm loci rd.Biocaml_unix.Fastq.name
                  (setup_one_read phmm rd)));
    printf "passed\n"
  with TestFailure m ->
    eprintf "%s" m

