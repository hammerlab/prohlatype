
open Prohlatype

let j =
  Post_analysis.of_json_file
    "FILL ME/res/2017_12_12_full_class1/006.json"

let j_bl =
  Post_analysis.reads_by_loci j

let j_bl_C =
  assoc_exn Nomenclature.C j_bl

let a1 = "C*05:111"
let a2 = "C*05:01:01:02"
let a3 = "C*08:02:01:01"

let ai = Alleles.Input.merge ~distance:Distances.WeightedPerSegment "../foreign/IMGTHLA/alignments/C"


let file1 = "FILL MEfastqs/PGV_006_1.fastq"
let file2 = "FILL MEfastqs/PGV_006_2.fastq"

module Pd = ParPHMM_drivers
module Pa = Post_analysis
module Ml = Pd.Multiple_loci

let char_of_two_ls l1 l2 =
  if l1 = l2 then      'E'
  else if l1 < l2 then 'L'
  else (* l1 > l2 *)   'G'

let pt = ParPHMM.construct ai |> unwrap
let prealigned_transition_model = true
let read_length = 125

let fpt1 =
  ParPHMM.setup_single_allele_forward_pass ~prealigned_transition_model
    read_length a1 pt
let alleles1 = [| a1 |]

let fpt2 =
  ParPHMM.setup_single_allele_forward_pass ~prealigned_transition_model
    read_length a2 pt
let alleles2 = [| a2 |]

let fpt3 =
  ParPHMM.setup_single_allele_forward_pass ~prealigned_transition_model
    read_length a3 pt
let alleles3 = [| a3 |]

open Versus_common

(* We're searching for a specific window. *)
let last_end_pos = 1372 + 12
let first_end_pos = last_end_pos - read_length

let specific_reads =
  List.map j_bl_C ~f:(fun (_, rn, _) -> rn)

let smap =
  List.map j_bl_C ~f:(fun (_, rn, r) -> rn, r)
  |> string_map_of_assoc

module Lp = ParPHMM.LogProbability

let describe lst =
  let sorted =
    (* descending. *)
    List.sort lst ~cmp:(fun (l1, _) (l2, _) ->
          Lp.compare l2 l1)
  in
  let paired =
      List.map2 (List.take sorted 2) (List.drop sorted 1)
        ~f:(fun (l1, a1) (l2, _) ->
              if Lp.close_enough l1 l2 then
                sprintf "%s=" a1
              else
                sprintf "%s>=" a1)
  in
  let str =
    sprintf "%s%s"
      (String.concat ~sep:"" paired)
      (snd (List.hd_exn (List.rev sorted)))
  in
  let values = string_of_list ~sep:";" ~f:(fun (l, _) -> Lp.to_string l) sorted in
  str, values

let do_it rn p rs re rc =
  let l1 = (callhd alleles1 fpt1 rs re rc).llhd in
  let l2 = (callhd alleles2 fpt2 rs re rc).llhd in
  let l3 = (callhd alleles3 fpt3 rs re rc).llhd in
  let s,v = describe [(l1, a1); (l2, a2); (l3, a3)] in
  printf "%s\t%s\t%s\t%s\n%!"
    rn p s v

let paired readname rs1 re1 rs2 re2 =
  let of_alp_list alp = (List.hd_exn alp).Pd.Alleles_and_positions.position in
  let of_aalp_pr = function
    | ParPHMM.Pass_result.Filtered _    -> invalid_argf "read was filtered ?!?"
    | ParPHMM.Pass_result.Completed alp -> of_alp_list alp
  in
  match StringMap.find readname smap with
  | Pa.Pr (Ml.FirstFiltered _) ->
      eprintf "%s was FirstFiltered ?" readname
  | Pa.Pr (Ml.FirstOrientedSecond fos ) ->
      let p1 = of_alp_list (fos.Ml.first) in
      if first_end_pos < p1 && p1 < last_end_pos then
        do_it readname "p 1" rs1 re1 fos.Ml.first_o;
      (*else
        printf "no %s %d %d %d\n%!" readname first_end_pos p1 last_end_pos; *)
      let p2 = of_aalp_pr (fos.Ml.second) in
      if first_end_pos < p2 && p2 < last_end_pos then
        do_it readname "p 2" rs2 re2 (not fos.Ml.first_o)
      (*else
        printf "no %s %d %d %d\n%!" readname first_end_pos p2 last_end_pos *)
  | Pa.Soi _ ->
      eprintf "%s supposed to be paired!" readname

let single rp readname rs re =
  let of_alp_list alp = (List.hd_exn alp).Pd.Alleles_and_positions.position in
  let take_regular r c = Pd.Alleles_and_positions.descending_cmp r c <= 0 in
  let mlo fp = Pd.Orientation.most_likely_between ~take_regular fp in
  match StringMap.find readname smap with
  | Pa.Soi (Ml.SingleRead or_) ->
      begin match mlo or_ with
      | ParPHMM.Pass_result.Filtered _ -> eprintf "%s filtered!" readname
      | ParPHMM.Pass_result.Completed (rc, rplst) ->
          let p = of_alp_list rplst in
          if first_end_pos < p && p < last_end_pos then
            do_it readname (sprintf "S %s" rp) rs re rc
          (*else
            printf "no %s %d %d %d\n%!" readname first_end_pos p last_end_pos *)
      end
  | Pa.Soi (Ml.PairedDependent _) ->
      eprintf "%s not paired dependent!" readname
  | Pa.Pr _ ->
      eprintf "%s not paired!" readname

let () =
  printf "read\ts_o_p\torder\tlikelihoods\n%!";
  Fastq.fold_paired_both
    ~specific_reads
    ~init:()
    ~f:(fun () r1 r2 -> Pd.Fastq_items.paired_untimed r1 r2 ~k:paired)
    ~ff:(fun () r -> Pd.Fastq_items.single_utimed r ~k:(single "1"))
    ~fs:(fun () r -> Pd.Fastq_items.single_utimed r ~k:(single "2"))
    file1
    file2
  |> ignore
