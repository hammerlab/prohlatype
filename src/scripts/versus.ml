
open Util

let j4 =
  Post_analysis.of_json_file
    "FILL ME/res/2017_12_12_full_class1/004.json" ;;
let j4_bl =
  Post_analysis.reads_by_loci j4  ;;
let j4_bl_A =
  List.Assoc.get Nomenclature.A j4_bl |> Option.value_exn ~msg:"" ;;
let j4_bl_A1, j4_bl_A2 =
  List.partition ~f:(fun (a, _, _) -> a = "A*24:53") j4_bl_A  ;;

let specific_reads =
  List.map j4_bl_A1 ~f:(fun (_, rn, _) -> rn) ;;

let smap =
  List.map j4_bl_A1 ~f:(fun (_, rn, r) -> rn, r) 
  |> string_map_of_assoc

  
let a1 = "A*24:02:01:01"
let a2 = "A*68:01:02:01"

let ai = Alleles.Input.merge ~distance:Distances.WeightedPerSegment "../foreign/IMGTHLA/alignments/A"

let file1 = "FILL MEfastqs/PGV_004_1.fastq"
let file2 = "FILL MEfastqs/PGV_004_2.fastq"

module Pd = ParPHMM_drivers
module Pa = Post_analysis
module Ml = Pd.Multiple_loci

let char_of_two_ls l1 l2 =
  if l1 = l2 then      'E'
  else if l1 < l2 then 'L'
  else (* l1 > l2 *)   'G'

let pt = ParPHMM.construct ai |> unwrap_ok
let prealigned_transition_model = true
let read_length = 125

let fpt1 =
  ParPHMM.setup_single_allele_forward_pass ~prealigned_transition_model
    read_length a1 pt 

let fpt2 =
  ParPHMM.setup_single_allele_forward_pass ~prealigned_transition_model
    read_length a2 pt

let callhd pt rs re rc =
  let _ = pt.ParPHMM.single rs re rc in
  List.hd_exn (pt.ParPHMM.best_allele_pos 1)

let paired readname rs1 re1 rs2 re2 =
  match StringMap.find readname smap with
  | Pa.Pr (Ml.FirstFiltered _) ->
      eprintf "%s was FirstFiltered ?" readname
  | Pa.Pr (Ml.FirstOrientedSecond { Ml.first_o; _}) ->
      let l11 = callhd fpt1 rs1 re1 first_o in
      let l12 = callhd fpt1 rs2 re2 (not first_o) in
      let l21 = callhd fpt2 rs1 re1 first_o in
      let l22 = callhd fpt2 rs2 re2 (not first_o) in
      printf "%s\tp 1\t%c\t%s\t%d\t%s\t%d\n%!"
        readname (char_of_two_ls l11.ParPHMM.llhd l21.ParPHMM.llhd) 
          (ParPHMM.Lp.to_string l11.ParPHMM.llhd)
          l11.ParPHMM.position
          (ParPHMM.Lp.to_string l21.ParPHMM.llhd) 
          l21.ParPHMM.position;
      printf "%s\tp 2\t%c\t%s\t%d\t%s\t%d\n%!"
        readname (char_of_two_ls l12.ParPHMM.llhd l22.ParPHMM.llhd) 
          (ParPHMM.Lp.to_string l12.ParPHMM.llhd)
          l12.ParPHMM.position
          (ParPHMM.Lp.to_string l22.ParPHMM.llhd) 
          l22.ParPHMM.position
  | Pa.Soi _ ->
      eprintf "%s supposed to be paired!" readname

let single rp readname rs re =
  let take_regular r c = Pd.Alleles_and_positions.descending_cmp r c <= 0 in
  let mlo fp = Pd.Orientation.most_likely_between ~take_regular fp in
  match StringMap.find readname smap with
  | Pa.Soi (Ml.SingleRead or_) ->
      begin match mlo or_ with
      | ParPHMM.Filtered _ -> eprintf "%s filtered!" readname
      | ParPHMM.Completed (rc, _) ->
        let l1 = callhd fpt1 rs re rc in
        let l2 = callhd fpt2 rs re rc in
        printf "%s\ts %s\t%c\t%s\t%d\t%s\t%d\n%!"
        readname
          rp
          (char_of_two_ls l1.ParPHMM.llhd l2.ParPHMM.llhd) 
          (ParPHMM.Lp.to_string l1.ParPHMM.llhd)
          l1.ParPHMM.position
          (ParPHMM.Lp.to_string l2.ParPHMM.llhd) 
          l2.ParPHMM.position;
      end
  | Pa.Soi (Ml.PairedDependent _) ->
      eprintf "%s not paired dependent!" readname
  | Pa.Pr _ ->
      eprintf "%s not paired!" readname

let () =
  printf "read\ts_o_p\tstate\tl: %s\tp: %s\tl %s\tp: %s\n%!"
    a1 a1 a2 a2;
  Fastq.fold_paired_both
    ~specific_reads
    ~init:()
    ~f:(fun () r1 r2 -> Pd.Fastq_items.paired_untimed r1 r2 ~k:paired) 
    ~ff:(fun () r -> Pd.Fastq_items.single_utimed r ~k:(single "1"))
    ~fs:(fun () r -> Pd.Fastq_items.single_utimed r ~k:(single "2"))
    file1
    file2
  |> ignore
