
open Util

let cargs ?(file="../foreign/IMGTHLA/alignments/A_nuc.txt") n =
  { Cache.alignment_file = file
  ; Cache.which = Some (Ref_graph.NumberOfAlts n)
  ; Cache.join_same_sequence = true
  } ;;

let all_args ?(file="../foreign/IMGTHLA/alignments/A_nuc.txt") () =
  { Cache.alignment_file = file
  ; Cache.which = None
  ; Cache.join_same_sequence = true
  } ;;

let graph ?file ?n () =
  let g =
    match n with
    | None   -> Cache.graph (all_args ?file ())
    | Some n -> Cache.graph (cargs ?file n)
  in
  n, g

let all_files =
  [ "../foreign/IMGTHLA/alignments/A_gen.txt"
  ; "../foreign/IMGTHLA/alignments/A_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/B_gen.txt"
  ; "../foreign/IMGTHLA/alignments/B_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/C_gen.txt"
  ; "../foreign/IMGTHLA/alignments/C_nuc.txt"
(*  ; "../foreign/IMGTHLA/alignments/ClassI_nuc.txt" screw your broken alignment! *)
  ; "../foreign/IMGTHLA/alignments/DMA_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DMA_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/DMB_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DMB_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/DOA_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DOA_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/DOB_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DOB_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/DPA1_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DPA1_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/DPB1_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DPB1_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/DPB2_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DPB2_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/DQA1_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DQA1_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/DQB1_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DQB1_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/DRA_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DRA_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/DRB1_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DRB3_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DRB4_gen.txt"
  ; "../foreign/IMGTHLA/alignments/DRB_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/E_gen.txt"
  ; "../foreign/IMGTHLA/alignments/E_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/F_gen.txt"
  ; "../foreign/IMGTHLA/alignments/F_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/G_gen.txt"
  ; "../foreign/IMGTHLA/alignments/G_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/HFE_gen.txt"
  ; "../foreign/IMGTHLA/alignments/HFE_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/H_gen.txt"
  ; "../foreign/IMGTHLA/alignments/H_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/J_gen.txt"
  ; "../foreign/IMGTHLA/alignments/J_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/K_gen.txt"
  ; "../foreign/IMGTHLA/alignments/K_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/L_gen.txt"
  ; "../foreign/IMGTHLA/alignments/L_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/MICA_gen.txt"
  ; "../foreign/IMGTHLA/alignments/MICA_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/MICB_gen.txt"
  ; "../foreign/IMGTHLA/alignments/MICB_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/P_gen.txt"
  ; "../foreign/IMGTHLA/alignments/TAP1_gen.txt"
  ; "../foreign/IMGTHLA/alignments/TAP1_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/TAP2_gen.txt"
  ; "../foreign/IMGTHLA/alignments/TAP2_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/V_gen.txt"
  ; "../foreign/IMGTHLA/alignments/V_nuc.txt"
  ; "../foreign/IMGTHLA/alignments/Y_gen.txt"
  ; "../foreign/IMGTHLA/alignments/Y_nuc.txt"
  ]

let more_than_one_sb g =
  Alleles.Map.fold g.Ref_graph.aindex ~f:(fun acc l all ->
    if List.length l > 1 then
      (all, l) :: acc
    else acc)
    ~init:[] g.Ref_graph.bounds

let () =
  all_files
  |> List.filter_map ~f:(fun file ->
      printf "--%s--\n%!" file;
      let _, g = graph ~file () in
      match more_than_one_sb g with
      | [] -> None
      | l -> Some (file, l))
  |> List.iter ~f:(fun (f, l) ->
      printf "in file %s: \n" f;
      List.iter l ~f:(fun (a, slst) ->
        printf "allele\t%s: %s\n"
          a (String.concat ~sep:", "
              (List.map slst ~f:(fun sep ->
                sprintf "%d -> %d" (fst sep.Ref_graph.start) sep.Ref_graph.end_)))))
