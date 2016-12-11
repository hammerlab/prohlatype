
open Util
#mod_use "src/scripts/common.ml";;
#use "src/scripts/merged_sensible_test.ml";;

let prefix = "B"
let gen_mp, nuc_mp, instr =
  Merge_mas.align_from_prefix ("../foreign/IMGTHLA/alignments/" ^ prefix)
  |> unwrap_ok ;;

let ref_instr =
  Merge_mas.map_instr_to_alignments nuc_mp.Mas_parser.reference
    ~gen:gen_mp.Mas_parser.ref_elems
    ~nuc:nuc_mp.Mas_parser.ref_elems instr |> unwrap_ok ;;

let new_ref_elems =
  Merge_mas.align_reference ref_instr ;;

let ref_position_check =
  Merge_mas.reference_positions_align ~seq:("reference:" ^ gen_mp.Mas_parser.reference) ;;

let gen_assoc = gen_mp.Mas_parser.alt_elems ;;
let nuc_assoc = nuc_mp.Mas_parser.alt_elems ;;
let same, just_gen, just_nuc = Merge_mas.same_and_diff ~gen_assoc ~nuc_assoc ;;

let nallele = "B*15:128";;

let (alt_inst, alt_als) = Merge_mas.same_alts instr same |> unwrap_ok ;;

let rdiff = Merge_mas.reference_as_diff ref_instr ;;
let instr_assoc = (gen_mp.Mas_parser.reference, rdiff) :: alt_inst ;;
let (trie, rmap) = Merge_mas.init_trie_and_map instr_assoc |> unwrap_ok ;;

                (* Add the alleles with just nucleic data. *)
let (diff_alt_lst, diff_map_lst) =
  Merge_mas.diff_alts ~verbose:true ~trie ~rmap ~na:just_nuc |> unwrap_ok ;;

let map_lst = List.map same ~f:(fun (a, _, _) -> (a,a))  ;;

let nuc = List.assoc nallele just_nuc ;;
let gene, (allele_resolution, _) = Nomenclature.parse nallele |> unwrap_ok ;;
let closest_allele_res = Nomenclature.Trie.nearest allele_resolution trie ;;
let alg = Merge_mas.RMap.find closest_allele_res rmap  ;;

let closest_allele_str =
  Nomenclature.resolution_and_suffix_opt_to_string ~gene closest_allele_res ;;
let name = sprintf "%s (nuc) -> %s (gen)" nallele closest_allele_str  ;;
let minstr =
  Merge_mas.merge_different_nuc_into_alignment nallele nuc alg |> unwrap_ok ;;

let malgn = Merge_mas.align_different ~verbose:true nallele minstr ;;

let t ~s1 ~s2 n =
  (List.nth (split_into_xons s1) n) = (List.nth (split_into_xons s2) n) ;;

let d ~s1 ~s2 n =
  printf "%s\n"
    (manual_comp_display
      (List.nth_exn (split_into_xons s1) n)
      (List.nth_exn (split_into_xons s2) n)) ;;

(*
let merged_graph = load prefix `Merge  ;;
let genetic_graph = load prefix `Genetic ;;
let nuclear_graph = load prefix `Nuclear ;;

let merged_seq = Ref_graph.sequence ~boundaries:true merged_graph nallele |> unwrap_ok ;;
let nuclear_seq = Ref_graph.sequence ~boundaries:true nuclear_graph nallele |> unwrap_ok ;;
let gallele = closest_allele_str
let genetic_seq = Ref_graph.sequence ~boundaries:true genetic_graph gallele |> unwrap_ok ;;

*)
