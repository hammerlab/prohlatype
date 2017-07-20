
#mod_use "src/scripts/common.ml";;
open Util

let agen = MSA_parser.from_file (Common.to_alignment_file "A_gen") ;;
let anuc = MSA_parser.from_file (Common.to_alignment_file "A_nuc") ;;

let dt mp n =
  let open MSA_parser in
  let al, allele = List.nth_exn mp.alt_elems n in
  allele_distances ~reference:mp.ref_elems ~allele, al

(*let bt_zip mp n m =
  let open MSA_parser in
  let al1, allele1 = List.nth_exn mp.alt_elems n in
  let al2, allele2 = List.nth_exn mp.alt_elems m in
  al1, al2, Zip2.zip2 ~reference:mp.ref_elems ~allele1 ~allele2
  *)

let bt mp n m =
  let open MSA_parser in
  let al1, allele1 = List.nth_exn mp.alt_elems n in
  let al2, allele2 = List.nth_exn mp.alt_elems m in
  al1, al2, allele_distances_between ~reference:mp.ref_elems ~allele1 ~allele2

