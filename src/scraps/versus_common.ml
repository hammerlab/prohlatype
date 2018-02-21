
open Prohlatype

let callhd alleles pt rs re rc =
  let _ = pt.ParPHMM.single rs re rc in
  let mt = pt.ParPHMM.per_allele_llhd_and_pos () in
  let la i = alleles.(i) in
  List.hd_exn (ParPHMM_drivers.Alleles_and_positions.of_mt la 1 mt)


