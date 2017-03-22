(* Typing via a Parametric PHMM. *)

let mp = Mas_parser.from_file (Common.to_alignment_file "A_gen") ;;

let spec =
  [ "A*01:01:01:01"
  ; "A*01:01:01:02N"
  ; "A*01:01:38L"
  ; "A*11:01:01:01"
  ; "A*02:01:01:01"
  ; "A*02:264"
  ]

let build ?read_size ?len ?n ?spec mp =
  let alleles, rlarr = ParPHMM.build_allele_and_rls ?n ?spec mp in
  let rlarr = match len with | None -> rlarr | Some len -> Array.sub rlarr ~pos:0 ~len in
  let conf, tlr = ParPHMM.build_matrix ?read_size rlarr in
  alleles,
  rlarr,
  conf,
  tlr,
  ParPHMM.forward_pass conf
