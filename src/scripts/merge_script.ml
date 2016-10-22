
open Util
open Common

let bounded lst =
  let open Mas_parser in
  let rec loop curb curs acc = function
    | []               -> List.rev ((curb, curs) :: acc)
    | Start _ :: t     -> loop curb             curs          acc                   t
    | End _ :: t       -> loop curb             curs          acc                   t
    | Boundary bp :: t -> loop (bp.idx, bp.pos) ""            ((curb, curs) :: acc) t
    | Sequence s :: t  -> loop curb             (curs ^ s.s)  acc                   t
    | Gap _ :: t       -> loop curb             curs          acc                   t
  in
  loop (-1, -1) "" [] lst

let boundary_info = fst
let sequence_info = snd

(* Nuc goes into gen! *)
let zip_align ~gen ~nuc =
  let rec loop acc g n =
    match n with
    | []        -> Ok ((List.rev acc) @
                       (List.map g ~f:(fun gs -> `Fill gs)))
    | ns :: nt  ->
        begin
          match g with
          | []       -> Error n
          | gs :: gt -> if (sequence_info gs) = (sequence_info ns) then
                          loop (`Merge (gs, boundary_info ns) :: acc) gt nt
                        else
                          loop (`Fill gs :: acc) gt n
        end
  in
  loop [] gen nuc

let prefix_from_f s =
  match String.split s ~on:(`Character '_') with
  | [p; _] -> p
  | _      -> invalid_argf "Missing '_' in %s" s

let setup_test (genf, nucf) =
  let prefix = prefix_from_f genf in
  (*assert ((prefix_from_f nucf) = prefix); *)
  let gen_mp = Mas_parser.from_file (to_alignment_file genf) in
  let nuc_mp = Mas_parser.from_file (to_alignment_file nucf) in
  let gen_bp = bounded gen_mp.Mas_parser.ref_elems in
  let nuc_bp = bounded nuc_mp.Mas_parser.ref_elems in
  prefix, gen_bp, nuc_bp

let test () =
  List.map paired_files ~f:(fun p ->
    let prefix, gen_bp, nuc_bp = setup_test p in
    prefix, zip_align ~gen:gen_bp ~nuc:nuc_bp)

let ok_prefix =
  [ "A"; "B"; "C"; "DMA"; "DMB"; "DOA"; "DOB"; "DPA1"; "DPB2"; "DQA1"
  ; "DRA" ; "F"; "G"; "HFE"; "H"; "J"; "K"; "L"; "MICA"; "TAP1"; "V"; "Y"
  ]

let err_prefix =
  [ "DRB", (* Comparing DRB1_gen vs DRB_nuc. 2nd exon. *)
  "CACGTTTCTTGTGGCAGCTTAAGTTTGAATGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTTGCTGGAAAGATGCATCTATAACCAAGAGGAGTCCGTGCGCTTCGACAGCGACGTGGGGGAGTACCGGGCGGTGACGGAGCTGGGGCGGCCTGATGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAGGCGGGCCGCGGTGGACACCTACTGCAGACACAACTACGGGGTTGGTGAGAGCTTCACAGTGCAGCGGCGAG
 4                                                                                                                                                                                                                             |                   |  ||
   CACGTTTCTTGTGGCAGCTTAAGTTTGAATGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTTGCTGGAAAGATGCATCTATAACCAAGAGGAGTCCGTGCGCTTCGACAGCGACGTGGGGGAGTACCGGGCGGTGACGGAGCTGGGGCGGCCTGATGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAGGCGGGCCGCGGTGGACACCTATTGCAGACACAACTACGGGGCTGTGGAGAGCTTCACAGTGCAGCGGCGAG"
  ; "DPB1", (* Mismatch in the 2nd exon *)
    "AGAATTACGTGTACCAGGGACGGCAGGAATGCTACGCGTTTAATGGGACACAGCGCTTCCTGGAGAGATACATCTACAACCGGGAGGAGTACGCGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGCTGCGGAGTACTGGAACAGCCAGAAGGACATCCTGGAGGAGAAGCGGGCAGTGCCGGACAGGGTATGCAGACACAACTACGAGCTGGACGAGGCCGTGACCCTGCAGCGCCGAG
  14         | | |                                                                             |  |                                                        |  |                                     |                    | |                      |  | |  |
     AGAATTACCTTTTCCAGGGACGGCAGGAATGCTACGCGTTTAATGGGACACAGCGCTTCCTGGAGAGATACATCTACAACCGGGAGGAGTTCGTGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGATGAGGAGTACTGGAACAGCCAGAAGGACATCCTGGAGGAGGAGCGGGCAGTGCCGGACAGGATGTGCAGACACAACTACGAGCTGGGCGGGCCCATGACCCTGCAGCGCCGAG"
  ; "DQB1", (* Mismatch in the 1st exon, didn't check the rest.*)
    "ATGTCTTGGAAGAAGTCTTTGCGGATCCCCGGAGACCTTCGGGTAGCAACTGTCACCTTGATGCTGGCGATCCTGAGCTCCTCACTGGCTGAGGGCAGAGACTCTCCCG
  10            |   |                  |        |         |            |    |      |  |  |
     ATGTCTTGGAAAAAGGCTTTGCGGATCCCCGGAGGCCTTCGGGCAGCAACTGTGACCTTGATGCTGTCGATGCTGAGCACCCCAGTGGCTGAGGGCAGAGACTCTCCCG"
  ; "E", (* 6th, last, exon. *)
    "GGAGCGACAGTGCCCAGGGGTCTGAGTCTCACAGCTTGTAA
   0
     GGAGCGACAGTGCCCAGGGGTCTGAGTCTCACAGCTTGTAAAG"
  ; "MICB", (* 2nd exon different. *)
   "AGCCCCACAGTCTTCGTTACAACCTCATGGTGCTGTCCCAGGATGAATCTGTGCAGTCAGGGTTTCTCGCTGAGGGACATCTGGATGGTCAGCCCTTCCTGCGCTATGACAGGCAGAAACGCAGGGCAAAGCCCCAGGGACAGTGGGCAGAAGATGTCCTGGGAGCTAAGACCTGGGACACAGAGACCGAGGACTTGACAGAGAATGGGCAAGACCTCAGGAGGACCCTGACTCATATCAAGGACCAGAAAGGAG
  2                                              |                                                                                                                         |
    AGCCCCACAGTCTTCGTTACAACCTCATGGTGCTGTCCCAGGATGGATCTGTGCAGTCAGGGTTTCTCGCTGAGGGACATCTGGATGGTCAGCCCTTCCTGCGCTATGACAGGCAGAAACGCAGGGCAAAGCCCCAGGGACAGTGGGCAGAAGATGTCCTGGGAGCTGAGACCTGGGACACAGAGACCGAGGACTTGACAGAGAATGGGCAAGACCTCAGGAGGACCCTGACTCATATCAAGGACCAGAAAGGAG"
  ; "TAP2", (* 1st exon, didn't check the rest: the C vs T 5th from end. *)
"ATGCGGCTCCCTGACCTGAGACCCTGGACCTCCCTGCTGCTGGTGGACGCGGCTTTACTGTGGCTGCTTCAGGGCCCTCTGGGGACTTTGCTTCCTCAAGGGCTGCCAGGACTATGGCTGGAGGGGACCCTGCGGCTGGGAGGGCTGTGGGGGCTGCTAAAGCTAAGAGGGCTGCTGGGATTTGTGGGGACACTGCTGCTCCCGCTCTGTCTGGCCACCCCCCTGACTGTCTCCCTGAGAGCCCTGGTCGCGGGGGCCTCACGTGCTCCCCCAGCCAGAGTCGCTTCAGCCCCTTGGAGCTGGCTGCTGGTGGGGTACGGGGCTGCGGGGCTCAGCTGGTCACTGTGGGCTGTTCTGAGCCCTCCTGGAGCCCAGGAGAAGGAGCAGGACCAGGTGAACAACAAAGTCTTGATGTGGAGGCTGCTGAAGCTCTCCAGGCCGGACCTGCCTCTCCTCGTTGCCGCCTTCTTCTTCCTTGTCCTTGCTGTCTTGG
1                                                                                                                            |
ATGCGGCTCCCTGACCTGAGACCCTGGACCTCCCTGCTGCTGGTGGACGCGGCTTTACTGTGGCTGCTTCAGGGCCCTCTGGGGACTTTGCTTCCTCAAGGGCTGCCAGGACTATGGCTGGAGGGGACCCTGCGGCTGGGAGGGCTGTGGGGGCTGCTAAAGCTAAGAGGGCTGCTGGGATTTGTGGGGACACTGCTGCTCCCGCTCTGTCTGGCCACCCCCCTGACTGTCTCCCTGAGAGCCCTGGTCGCGGGGGCCTCACGTGCTCCCCCAGCCAGAGTCGCTTCAGCCCCTTGGAGCTGGCTGCTGGTGGGGTACGGGGCTGCGGGGCTCAGCTGGTCACTGTGGGCTGTTCTGAGCCCTCCTGGAGCCCAGGAGAAGGAGCAGGACCAGGTGAACAACAAAGTCTTGATGTGGAGGCTGCTGAAGCTCTCCAGGCCGGACCTGCCTCTCCTCGTTGCCGCCTTCTTCTTCCTTGTCCTTGCTGTTTTGG"
  ; "P", (* No nuc *) ""
  ]

let align_from_prefix p =
  let gen_mp = Mas_parser.from_file (to_alignment_file (p ^ "_gen")) in
  let nuc_mp = Mas_parser.from_file (to_alignment_file (p ^ "_nuc")) in
  gen_mp, nuc_mp
  (*
  let gen = bounded gen_mp.Mas_parser.ref_elems in
  let nuc = bounded nuc_mp.Mas_parser.ref_elems in
  zip_align ~gen ~nuc >>= fun instr ->
    let all_alleles =
      *)

  (*
let until_boundary lst =
  let open Mas_parser in
  let rec loop acc = function
    | Boundary bp :: t -> Some (bp.idx, bp.pos, acc, t)
    | []               -> invalid_arg "Empty sequence end"
    | Start _ :: _     -> invalid_arg "Start when boundary searching"
    | End _ :: _       -> None
    | Sequence s :: t  -> loop (acc ^ s.s) t
    | Gap _ :: t       -> loop acc t
  in
  loop "" lst


let m ~nuc ~gen =
  let open Mas_parser in
  let rec sync_start n g =
    match n, g with
    | Start sn :: tn
    , Start sg :: tg -> sg - sn, tn, tg
    | _ -> invalid_argf "Reference sequences did not start with starts!"
  (*and sync_rest diff n g =  *)
  in
  sync_start nuc gen
  *)
