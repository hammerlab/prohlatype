(** Common script declarations *)
open Prohlatype

let time s f =
  let n = Sys.time () in
  let r = f () in
  Printf.printf "%-10s:%f\n" s (Sys.time () -. n);
  r

let (//) = Filename.concat
let imgthla_dir =
  try Sys.getenv "IMGTHLA_DIR"
  with _ -> "../foreign/IMGTHLA"

let to_alignment_file f = imgthla_dir // "alignments" // (f ^ ".txt")
let to_fasta_file f = imgthla_dir // "fasta" // (f ^ ".fasta")
let default_file = to_alignment_file "A_nuc"

let to_merge_prefix p = imgthla_dir // "alignments" // p

let common_alignment_files =
  [ "A_nuc"
  ; "B_nuc"
  ; "C_nuc"
  ; "DMA_nuc"
  ; "DMB_nuc"
  ; "DOA_nuc"
  ; "DOB_nuc"
  ; "DPA1_nuc"
  ; "DPB1_nuc"
  ; "DPB2_nuc"
  ; "DQA1_nuc"
  ; "DQB1_nuc"
  ; "DRA_nuc"
  ; "DRB_nuc"
  ]

let all_alignment_files =
  [ "A_gen"
  ; "A_nuc"
  ; "B_gen"
  ; "B_nuc"
  ; "C_gen"
  ; "C_nuc"
(*  ; "../foreign/IMGTHLA/alignments/ClassI_nuc.txt" screw your broken alignment! *)
  ; "DMA_gen"
  ; "DMA_nuc"
  ; "DMB_gen"
  ; "DMB_nuc"
  ; "DOA_gen"
  ; "DOA_nuc"
  ; "DOB_gen"
  ; "DOB_nuc"
  ; "DPA1_gen"
  ; "DPA1_nuc"
  ; "DPB1_gen"
  ; "DPB1_nuc"
  ; "DPB2_gen"
  ; "DPB2_nuc"
  ; "DQA1_gen"
  ; "DQA1_nuc"
  ; "DQB1_gen"
  ; "DQB1_nuc"
  ; "DRA_gen"
  ; "DRA_nuc"
  ; "DRB1_gen"
  ; "DRB3_gen"
  ; "DRB4_gen"
  ; "DRB_nuc"
  ; "E_gen"
  ; "E_nuc"
  ; "F_gen"
  ; "F_nuc"
  ; "G_gen"
  ; "G_nuc"
  ; "HFE_gen"
  ; "HFE_nuc"
  ; "H_gen"
  ; "H_nuc"
  ; "J_gen"
  ; "J_nuc"
  ; "K_gen"
  ; "K_nuc"
  ; "L_gen"
  ; "L_nuc"
  ; "MICA_gen"
  ; "MICA_nuc"
  ; "MICB_gen"
  ; "MICB_nuc"
  ; "P_gen"
  ; "TAP1_gen"
  ; "TAP1_nuc"
  ; "TAP2_gen"
  ; "TAP2_nuc"
  ; "V_gen"
  ; "V_nuc"
  ; "Y_gen"
  ; "Y_nuc"
  ]

module Test_sample = struct

  let fastq = "tools/test_reads.fastq"

  let all_reads ?(fastq=fastq) () =
    Fastq.all fastq

  let filter_spec ?n ?fastq lst =
    let s = string_set_of_list lst in
    let all =
      List.filter (all_reads ?fastq ()) ~f:(fun r ->
        StringSet.mem r.Biocaml_unix.Fastq.name s)
    in
    Option.value_map n ~default:all ~f:(List.take all)

  let a_reads ?n ?fastq () =
    filter_spec ?n ?fastq
      [ "HWI-ST1027:158:C13N2ACXX:4:1101:15995:139901"
      ; "HWI-ST1027:158:C13N2ACXX:5:1102:6092:9466"
      ; "HWI-ST1027:158:C13N2ACXX:5:1102:5678:41173"
      ; "HWI-ST1027:158:C13N2ACXX:4:1103:17479:35211"
      ; "HWI-ST1027:158:C13N2ACXX:4:1104:14960:83084"
      ; "HWI-ST1027:158:C13N2ACXX:4:1106:14619:108126"
      ; "HWI-ST1027:158:C13N2ACXX:4:1106:1479:135949"
      ; "HWI-ST1027:158:C13N2ACXX:4:1107:9360:124927"
      ; "HWI-ST1027:158:C13N2ACXX:4:1107:13096:132310"
      ; "HWI-ST1027:158:C13N2ACXX:4:1107:21263:177080"
      ; "HWI-ST1027:158:C13N2ACXX:4:1203:11689:53138"
      ; "HWI-ST1027:158:C13N2ACXX:4:1205:16853:196986"
      ; "HWI-ST1027:158:C13N2ACXX:4:1206:15603:157097"
      ; "HWI-ST1027:158:C13N2ACXX:4:1207:17079:65295"
      ; "HWI-ST1027:158:C13N2ACXX:4:1207:17763:187731"
      ; "HWI-ST1027:158:C13N2ACXX:4:1301:10531:23561"
      ; "HWI-ST1027:158:C13N2ACXX:5:1102:19242:109752"
      ; "HWI-ST1027:158:C13N2ACXX:5:1102:12920:133090"
      ; "HWI-ST1027:158:C13N2ACXX:5:1103:3198:132189"
      ; "HWI-ST1027:158:C13N2ACXX:5:1103:19164:172315"
      ; "HWI-ST1027:158:C13N2ACXX:5:1103:14508:191061"
      ; "HWI-ST1027:158:C13N2ACXX:5:1104:7066:95697"
      ; "HWI-ST1027:158:C13N2ACXX:4:1306:12781:21419"
      ; "HWI-ST1027:158:C13N2ACXX:5:1106:9020:74168"
      ; "HWI-ST1027:158:C13N2ACXX:4:1307:19028:72846"
      ; "HWI-ST1027:158:C13N2ACXX:4:1307:3693:156172"
      ; "HWI-ST1027:158:C13N2ACXX:4:1307:6855:157232"
      ; "HWI-ST1027:158:C13N2ACXX:4:1308:4454:19821"
      ; "HWI-ST1027:158:C13N2ACXX:5:1108:17624:18419"
      ; "HWI-ST1027:158:C13N2ACXX:5:1108:3646:55403"
      ]

  let b_reads ?n ?fastq () =
    filter_spec ?n ?fastq
      [ "HWI-ST1027:158:C13N2ACXX:5:1101:7937:10753"
      ; "HWI-ST1027:158:C13N2ACXX:4:1102:7125:25298"
      ; "HWI-ST1027:158:C13N2ACXX:4:1102:6901:130524"
      ; "HWI-ST1027:158:C13N2ACXX:4:1103:12568:102485"
      ; "HWI-ST1027:158:C13N2ACXX:4:1104:8149:5506"
      ; "HWI-ST1027:158:C13N2ACXX:4:1104:16705:25876"
      ; "HWI-ST1027:158:C13N2ACXX:4:1104:10339:94795"
      ; "HWI-ST1027:158:C13N2ACXX:4:1105:13021:13881"
      ; "HWI-ST1027:158:C13N2ACXX:4:1106:20858:45055"
      ; "HWI-ST1027:158:C13N2ACXX:4:1106:5434:144306"
      ; "HWI-ST1027:158:C13N2ACXX:4:1107:10571:36948"
      ; "HWI-ST1027:158:C13N2ACXX:4:1107:8746:127691"
      ; "HWI-ST1027:158:C13N2ACXX:4:1107:19615:142678"
      ; "HWI-ST1027:158:C13N2ACXX:4:1108:6312:59720"
      ; "HWI-ST1027:158:C13N2ACXX:4:1108:20782:96141"
      ; "HWI-ST1027:158:C13N2ACXX:4:1108:1479:116695"
      ; "HWI-ST1027:158:C13N2ACXX:4:1108:9611:154721"
      ; "HWI-ST1027:158:C13N2ACXX:4:1108:7022:179847"
      ; "HWI-ST1027:158:C13N2ACXX:4:1201:4223:60474"
      ; "HWI-ST1027:158:C13N2ACXX:4:1206:19230:115278"
      ; "HWI-ST1027:158:C13N2ACXX:4:1206:6062:175930"
      ; "HWI-ST1027:158:C13N2ACXX:4:1208:5087:139711"
      ; "HWI-ST1027:158:C13N2ACXX:4:1301:7242:65231"
      ; "HWI-ST1027:158:C13N2ACXX:4:1302:15441:87024"
      ; "HWI-ST1027:158:C13N2ACXX:4:1303:2290:78668"
      ; "HWI-ST1027:158:C13N2ACXX:4:1303:6352:196096"
      ; "HWI-ST1027:158:C13N2ACXX:5:1104:8547:51577"
      ; "HWI-ST1027:158:C13N2ACXX:4:1304:14330:81218"
      ; "HWI-ST1027:158:C13N2ACXX:4:1304:9581:124107"
      ; "HWI-ST1027:158:C13N2ACXX:5:1104:1713:95950"
      ]

  let c_reads ?n ?fastq () =
    filter_spec ?n ?fastq
      [ "HWI-ST1027:158:C13N2ACXX:4:1101:4388:50134"
      ; "HWI-ST1027:158:C13N2ACXX:4:1101:9239:124375"
      ; "HWI-ST1027:158:C13N2ACXX:4:1101:4015:192279"
      ; "HWI-ST1027:158:C13N2ACXX:5:1101:16349:89427"
      ; "HWI-ST1027:158:C13N2ACXX:4:1101:6237:195040"
      ; "HWI-ST1027:158:C13N2ACXX:4:1101:13938:198850"
      ; "HWI-ST1027:158:C13N2ACXX:4:1103:18201:180078"
      ; "HWI-ST1027:158:C13N2ACXX:4:1104:6594:68622"
      ; "HWI-ST1027:158:C13N2ACXX:4:1105:4076:41865"
      ; "HWI-ST1027:158:C13N2ACXX:4:1106:18100:2226"
      ; "HWI-ST1027:158:C13N2ACXX:4:1107:18092:60452"
      ; "HWI-ST1027:158:C13N2ACXX:4:1107:8655:185906"
      ; "HWI-ST1027:158:C13N2ACXX:4:1108:18968:154819"
      ; "HWI-ST1027:158:C13N2ACXX:4:1201:19363:125146"
      ; "HWI-ST1027:158:C13N2ACXX:4:1201:20686:157860"
      ; "HWI-ST1027:158:C13N2ACXX:4:1203:15502:106406"
      ; "HWI-ST1027:158:C13N2ACXX:4:1204:6841:105946"
      ; "HWI-ST1027:158:C13N2ACXX:4:1206:3265:170163"
      ; "HWI-ST1027:158:C13N2ACXX:4:1207:14600:3171"
      ; "HWI-ST1027:158:C13N2ACXX:4:1208:5092:79650"
      ; "HWI-ST1027:158:C13N2ACXX:4:1301:20015:35846"
      ; "HWI-ST1027:158:C13N2ACXX:4:1303:6209:88459"
      ; "HWI-ST1027:158:C13N2ACXX:4:1304:5363:29560"
      ; "HWI-ST1027:158:C13N2ACXX:5:1104:2289:62810"
      ; "HWI-ST1027:158:C13N2ACXX:4:1305:18220:25626"
      ; "HWI-ST1027:158:C13N2ACXX:4:1305:13891:74526"
      ; "HWI-ST1027:158:C13N2ACXX:4:1306:14850:164292"
      ; "HWI-ST1027:158:C13N2ACXX:5:1106:13650:65412"
      ; "HWI-ST1027:158:C13N2ACXX:5:1106:13165:84883"
      ; "HWI-ST1027:158:C13N2ACXX:5:1106:16542:125493"
      ]

end (* Test_sample *)
