(* ./test.native ../foreign/IMGTHLA/ *)
open Printf 

let to_suffix = function
  | `gen  -> "gen"
  | `nuc  -> "nuc"
  | `prot -> "prot"

let to_fnames ?(suffix=`gen) dir =
  let ss = to_suffix suffix in
  List.map (fun s -> Filename.concat dir (sprintf "%s_%s.txt" s ss)) 
    [ "alignments/A"
    ; "alignments/B"
    ; "alignments/C"
    ; "alignments/DMA"
    ; "alignments/DMB"
    ; "alignments/DOA"
    ; "alignments/DOB"
    ; "alignments/DPA1"
    ; "alignments/DPB1"
    ; "alignments/DPB2"
    ; "alignments/DQA1"
    ; "alignments/DQB1"
    ; "alignments/DRA"
    ; "alignments/DRB1"
    ; "alignments/DRB3"
    ; "alignments/DRB4"
    ; "alignments/E"
    ; "alignments/F"
    ; "alignments/G"
    ; "alignments/HFE"
    ; "alignments/H"
    ; "alignments/J"
    ; "alignments/K"
    ; "alignments/L"
    ; "alignments/MICA"
    ; "alignments/MICB"
    ; "alignments/P"
    ; "alignments/TAP1"
    ; "alignments/TAP2"
    ; "alignments/V"
    ; "alignments/Y"
    ] 


let () =
  if !Sys.interactive || Array.length Sys.argv = 1 then
    ()
  else
    to_fnames Sys.argv.(1)
    |> List.iter (fun f ->
        let _p = Mas_parser.from_file f in
        Printf.printf "parsed %s\n" f)
