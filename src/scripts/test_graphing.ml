
open Printf

let () =
  let arg_length = Array.length Sys.argv in
  if !Sys.interactive then
    ()
  else if arg_length <= 1 then
    print_endline "Please specify IMGT alignments directory."
  else
    let dir = Sys.argv.(1) in
    let num_alt_to_add, fname =
      if arg_length > 2 then
        let n = int_of_string Sys.argv.(2) in
        Some n, sprintf "A_gen_w%d" n
      else
        None, "A_gen_all"
    in
    let f = Filename.concat dir "A_gen.txt" in
    let g = To_graph.construct ?num_alt_to_add [] (Mas_parser.from_file f) in
    exit (Ref_graph.output ~short:true fname g)
