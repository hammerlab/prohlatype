

let () =
  let fname = Sys.argv.(1) in
  let outfn = Filename.basename fname |> Filename.chop_extension in
  let r = Parser.parse_f fname in
  let g = To_graph.construct r in
  To_graph.output ~dot:true ~_open:false ~short:true outfn g
  |> ignore
