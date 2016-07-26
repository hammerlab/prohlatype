
open Util

let (//) = Filename.concat
let root_dir = "../foreign/IMGTHLA/"

let all_args ?(file="../foreign/IMGTHLA/alignments/A_nuc.txt") () =
  { Cache.alignment_file = file
  ; Cache.which = None
  ; Cache.join_same_sequence = true
  }

let test_graph file g pos =
  let open Ref_graph in
  match adjacents_at g ~pos with
  | Error e ->
      printf "Couldn't find adjacents at %d for %s because of %s" pos file e
  | Ok (ens, es, stack) ->
      print_endline (edge_node_set_to_table g.aindex ens)

let test_file file =
  let open Ref_graph in
  let all_args = all_args ~file:(root_dir // "alignments" // (file ^ ".txt")) () in
  let g = Cache.graph all_args in
  let st, en = Ref_graph.range g in
  for i = st to en - 1 do
    printf "------------testing %d -------------\n" i;
    test_graph file g i
  done

let () =
  let n = Array.length Sys.argv in
  if !Sys.interactive then
    ()
  else
    let file = if n <= 1 then "A_nuc" else Sys.argv.(1) in
    test_file file
