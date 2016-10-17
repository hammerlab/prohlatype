(* Test that we can detect adjacents at all positions along the graph. *)

open Util
open Common

let test_graph file g pos =
  let open Ref_graph in
  match adjacents_at g ~pos with
  | Error e ->
      printf "Couldn't find adjacents at %d for %s because of %s" pos file e
  | Ok (ens, es, stack) ->
      print_endline (edge_node_set_to_table g.aindex ens)

let test_file file =
  let open Ref_graph in
  let all_args = Cache.graph_args ~file:(to_alignment_file file) () in
  let g = Cache.graph all_args in
  let st, en = Ref_graph.range g in
  for i = st to en - 1 do
    printf "------------testing %d -------------\n" i;
    test_graph file g i
  done

(* Optionally take as input the without-extension name of an alignment
   file to run as the test against.
   ex. A_nuc, F_gen
   default: A_nuc *)
let () =
  let n = Array.length Sys.argv in
  if !Sys.interactive then
    ()
  else
    let file = if n <= 1 then "A_nuc" else Sys.argv.(1) in
    test_file file
