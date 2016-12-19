(* Test that we can detect adjacents at all positions along the graph. *)

open Util
open Common

let adjacents_at_test { Ref_graph.g; aindex; offset; posarr; bounds} ~pos =
  Ref_graph.adjacents_at_private g aindex offset posarr bounds ~pos

let test_graph file g pos =
  let open Ref_graph in
  match adjacents_at_test g ~pos with
  | Error e ->
      printf "Couldn't find adjacents at %d for %s because of %s" pos file e
  | Ok (edge_node_set, _, _) ->
      print_endline (EdgeNodeSet.to_table g.aindex edge_node_set)

let test_graph_fail file g pos =
  let open Ref_graph in
  match adjacents_at_test g ~pos with
  | Error e ->
      invalid_argf "Couldn't find adjacents at %d for %s because of %s" pos file e
  | Ok _ ->
      ()

let test_file file =
  let open Ref_graph in
  let input = AlignmentFile (to_alignment_file file)  in
  let all_args = Cache.graph_args ~input () in
  let g = Cache.graph all_args in
  let st, en = Ref_graph.range g in
  for i = st to en - 1 do
    if i mod 10 = 0 then
      printf "------------testing %d -------------\n" i;
    test_graph_fail file g i
  done

(* Optionally take as input the without-extension name of an alignment
   file to run as the test against.
   ex. A_nuc, F_gen
   default: A_nuc *)
let () =
  let n = Array.length Sys.argv in
  if !Sys.interactive then () else begin
    Ref_graph.adjacents_debug_ref := true;
    let file = if n <= 1 then "A_nuc" else Sys.argv.(1) in
    test_file file
  end
