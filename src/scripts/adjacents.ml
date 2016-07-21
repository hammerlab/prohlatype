

#use "src/scripts/utils.ml"
#use "src/scripts/create_index.ml"


let test pos =
  let open Ref_graph in
  let ens, es, stack = sequence_nodes_at gall ~pos |> unwrap_ok in
  print_endline (edge_node_set_to_table gall.aindex ens) ;;

let do_test () =
  for i = 1 to 1200 do
    printf "------------testing %d -------------\n" i;
    test i
  done

