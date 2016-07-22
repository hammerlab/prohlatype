
open Util

let (//) = Filename.concat
let root_dir = "../foreign/IMGTHLA/"

let all_args ?(file="../foreign/IMGTHLA/alignments/A_nuc.txt") () =
  { Cache.alignment_file = file
  ; Cache.which = None
  ; Cache.normalize = true
  }

let test_graph g pos =
  let open Ref_graph in
  let ens, es, stack = sequence_nodes_at g ~pos |> unwrap_ok in
  print_endline (edge_node_set_to_table g.aindex ens) ;;

let test_file file =
  let open Ref_graph in
  let all_args = all_args ~file:(root_dir // "alignments" // (file ^ ".txt")) () in
  let g = Cache.graph all_args in
  let st, en =
     g.bounds
     |> Alleles.Map.fold g.aindex ~init:(max_int, min_int)
        ~f:(fun p sep_lst _allele ->
            List.fold_left sep_lst ~init:p ~f:(fun (st, en) sep ->
              (min st (fst sep.Ref_graph.start)), (max en sep.Ref_graph.end_)))
  in
  for i = st to en - 1 do
    printf "------------testing %d -------------\n" i;
    test_graph g i
  done

let () =
  let n = Array.length Sys.argv in
  if !Sys.interactive then
    ()
  else
    let file = if n <= 1 then "A_nuc" else Sys.argv.(1) in
    test_file file
