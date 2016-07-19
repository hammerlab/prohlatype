
open Util

let e_n_lst_to_out aindex lst =
  let open Ref_graph in
  List.iter lst ~f:(fun (e, n) ->
    printf "%s -- %s\n"
      (Nodes.vertex_name n)
      (Alleles.Set.to_human_readable ~compress:true aindex e))

let nodesset_to_string s =
  let open Ref_graph in
  NodesSet.fold ~f:(fun n a -> Nodes.vertex_name n :: a) s ~init:[]
  |> String.concat ~sep:"; "
  |> sprintf "{%s}"


