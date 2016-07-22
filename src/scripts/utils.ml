
open Util

let e_n_lst_to_out aindex lst =
  let open Ref_graph in
  List.iter lst ~f:(fun (e, n) ->
    printf "%s -- %s\n"
      (Nodes.vertex_name n)
      (Alleles.Set.to_human_readable ~compress:true aindex e))

