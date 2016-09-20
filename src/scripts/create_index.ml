(* Load a graph and index. *)
open Common

let k = 10 ;;

let gall, idxall =
  Cache.graph_and_two_index
    { Cache.k = k ; Cache.g = cache_arg () } ;;

let g_and_idx ?(k=10) ?file ?n () =
  let g2 =
    Cache.graph_and_two_index
      { Cache.k = k; Cache.g = cache_arg ?n ?file () }
  in
  n, g2

let lst = [3; 5; 10; 15; 20; 50; 100; 150; 200; 250]
let idxs = List.map lst ~f:(fun n -> g_and_idx ~n:(Ref_graph.NumberOfAlts n) ())

