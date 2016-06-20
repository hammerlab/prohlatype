
open Util

type graph_args =
  { alignment_file  : string       (* this is really a path, basename below *)
  ; which           : To_graph.construct_which_args option
  }

let graph_args_to_string { alignment_file; which } =
  sprintf "%s_%s"
    (Filename.basename alignment_file)
    (match which with
      | None -> "All"
      | Some w -> To_graph.construct_which_args_to_string w)

let dir = ".cache"

let graph_cache_dir = Filename.concat dir "graphs"

let graph =
  let dir = Filename.concat (Sys.getcwd ()) graph_cache_dir in
  disk_memoize ~dir graph_args_to_string
    (fun { alignment_file; which } ->
       To_graph.construct_from_file ?which alignment_file)

type index_args =
  { k : int
  ; g : graph_args
  }

let index_args_to_string {k; g} =
  sprintf "%d_%s" k (graph_args_to_string g)

let index_cache_dir = Filename.concat dir "indices"

let graph_and_two_index =
  let dir = Filename.concat (Sys.getcwd()) index_cache_dir in
  disk_memoize ~dir index_args_to_string
    (fun {k; g} ->
        let ai, gr = graph g in
        let id = Graph_index.create ~k gr in
        gr, ai, id)
