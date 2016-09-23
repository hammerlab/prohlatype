
open Util

type graph_args =
  { alignment_file      : string       (* this is really a path, basename below *)
  ; which               : Ref_graph.construct_which_args option
  ; join_same_sequence  : bool
  }

let graph_arg ?n ?(join_same_sequence=true) ~file () =
  { alignment_file = file
  ; which = n
  ; join_same_sequence
  }

let graph_args_to_string { alignment_file; which; join_same_sequence } =
  sprintf "%s_%s_%b"
    (Filename.basename alignment_file)
    (match which with
      | None -> "All"
      | Some w -> Ref_graph.construct_which_args_to_string w)
    join_same_sequence

let dir = ".cache"

let graph_cache_dir = Filename.concat dir "graphs"

let graph_no_cache { alignment_file; which; join_same_sequence } =
  Ref_graph.construct_from_file ~join_same_sequence ?which alignment_file

let graph =
  let dir = Filename.concat (Sys.getcwd ()) graph_cache_dir in
  disk_memoize ~dir graph_args_to_string graph_no_cache

type index_args =
  { k : int
  ; g : graph_args
  }

let index_args_to_string {k; g} =
  sprintf "%d_%s" k (graph_args_to_string g)

let index_cache_dir = Filename.concat dir "indices"

let graph_and_two_index_no_cache {k; g} =
  let gr = graph_no_cache g in
  let id = Index.create ~k gr in
  gr, id

let graph_and_two_index =
  let dir = Filename.concat (Sys.getcwd ()) index_cache_dir in
  disk_memoize ~dir index_args_to_string graph_and_two_index_no_cache

