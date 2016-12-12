
open Util

(* File system logic *)
let make_full_path ?(perm=0o777) path =
  let rec loop p = function
    | []     -> ()
    | h :: t ->
        let pp = Filename.concat p h in
        if not (Sys.file_exists pp) then Unix.mkdir pp perm;
        loop pp t
  in
  (* TODO: Probably need to do something different on Windows. *)
  let start = if Filename.is_implicit path then "" else Filename.dir_sep in
  String.split ~on:(`String Filename.dir_sep) path
  |> List.filter ~f:(fun s -> not (String.is_empty s))
  |> loop start

(*  Desired behavior:
If cache flag override -> compute it
If cache exists -> load it.
  If file exists ->
     read up to date
     if matches, return cached.
     if not
       warn
       create new,
       cache
  file doesn't exist.
    warn and return cached.
read, cache and return. *)

let disk_memoize ?dir ?up_to_date arg_to_string f =
  let dir = Option.value dir ~default:(Filename.get_temp_dir_name ()) in
  fun ?(skip_disk_cache=false) arg ->
    if skip_disk_cache then
      f arg
    else
      let file = Filename.concat dir (arg_to_string arg) in
      let save r =
        if not (Sys.file_exists dir) then make_full_path dir;
        let o = open_out file in
        Marshal.to_channel o r [];
        close_out o;
        r
      in
      let load () =
        let i = open_in file in
        let r = Marshal.from_channel i in
        close_in i;
        r
      in
      if Sys.file_exists file then begin
        let r = load () in
        match up_to_date with
        | None       -> r
        | Some check -> if check arg r then r else save (f arg)
      end else begin
        save (f arg)
      end

type graph_args =
  { input               : Ref_graph.input
  ; which               : Ref_graph.construct_which_args option
  ; join_same_sequence  : bool
  ; remove_reference    : bool
  }

let graph_args ?n ?(join_same_sequence=true) ?(remove_reference=false) ~input () =
  { input
  ; which = n
  ; join_same_sequence
  ; remove_reference
  }

let graph_args_to_string ga =
  sprintf "%s_%s_%b_%b"
    (Ref_graph.input_to_string ga.input)
    (match ga.which with
      | None -> "All"
      | Some w -> Ref_graph.construct_which_args_to_string w)
    ga.join_same_sequence
    ga.remove_reference

let dir = ".cache"

let graph_cache_dir = Filename.concat dir "graphs"

let graph_no_cache
  { input; which; join_same_sequence; remove_reference; _ } =
  Ref_graph.construct_from_file ~join_same_sequence ~remove_reference ?which
    input

(* TODO: We should add date parsing, so that we can distinguish different
   graphs, and make sure that we can type against the most recent. Or compare
   typing across IMGT release dates. *)
let recent_check arg graph =
  let file =
    match arg.input with
    | Ref_graph.AlignmentFile f -> f
    | Ref_graph.MergeFromPrefix (p,_) -> p ^ "_nuc.txt"
  in
  if Sys.file_exists file then begin
    let ic = open_in file in
    match Mas_parser.parse_align_date ic with
    | None    ->
        eprintf "Did not find a Sequence alignment date for %s, will use current.\n" file;
        close_in ic;
        true
    | Some ad ->
        close_in ic;
        let recent = ad = graph.Ref_graph.align_date in
        if recent then recent else begin
          eprintf "Cached sequence alignment graph is not recent, rebuilding.\n";
          false
        end
  end else begin
    eprintf "Alignment file %s is missing, using cached graph.\n" file;
    true
  end

let invalid_arg_on_error action f =
  (fun arg ->
    (* Raise exception on error to avoid caching! *)
    match f arg with
    | Error e -> invalid_argf "Failed to %s: %s" action e;
    | Ok o    -> o)

let graph =
  let dir = Filename.concat (Sys.getcwd ()) graph_cache_dir in
  disk_memoize ~dir ~up_to_date:recent_check graph_args_to_string
    (invalid_arg_on_error "construct graph" graph_no_cache)

type index_args =
  { k          : int
  ; graph_args : graph_args
  }

let index_args_to_string {k; graph_args} =
  sprintf "%d_%s" k (graph_args_to_string graph_args)

let index_cache_dir = Filename.concat dir "indices"

let graph_and_two_index_no_cache {k; graph_args} =
  graph_no_cache graph_args >>= fun gr ->
    let id = Index.create ~k gr in
    Ok (gr, id)

let graph_and_two_index =
  let dir = Filename.concat (Sys.getcwd ()) index_cache_dir in
  disk_memoize ~dir index_args_to_string
    (invalid_arg_on_error "construct graph and index" graph_and_two_index_no_cache)
