
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
  { input   : Alleles.Input.t
  ; arg     : Ref_graph.construction_arg
  }

let graph_args ~arg ~input =
  { input; arg }

let graph_args_to_string {arg; input} =
  sprintf "%s_%s"
    (Alleles.Input.to_string input)
    (Ref_graph.construction_arg_to_string arg)

let dir = ".cache"

let graph_cache_dir = Filename.concat dir "graphs"

let graph_no_cache { input; arg } =
  Ref_graph.construct ~arg input

(* TODO: We should add date parsing, so that we can distinguish different
   graphs, and make sure that we can type against the most recent. Or compare
   typing across IMGT release dates. *)
let recent_check to_input to_date arg dateable =
  let file =
    match to_input arg with
    | Alleles.Input.AlignmentFile (f,_)     -> f
    | Alleles.Input.MergeFromPrefix (p,_,_) -> p ^ "_nuc.txt"
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
        let recent = (ad = to_date dateable) in
        if recent then recent else begin
          eprintf "Cached sequence alignment graph is not recent, rebuilding.\n";
          false
        end
  end else begin
    eprintf "Alignment file %s is missing, using cached object.\n" file;
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
  let up_to_date =
    recent_check (fun i -> i.input) (fun g -> g.Ref_graph.align_date)
  in
  disk_memoize ~dir ~up_to_date graph_args_to_string
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

type par_phmm_args =
  { pinput    : Alleles.Input.t
  ; selectors : Alleles.Selection.t list
  ; read_size : int
  }

let par_phmm_args ~input ~selectors ~read_size =
  { pinput = input; selectors ; read_size }

let par_phmm_args_to_string {pinput; selectors; read_size} =
  sprintf "%s_%s_%d"
    (Alleles.Input.to_string pinput)
    (Alleles.Selection.list_to_string selectors)
    read_size

let par_phmm_cache_dir =
  Filename.concat dir "parhmm"

let par_phmm_no_cache { pinput; selectors ; read_size } =
  ParPHMM2.construct pinput selectors

let par_phmm =
  let dir = Filename.concat (Sys.getcwd ()) par_phmm_cache_dir in
  let up_to_date =
    recent_check (fun i -> i.pinput) (fun p -> p.ParPHMM2.align_date)
  in
  disk_memoize ~dir ~up_to_date par_phmm_args_to_string
    (invalid_arg_on_error "construct parphmm" par_phmm_no_cache)
