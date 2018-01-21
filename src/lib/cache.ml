
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

let disk_memoize ?dir ?up_to_date ?after_load arg_to_string f =
  let dir = Option.value dir ~default:(Filename.get_temp_dir_name ()) in
  fun ?(skip_disk_cache=false) arg ->
    if skip_disk_cache then
      f arg
    else begin
      let file = Filename.concat dir (arg_to_string arg) in
      let save r = match r with
        | Error e -> r
        | Ok ro   ->
            if not (Sys.file_exists dir) then make_full_path dir;
            let o = open_out file in
            Marshal.to_channel o ro [];
            close_out o;
            r
      in
      let load () =
        let i = open_in file in
        let r = Marshal.from_channel i in
        close_in i;
        Option.iter after_load ~f:(fun f -> f r);
        r
      in
      if Sys.file_exists file then begin
        let r = load () in
        match up_to_date with
        | None       -> r
        | Some check -> if check arg r then Ok r else save (f arg)
      end else begin
        save (f arg)
      end
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
    | Alleles.Input.AlignmentFile { path; _}              -> path
    | Alleles.Input.MergeFromPrefix { prefix_path = p; _} -> p ^ "_nuc.txt"
  in
  if Sys.file_exists file then begin
    let ic = open_in file in
    match MSA.Parser.find_header_lines ic with
    | None    ->
        eprintf "Did not find a Sequence alignment date for %s, will use current.\n" file;
        close_in ic;
        true
    | Some (_, ad) ->
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

let graph =
  let dir = Filename.concat (Sys.getcwd ()) graph_cache_dir in
  let up_to_date =
    recent_check (fun i -> i.input) (fun g -> g.Ref_graph.align_date)
  in
  disk_memoize ~dir ~up_to_date graph_args_to_string
    graph_no_cache

type par_phmm_args =
  { pinput    : Alleles.Input.t
  ; read_size : int
  }

let par_phmm_args ~input ~read_size =
  { pinput = input; read_size }

let par_phmm_args_to_string {pinput; read_size} =
  sprintf "%s_%d"
    (Alleles.Input.to_string pinput)
    read_size

let par_phmm_cache_dir =
  Filename.concat dir "parhmm"

let par_phmm_no_cache { pinput; read_size } =
  ParPHMM.construct pinput

let par_phmm =
  let dir = Filename.concat (Sys.getcwd ()) par_phmm_cache_dir in
  let up_to_date =
    recent_check (fun i -> i.pinput) (fun p -> p.ParPHMM.align_date)
  in
  disk_memoize ~dir ~up_to_date par_phmm_args_to_string par_phmm_no_cache
