(* ./test.native ../foreign/IMGTHLA/ *)
open Printf
module String = Sosa.Native_string

let to_suffix = function
  | `gen  -> "gen.txt"
  | `nuc  -> "nuc.txt"
  | `prot -> "prot.txt"

let to_fnames ?fname ?suffix dir =
  let fname_f = match fname with | None -> fun _ -> true | Some s -> (=) s in
  let suffix_f =
    match suffix with
    | None -> fun _ -> true
    | Some t ->
        let suffix = to_suffix t in
        String.is_suffix ~suffix
  in
  let not_swap_f s = not (String.is_prefix ~prefix:"." s) in
  Sys.readdir dir
  |> Array.to_list
  |> List.filter (fun f -> fname_f f && suffix_f f && not_swap_f f)
  |> List.map (Filename.concat dir)

let () =
  let n = Array.length Sys.argv in
  if !Sys.interactive then
    ()
  else if n <= 1 then
    print_endline "Please specify IMGT alignments directory."
  else
    let fname = if n <= 2 then None else Some (Sys.argv.(2)) in
    to_fnames ?fname ~suffix:`nuc Sys.argv.(1)
    |> List.iter (fun f ->
        try
          let _p = Mas_parser.from_file f in
          printf "parsed %s\n" f
        with e ->
          eprintf "failed to parse %s\n" f;
          raise e)
