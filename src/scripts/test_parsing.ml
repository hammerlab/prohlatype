(* ./test.native ../foreign/IMGTHLA/ *)
open Printf 
module String = Sosa.Native_string

let to_suffix = function
  | `gen  -> "gen.txt"
  | `nuc  -> "nuc.txt"
  | `prot -> "prot.txt"

let to_fnames ?suffix dir =
  let filter =
      match suffix with
      | None ->
        fun l -> l
      | Some t ->
        let suffix = to_suffix t in
        List.filter (String.is_suffix ~suffix)
  in
  Sys.readdir dir
  |> Array.to_list
  |> filter
  |> List.map (Filename.concat dir)

let () =
  if !Sys.interactive || Array.length Sys.argv = 1 then
    ()
  else
    to_fnames Sys.argv.(1)
    |> List.iter (fun f ->
        try
        let _p = Mas_parser.from_file f in
        Printf.printf "parsed %s\n" f
        with e ->
          Printf.eprintf "failed to parse %s\n" f;
          raise e)
