
open Util

#require "csv" ;;

module M = Map.Make (struct type t = int * int let compare = compare end)

let add expected_gene mp e =
  let open Nomenclature in
  match e with
  | (_, (_, One d))     -> invalid_argf "found only one digit resolution in result: %d" d
  | (v, (al, Two (a,b)))
  | (v, (al, Three (a, b, _)))
  | (v, (al, Four (a, b, _, _))) ->
      if al <> expected_gene then
        invalid_argf "Found different gene %s than expected: %s" al expected_gene
      else
        match M.find (a, b) mp with
        | exception Not_found -> M.add (a,b) (float_of_string v) mp
        | _current_value      -> mp

let compress expected_gene fname =
  let ic = Scanf.Scanning.from_file fname in
  let add = add expected_gene in
  let rec loop mp =
    try
      let p =
        Scanf.bscanf ic "%s\t%s\t%s\n" (fun v _ a -> v, Nomenclature.parse a |> unwrap_ok)
      in
      loop (add mp p)
    with End_of_file -> Scanf.Scanning.close_in ic; mp
  in
  let oc = open_out (fname ^ ".c") in
  loop M.empty
  |> M.bindings
  |> List.sort ~cmp:(fun (_,v1) (_,v2) -> compare v2 v1)  (* lower is better *)
  |> List.iter ~f:(fun ((a,b), v) -> fprintf oc "%f\t%s*%.2d:%.2d\n" v expected_gene a b);
  close_out oc

