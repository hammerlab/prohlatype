(* ./test.native ../foreign/IMGTHLA/ *)
open Mas_parser
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

let id x = x

let starts_with_start = 
  "Sequence elements different from previous",
  function
  | []     -> true    (* A sequence could be identical to the reference
                         and therefore the parsed alt will be [] *)
  | h :: t ->
      List.fold_left (fun (s, p) n ->
        match p, n with
        | Sequence s1, Sequence s2 when s1.start = s2.start -> 
          begin
            Printf.printf "p %s n %s\n" (to_string id p) (to_string id n);
            (false && s, n)
          end
        | Gap g1, Gap g2 when g1.start = g2.start           ->
          begin
            Printf.printf "p %s n %s\n" (to_string id p) (to_string id n);
            (false && s, n)
          end
        | Unknown u1, Unknown u2 when u1.start = u2.start   ->
          begin
            Printf.printf "p %s n %s\n" (to_string id p) (to_string id n);
            (false && s, n)
          end
        | _ -> (true && s, n))
        (true, h) t
      |> fst

exception TestFailed of string

let check (desc, pred) allele lst = 
  if pred lst then
    () 
  else
    raise (TestFailed (sprintf "%s failed for %s" desc allele))

let all_sequences_in_result f r =
  check f r.reference r.ref_elems;
  List.iter (fun (al, el) -> check f al el) r.alt_elems

let test_result r =
  all_sequences_in_result starts_with_start r

let () =
  let n = Array.length Sys.argv in
  if !Sys.interactive then
    ()
  else if n <= 1 then
    print_endline "Please specify IMGT alignments directory."
  else
    let fname = if n <= 2 then None else Some (Sys.argv.(2)) in
    to_fnames ?fname Sys.argv.(1)
    |> List.iter (fun f ->
        try
          let p = from_file f in
          test_result p;
          printf "parsed and checed %s\n" f
        with e ->
          eprintf "failed to parse %s\n" f;
          raise e)
