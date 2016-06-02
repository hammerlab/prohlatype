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

let starts_with_start = 
  "All sequences start with a Start",
  fun lst ->
    let rec test_loop = function
      | []              -> false
      | Start _ :: _t   -> true
      (* We can have Boundaries before the sequence actually starts. *)
      | Boundary _ :: t -> test_loop t
      | _               -> false
    in
    test_loop lst

let ends_with_end =
  "All sequences end with an End",
  fun lst ->
    let rec test_loop = function
      | []              -> false
      | End _ :: _t     -> true
      | Boundary _ :: t -> test_loop t
      | _               -> false
    in
    test_loop (List.rev lst)

exception Double of string

let theres_an_end_for_every_start =
  "There is an end for every start",
  fun lst ->
    try
      let c =
        List.fold_left (fun in_data a ->
            match a with
            | Start _ -> if in_data then raise (Double "start") else true
            | End _   -> if not in_data then raise (Double "end") else false
            | _       -> in_data)
          false lst
      in
      not c
    with Double s ->
      eprintf "Found double %s" s;
      false

let sequence_have_diff_elemns =
  "Sequence elements different from previous",
  function
  | []     -> true    (* A sequence could be identical to the reference
                         and therefore the parsed alt will be [] *)
  | h :: t ->
      List.fold_left (fun (s, p) n ->
        let false_ () = 
          Printf.printf "p %s n %s\n" (al_el_to_string p) (al_el_to_string n);
          (false && s, n)
        in
        match p, n with
        | Sequence s1, Sequence s2 when s1.start = s2.start -> false_ ()
        | Gap g1, Gap g2 when g1.start = g2.start           -> false_ ()
        | _                                                 -> (true && s, n))
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
  [ starts_with_start
  ; ends_with_end
  ; sequence_have_diff_elemns 
  ; theres_an_end_for_every_start ]
  |> List.iter (fun check -> all_sequences_in_result check r)

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
