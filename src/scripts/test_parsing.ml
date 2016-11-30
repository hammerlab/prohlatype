(* Test the Alignment file parsing logic and invariants about their contained
   sequences:
    - All sequences have a Start, with only Gaps and Boundaries before the start
    - All sequences have an End that may be followed only by Boundaries and Gaps
    - There is an end for every start
    - Sequence elements are different from previous

ex:
  $ ./test_parsing.native [option alignment file]
  
  defaults to using all alignment files in IMGTHLA directory
*)

open Util
open Common
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
  let not_swap s = not (String.is_prefix ~prefix:"." s) in
  let not_readme s =
    (* Get the original, not Sosa strings *)
    not ((StringLabels.lowercase (Filename.chop_extension (Filename.basename s)))
          = "readme")
  in
  Sys.readdir dir
  |> Array.to_list
  |> List.filter ~f:(fun f -> fname_f f && suffix_f f && not_swap f && not_readme f)
  |> List.map ~f:(Filename.concat dir)

let starts_with_start =
  let open Mas_parser in
  "All sequences have a Start, with only Gaps and Boundaries before the start.",
  fun _allele lst ->
    let rec test_loop = function
      | []              -> false
      | Start _ :: _t   -> true
      (* We can have Boundaries and Gaps before the sequence actually starts. *)
      | Gap _ :: t      -> test_loop t
      | Boundary _ :: t -> test_loop t
      | _               -> false
    in
    test_loop lst

let ends_with_end =
  let open Mas_parser in
  "All sequences have an End that may be followed only by Boundaries and Gaps.",
  fun _allele lst ->
    let rec test_loop = function
      | []              -> false
      | End _ :: _t     -> true
      | Gap _ :: t      -> test_loop t
      | Boundary _ :: t -> test_loop t
      | _               -> false
    in
    test_loop (List.rev lst)

exception Double of string

let theres_an_end_for_every_start =
  let open Mas_parser in
  "There is an end for every start",
  fun _allele lst ->
    try
      let c =
        List.fold_left ~f:(fun in_data a ->
            match a with
            | Start _ -> if in_data then raise (Double "start") else true
            | End _   -> if not in_data then raise (Double "end") else false
            | _       -> in_data)
          ~init:false lst
      in
      not c
    with Double s ->
      eprintf "Found double %s" s;
      false

let sequence_have_diff_elemns =
  let open Mas_parser in
  "Sequence elements are different from previous",
  fun _allele lst ->
    match lst with
    | []     -> true    (* A sequence could be identical to the reference
                          and therefore the parsed alt will be [] *)
    | h :: t ->
        List.fold_left ~f:(fun (s, p) n ->
          let false_ () =
            Printf.printf "p %s n %s\n" (al_el_to_string p) (al_el_to_string n);
            (false && s, n)
          in
          match p, n with
          | Sequence s1, Sequence s2 when s1.start = s2.start -> false_ ()
          | Gap g1, Gap g2 when g1.gstart = g2.gstart         -> false_ ()
          | _                                                 -> (true && s, n))
          ~init:(true, h) t
        |> fst

let we_can_parse_the_allele_name =
  let open Mas_parser in
  "We can parse the allele name",
  fun allele _lst ->
    match Nomenclature.parse allele with
    | Error _ -> false
    | Ok _    -> true

exception TestFailed of string

let check (desc, pred) allele lst =
  if pred allele lst then
    ()
  else
    raise (TestFailed (sprintf "%s failed for %s" desc allele))

let all_sequences_in_result f r =
  let open Mas_parser in
  check f r.reference r.ref_elems;
  List.iter ~f:(fun (al, el) -> check f al el) r.alt_elems

let test_result r =
  [ starts_with_start
  ; ends_with_end
  ; sequence_have_diff_elemns
  ; theres_an_end_for_every_start
  ; we_can_parse_the_allele_name
  ]
  |> List.iter ~f:(fun check -> all_sequences_in_result check r)


let () =
  let n = Array.length Sys.argv in
  if !Sys.interactive then
    ()
  else
    let fname = if n <= 1 then None else Some (Sys.argv.(2)) in
    to_fnames ?fname (imgthla_dir // "alignments")
    |> List.iter ~f:(fun f ->
        try
          let p = Mas_parser.from_file f in
          test_result p;
          printf "parsed and checked %s\n" f
        with e ->
          eprintf "failed to parse %s\n" f;
          raise e)
