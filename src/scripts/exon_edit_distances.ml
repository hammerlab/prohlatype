
open Util
(* #mod_use "src/scripts/common.ml";; *)
open Common

let time f =
  let start = Sys.time () in
  let x = f () in
  let stop = Sys.time () in
  Printf.printf "Time: %f s\n" (stop -. start);
  x

let prefix = "B"
let gen_mp, nuc_mp, instr =
  Merge_mas.align_from_prefix ("../foreign/IMGTHLA/alignments/" ^ prefix)
  |> unwrap_ok ;;

let boundary_char = '|'

let split_on_boundaries = String.split ~on:(`Character boundary_char)

let all_seq_by_boundaries ?(verbose=false) mp =
  let open MSA_parser in
  let rseq =
    reference_sequence ~boundary_char mp
    |> split_on_boundaries
  in
  let amp =
    List.map mp.alt_elems ~f:(fun (al, allele) ->
      if verbose then printf "at %s%!\n" al;
      al,
      allele_sequence ~boundary_char ~reference:mp.ref_elems ~allele ()
      |> split_on_boundaries)
  in
  (mp.reference, rseq) :: amp

let dist s1 s2 =
  if String.is_empty s1 then
    String.length s2
  else if String.is_empty s2 then
    String.length s1
  else
    Edist.Needleman_wunsch.distance s1 s2

module DMap = Map.Make (struct
  type t = string * string
  let compare = compare
end)

let canonical a1 a2 =
  min a1 a2 , max a1 a2

let compute_distance_map ?(verbose=false) ?versus dist mp =
  let open MSA_parser in
  let against_allele allele exons map =
    List.fold_left ~init:map ~f:(fun map (allele2, exons2) ->
      if verbose then printf "distance of %s vs %s\n%!" allele allele2;
      let key   = canonical allele allele2 in
      let data  = dist exons exons2 in
      DMap.add key data map)
  in
  let rec loop map = function
    | []                    -> map
    | (allele, exons) :: tl -> loop (against_allele allele exons map tl) tl
  in
  let nuc_seqs = all_seq_by_boundaries ~verbose mp in
  match versus with
  | None    ->
      loop DMap.empty nuc_seqs
  | Some vs ->
      let targets = List.filter nuc_seqs ~f:(fun (a, _) -> StringSet.mem a vs) in
      List.fold_left nuc_seqs ~init:DMap.empty ~f:(fun map (allele, exon) ->
        if StringSet.mem allele vs then (* Assume self-similarity *)
          map
        else
          against_allele allele exon map targets)

let output fname dmap =
  let oc = open_out fname in
  fprintf oc "allele vs, distances \n";
  DMap.iter (fun (a1, a2) lst ->
    fprintf oc "%s:%s,%s\n" a1 a2
      (String.concat ~sep:"," (List.map ~f:(sprintf "%d") lst))) dmap;
  close_out oc

let () =
  if !Sys.interactive then () else
    let n = Array.length Sys.argv in
    if n <= 1 then
      let dmap = time (fun () ->
          compute_distance_map ~verbose:true (List.map2 ~f:dist) nuc_mp)
      in
      output (sprintf "%s_total_similarity.csv" prefix) dmap
    else
      let open MSA_parser in
      let versus =
        List.fold_left gen_mp.alt_elems ~init:(StringSet.singleton gen_mp.reference)
          ~f:(fun s (allele, _) -> StringSet.add allele s)
      in
      let dmap = time (fun () ->
          compute_distance_map ~verbose:true (List.map2 ~f:dist) ~versus
            nuc_mp)
      in
      output (sprintf "%s_gen_similarity.csv" prefix) dmap
