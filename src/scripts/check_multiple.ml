(* Find alleles that have discontinuities in their sequence data:
  ex ***ACA***ACA***, and report the alignment start and stop positions. *)
open Common
open Util

let graph ~file ?n () =
  let g = Cache.(graph (graph_args ?n ~file ())) in
  n, g

let more_than_one_sb g =
  Alleles.Map.fold g.Ref_graph.aindex ~f:(fun acc l all ->
    if List.length l > 1 then
      (all, l) :: acc
    else acc)
    ~init:[] g.Ref_graph.bounds

let () =
  all_alignment_files
  |> List.filter_map ~f:(fun file ->
      printf "--%s--\n%!" file;
      let _, g = graph ~file:(to_alignment_file file) () in
      match more_than_one_sb g with
      | [] -> None
      | l -> Some (file, l))
  |> List.iter ~f:(fun (f, l) ->
      printf "in file %s: \n" f;
      List.iter l ~f:(fun (a, slst) ->
        printf "allele\t%s: %s\n"
          a (String.concat ~sep:", "
              (List.map slst ~f:(fun sep ->
                sprintf "%d -> %d" (fst sep.Ref_graph.start) sep.Ref_graph.end_)))))
