(* Tests that the merging produces sensible results. *)
open Common
open Util

let to_input prefix = function
  | `Merge    -> Ref_graph.MergeFromPrefix (to_merge_prefix prefix)
  | `Genetic  -> Ref_graph.AlignmentFile (to_alignment_file (prefix ^ "_gen"))
  | `Nuclear  -> Ref_graph.AlignmentFile (to_alignment_file (prefix ^ "_nuc"))

let load prefix t =
  Cache.(graph (graph_args ~input:(to_input prefix t) ()))

let split_into_xons = String.split ~on:(`Character '|')

let comp ~merged_seq ~genetic_seq ~nuclear_seq (nuc, gen) =
  let labels = sprintf "nuc %s: " nuc , sprintf "gen %s: " gen in
  if nuc = gen then begin
    if merged_seq <> genetic_seq then begin
      let mxs = split_into_xons merged_seq in
      let gxs = split_into_xons genetic_seq in
      let mxs_n = List.length mxs in
      let gxs_n = List.length gxs in
      if mxs_n <> gxs_n then
        error "Merged list for %s doesn't have the same number %d of xon elements as genetic %d %s"
          nuc mxs_n gxs_n gen
      else begin
        List.iter2 mxs gxs ~f:(fun m g ->
          if m <> g then
            printf "while testing %s vs %s\n%s\n" nuc gen
                (manual_comp_display ~labels m g));
        error "%s vs %s don't match!" nuc gen
      end
    end else begin
      Printf.printf "equal!\n";
      Ok ()
    end
  end else
    let mxs = split_into_xons merged_seq in
    let gxs = split_into_xons genetic_seq in
    let nxs = split_into_xons nuclear_seq in
    let mxs_n = List.length mxs in
    let gxs_n = List.length gxs in
    if mxs_n <> gxs_n then
      error "Merged list for %s doesn't have the same number %d of xon elements as genetic %d %s"
        nuc mxs_n gxs_n gen
    else
    (*let nxs_n = List.length nxs in
      Not necessarily a valid test since some sequences (ex 'N') may be missing exons.
    if mxs_n / 2 <> nxs_n then
      error "Merged list for %s doesn't have %d 2x+1 of xon elements as nuclear %d"
        nuc mxs_n nxs_n
    else *)
      let to_type i = if i mod 2 = 0 then "intron" else "exon" in
      List.map2 mxs gxs ~f:(fun m g -> (m, g))
      |> list_fold_ok ~init:0 ~f:(fun i (m, g) ->
          if i mod 2 = 0 then (* Intron, compare m to g *)
            if m = g then Ok (i + 1) else begin
              printf "while testing %s vs %s at %d %s\n%s\n" nuc gen
                i (to_type i)
                (manual_comp_display ~labels m g);
              error "%d element (%s) not equal between %s and %s"
                i (to_type i) nuc gen
            end
          else (* Exon, compare to nuclear. *)
            match List.nth nxs (i / 2) with
            | None    -> error "Couldn't find %dth exon for %s" (i / 2) nuc
            | Some ex ->
                if ex = String.empty then
                  if m = g then Ok (i + 1) else begin
                    printf "while testing %s vs %s at %d %s\n%s\n" nuc gen
                      i (to_type i)
                      (manual_comp_display ~labels m g);
                    error "%d element (%s) not equal between %s and %s"
                      i (to_type i) nuc gen
                  end
                else if ex = m then
                  Ok (i + 1)
                else
                  error "%dth %d exon for %s, doesn't match nuc:\n%s\n" (i / 2) i nuc
                    (manual_comp_display ~labels ex m))
      >>= fun _n -> Ok ()

let () =
  if !Sys.interactive then () else begin
    let n = Array.length Sys.argv in
    let prefix = if n < 2 then "A" else Sys.argv.(1) in
    let merged_graph = load prefix `Merge in
    let genetic_graph = load prefix `Genetic in
    let nuclear_graph = load prefix `Nuclear in
    List.iter merged_graph.Ref_graph.merge_map ~f:(fun ((nuc, gen) as p) ->
      printf "comparing %s vs %s %!" nuc gen;
      Ref_graph.sequence ~boundaries:true merged_graph nuc >>= begin fun merged_seq ->
        Ref_graph.sequence ~boundaries:true genetic_graph gen >>= fun genetic_seq ->
          Ref_graph.sequence ~boundaries:true nuclear_graph nuc >>= fun nuclear_seq ->
            comp ~merged_seq ~genetic_seq ~nuclear_seq p
      end
      |> function
        | Ok () -> ()
        | Error e -> eprintf "%s\n" e; exit 1);
    printf "All tests passed for %s\n" prefix
  end
