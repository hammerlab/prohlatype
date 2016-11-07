(* Tests that the merging produces sensible results. *)
open Common
open Util

let () =
  if !Sys.interactive then () else begin
    let n = Array.length Sys.argv in
    let prefix = if n < 2 then "A" else Sys.argv.(1) in
    let merged_graph =
      Cache.(graph (graph_args ~input:(Ref_graph.MergeFromPrefix (to_merge_prefix prefix)) ()))
    in
    let genetic_graph =
      Cache.(graph (graph_args ~input:(Ref_graph.AlignmentFile (to_alignment_file (prefix ^ "_gen"))) ()))
    in
    List.iter merged_graph.Ref_graph.merge_map ~f:(fun (nuc, gen) ->
      if nuc <> gen then () else begin
        let oe =
          Ref_graph.sequence merged_graph nuc >>= fun merged_seq ->
            Ref_graph.sequence genetic_graph gen >>= fun genetic_seq ->
              if merged_seq <> genetic_seq then begin
                printf "while testing %s vs %s\n%s\n" nuc gen
                  (manual_comp_display merged_seq genetic_seq);
                error "%s vs %s don't match!" nuc gen
              end else
                Ok ()
        in
        begin match oe with
        | Ok () -> ()
        | Error e ->
            printf "%s\n" e;
            exit 1
        end
      end);
    printf "All tests passed for %s\n" prefix
  end
