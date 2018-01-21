(* Tests that the merging produces sensible results. *)
open Common
open Prohlatype

let to_input prefix distance =
  let open Alleles.Input in function
  | `Merge    -> merge (to_merge_prefix prefix) ~distance
  | `Genetic  -> alignment (to_alignment_file (prefix ^ "_gen")) ~distance
  | `Nuclear  -> alignment (to_alignment_file (prefix ^ "_nuc")) ~distance

let load prefix d t =
  let arg = Ref_graph.default_construction_arg in
  match Cache.(graph (graph_args ~input:(to_input prefix d t) ~arg)) with
  | Error e -> failwithf "failed to load %s : %s" prefix e
  | Ok g -> g

let split_into_xons s =
  String.split s ~on:(`Character '|')
  |> List.filter ~f:((<>) "")

let list_zip = List.map2 ~f:(fun a b -> (a,b))

let test_same ~merged_seq ~genetic_seq allele =
  let labels = sprintf "merged %s: " allele , sprintf "genetic %s: " allele in
  if merged_seq <> genetic_seq then begin
    let mxs = split_into_xons merged_seq in
    let gxs = split_into_xons genetic_seq in
    let mxs_n = List.length mxs in
    let gxs_n = List.length gxs in
    if mxs_n <> gxs_n then
      error "Merged list for %s doesn't have the same number %d of xon elements as genetic %d\n%s\n%s"
        allele mxs_n gxs_n
          (String.concat ~sep:" | " mxs)
          (String.concat ~sep:" | " gxs)
    else begin
      list_zip mxs gxs
      |> list_fold_ok ~init:() ~f:(fun () (m, g) ->
          if m <> g then
            error "while testing inserting same nucleic sequence in genetic sequence: %s\n%s\n" allele
              (manual_comp_display ~labels m g)
          else
              Ok ())
    end
  end else
    Ok ()

let compare_different_lengths ~s ~b =
  let m = String.length s in
  let n = String.length b in
  if String.compare_substring (s, 0, m) (b, 0, m) = 0 then
    `Left
  else if String.compare_substring (s, 0, m) (b, (n - m), m) = 0 then
    `Right
  else
    `NotEqual

(* Since this is a possibility of choosing the inner segment,
   we have no power to determine a difference. *)
let difference_in_non_coding_region g1 g2 =
  let open Nomenclature in
  match parse g1 with
  | Error _               -> false
  | Ok (_, (One _, _))    -> false
  | Ok (_, (Two _, _))    -> false
  | Ok (_, (Three _, _))  -> false
  | Ok (_, (Four (a1, b1, c1, d1), _))  ->
      match parse g2 with
      | Error _               -> false
      | Ok (_, (One _, _))    -> false
      | Ok (_, (Two _, _))    -> false
      | Ok (_, (Three _, _))  -> false
      | Ok (_, (Four (a2, b2, c2, d2), _))  ->
          a1 = a1 && b1 = b2 && c1 = c2 && d1 <> d2

let test_diff ~merged_seq ~genetic_seq ~nuclear_seq (nuc, gen) =
  let desc = sprintf "%s -> %s" nuc gen in
  let labels s = sprintf "%s %s: " s nuc , sprintf "gen %s: " gen in
  if nuc = gen then begin
    if merged_seq = genetic_seq then
      Ok ()
    else
      error "%s merged_seq <> genetic_seq" desc
  end else if merged_seq = genetic_seq
            && not (difference_in_non_coding_region nuc gen) then
    error "%s merged_seq = genetic_seq" desc
  else if merged_seq = nuclear_seq then
    error "%s merged_seq = nuclear_seq" desc
  else
    let mxs = split_into_xons merged_seq in
    let gxs = split_into_xons genetic_seq in
    let nxs = split_into_xons nuclear_seq in
    let mxs_n = List.length mxs in
    let gxs_n = List.length gxs in
    if mxs_n <> gxs_n then
      error "Merged list for %s doesn't have the same number %d of xon \
        elements as imputed genetic %d %s:\n%s\n%s"
        nuc mxs_n gxs_n gen
        (String.concat ~sep:" | " mxs)
        (String.concat ~sep:" | " gxs)
    else
      let to_type i = if i mod 2 = 0 then "intron" else "exon" in
      list_zip mxs gxs
      |> list_fold_ok ~init:0 ~f:(fun i (m, g) ->
        if i mod 2 = 0 then (* Intron, compare m to g *)
          if m = g then
            Ok (i + 1)
          else begin
            error "while testing %s at %d %s\n%s\n" desc
              i (to_type i) (manual_comp_display ~labels:(labels "mgd") m g)
          end
        else (* Exon, compare to nuclear. *)
          let ex = Option.value (List.nth nxs (i / 2)) ~default:"" in
          if ex = String.empty then
            if m = g then
              Ok (i + 1)
            else begin
              error "while testing %s at %d %s with empty nucleic exon \
                comparing vs genetic:\n%s\n" desc i (to_type i)
                (manual_comp_display ~labels:(labels "mgd") m g)
            end
          else if ex = m then
            Ok (i + 1)
          else begin
            match compare_different_lengths ~s:ex ~b:m with
            | `Left     -> printf "%s %d exon matches to the left.\n" desc (i/2 + 1);
                           Ok (i+1)
            | `Right    -> printf "%s %d exon matches to the right.\n" desc (i/2 + 1);
                           Ok (i+1)
            | `NotEqual -> error "i: %d %d exon for %s, doesn't match nuc:\n%s\n"
                             i (i/2 + 1) nuc (manual_comp_display ~labels:(labels "nuc") ex m)
          end)
      >>= fun _n -> Ok ()

let alters_to_string = function
  | []  -> ""
  | lst -> (string_of_list lst ~sep:"," ~f:MSA.Alteration.to_string)

let is_altered allele graph =
  match List.Assoc.get allele graph.Ref_graph.merge_map with
  | None | Some []  -> None
  | Some (alt :: _) -> Some alt.MSA.Alteration.allele

let () =
  if !Sys.interactive then
    ()
  else begin
    let n = Array.length Sys.argv in
    let prefix = if n < 2 then "A" else Sys.argv.(1) in
    Alter_MSA.Impute.debug := true;
    Alter_MSA.Merge.debug := true;
    Ref_graph.debug := true;
    let open Alleles.Input in
    printf "constructing genetic graph\n";
    let genetic_graph_trie = load prefix Distances.Trie `Genetic in
    printf "constructing nuclear graph\n";
    let nuclear_graph_trie = load prefix Distances.Trie `Nuclear in
    printf "constructing merged graph\n";
    let merged_graph_trie = load prefix Distances.Trie `Merge in
    let comp merged_graph genetic_graph nuclear_graph =
      printf "testing\n";
      List.iter merged_graph.Ref_graph.merge_map ~f:(fun (nuc, alters) ->
        printf "checking: %s: %s\n%!" nuc (alters_to_string alters);
        Ref_graph.sequence ~boundaries:true merged_graph nuc >>= begin fun merged_seq ->
          match alters with
          | [] -> assert false
          | { MSA.Alteration.allele ; _} :: _ ->
              begin match is_altered nuc nuclear_graph with
                | None ->
                    Ref_graph.sequence ~boundaries:true genetic_graph allele >>= fun genetic_seq ->
                      Ref_graph.sequence ~boundaries:true nuclear_graph nuc >>= fun nuclear_seq ->
                        test_diff ~merged_seq ~genetic_seq ~nuclear_seq (nuc, allele)
                | Some other ->
                    printf "\tskipping because %s has an alteration: %s\n" allele other;
                    Ok ()
              end
        end
        |> function
          | Ok () -> printf "equal\n"
          | Error e -> printf ": ERROR %s\n%!" e; exit 1);
      printf "All tests passed for %s\n" prefix
    in
    comp merged_graph_trie genetic_graph_trie nuclear_graph_trie;
    printf "constructing genetic graph\n";
    let genetic_graph_ws = load prefix Distances.WeightedPerSegment `Genetic in
    printf "constructing nuclear graph\n";
    let nuclear_graph_ws = load prefix Distances.WeightedPerSegment `Nuclear in
    printf "constructing merged graph\n";
    let merged_graph_ws = load prefix Distances.WeightedPerSegment `Merge in
    comp merged_graph_ws genetic_graph_ws nuclear_graph_ws
  end
