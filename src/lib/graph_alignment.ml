
open Util
open Ref_graph

(** Alignment. *)
let align_sequences ~s1 ~o1 ~s2 ~o2 =
  let l1 = String.length s1 in
  let l2 = String.length s2 in
  let rec loop i m =
    let i1 = i + o1 in
    let i2 = i + o2 in
    if i1 >= l1 then
      if i2 >= l2 then
        `Both m
      else
        `First (m, i2)
    else if i2 >= l2 then
      `Second (m, i1)
    else
      loop (i + 1) (if s1.[i1] = s2.[i2] then m else m + 1)
  in
  loop 0 0

(* TODO:
   There are 2 pieces of information about a graph that are relevant
   for alignment:

   1. Which edges a path through the graph traverses, since edges give us alleles
   2. If a path is entirely on a vertex (ex. "AC" in "GGACTT") then we need to
      know if the vertex ("GGACTT") is only found in specific alleles
   *)

type 'e pa =
  { mismatches : int
  ; edge       : 'e             (*BitSet.t        too big to carry around? *)
  ; alignment  : 'e alignment
  }
and 'e alignment = Leaf of int
                 | Partial of 'e pa list

let rec map_alignment f = function
  | Leaf n     -> Leaf n
  | Partial pl -> Partial (List.map pl ~f:(fun p ->
      { p with edge = f p.edge ; alignment = map_alignment f p.alignment }))

(* Add max number of mismatches as argument  *)
let align ?(mub=max_int) g index search_seq =
  let open Nodes in
  let align_against_seq ~search_pos ~node_seq ~node_offset =
    (*let () =
      Printf.printf "align_against_seq search_pos %d \t node_seq %s \t node_offset %d\n"
              search_pos node_seq node_offset
    in*)
    let reached_end_of =
      align_sequences ~s1:search_seq ~o1:search_pos ~s2:node_seq ~o2:node_offset
    in
    match reached_end_of with
    | `Both m
    | `First (m, _)   -> `Finished m      (* end of the search string *)
    | `Second (m, so) -> `GoOn (m, so)
  in
  let over_all_sequence_nodes node m0 ~search_pos =
    let rec descend node m0 ~search_pos =
      (*let () = Printf.printf "descend m0 %d search_pos %d\n" m0 search_pos in *)
      G.fold_succ_e (fun (_, edge, vs) acc ->
          match vs with
          | S _             -> invalid_argf "Another start?"
          | E _             -> acc
          | B _             -> (* This is a weird case since all prev nodes point at a boundary
                                  the acc in this case should be [] *)
                               assert (acc = []);
                               descend vs m0 ~search_pos
          | N (_, node_seq) ->
              match align_against_seq ~search_pos ~node_seq ~node_offset:0 with
              | `Finished mismatches            ->
                  let tm = m0 + mismatches in
                  if tm > mub then acc else
                    { mismatches ; edge; alignment = Leaf tm } :: acc
              | `GoOn (mismatches, search_pos)  ->
                  let tm = m0 + mismatches in
                  if tm > mub then acc else
                    let alignment = Partial (descend vs tm ~search_pos) in
                    { mismatches ; edge; alignment } :: acc)
        g node []
    in
    descend node m0 ~search_pos
  in
  let open Graph_index in
  let km1 = k index - 1 in
  starting_with index search_seq >>= (fun lst ->
      List.filter_map lst ~f:(fun { alignment; sequence; offset} ->
          match align_against_seq ~search_pos:km1 ~node_seq:sequence ~node_offset:offset with
          | `Finished mismatches  ->
              if mismatches > mub then None else
                Some (Leaf mismatches)
          | `GoOn (mismatches, search_pos) ->
              if mismatches > mub then None else
                Some (Partial (over_all_sequence_nodes (N (alignment, sequence)) mismatches ~search_pos)))
      |> fun als -> Ok als)

let name_edges_in_alignment aindex =
  map_alignment (Alleles.Set.to_human_readable aindex)

let rec best_of_paths = function
  | Leaf n      -> n
  | Partial pl  ->
      let bps = List.map ~f:(fun pa -> best_of_paths pa.alignment) pl in
      List.fold_left ~f:min ~init:max_int bps
      (* if pl = [] then there are no alignments that satisfy the mub *)

let to_weights lst =
  let flst = List.map ~f:float_of_int lst in
  let ilst = List.map ~f:(fun x -> 1. /. (1. +. x)) flst in
  let s = List.fold_left ~f:(+.) ~init:0. ilst in
  List.map ~f:(fun x -> x /. s) ilst

let init_alignment_map aindex =
  Alleles.Map.make aindex 0.

(* Weighing Alignments ... inference *)
let alignments_to_weights ?(base_weight=1.) amap al =
  let rec descend w0 = function
  | Leaf _     -> ()
  | Partial pl -> let weights = to_weights (List.map pl ~f:(fun pa -> pa.mismatches)) in
                  List.iter2 pl weights ~f:(fun pa w ->
                    let w_i = w *. w0 in
                    Alleles.Map.update_from pa.edge amap ((+.) w_i);
                    descend w_i pa.alignment)
  in
  descend base_weight al

let most_likely aindex amap =
  Alleles.Map.fold aindex ~init:[] ~f:(fun acc v allele ->
    if v > 0. then (v,allele) :: acc else acc) amap
  |> List.sort ~cmp:(fun ((v1 : float), _) (v2,_) -> compare v2 v1)

