
open Util

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
      let c1 = String.get_exn s1 i1 in
      let c2 = String.get_exn s2 i2 in
      loop (i + 1) (if c1 = c2 then m else m + 1)
  in
  loop 0 0

let num_mismatches_against_seq s =
  fun ~search_pos ~node_seq ~node_offset ->
    match align_sequences ~s1:s ~o1:search_pos ~s2:node_seq ~o2:node_offset with
    | `Both m
    | `First (m, _)   -> `Finished m      (* end of the search string *)
    | `Second (m, so) -> `GoOn (m, so)

let align_sequences_backwards ~s1 ~o1 ~s2 ~o2 =
  let rec loop i m =
    let i1 = o1 - i in
    let i2 = o2 - i in
    if i1 < 0 then
      if i2 < 0 then
        `Both m
      else
        `First (m, i2)
    else if i2 < 0 then
      `Second (m, i1)
    else
      let c1 = String.get_exn s1 i1 in
      let c2 = String.get_exn s2 i2 in
      loop (i + 1) (if c1 = c2 then m else m + 1)
  in
  loop 0 0

let num_mismatches_against_seq_backwards s =
  fun ~search_pos ~node_seq ->
    (* Can't think of a good usecase for this method to export node_offset.
       If we wanted to run the entire alignment algorithm in reverse ?  *)
    let o2 = String.length node_seq - 1 in
    match align_sequences_backwards ~s1:s ~o1:search_pos ~s2:node_seq ~o2 with
    | `Both m
    | `First (m, _)   -> `Finished m      (* end of the search string *)
    | `Second (m, so) -> `GoOn (m, so)

let debug_ref = ref false
let fail_ref = ref false

module EdgeNodeOffsetQueue = MoreLabels.Set.Make
  (struct
    open Ref_graph

    (* TODO: I'm probably boxing and unboxing too much here. *)
    type t = alignment_position * sequence * int * Edges.t

    let compare (p1, s1, o1, e1) (p2, s2, o2, e2) =
      let r = compare_alignment_position p1 p2 in
      if r <> 0 then r else
        let r = compare_sequence s1 s2 in
        if r <> 0 then r else
          (* We should not have 2 of the same node/edges with different
             offsets, so ignore the offset in the comparator! This has
             the nice consequence that we won't add the same quadruple
             when we consider the same node but different edges into
             that node.*)
          Edges.compare e1 e2
  end)

let at_min_position q =
  let rec loop p q acc =
    if q = EdgeNodeOffsetQueue.empty then
      q, acc
    else
      let (pn, _, _, _) as me = EdgeNodeOffsetQueue.min_elt q in
      match acc with
      | [] -> loop pn (EdgeNodeOffsetQueue.remove me q) [me]
      | _  -> if pn = p then
                loop p (EdgeNodeOffsetQueue.remove me q) (me :: acc)
              else
                q, acc
  in
  loop min_int q []

let compute_mismatches gt search_seq pos =
  let open Ref_graph in
  let open Nodes in
  let open Index in
  let search_str_length = String.length search_seq in
  let mis_map = Alleles.Map.make gt.aindex 0 in
  let nmas = num_mismatches_against_seq search_seq in
  (* We need a globally valid position. *)
  let pos = pos.alignment + pos.offset in
  let rec match_and_add_succ queue (p, node_seq, search_pos, edge) =
    match nmas ~search_pos ~node_seq ~node_offset:0 with
    | `Finished mismatches            ->
        Alleles.Map.update_from edge mis_map ((+) mismatches);
        queue
    | `GoOn (mismatches, search_pos)  ->
        Alleles.Map.update_from edge mis_map ((+) mismatches);
        add_next_sequence_nodes search_pos (N (p, node_seq)) queue
  and add_if_sequence search_pos (_, edge, node) queue =
    if !debug_ref then
      eprintf "Adding to queue %s -> %s\n"
        (Alleles.Set.to_human_readable gt.aindex edge)
        (vertex_name node);
    match node with
    | S _       -> invalid_argf "Asked to add a start node: %s" (vertex_name node)
    | E _       -> queue
    | B _       -> add_next_sequence_nodes search_pos node queue
    | N (p, ns) -> EdgeNodeOffsetQueue.add (p, ns, search_pos, edge) queue
  and add_next_sequence_nodes search_pos node queue =
    G.fold_succ_e (add_if_sequence search_pos) gt.g node queue
  in
  let rec assign_loop q =
    if EdgeNodeOffsetQueue.is_empty q then
      mis_map
    else
      let nq, elst = at_min_position q in
      assign_loop (List.fold_left elst ~init:nq ~f:match_and_add_succ)
  in
  Ref_graph.adjacents_at gt ~pos >>= (fun (edge_node_set, seen_alleles, _) ->
    (* TODO. For now assume that everything that isn't seen has a full mismatch,
       this isn't strictly true since the Start of that allele could be within
       the range of the search str.
       - One approach would be to add the other starts, to the adjacents results.
       *)
    let not_seen = Alleles.Set.complement gt.aindex seen_alleles in
    Alleles.Map.update_from not_seen mis_map (fun _ -> search_str_length);
    let startq =
      EdgeNodeSet.fold edge_node_set ~init:EdgeNodeOffsetQueue.empty
        (* Since the adjacents aren't necessarily at pos we have extra
           bookkeeping at the start of the recursion. *)
        ~f:(fun (edge, node) queue ->
              match node with
              | S _ | E _ | B _ ->
                  invalid_argf "Asked to compute mismatches at %s, not a sequence node"
                    (vertex_name node)
              | N (p, node_seq)  ->
                  let node_offset, start_mismatches =
                    if p <= pos then
                      pos - p, 0
                    else
                      0, p - pos
                  in
                  match nmas ~search_pos:0 ~node_seq ~node_offset with
                  | `Finished mismatches            ->
                      let mm = mismatches + start_mismatches in
                      Alleles.Map.update_from edge mis_map ((+) mm);
                      queue
                  | `GoOn (mismatches, search_pos)  ->
                      let mm = mismatches + start_mismatches in
                      Alleles.Map.update_from edge mis_map ((+) mm);
                      add_next_sequence_nodes search_pos node queue)
    in
    Ok (assign_loop startq))

module Ms = Set.Make (
  struct
    (* offset in search-string/read, Node *)
    type t = int * Ref_graph.Edges.t option
    let compare = compare
  end)

(** The propogate and back fill algorithm. *)
let compute_mismatches_old ?(search_pos_start=0) gt search_seq pos =
  let open Ref_graph in
  let open Nodes in
  let open Index in
  let {g; aindex; bounds; _} = gt in
  let search_str_length = String.length search_seq in
  let mis_map = Alleles.Map.make aindex 0 in
  let nmas = num_mismatches_against_seq search_seq in
  let nmasb = num_mismatches_against_seq_backwards search_seq in
  (* Return the edges that we have seen -> we have to backfill those we haven't seen! *)
  let assign_mismatches_for_node mismatches node =
    let init = Alleles.Set.empty () in
    G.fold_succ_e (fun (_, e, _) ea ->
      if mismatches <> 0 then begin
        if !debug_ref then
          printf "Asgn_msm NODE adding %d mismatches to: %s\n"
            mismatches (Alleles.Set.to_human_readable aindex e);
        Alleles.Map.update_from e mis_map ((+) mismatches)
      end;
      Alleles.Set.union ea e)
      g node init
  in
  let assign_mismatches_to_edge mismatches edge =
    if mismatches <> 0 then begin
      if !debug_ref then
        printf "Asgn_msm EDGE adding %d mismatches to: %s\n"
          mismatches (Alleles.Set.to_human_readable aindex edge);
      Alleles.Map.update_from edge mis_map ((+) mismatches)
    end
  in
  let if_not_visited ?edge search_pos st t e =
    (* TODO: Is this Ms set the most efficient way to perform this Marking? *)
    let k = search_pos, edge in
    if Ms.mem k st then begin
      if !debug_ref then
        printf "already visited! %d %s\n" search_pos
          (match edge with | None -> "no edge" | Some e ->
              Alleles.Set.to_human_readable aindex e);
      t st end else  begin
        if !debug_ref then
          printf "visiting %d %s\n" search_pos
            (match edge with | None -> "no edge" | Some e ->
                Alleles.Set.to_human_readable aindex e);
          e (Ms.add k st)
          end
  in
  let outside_start_end pos allele =
    let sep_lst = Alleles.Map.get aindex bounds allele in
    let inside_start_stop =
      List.exists sep_lst ~f:(fun {start; end_} ->
        (fst start) <= pos && pos < end_)
    in
    not inside_start_stop
  in
  let rec assign_ending_edge_mismatch search_pos edge st =
    if_not_visited ~edge search_pos st (fun x -> x) (fun st ->
      let mismatches = search_str_length - search_pos in
      if !debug_ref then
        printf "Adding %d mismatches to %s for ending edge\n"
          mismatches (Alleles.Set.to_human_readable aindex edge);
      assign_mismatches_to_edge mismatches edge;
      st)
  and assign_starting_backfill_mismatch search_pos edge ea st =
    if_not_visited ~edge search_pos st (fun st -> ea, st) (fun st ->
      if !debug_ref then
        printf "Adding %d mismatches to %s for start backfill\n"
          search_pos (Alleles.Set.to_human_readable aindex edge);
      assign_mismatches_to_edge search_pos edge;
      Alleles.Set.union edge ea, st)
  and next_sequence_node search_pos (_, edge, node) st =
    match node with
    | S _             -> invalid_argf "Another Start?"
    | E _             -> assign_ending_edge_mismatch search_pos edge st
    | B _ as vs       -> G.fold_succ_e (next_sequence_node search_pos) g vs st
    | N (_, node_seq) -> compute_mismatch ~edge search_pos ~node_seq ~node st
  and prev_sequence_node search_pos (node, edge, _) ((ea, st) as stp) =
    match node with
    | S _             -> assign_starting_backfill_mismatch search_pos edge ea st
    | E _             -> invalid_argf "Found an End working backwards!"
    | B _ as vs       -> G.fold_pred_e (prev_sequence_node search_pos) g vs stp
    | N (_, node_seq) -> compute_backfill search_pos ~node_seq ~node edge stp
  and compute_mismatch ?edge ?(node_offset=0) search_pos ~node_seq ~node st =
    if_not_visited ?edge search_pos st (fun x -> x)
      (fun st ->
        let () = if !debug_ref then printf "at %d\t %s\t %d\t" search_pos node_seq node_offset in
        let after_align edge_opt mismatches st f =
          match edge_opt with
          | Some edge ->
              (*assign_mismatches_to_edge mismatches edge; *)
              G.iter_pred_e (fun (_, e, _) ->
                if mismatches <> 0 then begin
                  if !debug_ref then
                    printf "Asgn_msm NODE Edge adding %d mismatches to: %s\n"
                      mismatches (Alleles.Set.to_human_readable aindex e);
                  Alleles.Map.update_from e mis_map ((+) mismatches)
                end)
                g node;
              f st
          | None      ->
              let seen_edges = assign_mismatches_for_node mismatches node in
              if Alleles.Set.all aindex seen_edges then
                f st
              else begin
                let node_seq_length = String.length node_seq in
                let prev_seq_pos = search_pos + node_seq_length in
                let seen_edges_after_backfill, nst =
                  G.fold_succ_e (fun (_me, _e, p) acc ->
                    G.fold_pred_e (fun ((possibly_me, e, _) as v) acc ->
                      if possibly_me = node then acc else begin
                        if !debug_ref then
                          printf "prev seq along %s of %s\n"
                            (Alleles.Set.to_human_readable aindex e)
                            (Nodes.vertex_name possibly_me);
                        prev_sequence_node prev_seq_pos v acc end)
                      g p acc) g node (seen_edges, st)
                in
                let not_seen = Alleles.Set.complement aindex seen_edges_after_backfill in
                if Alleles.Set.is_empty not_seen then
                  f st
                else begin
                  let node_pos = Nodes.position node in
                  Alleles.Set.iter aindex not_seen ~f:(fun allele ->
                    if outside_start_end node_pos allele then begin
                      if !debug_ref then
                        printf "Adding %d mismatches to %s for outside start/stop\n"
                          node_seq_length allele;
                      Alleles.Map.update_spec aindex mis_map allele ((+) node_seq_length);
                      Alleles.Set.clear aindex not_seen allele  (* now we've seen it! *)
                    end);
                  if Alleles.Set.is_empty not_seen then
                    f st
                  else begin
                    let msg =
                      sprintf "Still haven't found all of them! at %s\n only saw: %s\n missing %s"
                        (Nodes.vertex_name node) (Alleles.Set.to_human_readable aindex seen_edges_after_backfill)
                        (Alleles.Set.to_human_readable aindex (Alleles.Set.complement aindex seen_edges_after_backfill))
                    in
                    if !fail_ref then invalid_arg msg;
                    if !debug_ref then print_endline msg;
                    f nst
                  end
                end
              end
        in
      match nmas ~search_pos ~node_seq ~node_offset with
      | `Finished mismatches            ->
          let () = if !debug_ref then printf "Finished %d\n" mismatches in
          after_align edge mismatches st (fun st -> st)
      | `GoOn (mismatches, search_pos)  ->
          let () = if !debug_ref then
            printf "GoOn %d (mismatches) %d (search_pos)\n"
              mismatches search_pos
          in
          after_align edge mismatches st (G.fold_succ_e (next_sequence_node search_pos) g node))
  and compute_backfill search_pos ~node_seq ~node edge (es, st) =
    if_not_visited ~edge search_pos st (fun st -> es, st)
      (fun st ->
        let nes = Alleles.Set.union es edge in
        let () = if !debug_ref then printf "working backwards at %d\t %s\t" search_pos node_seq in
        match nmasb ~search_pos:(search_pos - 1) ~node_seq with
        | `Finished mismatches            ->
            let () = if !debug_ref then
              printf "Finished backwards %d\n" mismatches
            in
            assign_mismatches_to_edge mismatches edge;
            nes, st
        | `GoOn (mismatches, search_pos)  ->
            let () = if !debug_ref then
              printf "GoOn backwards %d (mismatches) %d (search_pos)\n"
                mismatches search_pos
            in
            assign_mismatches_to_edge mismatches edge;
            G.fold_pred_e (prev_sequence_node search_pos) g node (nes, st))
  in
  let { sequence; alignment; offset } = pos in
  let _st =
    let node = N (alignment, sequence) in
    compute_mismatch ~node_offset:offset search_pos_start ~node_seq:sequence
      ~node Ms.empty
  in
  mis_map

let to_weights lst =
  let flst = List.map ~f:float_of_int lst in
  let ilst = List.map ~f:(fun x -> 1. /. (1. +. x)) flst in
  let s = List.fold_left ~f:(+.) ~init:0. ilst in
  List.map ~f:(fun x -> x /. s) ilst

let init_alignment_map aindex =
  Alleles.Map.make aindex 0.

(* Weighing Alignments ... inference *)

let most_likely aindex amap =
  Alleles.Map.fold aindex ~init:[] ~f:(fun acc v allele ->
    if v > 0. then (v,allele) :: acc else acc) amap
  |> List.sort ~cmp:(fun ((v1 : float), _) (v2,_) -> compare v2 v1)

