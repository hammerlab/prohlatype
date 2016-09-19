
open Util

let rec merge_splst l1 l2 =
  match l1, l2 with
  | [], [] -> []
  | [], ls -> ls
  | ls, [] -> ls
  | h1 :: t1, h2 :: t2 ->
      let o1 = fst h1 in
      let o2 = fst h2 in
      if o1 = o2 then
        (o1, Alleles.Set.union (snd h1) (snd h2)) :: merge_splst t1 t2
      else if o1 < o2 then
        h1 :: merge_splst t1 l2
      else
        h2 :: merge_splst l1 t2

module NodeMapQueue = struct
  include MoreLabels.Map.Make
    (struct

      open Ref_graph
      type t =  Nodes.t
      let compare = Nodes.compare_by_position_first
    end)

  let at_min_position q =
    let open Ref_graph in
    let rec loop p q acc =
      if is_empty q then
        q, acc
      else
        let n, _ as me = min_binding q in
        match acc with
        | [] -> loop (Nodes.position n) (remove n q) [me]
        | _  -> if Nodes.position n = p then
                  loop p (remove n q) (me :: acc)
                else
                  q, acc
    in
    loop min_int q []

  let add_to_queue q key data =
    match find key q with
    | exception Not_found -> add ~key ~data q
    | offlst -> add ~key ~data:(merge_splst offlst data) q

  end

(* A semi_group? *)
module type Alignment_group = sig
  type a
  val zero : a
  val incr : pos:int -> v:int -> a -> a
  val merge : a -> a -> a
  val acc_to_string : a -> string
end

let debug_ref = ref false
let fail_ref = ref false

let search_pos_edge_lst_to_string aindex l =
  List.map l ~f:(fun (o, e) -> sprintf "%d, %s" o
    (Alleles.Set.to_human_readable ~compress:true aindex e))
  |> String.concat ~sep:";"
  |> sprintf "[%s]"

module Align (Ag : Alignment_group) = struct

  let align ~s1 ~o1 ~s2 ~o2 =
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
        loop (i + 1) (if c1 = c2 then m else Ag.incr ~pos:i ~v:1 m)
    in
    loop 0 Ag.zero

  let align_against s =
    fun ~search_pos ~node_seq ~node_offset ->
      match align ~s1:s ~o1:search_pos ~s2:node_seq ~o2:node_offset with
      | `Both m
      | `First (m, _)   -> `Finished m      (* end of the search string *)
      | `Second (m, so) -> `GoOn (m, so)

  let compute gt search_seq pos =
    let open Ref_graph in
    let open Nodes in
    let open Index in
    let search_str_length = String.length search_seq in
    let mis_map = Alleles.Map.make gt.aindex Ag.zero in
    let nmas = align_against search_seq in
    (* We need a globally valid position, maybe just pass this integer as
      argument instead? *)
    let pos = pos.alignment + pos.offset in
    let assign ?node edge_set ~position mismatches =
      if !debug_ref then
        eprintf "Assigning %d to %s because of:%s\n%!" mismatches
          (Alleles.Set.to_human_readable ~max_length:20000 gt.aindex edge_set)
          (Option.value ~default:"---" (Option.map node ~f:Nodes.vertex_name));
      Alleles.Map.update_from edge_set mis_map (Ag.incr ~pos:position ~v:mismatches)
    in
    let merge_assign ?node edge_set ag =
      if !debug_ref then
        eprintf "Merging %s to %s because of:%s\n%!" (Ag.acc_to_string ag)
          (Alleles.Set.to_human_readable ~max_length:20000 gt.aindex edge_set)
          (Option.value ~default:"---" (Option.map node ~f:Nodes.vertex_name));
      Alleles.Map.update_from edge_set mis_map (Ag.merge ag)
    in
    let rec match_and_add_succ queue ((node, splst) as ns) =
      match node with
      | S _               -> invalid_argf "How did a Start get here %s!" (vertex_name node)
      | B _               -> add_successors queue ns
      | E _               ->
          List.iter splst ~f:(fun (search_pos, edge) ->
              assign ~node edge ~position:search_pos (search_str_length - search_pos));
          queue
      | N (_p, node_seq)  ->
          let nsplst =
            List.filter_map splst ~f:(fun (search_pos, edge) ->
              match nmas ~search_pos ~node_seq ~node_offset:0 with
              | `Finished local_mismatches            ->
                  merge_assign ~node edge local_mismatches;
                  None
              | `GoOn (local_mismatches, search_pos)  ->
                  merge_assign ~node edge local_mismatches;
                  Some (search_pos, edge))
          in
          match nsplst with
          | [] -> queue
          | _  -> add_successors queue (node, nsplst)
    and add_edge_node splst (_, edge, node) queue =
      let nsplst =
        List.filter_map splst ~f:(fun (sp, ep) ->
          let i = Alleles.Set.inter edge ep in
          if Alleles.Set.is_empty i then None else Some (sp, i))
      in
      if !debug_ref then begin
        eprintf "Considering adding to queue %s -> %s -> %s\n%!"
          (search_pos_edge_lst_to_string gt.aindex splst)
          (vertex_name node)
          (search_pos_edge_lst_to_string gt.aindex nsplst)
      end;
      NodeMapQueue.add_to_queue queue node nsplst
    and add_successors queue (node, splst) =
      G.fold_succ_e (add_edge_node splst) gt.g node queue
    in
    let rec assign_loop q =
      if NodeMapQueue.is_empty q then
        mis_map
      else
        let nq, elst = NodeMapQueue.at_min_position q in
        assign_loop (List.fold_left elst ~init:nq ~f:match_and_add_succ)
    in
    Ref_graph.adjacents_at gt ~pos >>= (fun (edge_node_set, seen_alleles, _) ->
      (* TODO. For now assume that everything that isn't seen has a full mismatch,
        this isn't strictly true since the Start of that allele could be within
        the range of the search str.
        - One approach would be to add the other starts, to the adjacents results.
        *)
      let not_seen = Alleles.Set.complement gt.aindex seen_alleles in
      (* This is a no-op on empty set but for readability, adding if *)
      if not (Alleles.Set.is_empty not_seen) then
        assign not_seen ~position:0 search_str_length;
      let startq =
        EdgeNodeSet.fold edge_node_set ~init:NodeMapQueue.empty
          (* Since the adjacents aren't necessarily at pos we have extra
            bookkeeping at the start of the recursion. *)
          ~f:(fun (edge, node) queue ->
                match node with
                | S _              ->
                    invalid_argf "Asked to compute mismatches at %s, not a sequence node"
                      (vertex_name node)
                | E _              ->
                    assign ~node edge ~position:0 search_str_length;
                    queue
                | B (p, _)         ->
                    let dist = p - pos in
                    if dist <= 0 then
                      add_successors queue (node, [0, edge])
                    else if dist < search_str_length then begin
                      assign ~node edge ~position:0 dist;
                      add_successors queue (node, [dist, edge])
                    end else begin
                      assign ~node edge ~position:0 search_str_length;
                      queue (* Nothing left to match. *)
                    end
                | N (p, node_seq)  ->
                    let nmas_and_assign ~node_offset ~start_mismatches =
                      let start_ag = Ag.incr ~pos:0 ~v:start_mismatches Ag.zero in
                      match nmas ~search_pos:start_mismatches ~node_seq ~node_offset with
                      | `Finished mismatches            ->
                          merge_assign ~node edge (Ag.merge mismatches start_ag);
                          queue
                      | `GoOn (mismatches, search_pos)  ->
                          merge_assign ~node edge (Ag.merge mismatches start_ag);
                          add_successors queue (node, [search_pos, edge])
                    in
                    let dist = p - pos in
                    if dist <= 0 then
                      nmas_and_assign ~node_offset:(-dist) ~start_mismatches:0
                    else if dist < search_str_length then
                      nmas_and_assign ~node_offset:0 ~start_mismatches:dist
                    else begin
                      assign ~node edge ~position:0 search_str_length;
                      queue
                    end)
      in
      Ok (assign_loop startq))

end (* Align *)

module Mismatches = Align (struct
  type a = int
  let zero = 0
  let incr ~pos ~v m = v + m
  let merge m1 m2 = m1 + m2
  let acc_to_string = sprintf "%d"
end)

let num_mismatches_against_seq = Mismatches.align_against
let compute_mismatches = Mismatches.compute

module PositionMismatches = Align (struct
  type a = (int * int) list
  let zero = []
  let incr ~pos ~v m = (pos, v) :: m
  let merge m1 m2 = List.sort ~cmp:compare (m1 @ m2)
  let acc_to_string l =
    List.map l ~f:(fun (p,m) -> sprintf "(%d,%d)" p m)
    |> String.concat ~sep:"; "

end)

let align_sequences_lst = PositionMismatches.align_against
let compute_mismatches_lst = PositionMismatches.compute

(* This method is a bad strawman... Would probably be much faster to
   go back to the original file and apply the Mas_parser changes to the
   reference. Only for comparison purposes... It probably makes sense to
   parameterize Ref_graph.sequence into a more general 'fold'? *)
let manual_mismatches gt search_seq pos =
  let open Ref_graph in
  let p = pos.Index.alignment + pos.Index.offset in
  let contains_p sep_lst =
    List.exists sep_lst ~f:(fun sep -> (fst sep.start) <= p && p <= sep.end_)
  in
  let n = String.length search_seq in
  let s = sequence ~start:(`Pad p) ~stop:(`Pad n) gt in
  let nmas = num_mismatches_against_seq search_seq ~search_pos:0 in
  Alleles.Map.map gt.aindex gt.bounds ~f:(fun sep_lst allele ->
    if contains_p sep_lst then
      s allele >>= fun graph_seq ->
          Ok (nmas ~node_seq:graph_seq ~node_offset:0)
    else
      Ok (`Finished n))

let to_weights lst =
  let flst = List.map ~f:float_of_int lst in
  let ilst = List.map ~f:(fun x -> 1. /. (1. +. x)) flst in
  let s = List.fold_left ~f:(+.) ~init:0. ilst in
  List.map ~f:(fun x -> x /. s) ilst

(* Weighing Alignments ... inference *)

let most_likely aindex amap =
  Alleles.Map.fold aindex ~init:[] ~f:(fun acc v allele ->
    if v > 0. then (v,allele) :: acc else acc) amap
  |> List.sort ~cmp:(fun ((v1 : float), _) (v2,_) -> compare v2 v1)

