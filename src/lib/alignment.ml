
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

module Segment = struct

  type t =
    { matches     : int
    ; mismatches  : int
    }

  let zero = { matches = 0; mismatches = 0 }
  let mismatches n = { matches = 0; mismatches = n }
  let matches    n = { matches = n; mismatches = 0 }

  let matched ?(n=1) m = { m with matches = m.matches + n }
  let mismatched ?(n=1) m = { m with mismatches = m.mismatches + n }
  let to_string m =
    sprintf "{m : %d, ms: %d}" m.matches m.mismatches

end (* Segment *)

let alignment_over_segments ~segment1 ~offset1 ~segment2 ~offset2 =
  let open Segment in
  let l1 = String.length segment1 in
  let l2 = String.length segment2 in
  let rec loop i m =
    let i1 = i + offset1 in
    let i2 = i + offset2 in
    if i1 >= l1 then
      if i2 >= l2 then
        `Both m
      else
        `First (m, i2)
    else if i2 >= l2 then
      `Second (m, i1)
    else
      let c1 = String.get_exn segment1 i1 in
      let c2 = String.get_exn segment2 i2 in
      if c1 = c2 then
        loop (i + 1) (matched m)
      else
        loop (i + 1) (mismatched m)
  in
  loop 0 zero

module type Alignment_config = sig

  type a                (* Value to track for each allele /edge. *)
  type stop             (* Global (across graph) value to track. *)
  type stop_parameter   (* Parameter to stopping logic. *)

  val zero : a
  val to_string : a -> string

  val stop_init : stop
  val stop_to_string : stop -> string

  val from_alignment : position:int -> Segment.t -> a

  val update : stop_parameter -> stop -> current:a -> aligned:a -> (a * stop) option

end

let debug_ref = ref false

let search_pos_edge_lst_to_string aindex l =
  List.map l ~f:(fun (o, e) -> sprintf "%d, %s" o
    (Alleles.Set.to_human_readable ~compress:true aindex e))
  |> String.concat ~sep:";"
  |> sprintf "[%s]"

module Align (Ag : Alignment_config) = struct

  let align_against s =
    let align = alignment_over_segments ~segment1:s in
    fun ~search_pos ~node_seq ~node_offset ->
      match align ~offset1:search_pos ~segment2:node_seq ~offset2:node_offset with
      | `Both m
      | `First (m, _)   -> `Finished m      (* end of the search string *)
      | `Second (m, so) -> `GoOn (m, so)

  let compute gt stop_parameter search_seq pos =
    let open Ref_graph in
    let open Nodes in
    let open Index in

    let nmas = align_against search_seq in
    let search_str_length = String.length search_seq in
    (* We need a globally valid position, maybe just pass this integer as
      argument instead? *)
    let pos = pos.alignment + pos.offset in

    (* State and updating. *)
    let mis_map = Alleles.Map.make gt.aindex Ag.zero in

    let ag_update = Ag.update stop_parameter in
    let assign ?node edge_set ~search_pos stop sa =
      let aligned = Ag.from_alignment ~position:search_pos sa in
      if !debug_ref then
        eprintf "Assigning %s at %d -> %s to %s because of:%s\n%!"
          (Segment.to_string sa)
          search_pos
          (Ag.to_string aligned)
          (Alleles.Set.to_human_readable ~max_length:20000 gt.aindex edge_set)
          (Option.value ~default:"---" (Option.map node ~f:Nodes.vertex_name));

      Alleles.Map.update_from_and_fold edge_set mis_map ~init:(Some stop)
        ~f:(fun stop_opt current ->
              match stop_opt with
              | None      -> current, None
              | Some stop ->
                  match ag_update ~current ~aligned stop with
                  | None        -> current, None
                  | Some (n, s) -> n, (Some s))
        (* NOTE: One could move the optionality into an 'update_from_and_fold_opt'
          if the underlying fold (in Alleles.Map) had the appropriate
          short-circuiting logic. I think that would be cleaner, because then we
          would stop updating asap, and avoid this double option matching. In
          practice, computing likelihood probably isn't that taxing, so we'll
          fold over the entire edge set. *)
    in

    let rec add_successors queue (node, splst) =
      G.fold_succ_e (add_edge_node splst) gt.g node queue
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
    (* Match a sequence, and add successor nodes to the queue for processing.
       XXX Return an optional current queue, where None mean stop.  *)
    and match_and_add_succ state_opt ((node, splst) as ns) =
      Option.bind state_opt ~f:(fun (stop, queue) ->
        match node with
        | S _               -> invalid_argf "How did a Start get here %s!" (vertex_name node)
        | B _               -> Some (stop, add_successors queue ns)
        | E _               ->
            List.fold_left splst ~init:(Some stop) ~f:(fun stop_opt (search_pos, edge) ->
              Option.bind stop_opt ~f:(fun stop ->
                let sa = Segment.mismatches (search_str_length - search_pos) in
                assign ~node edge ~search_pos stop sa))
            |> Option.map ~f:(fun s -> s, queue)

        | N (_p, node_seq)  ->
            List.fold_left splst ~init:(Some (stop, [])) ~f:(fun stop_opt (search_pos, edge) ->
              Option.bind stop_opt ~f:(fun (stop, acc) ->
                match nmas ~search_pos ~node_seq ~node_offset:0 with
                | `Finished local_mismatches            ->
                    assign ~node edge ~search_pos stop local_mismatches
                    |> Option.map ~f:(fun s -> s, acc)
                | `GoOn (local_mismatches, search_pos)  ->
                    assign ~node edge ~search_pos stop local_mismatches
                    |> Option.map ~f:(fun s -> s, (search_pos, edge) :: acc)))
            |> Option.map ~f:(fun (s, acc) ->
                match acc with
                | [] -> (s, queue)
                | _  -> (s, add_successors queue (node, acc))))
    in
    let rec assign_loop (stop, q) =
      if NodeMapQueue.is_empty q then
        `Finished mis_map
      else
        let nq, elst = NodeMapQueue.at_min_position q in
        match List.fold_left elst ~init:(Some (stop, nq)) ~f:match_and_add_succ with
        | Some sq -> assign_loop sq
        | None    -> `Stopped mis_map
    in
    Ref_graph.adjacents_at gt ~pos >>= begin fun (edge_node_set, seen_alleles, _) ->
      (* TODO. For now assume that everything that isn't seen has a full mismatch,
        this isn't strictly true since the Start of that allele could be within
        the range of the search str.
        - One approach would be to add the other starts, to the adjacents results.

        *)
      let not_seen = Alleles.Set.complement gt.aindex seen_alleles in
      let mismatch_whole_read = Segment.mismatches search_str_length in
      let assigned_not_seen = assign not_seen ~search_pos:0 Ag.stop_init mismatch_whole_read in
      match assigned_not_seen with
      | None ->
          (* This is also a weird case where we may stop aligning because of the
             alleles that we haven't seen. Should we communicate this explicitly to
             the 'Ag' logic? At least, we're communicating this condition via the
             `Stopped | `Finished distinction.  *)
          Ok (`Stopped mis_map)
      | Some stop ->
          let startq_opt =
            EdgeNodeSet.fold edge_node_set ~init:(Some (stop, NodeMapQueue.empty))
              (* Since the adjacents aren't necessarily at pos we have extra
                bookkeeping at the start of the recursion. *)
              ~f:(fun (edge, node) state_opt ->
                    Option.bind state_opt ~f:(fun (stop, queue) ->
                      let assign_start = assign ~node edge ~search_pos:0 stop in
                      match node with
                      | S _              ->
                          invalid_argf "Asked to compute mismatches at %s, not a sequence node"
                            (vertex_name node)
                      | E _              ->
                          assign_start (Segment.mismatches search_str_length)
                          |> Option.map ~f:(fun s -> (s, queue))
                      | B (p, _)         ->
                          let dist = p - pos in
                          if dist <= 0 then
                            Some (stop, add_successors queue (node, [0, edge]))
                          else if dist < search_str_length then begin
                            assign_start (Segment.mismatches dist)
                            |> Option.map ~f:(fun s -> (s, add_successors queue (node, [dist, edge])))
                          end else begin
                            assign_start (Segment.mismatches search_str_length)
                            |> Option.map ~f:(fun s -> (s, queue (* Nothing left to match. *)))
                          end
                      | N (p, node_seq)  ->
                          let nmas_and_assign ~node_offset ~start_mismatches =
                            match nmas ~search_pos:start_mismatches ~node_seq ~node_offset with
                            | `Finished mismatches            ->
                                assign_start (Segment.mismatched ~n:start_mismatches mismatches)
                                |> Option.map ~f:(fun s -> (s, queue))
                            | `GoOn (mismatches, search_pos)  ->
                                assign_start (Segment.mismatched ~n:start_mismatches mismatches)
                                |> Option.map ~f:(fun s -> (s, add_successors queue (node, [search_pos, edge])))
                          in
                          let dist = p - pos in
                          if dist <= 0 then
                            nmas_and_assign ~node_offset:(-dist) ~start_mismatches:0
                          else if dist < search_str_length then
                            nmas_and_assign ~node_offset:0 ~start_mismatches:dist
                          else begin
                            assign_start (Segment.mismatches search_str_length)
                            |> Option.map ~f:(fun s -> (s, queue))
                          end))
          in
          match startq_opt with
          | None        -> Ok (`Stopped mis_map)
          | Some startq -> Ok (assign_loop startq)
    end

end (* Align *)

module Mismatches = Align (struct

  type a = int
  type stop = float                   (* Current sum. *)
  type stop_parameter = int * float   (* Number of edges, max mean. *)

  let zero = 0
  let to_string = sprintf "%d"

  let stop_init = 0.
  let stop_to_string = sprintf "sum %f"

  let from_alignment ~position sa = sa.Segment.mismatches

  let update (number_elements, max_mean) stop ~current ~aligned =
    let newal = current + aligned in
    let nstop = stop +. (float aligned) in
    if nstop /. (float number_elements) > max_mean then
      None
    else
      Some (newal, nstop)

end)

let num_mismatches_against_seq = Mismatches.align_against
let compute_mismatches = Mismatches.compute

module PositionMismatches = Align (struct

  type a = (int * int) list
  type stop = float                   (* Current sum. *)
  type stop_parameter = int * float   (* Number of edges, max mean. *)

  let zero = []
  let to_string l =
    List.map l ~f:(fun (p,m) -> sprintf "(%d,%d)" p m)
    |> String.concat ~sep:"; "

  let stop_init = 0.
  let stop_to_string = sprintf "sum %f"

  let from_alignment ~position sa = [position, sa.Segment.mismatches]

  let update (number_elements, max_mean) stop ~current ~aligned =
    let newal = current @ aligned in
    let nstop = stop +. (float (List.length newal)) in
    if nstop /. (float number_elements) > max_mean then
      None
    else
      Some (newal, nstop)

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
      Ok (`Finished (Segment.mismatches n)))

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

