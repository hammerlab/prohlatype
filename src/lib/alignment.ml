
open Util

module NodeMapQueue(V : sig
  type a
  val merge : a -> a -> a
end) : sig

  type 'a t

  type key = Ref_graph.Nodes.t

  val at_min_position : V.a t -> V.a t * (key * V.a) list

  val add_to_queue : V.a t -> key -> V.a -> V.a t

  val empty : V.a t

  val is_empty : V.a t -> bool

end = struct

  include Map.Make
    (struct
      type t = Ref_graph.Nodes.t
      let compare = Ref_graph.Nodes.compare_by_position_first
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
    | offlst -> add ~key ~data:(V.merge offlst data) q

  end (* NodeMapQueue *)

module Nmq = NodeMapQueue (struct

  type a = (int * Alleles.Set.t) list

  let rec merge l1 l2 =
    match l1, l2 with
    | [], [] -> []
    | [], ls -> ls
    | ls, [] -> ls
    | h1 :: t1, h2 :: t2 ->
        let o1 = fst h1 in
        let o2 = fst h2 in
        if o1 = o2 then
          (o1, Alleles.Set.union (snd h1) (snd h2)) :: merge t1 t2
        else if o1 < o2 then
          h1 :: merge t1 l2
        else
          h2 :: merge l1 t2
end)


(* How we measure the quality of alignment over a segment, a segment refers to
   a contiguous sequence, whether in the node of our graph, or to other
   potential paths due to start/end offsets.  *)
module type Segment_alignment_measure = sig

  type t
  val mismatches : int -> t
  val to_string : t -> string

  type thread
  val thread_length : thread -> int
  val thread_update : thread -> int -> char -> t -> t

end

module MakeThreadAlignment (Am : Segment_alignment_measure) : sig

  type segment_result =
    | Both of Am.t
    | First of Am.t * int
    | Second of Am.t * int

  val alignment : ?start:Am.t -> thread:Am.thread -> thread_offset:int -> segment:string ->
    offset:int -> segment_result

end = struct

  type segment_result =
    | Both of Am.t
    | First of Am.t * int
    | Second of Am.t * int

  let alignment ?start ~thread ~thread_offset ~segment ~offset =
    let l1 = Am.thread_length thread in
    let l2 = String.length segment in
    let rec loop i m =
      let i1 = i + thread_offset in
      let i2 = i + offset in
      if i1 >= l1 then
        if i2 >= l2 then
          Both m
        else
          First (m, i2)
      else if i2 >= l2 then
        Second (m, i1)
      else
        let c2 = String.get_exn segment i2 in
        loop (i + 1) (Am.thread_update thread i1 c2 m)
    in
    let start = Option.value start ~default:(Am.mismatches 0) in
    loop 0 start

end (* MakeThreadAlignment *)

module type Alignment_config = sig

  include Segment_alignment_measure

  (* Value to track for each allele /edge. *)
  type a
  val init : a
  val aligned_to_string : a -> string

  (* Global (across graph) value to track, use this to terminate alignment
     early, by returning None in [update]. *)
  type stop
  type stop_parameter   (* Parameter to stopping logic. *)

  val stop_init : stop
  val stop_to_string : stop -> string

  val from_segment : position:int -> t -> a

  val update_state : ?early_stop:stop_parameter -> stop -> current:a -> aligned:a -> (a * stop) option

end

let debug_ref = ref false

let search_pos_edge_lst_to_string aindex l =
  List.map l ~f:(fun (o, e) -> sprintf "%d, %s" o
    (Alleles.Set.to_human_readable ~compress:true aindex e))
  |> String.concat ~sep:";"
  |> sprintf "[%s]"

module Align (Ag : Alignment_config) = struct

  module Ta = MakeThreadAlignment(Ag)

  type node_result =
    | Finished of Ag.t
    | GoOn of Ag.t * int

  let align_against s =
    let open Ta in
    let align = alignment ~thread:s in
    fun ?start ~search_pos ~node_seq ~node_offset () ->
      match align ?start ~thread_offset:search_pos ~segment:node_seq ~offset:node_offset with
      | Both m
      | First (m, _)   -> Finished m      (* end of the search string *)
      | Second (m, so) -> GoOn (m, so)

  let compute ?early_stop gt search_seq pos =
    let open Ref_graph in
    let open Nodes in
    let open Index in

    let nmas = align_against search_seq in
    let search_str_length = Ag.thread_length search_seq in
    (* We need a globally valid position, maybe just pass this integer as
      argument instead? *)
    let pos = pos.alignment + pos.offset in

    (* State and updating. *)
    let mis_map = Alleles.Map.make gt.aindex Ag.init in

    let ag_update = Ag.update_state ?early_stop in
    let assign ?node edge_set ~search_pos stop sa =
      let aligned = Ag.from_segment ~position:search_pos sa in
      if !debug_ref then
        eprintf "Assigning %s at %d -> %s to %s because of:%s\n%!"
          (Ag.to_string sa)
          search_pos
          (Ag.aligned_to_string aligned)
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
      Nmq.add_to_queue queue node nsplst
    (* Match a sequence, and add successor nodes to the queue for processing.
       Return an optional current queue and stop-state, where None mean stop. *)
    and match_and_add_succ state_opt ((node, splst) as ns) =
      Option.bind state_opt ~f:(fun (stop, queue) ->
        match node with
        | S _               ->
            invalid_argf "How did a Start get here %s!" (vertex_name node)
        | B _               ->
            Some (stop, add_successors queue ns)
        | E _               ->
            List.fold_left splst ~init:(Some stop) ~f:(fun stop_opt (search_pos, edge) ->
              Option.bind stop_opt ~f:(fun stop ->
                let sa = Ag.mismatches (search_str_length - search_pos) in
                assign ~node edge ~search_pos stop sa))
            |> Option.map ~f:(fun s -> s, queue)

        | N (_p, node_seq)  ->
            List.fold_left splst ~init:(Some (stop, [])) ~f:(fun stop_opt (search_pos, edge) ->
              Option.bind stop_opt ~f:(fun (stop, acc) ->
                match nmas ~search_pos ~node_seq ~node_offset:0 () with
                | Finished local_mismatches            ->
                    assign ~node edge ~search_pos stop local_mismatches
                    |> Option.map ~f:(fun s -> s, acc)
                | GoOn (local_mismatches, search_pos)  ->
                    assign ~node edge ~search_pos stop local_mismatches
                    |> Option.map ~f:(fun s -> s, (search_pos, edge) :: acc)))
            |> Option.map ~f:(fun (s, acc) ->
                match acc with
                | [] -> (s, queue)
                | _  -> (s, add_successors queue (node, acc))))
    in
    let rec assign_loop (stop, q) =
      if Nmq.is_empty q then
        `Finished mis_map
      else
        let nq, elst = Nmq.at_min_position q in
        match List.fold_left elst ~init:(Some (stop, nq)) ~f:match_and_add_succ with
        | Some sq -> assign_loop sq
        | None    -> `Stopped mis_map
    in
    Ref_graph.adjacents_at gt ~pos >>= begin fun {edge_node_set; seen_alleles} ->
      (* TODO. For now assume that everything that isn't seen has a full mismatch,
        this isn't strictly true since the Start of that allele could be within
        the range of the search str.
        - One approach would be to add the other starts, to the adjacents results.  *)
      let not_seen = Alleles.Set.complement gt.aindex seen_alleles in
      let mismatch_whole_read = Ag.mismatches search_str_length in
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
            EdgeNodeSet.fold edge_node_set ~init:(Some (stop, Nmq.empty))
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
                          assign_start (Ag.mismatches search_str_length)
                          |> Option.map ~f:(fun s -> (s, queue))
                      | B (p, _)         ->
                          let dist = p - pos in
                          if dist <= 0 then
                            Some (stop, add_successors queue (node, [0, edge]))
                          else if dist < search_str_length then begin
                            assign_start (Ag.mismatches dist)
                            |> Option.map ~f:(fun s -> (s, add_successors queue (node, [dist, edge])))
                          end else begin
                            assign_start (Ag.mismatches search_str_length)
                            |> Option.map ~f:(fun s -> (s, queue (* Nothing left to match. *)))
                          end
                      | N (p, node_seq)  ->
                          let nmas_and_assign ~node_offset ~start_mismatches =
                            match nmas ~search_pos:start_mismatches ~node_seq ~node_offset
                                    ~start:(Ag.mismatches start_mismatches) () with
                            | Finished mismatches            ->
                                assign_start mismatches
                                |> Option.map ~f:(fun s -> (s, queue))
                            | GoOn (mismatches, search_pos)  ->
                                assign_start mismatches
                                |> Option.map ~f:(fun s -> (s, add_successors queue (node, [search_pos, edge])))
                          in
                          let dist = p - pos in
                          if dist <= 0 then
                            nmas_and_assign ~node_offset:(-dist) ~start_mismatches:0
                          else if dist < search_str_length then
                            nmas_and_assign ~node_offset:0 ~start_mismatches:dist
                          else begin
                            assign_start (Ag.mismatches search_str_length)
                            |> Option.map ~f:(fun s -> (s, queue))
                          end))
          in
          match startq_opt with
          | None        -> Ok (`Stopped mis_map)
          | Some startq -> Ok (assign_loop startq)
    end

  let mismatches = Ag.mismatches

end (* Align *)

module MismatchesCounts = struct

  type t =
    { matches     : int
    ; mismatches  : int
    }

  let zero         = { matches = 0; mismatches = 0 }
  let mismatches n = { matches = 0; mismatches = n }
  let matches    n = { matches = n; mismatches = 0 }

  let matched ?(n=1) m = { m with matches = m.matches + n }
  let mismatched ?(n=1) m = { m with mismatches = m.mismatches + n }
  let to_string m =
    sprintf "{m : %d, ms: %d}" m.matches m.mismatches

  type thread = string
  let thread_length = String.length
  let thread_update t p c m =
    let cs = String.get_exn t p in
    if cs = c then matched m else mismatched m

end


module Mismatches_config = struct

  include MismatchesCounts

  type a = int
  let init = 0
  let aligned_to_string = sprintf "%d"

  type stop = float                   (* Current sum. *)
  type stop_parameter = int * float   (* Number of edges, max mean. *)
  let stop_init = 0.
  let stop_to_string = sprintf "sum %f"

  let from_segment ~position sa = sa.mismatches

  let update_state ?early_stop stop ~current ~aligned =
    let newal = current + aligned in
    match early_stop with
    | None ->
        Some (newal, stop)
    | Some (number_elements, max_mean) ->
        let nstop = stop +. (float aligned) in
        if nstop /. (float number_elements) > max_mean then
          None
        else
          Some (newal, nstop)

end  (* Mismatches_config *)

module Mismatches = Align (Mismatches_config)

let num_mismatches_against_seq = Mismatches.align_against
let compute_mismatches = Mismatches.compute

module PositionMismatches_config = struct

  include MismatchesCounts

  type a = (int * int) list
  let init = []
  let aligned_to_string l =
    List.map l ~f:(fun (p,m) -> sprintf "(%d,%d)" p m)
    |> String.concat ~sep:"; "

  type stop = float                   (* Current sum. *)
  type stop_parameter = int * float   (* Number of edges, max mean. *)
  let stop_init = 0.
  let stop_to_string = sprintf "sum %f"

  let from_segment ~position sa = [position, sa.mismatches]

  let update_state ?early_stop stop ~current ~aligned =
    let newal = current @ aligned in
    match early_stop with
    | None -> Some (newal, stop)
    | Some (number_elements, max_mean) ->
        let nstop = stop +. (float (List.length newal)) in
        if nstop /. (float number_elements) > max_mean then
          None
        else
          Some (newal, nstop)

end (* PositionMismatches_config *)

module PositionMismatches = Align (PositionMismatches_config)

let align_sequences_lst = PositionMismatches.align_against
let compute_mismatches_lst = PositionMismatches.compute

module PhredLikelihood_config = struct

  type t =
    { mismatches : float  (* number of mismatches, for stop logic. *)
    ; sum_llhd   : float  (* log likelihood *)
    }

  let mismatches n =
    let default_error = log 0.01 in
    let mismatches = float n in
    { mismatches
    ; sum_llhd = mismatches *. default_error}

  let to_string { mismatches; sum_llhd } =
    sprintf "{mismatches: %f; sum_llhd: %f}" mismatches sum_llhd

  type thread = string * float array
  let thread_length (s, _) = String.length s
  let thread_update (s, a) p c m =
    let cs = String.get_exn s p in
    let er = Array.get a p in
    if cs = c then
      { m with sum_llhd = m.sum_llhd +. log1p (-.er) }
    else
      { mismatches = m.mismatches +. 1.
      ; sum_llhd = m.sum_llhd +. log er
      }

  type a = t
  let init = mismatches 0
  let aligned_to_string = to_string

  type stop = float                   (* Current sum. *)
  type stop_parameter = int * float   (* Number of edges, max mean. *)
  let stop_init = 0.
  let stop_to_string = sprintf "sum %f"

  let from_segment ~position v = v
  let update_state ?early_stop stop ~current ~aligned =
    let newal =
      { mismatches = current.mismatches +. aligned.mismatches
      ; sum_llhd   = current.sum_llhd +. aligned.sum_llhd
      }
    in
    match early_stop with
    | None -> Some (newal, stop)
    | Some (number_elements, max_mean) ->
        let nstop = stop +. aligned.mismatches in
        if nstop /. (float number_elements) > max_mean then
          None
        else
          Some (newal, nstop)

end (* PhredLikelihood_config *)

module PhredLikelihood = Align (PhredLikelihood_config)

let align_sequences_plhd = PhredLikelihood.align_against
let compute_plhd = PhredLikelihood.compute

(* This method is a bad strawman... Would probably be much faster to
   go back to the original file and apply the Mas_parser changes to the
   reference. Only for comparison purposes... It probably makes sense to
   parameterize Ref_graph.sequence into a more general 'fold'?  *)
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
          Ok (nmas ~node_seq:graph_seq ~node_offset:0 ())
    else
      let mn = Mismatches.mismatches n in
      Ok (Mismatches.Finished mn))

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

