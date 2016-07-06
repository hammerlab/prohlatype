(** Simple construction of alignment graphs. *)
open Nonstd

type allele = string
  [@@deriving eq, ord]
type uuid = int
  [@@deriving eq, ord]
module Node = struct
  type t =
    | Start of allele list * uuid
    | End of uuid
    | Boundary of uuid
    | Gap of int * uuid
    | Sequence of string * uuid
    [@@deriving eq, ord]

  let unique_name v =
    let q s = sprintf "\"%s\"" s in
    let gap s = String.make s '.' in
    (match v with
    | Start (alleles, uuid) -> sprintf "Start %d" uuid
    | End uuid -> sprintf "End %d" uuid
    | Boundary uuid -> sprintf "Boundary %d" uuid
    | Sequence (seq, uuid) -> sprintf "Seq %s %d" seq uuid
    | Gap (s, uuid) -> sprintf "Gap %s %d" (gap s) uuid
    ) |> q

  let name v =
    let gap s = String.make s '.' in
    (match v with
    | Start (alleles, uuid) -> sprintf "Start"
    | End uuid -> sprintf "End"
    | Boundary uuid -> sprintf "Boundary"
    | Sequence (seq, uuid) -> sprintf "%s" seq
    | Gap (s, uuid) -> sprintf "%s" (gap s)
    )

  let color v =
    match v with
    | Start (alleles, uuid) -> 0x00FF00EEl
    | End uuid -> 0x0000FF99l
    | Boundary uuid
    | Sequence (_, uuid) | Gap (_, uuid)
      -> 0x00000000l

  let shape v =
    match v with
    | Start (_, uuid) | End uuid -> `Oval
    | Boundary uuid -> `Diamond
    | Sequence (_, uuid) | Gap (_, uuid)
      -> `Box

  let style v =
    match v with
    | Start (_, uuid) | End uuid -> `Filled
    | Boundary uuid -> `Bold
    | Sequence (_, uuid) -> `Solid
    | Gap (_, uuid) -> `Solid

  let get_uuid v =
    match v with
    | Start (alleles, uuid) -> uuid
    | End uuid -> uuid
    | Boundary uuid -> uuid
    | Sequence (seq, uuid) -> uuid
    | Gap (s, uuid) -> uuid

  let hash = Hashtbl.hash
end

module ASet = Set.Make(struct type t = allele let compare = String.compare end)
module AlleleSet = struct
  include ASet
  let to_list t = ASet.fold (fun e acc -> e::acc) t []
end

module Edge = struct
  type t = AlleleSet.t
  let compare = AlleleSet.compare
  let equal = AlleleSet.equal
  let default = AlleleSet.empty

  let to_allele_list_string (_, allele_set, _) =
    AlleleSet.fold (fun e acc -> if acc = "" then e else acc ^ "," ^ e) allele_set ""
end

module G = Graph.Imperative.Digraph.ConcreteLabeled(Node)(Edge)

let uuid = ref 0
let create_node t =
  let v = G.V.create t in
  uuid := !uuid + 1;
  v

(* This type is used to hold a reference to any node which doesn't yet have all
   of the alleles coming into it, going out of it to another node. We track which
   alleles still haven't been accounted for (point to another node from here) so
   that we can 'garbage collect' the node when this list is empty. *)
type dangling_end = {
  incoming_only_alleles: AlleleSet.t;
  node: G.V.t;
}

let new_end node alleles =
  { incoming_only_alleles = AlleleSet.of_list alleles; node  }

(* This is the reference graph we'll be creating. ends is a list of unconnection
    nodes that we could possibly be extending the graph from. A "start" node is
    implicit in this list at all times for any yet-unseed allele. *)
type graph = { graph : G.t; ends : dangling_end list ; }


(* Utilities *)

let output_dot file_name {graph; _} =
  let oc = open_out file_name in
  let module Dot = Graph.Graphviz.Dot (
    struct
      include G
      let graph_attributes g = []
      let default_vertex_attributes v = []
      let vertex_attributes v = [`Shape (Node.shape v); `Label (Node.name v);
                                 `FillcolorWithTransparency (Node.color v);
                                 `Style (Node.style v)]
      let get_subgraph _v = None
      let vertex_name n = Node.unique_name n
      let default_edge_attributes _t = [`Color 4711]
      let edge_attributes e = [`Label (Edge.to_allele_list_string e)]
    end)
  in
  Dot.output_graph oc graph;
  close_out oc

(** Return AlleleSet.t of alleles coming into a given node. *)
let incoming_alleles ~graph node  =
    G.fold_pred_e (fun (p, allele_set, n) acc ->
      AlleleSet.union allele_set acc)
      graph node AlleleSet.empty

(** Return true if the given node has the specified incoming allele. *)
let is_node_of_allele ~allele ~graph node =
  G.fold_pred_e (fun (p, allele_set, n) acc ->
      if acc
      then acc
      else AlleleSet.mem allele allele_set)
    graph node false

let to_graph {graph; _} = graph

let update_ends_for_alleles alleles (ends:dangling_end list) =
  let module S = AlleleSet in
  let alleles = S.of_list alleles in
  List.fold ends
    ~f:(fun (froms, ends) ({node; incoming_only_alleles} as e) ->
        let common_alleles = S.inter alleles incoming_only_alleles in
        match S.is_empty common_alleles with
        | true -> (froms, e::ends)
        | false ->
          let new_ends =
            let remaining_alleles = S.diff incoming_only_alleles alleles in
            let new_end = {node; incoming_only_alleles = remaining_alleles} in
            if S.is_empty remaining_alleles
            then ends
            else new_end::ends
          in
          let new_from = new_end node (AlleleSet.to_list common_alleles) in
          (new_from::froms, new_ends)
      ) ~init:([], [])

let all_alleles_from_ends ends =
  List.fold ends ~init:[] ~f:(fun acc {incoming_only_alleles; _} ->
      (AlleleSet.to_list incoming_only_alleles) @ acc)
  |> List.dedup

let insert_node_or_start new_node alleles { graph; ends; } =
  assert ((List.length alleles) > 0);
  let open Node in
  let () = G.add_vertex graph new_node in
  let new_end = new_end new_node alleles in
  let from_ends, new_ends = update_ends_for_alleles alleles ends in
  if (List.length from_ends) = 0 then
    (* We're starting a new sequence, as we're not continuing from a previous
       sequence. *)
    let start = create_node (Start (alleles, !uuid)) in
    G.add_vertex graph start;
    G.add_edge_e graph (start, AlleleSet.of_list alleles, new_node);
    { graph; ends = new_end::new_ends; }
  else begin
    (* As we add edges from nodes which contain our alleles of interest, we make
       sure that we're adding edges for every allele we're expecting to: if
       neglected_alleles has alleles in it, we need to account for them (and do so
       in the match statement below). *)
    let neglected_alleles = List.fold from_ends ~init:(AlleleSet.of_list alleles)
        ~f:(fun neglected {node; incoming_only_alleles} ->
            let alleles_from_node =
              AlleleSet.(inter incoming_only_alleles (of_list alleles))
            in
            let () = G.add_edge_e graph (node, alleles_from_node, new_node) in
            (* Remove alleles from neglected that are now accounted for: *)
            AlleleSet.diff neglected incoming_only_alleles) in
    match AlleleSet.is_empty neglected_alleles with
    | true -> { graph; ends = new_end::new_ends; }
    | false ->
      (* Since we have alleles that didn't exist in previous nodes, We need to
         make a start node for these alleles and then connect them from that
         start to the new node *)
      let start = create_node
          (Start (AlleleSet.to_list neglected_alleles, !uuid)) in
      G.add_vertex graph start;
      G.add_edge_e graph (start, neglected_alleles, new_node);
      { graph; ends = new_end::new_ends; }
  end


let insert_nexus_node new_node { graph; ends } =
  let open Node in
  let all_alleles = all_alleles_from_ends ends in
  let new_end = new_end new_node all_alleles in
  G.add_vertex graph new_node;
  List.iter ends ~f:(fun {node; _} ->
      let alleles = incoming_alleles ~graph node in
      G.add_edge_e graph (node, alleles, new_node)
    );
  { graph; ends = [new_end] }


(******************************************************************************)
(*                               API                                          *)
(******************************************************************************)

let create () =
  { graph = G.create ();  ends = []; }

let seq sequence alleles g =
  let new_node = create_node Node.(Sequence (sequence, !uuid)) in
  insert_node_or_start new_node alleles g

let gap size alleles g =
  let new_node = create_node Node.(Gap (size, !uuid)) in
  insert_node_or_start new_node alleles g

let boundary g =
  let new_node = create_node Node.(Boundary !uuid) in
  insert_nexus_node new_node g

let finish g =
  let new_node = create_node Node.(End !uuid) in
  insert_nexus_node new_node g

let output ?(open_=true) file_name ({graph; _ } as g)=
  output_dot (file_name ^ ".dot") g;
  let r =
    Sys.command (sprintf "dot -Tpdf %s.dot -o %s.pdf" file_name file_name) in
  if r = 0 && open_
  then Sys.command (sprintf "open %s.pdf" file_name)
  else r
