
open Util
open Graph

(* TODO:
  - Hashcons the sequences
  - Turn fold_succ_e from O(n) into something better
*)

type allele_name = string [@@deriving eq, ord]
type alignment_position = int [@@deriving eq, ord]
type start = alignment_position * allele_name [@@deriving eq, ord]
type end_ = alignment_position [@@deriving eq, ord]
type sequence = string [@@deriving eq, ord]

(* start end pairs *)
type sep = { start : start ; end_ : end_ }

module Nodes = struct

  type t =
    | S of start
    | E of end_
    | B of alignment_position * int      (* Boundary of position and count *)
    | N of alignment_position * sequence   (* Sequences *)
    [@@deriving eq, ord]

  let vertex_name ?(short=true) = function
    | S (n, s)  -> sprintf "\"S%d-%s\"" n s
    | E n       -> sprintf "\"E%d\"" n
    | B (_, n)  -> sprintf "\"B%d\"" n
    | N (n, s)  -> sprintf "\"%d%s\"" n (if short then short_seq s else s)

  let hash = Hashtbl.hash
end

module Edges = struct
  type t = Alleles.Set.t
  let hash = Hashtbl.hash
  let compare = Alleles.Set.compare
  let equal = Alleles.Set.equals
  let default = Alleles.Set.empty ()
end

module G = Imperative.Digraph.ConcreteLabeled(Nodes)(Edges)

exception Found of Nodes.t

let next_node_along aindex allele g ~from =
  try
    G.fold_succ_e (fun (_, bt, vs) n ->
        if Alleles.Set.is_set aindex bt allele then raise (Found vs) else n)
      g from None
  with Found v ->
    Some v

let fold_along_allele aindex ~start allele g ~f ~init =
  let next = next_node_along aindex allele g in
  let rec loop from (acc, stop) =
    if stop then
      acc
    else
      match next ~from with
      | None    -> acc
      | Some vs -> loop vs (f acc vs)
  in
  loop start (f init start)

let between g start stop =
  let ng = G.create () in
  G.add_vertex ng start;
  let rec add_from node =
    G.iter_succ_e (fun ((_, _, n) as e) ->
      G.add_vertex ng n;
      G.add_edge_e ng e;
      if n <> stop then add_from n)
      g node
  in
  add_from start;
  ng

(** Output **)

(* TODO:
  - When constructing the dot files, it would be nice if the alleles (edges),
    were in some kind of consistent order. *)

let output_dot ?short fname (aindex, g) =
  let oc = open_out fname in
  let module Dot = Graphviz.Dot (
    struct
      include G
      let graph_attributes _g = []
      let default_vertex_attributes _g = []
      let vertex_name = Nodes.vertex_name ?short
      let vertex_attributes _v = [`Shape `Box]
      let get_subgraph _v = None

      let default_edge_attributes _t = [`Color 4711]
      let edge_attributes e = [`Label (Alleles.Set.to_human_readable aindex (G.E.label e))]

    end)
  in
  Dot.output_graph oc g;
  close_out oc

let output ?(pdf=true) ?(open_=true) ~short fname (aindex, g) =
  output_dot ~short (fname ^ ".dot") (aindex, g);
  let r =
    if pdf then
      Sys.command (sprintf "dot -Tpdf %s.dot -o %s.pdf" fname fname)
    else
      -1
  in
  if r = 0 && open_ then
    Sys.command (sprintf "open %s.pdf" fname)
  else
    r

let save fname g =
  let oc = open_out fname in
  Marshal.to_channel oc g [];
  close_out oc

let load fname =
  let ic = open_in fname in
  let g : G.t = (Marshal.from_channel ic) in
  close_in ic;
  g
