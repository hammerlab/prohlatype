
open Graph
open Nonstd
module String = Sosa.Native_string

let invalid_argf fmt = ksprintf invalid_arg fmt

(* TODO:
  - Hashcons the sequences
  - Turn fold_succ_e from O(n) into something better
*)

let short_seq s =
  let n = String.length s in
  if n > 10 then
    sprintf "%s...%s"
      (String.slice_exn ~finish:4 s) (String.slice_exn ~start:(n-3) s)
  else
    s

module Nodes = struct

  type t =
    | S of int * string
    | E of int
    | B of int * int      (* Boundary of position and count *)
    | N of int * string   (* Sequences *)

  let vertex_name ?(short=true) = function
    | S (n, s)  -> sprintf "\"S%d-%s\"" n s
    | E n       -> sprintf "\"E%d\"" n
    | B (_, n)  -> sprintf "\"B%d\"" n
    | N (n, s)  -> sprintf "\"%d%s\"" n (if short then short_seq s else s)

  let compare = Pervasives.compare
  let equal = (=)
  let hash = Hashtbl.hash
end

module Edges = struct
  type t = String.t           (* Which allele *)
  let compare = String.compare
  let hash = Hashtbl.hash
  let equal = (=)
  let default = ""
end

module G = Imperative.Digraph.ConcreteLabeled(Nodes)(Edges)

(** Traversing **)
exception Found of Nodes.t

let next_node_along allele g ~from =
  try
    G.fold_succ_e (fun (_, l, vs) n ->
        if l = allele then raise (Found vs) else n)
      g from None
  with Found v ->
    Some v

let previous_node_along allele g ~from =
  try
    G.fold_pred_e (fun (pv, l, _) n ->
        if l = allele then raise (Found pv) else n)
      g from None
  with Found v ->
    Some v

let fold_along_allele ~start allele g ~f ~init =
  let next = next_node_along allele g in
  let rec loop from (acc, stop) =
    if stop then
      acc
    else
      match next ~from with
      | None    -> acc
      | Some vs -> loop vs (f acc vs)
  in
  loop start (f init start)

(* If this was a GADT we could encode the p s v inside the
   N that we're matching! *)
let vsplit_at g allele ~pos p s v =
  let open Nodes in
  let index = pos - p in
  let fs, sn = String.split_at s ~index in
  let pr = G.pred_e g v in
  let su = G.succ_e g v in
  G.remove_vertex g v;
  let v1 = N (p, fs) in
  G.add_vertex g v1;
  List.iter ~f:(fun (p, l, _) -> G.add_edge_e g (G.E.create p l v1)) pr;
  let v2 = N (pos, sn) in
  G.add_vertex g v2;
  G.add_edge_e g (G.E.create v1 allele v2);
  List.iter ~f:(fun (_, l, s) -> G.add_edge_e g (G.E.create v2 l s)) su;
  (v1, v2)

(** Output **)

(* TODO:
  - When constructing the dot files, it would be nice if the alleles (edges),
    were in some kind of consistent order. *)

let output_dot ?short fname g =
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
      let edge_attributes e = [`Label (G.E.label e)]

    end)
  in
  Dot.output_graph oc g;
  close_out oc

let output ?(pdf=true) ?(open_=true) ~short fname g =
  output_dot ~short (fname ^ ".dot") g;
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
