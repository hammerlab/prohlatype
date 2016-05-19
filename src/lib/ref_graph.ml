
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

module Sequences = struct

  type t =
    | S of string             (* Start of allele *)
    | E                       (* End *)
    | B of int                (* Boundary *)
    | N of int * String.t     (* Sequences *)

  let vertex_name ?(short=true) = function
    | S a       -> sprintf "\"S%s\"" a
    | E         -> sprintf "\"E\""
    | B n       -> sprintf "\"B%d\"" n
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

module G = Imperative.Digraph.ConcreteLabeled(Sequences)(Edges)

(** Traversing **)
exception Found of Sequences.t

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

let fold_along_allele ?start allele g ~f ~init =
  let start = Option.value start ~default:(Sequences.S allele) in
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
  let open Sequences in
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

let find_or_create_split_at allele g pos start =
  let open Sequences in
  (*let module M = struct exception Done of G.vertex end in *)
  fold_along_allele ~start allele g ~init:(start, start)
    ~f:(fun (pv1, pv2) v ->
          match v with
          | S _       -> (pv2, v), false
          | E         -> (pv2, v), true
          | B _       -> (pv2, v), false
          | N (p, s)  ->
            if pos >= p + (String.length s) then  (* too long *)
              (pv2, v), false
            else if p = pos then                  (* found it exactly! *)
              (pv2, v), true
            else                                  (* in this nucleotide it! *)
              vsplit_at g allele ~pos p s v, true)

let vsplit_at_two g allele ~pos1 ~pos2 p s v =
  let open Sequences in
  let fs, sn = String.split_at s ~index:(pos1 - p) in
  let sn, th = String.split_at sn ~index:(pos2 - pos1) in
  let pr = G.pred_e g v in
  let su = G.succ_e g v in
  G.remove_vertex g v;
  let v1 = N (p, fs) in
  G.add_vertex g v1;
  List.iter ~f:(fun (p, l, _) -> G.add_edge_e g (G.E.create p l v1)) pr;
  let v2 = N (pos1, sn) in
  G.add_vertex g v2;
  G.add_edge_e g (G.E.create v1 allele v2);
  let v3 = N (pos2, th) in
  G.add_vertex g v3;
  G.add_edge_e g (G.E.create v2 allele v3);
  List.iter ~f:(fun (_, l, s) -> G.add_edge_e g (G.E.create v3 l s)) su;
  (v1, v2, v3)

let find_or_create_two_breaks allele g (pos1, pos2) start =
  let open Sequences in
  if pos1 >= pos2 then invalid_argf "pos must be increasing" else
  fold_along_allele allele g ~init:start
    ~f:(fun (pv1, pv2) v ->
          match v with
          | S _       -> (pv1, pv2), false
          | E         -> (pv1, pv2), true
          | B _       -> (pv2, v), false
          | N (p, s)  ->
            let len = String.length s in
            if pos1 >= p + len then       (* haven't found first position. *)
              (pv2, v), false
            else if p = pos1 then begin   (* first position matches !*)
              if pos2 < p + len then begin
                let v1, v2 = vsplit_at g allele ~pos:pos2 p s v in
                (v1, v2), true
              end else
                (* Won't change v *)
                let _v1, v2 = find_or_create_split_at allele g pos2 v in
                (v, v2), true
            end else begin
              if pos2 < p + len then begin
                let _v1, v2, v3 = vsplit_at_two g allele ~pos1 ~pos2 p s v in
                (v2, v3), true
              end else
                let _v1, v2 = vsplit_at g allele ~pos:pos1 p s v in
                let _v2, v3 = find_or_create_split_at allele g pos2 v2 in
                (v2, v3), true
            end)

(** Output **)

(* TODO:
  - When constructing the dot files, it would be nice if the alleles, edges,
    were in some kind of consistent order. *)

let output_dot ?short fname g =
  let oc = open_out fname in
  let module Dot = Graphviz.Dot (
    struct
      include G
      let graph_attributes _g = []
      let default_vertex_attributes _g = []
      let vertex_name = Sequences.vertex_name ?short
      let vertex_attributes _v = [`Shape `Box]
      let get_subgraph _v = None

      let default_edge_attributes _t = [`Color 4711]
      let edge_attributes e = [`Label (G.E.label e)]

    end)
  in
  Dot.output_graph oc g;
  close_out oc

let output ?(dot=true) ?(_open=true) ~short fname g =
  output_dot ~short (fname ^ ".dot") g;
  let r =
    if dot then
      Sys.command (sprintf "dot -Tpdf %s.dot -o %s.pdf" fname fname)
    else
      -1
  in
  if r = 0 && _open then
    Sys.command (sprintf "open %s.pdf" fname)
  else
    r
