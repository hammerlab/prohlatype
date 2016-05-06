
(* TODO:
  Hashcons the strings?
*)

(* OCamlgraph *)
open Graph
open Printf
open Mas_parser

module String = Sosa.Native_string

let invalid_argf fmt = ksprintf invalid_arg fmt

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

let seq_elem_to_vertex = function
  | Start (_, a)    -> G.V.create (Sequences.S a)
  | End _           -> G.V.create Sequences.E
  | Boundary (n, p) -> G.V.create (Sequences.B n)
  | Nuc (p, s)      -> G.V.create (Sequences.N (p, s))
  | Gap _           -> invalid_argf "Not a valid vertex"

let add_reference allele g pv se =
  match se with
  | Start _ | End _ | Boundary _ | Nuc _ ->
    let v = seq_elem_to_vertex se in
    G.add_vertex g v;
    G.add_edge_e g (G.E.create pv allele v);
    v
  | Gap _ -> pv (* Ignore gaps, since no nodes to create. *)

let next_along allele g v =
  let module M = struct exception Found of Sequences.t end in
  try
    G.fold_succ_e (fun (_, l, vs) n ->
        if l <> allele then n else raise (M.Found vs))
      g v None
  with M.Found v ->
    Some v

let prev_along allele g v =
  let module M = struct exception Found of Sequences.t end in
  try
    G.fold_pred_e (fun (pv, l, _) n ->
        if l <> allele then n else raise (M.Found pv))
      g v None
  with M.Found v ->
    Some v

let first_boundary = Sequences.B 0

let fold_along_allele ?(start=first_boundary) allele g ~f ~init =
  let next = next_along allele g in
  let rec loop v acc =
    match next v with
    | None    -> acc
    | Some vs -> loop vs (f acc vs)
  in
  loop start (f init start)

let vsplit_at g allele ~pos p s v =
  let open Sequences in
  let index = pos - p in
  let fs, sn = String.split_at s ~index in
  let pr = G.pred_e g v in
  let su = G.succ_e g v in
  G.remove_vertex g v;
  let v1 = N (p, fs) in
  G.add_vertex g v1;
  List.iter (fun (p, l, _) -> G.add_edge_e g (G.E.create p l v1)) pr;
  let v2 = N (pos, sn) in
  G.add_vertex g v2;
  G.add_edge_e g (G.E.create v1 allele v2);
  List.iter (fun (_, l, s) -> G.add_edge_e g (G.E.create v2 l s)) su;
  (v1, v2)

let find_or_create_break allele g pos start =
  let open Sequences in
  let module M = struct exception Done of G.vertex end in
  try
    fold_along_allele ~start allele g ~init:None ~f:(fun pv v ->
      match v with
      | S _ | E   -> pv
      | B _       -> Some v
      | N (p, s)  ->
        if p = pos then
          let v = match pv with | Some b -> b | None -> v in
          raise (M.Done v)
        else if pos >= p + (String.length s) then
          None
        else  (* found it! *)
          let _v1, v2 = vsplit_at g allele ~pos p s v in
          raise (M.Done v2))
    |> fun _ -> invalid_argf "reached end of find_or_create_break %d" pos
  with M.Done v -> v

let vsplit_at_two g allele ~pos1 ~pos2 p s v =
  let open Sequences in
  let fs, sn = String.split_at s ~index:(pos1 - p) in
  let sn, th = String.split_at sn ~index:(pos2 - pos1) in
(*   let () = Printf.printf "%s|\t%s|\t%s|\t%s\n" s fs sn th in *)
  let pr = G.pred_e g v in
  let su = G.succ_e g v in
  G.remove_vertex g v;
  let v1 = N (p, fs) in
  G.add_vertex g v1;
  List.iter (fun (p, l, _) -> G.add_edge_e g (G.E.create p l v1)) pr;
  let v2 = N (pos1, sn) in
  G.add_vertex g v2;
  G.add_edge_e g (G.E.create v1 allele v2);
  let v3 = N (pos2, th) in
  G.add_vertex g v3;
  G.add_edge_e g (G.E.create v2 allele v3);
  List.iter (fun (_, l, s) -> G.add_edge_e g (G.E.create v3 l s)) su;
  (v1, v2, v3)

let find_or_create_two_breaks allele g (pos1, pos2) =
  if pos1 >= pos2 then invalid_argf "pos must be increasing" else
  let open Sequences in
  let module M = struct exception Done of G.vertex * G.vertex end in
  try
    let start = S allele in
    fold_along_allele ~start allele g ~init:None ~f:(fun pv v ->
      match v with
      | S _ | E   -> pv
      | B _       -> Some v
      | N (p, s)  ->
        let len = String.length s in
        if pos1 >= p + len then
(*          let () = Printf.printf "case1\n" in *)
          None
        else if p = pos1 then begin
          if pos2 < p + len then begin
(*             let () = Printf.printf "case2: %d\n" p in  *)
            let v1, v2 = vsplit_at g allele ~pos:pos2 p s v in
            (*let vp = match pv with Some b -> b | None -> v1 in *)
            raise (M.Done (v1, v2))
          end else
(*             let () = Printf.printf "case3: %d\n" p in *)
            (* Won't change v *)
            (*let vp = match pv with Some b -> b | None -> v in *)
            raise (M.Done (v, find_or_create_break allele g pos2 v))
        end else begin
          if pos2 < p + len then begin
(*              let () = Printf.printf "case4: %d\n" p in  *)
            let _v1, v2, v3 = vsplit_at_two g allele ~pos1 ~pos2 p s v in
            raise (M.Done (v2, v3))
          end else
(*              let () = Printf.printf "case5: %d\n" p in  *)
            let _v1, v2 = vsplit_at g allele ~pos:pos1 p s v in
            raise (M.Done (v2, find_or_create_break allele g pos2 v2))
        end)
    |> (fun _ ->
        invalid_argf "reached end of find_or_create_two_breaks allele %s %d %d" allele pos1 pos2)
  with M.Done (v1, v2) -> v1, v2

let insert_gap ref_allele g (allele, pos, gap_length) =
  let v1, v2 = find_or_create_two_breaks ref_allele g (pos, pos + gap_length) in
  match prev_along ref_allele g v1 with
  | None    -> invalid_argf "No previous node for pos %d" pos
  | Some pv -> G.add_edge_e g (G.E.create pv allele v2); v2

let insert_nuc ref_allele g (allele, pos, seq) =
  let open Sequences in
  let v1, v2 = find_or_create_two_breaks ref_allele g (pos, pos + String.length seq) in
  match prev_along ref_allele g v1 with
  | None    -> invalid_argf "No previous node for pos %d" pos
  | Some pv ->
    let vm = N (pos, seq) in
    (* Will adding be a no-op if this seq is already present. *)
    G.add_vertex g vm;
    G.add_edge_e g (G.E.create pv allele vm);
    G.add_edge_e g (G.E.create vm allele v2);
    v2

let add_alt ref_allele g allele pv se =
  match se with
  | Start _ | End _   ->
    let v = seq_elem_to_vertex se in
    G.add_vertex g v;
    G.add_edge_e g (G.E.create pv allele v);
    v
  | Boundary (_, _) ->
    let b = seq_elem_to_vertex se in
    assert (G.mem_vertex g b);
    G.add_edge_e g (G.E.create pv allele b);
    b
  | Nuc (p, s) ->
    insert_nuc ref_allele g (allele, p, s)
  | Gap (p, gl) ->
    insert_gap ref_allele g (allele, p, gl)

let add_elems g add_seq_elem = function
  | h :: t ->
      let v0 = seq_elem_to_vertex h in
      G.add_vertex g v0;
      List.fold_left add_seq_elem v0 t |> ignore
  | _ -> invalid_argf "Empty sequence"

let add_non_ref ref_allele g (allele, allele_seq) =
  printf "adding %s \n" allele;
  add_elems g (add_alt ref_allele g allele) allele_seq

let construct ?(num_too_add=max_int) (ref_allele, ref_seq, alg_assoc) =
  let ref_length = List.length ref_seq in
  let num_alleles = List.length alg_assoc in
  let g = G.create ~size:(ref_length * num_alleles) () in
  add_elems g (add_reference ref_allele g) ref_seq;
  List.iteri (fun i p -> if i < num_too_add then
                 add_non_ref ref_allele g p) alg_assoc;
  g

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
