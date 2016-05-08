
open Graph        (* OCamlgraph *)
open Nonstd
open Mas_parser
open Ref_graph
module String = Sosa.Native_string

let invalid_argf fmt = ksprintf invalid_arg fmt

(** Construction from Mas_parser *)
let seq_elem_to_vertex se =
  let open Mas_parser in
  match se with
  | Start (_, a)    -> G.V.create (Sequences.S a)
  | End _           -> G.V.create Sequences.E
  | Boundary (n, _) -> G.V.create (Sequences.B n)
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

let insert_gap reference g (allele, pos, gap_length) =
  let from, v2 = find_or_create_two_breaks reference g (pos, pos + gap_length) in
  match previous_node_along reference g ~from with
  | None    -> invalid_argf "No previous node for pos %d" pos
  | Some pv -> G.add_edge_e g (G.E.create pv allele v2); v2

let insert_nuc reference g (allele, pos, seq) =
  let open Sequences in
  let from, v2 = find_or_create_two_breaks reference g (pos, pos + String.length seq) in
  match previous_node_along reference g ~from with
  | None    -> invalid_argf "No previous node for pos %d" pos
  | Some pv ->
    (* Will adding be a no-op if this seq is already present. *)
    let vm = N (pos, seq) in
    G.add_edge_e g (G.E.create pv allele vm);
    G.add_vertex g vm;
    G.add_edge_e g (G.E.create vm allele v2);
    v2

module SeqSet = Set.Make (struct
  type t = string Mas_parser.seq_elems
  let compare = compare
end)

let add_alt ~reference g ref_gaps_set allele pv se =
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
    insert_nuc reference g (allele, p, s)
  | Gap (p, gl) ->
    (* Gap is in reference, so the appropriate break already exists, just need
      to find the appropriate reference edge and replicated it. *)
    if SeqSet.mem se ref_gaps_set then
      (* Not create! *)
      let v = find_or_create_break reference g p pv in
      G.add_edge_e g (G.E.create pv allele v);
      v
    else
      insert_gap reference g (allele, p, gl)

let add_elems g add_seq_elem = function
  | h :: t ->
      let v0 = seq_elem_to_vertex h in
      G.add_vertex g v0;
      List.fold_left ~f:add_seq_elem ~init:v0 t |> ignore
  | _ -> invalid_argf "Empty sequence"

let add_non_ref reference g ref_gaps_set (allele, allele_seq) =
  printf "adding %s \n" allele;
  add_elems g (add_alt ~reference g ref_gaps_set allele) allele_seq

let construct ?(num_alt_to_add=max_int) {reference; ref_elems; alt_elems} =
  let ref_length = List.length ref_elems in
  let num_alleles = List.length alt_elems in
  let g = G.create ~size:(ref_length * num_alleles) () in
  add_elems g (add_reference reference g) ref_elems;
  let ref_gaps_set =
    List.fold_left ref_elems ~init:SeqSet.empty
      ~f:(fun s v->
        match v with | Gap _ -> SeqSet.add v s | _ -> s)
  in
  List.iteri ~f:(fun i p ->
    if i < num_alt_to_add then
      add_non_ref reference g ref_gaps_set p)
  alt_elems;
  g


