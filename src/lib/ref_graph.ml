
open Util
open Graph

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

type start = int * string
type end_ = int

(* start end pairs *)
type sep = { start : start ; end_ : end_ }

module Nodes = struct

  type t =
    | S of start
    | E of end_
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

(** Traversal **)
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

module Tg = Topological.Make_stable (G)

(* Fold over whole k-mers in s, but return an assoc of length remaining
   and suffixes that are not matched. *)
let only_over_whole_kmers ?(canonical=false) k s f init =
  let l = String.length s in
  let rec loop i (a, lst) =
    if i = l then (a, lst) else
      let sub_length = l - i in
      if sub_length >= k then
        loop (i + 1) (f a (String.sub_exn s i k), lst)
      else
        loop (i + 1) (a, (k - sub_length, String.sub_exn s i sub_length)::lst)
  in
  loop 0 (init, [])
  |> (fun (s, l) -> s, if canonical then List.rev l else l)

let fold_kmers ?(canonical=false) ~k g f init =
  let open Nodes in
  let rec over_succ state node = function
    | [] -> state
    | k_seq_lst ->
        (* Since we're folding over each edge, we'll get the successor node 'v'
          for every allele edge. *)
        G.fold_succ (fun node ((state, set) as st) ->
            if List.mem node ~set then st else
              (* This logic makes it DFS *)
              let nstate = proc_succ state k_seq_lst node in
              nstate, node :: set) g node (state, [])
        |> fst
  and proc_succ state k_seq_lst node =
    match node with
      | S _       -> invalid_argf "Start should not be a successor!"
      | E _       -> state
      | B _       -> over_succ state node k_seq_lst
      | N (_, s)  ->
          let l = String.length s in
          List.fold_left k_seq_lst ~init:(state,[])
            ~f:(fun (state, acc) (k, seq) ->
                  if k <= l then
                    let ns = String.concat [seq; String.sub_exn s 0 k] in
                    let nstate = f state ns in
                    nstate, acc
                  else
                    state, (k - l, String.concat [seq; s]) :: acc)
          |> fun (state, ksl) ->
                over_succ state node (if canonical then List.rev ksl else ksl)
  in
  let proc node state =
    match node with
    | S _ | E _ | B _ -> state
    | N (_, s)        ->
      let nstate, lst = only_over_whole_kmers ~canonical k s f state in
      over_succ nstate node lst
  in
  Tg.fold proc g init

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

let create_kmer_table ~k g f i =
  let t = Kmer_table.make k i in
  fold_kmers ~k g (Kmer_table.update f) t

let kmer_counts ~k g =
  create_kmer_table ~k g ((+) 1) 0

(*
let kmer_nodes ~k g =
  create_kmer_table ~k g (fun acc node -> node :: acc) [] 
*)
(*let lookup s (k, kmt) =
  let prefix = String.take ~index:k s in
  *)

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

let save fname g =
  let oc = open_out fname in
  Marshal.to_channel oc g [];
  close_out oc

let load fname =
  let ic = open_in fname in
  let g : G.t = (Marshal.from_channel ic) in
  close_in ic;
  g
