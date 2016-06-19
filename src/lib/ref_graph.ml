
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

type allele_name = string [@@deriving eq, ord]
type start = int * allele_name [@@deriving eq, ord]
type end_ = int [@@deriving eq, ord]

(* start end pairs *)
type sep = { start : start ; end_ : end_ }

module Nodes = struct

  type t =
    | S of start
    | E of end_
    | B of int * int      (* Boundary of position and count *)
    | N of int * string   (* Sequences *)
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

module Tg = Topological.Make_stable (G)

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


(* I want to avoid all of the calls to 'String.sub_ below so this call will
   just pass the (1) whole string, (2) an offset, and (3) k; putting the
   responsibility on 'f' to deal with not full length strings. *)

(* Fold over whole k-mers in s, but return an assoc of length remaining
   and suffixes that are not matched. *)
let over_kmers_of_string k p s ~f ~init =
  let l = String.length s in
  let rec loop index acc =
    if index = l then acc else
      let sub_length = l - index in
      if sub_length >= k then
        loop (index + 1) (f acc p s ~index ~length:`Whole)
      else
        loop (index + 1) (f acc p s ~index ~length:(`Part sub_length))
  in
  loop 0 init

type ('a, 'b) kmer_fold_state =
  { full    : 'a
  ; partial : 'b list
  }

let update_full f state = { state with full = f state }
let update_partial f state = { state with partial = f state }

let fold_kmers ~k g ~f ~fpartial ~init =
  let open Nodes in
  let rec fill_partial_matches node ({ partial ; _} as state) =
    match partial with
    | [] -> state
    | ps ->
        G.fold_succ (fun node state -> partial_state_on_successors state ps node)
          g node state
  and partial_state_on_successors state partial = function
    | S _             -> invalid_argf "Start should not be a successor!"
    | E _             -> state
    | B _ as n        -> fill_partial_matches n state (* skip node *)
    | (N (p, s) as n) -> fpartial state partial p s
                         |> fill_partial_matches n
  in
  let proc node state =
    match node with
    | S _ | E _ | B _ -> state
    | N (p, s)        -> over_kmers_of_string k p s ~f ~init:state
                         |> fill_partial_matches node
  in
  Tg.fold proc g init

let create_kmer_table ~k g f init =
  let init = { full = Kmer_table.make k init ; partial = [] } in
  let new_f state p s ~index ~length =
    let posit = p, s, index in
    match length with
    | `Whole  ->
        let i = Kmer_to_int.encode s ~pos:index ~len:k in     (* TODO: use consistent args *)
        Kmer_table.update (f posit) state.full i;
        state
    | `Part p ->
        let is = Kmer_to_int.encode s ~pos:index ~len:p in
        { state with partial = (k - p, is, posit) :: state.partial }
  in
  let fpartial state kseqlst p s =
    let l = String.length s in
    List.fold_left kseqlst ~init:(state, [])
      ~f:(fun (state, acc) (krem, curp, posit) ->
            if krem <= l then
              let pn = Kmer_to_int.encode s ~pos:0 ~len:krem ~ext:curp in
              Kmer_table.update (f posit) state.full pn;
              state, acc
            else
              let pn = Kmer_to_int.encode s ~pos:0 ~len:l ~ext:curp in
              let na = (krem - l, pn, posit) :: acc in
              state, na)
    |> (fun (s, plst) -> { s with partial = plst })
  in
  fold_kmers ~k g ~f:new_f ~fpartial ~init

let kmer_counts ~k g =
  let f _ c = c + 1 in
  create_kmer_table ~k g f 0

let kmer_list ~k g =
  let f p acc = p :: acc in
  create_kmer_table ~k g f []

(*
val cross_boundary : ('a * string * int) list t -> (string * ('a * string * int) list) list
*)
let cross_boundary kt =
  let k = Kmer_table.k kt in
  Kmer_table.fold kt ~init:(0, []) ~f:(fun (i, acc) lst ->
      let p = Kmer_to_int.decode ~k i in
      let ni = i + 1 in
      let cross = List.filter lst ~f:(fun (_, s, o) -> o > String.length s - k) in
      match cross with
      | []   -> (ni, acc)
      | glst -> (ni, (p, glst) :: acc))
  |> snd

let starting_with index s =
  let k = Kmer_table.k index in
  match String.sub s ~index:0 ~length:k with
  | None    -> error "Not long enough %s for index size %d" s k
  | Some ss -> Ok (Kmer_table.lookup index ss)

let align_sequences ~s1 ~o1 ~s2 ~o2 =
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
      loop (i + 1) (if s1.[i1] = s2.[i2] then m else m + 1)
  in
  loop 0 0

(* TODO:
   There are 2 pieces of information about a graph that are relevant
   for alignment:

   1. Which edges a path through the graph traverses, since edges give us alleles
   2. If a path is entirely on a vertex (ex. "AC" in "GGACTT") then we need to
      know if the vertex ("GGACTT") is only found in specific alleles
   *)

type 'e pa =
  { mismatches : int
  ; edge       : 'e             (*BitSet.t        too big to carry around? *)
  ; alignment  : 'e alignment
  }
and 'e alignment = Leaf of int
                 | Partial of 'e pa list

let rec map_alignment f = function
  | Leaf n     -> Leaf n
  | Partial pl -> Partial (List.map pl ~f:(fun p ->
      { p with edge = f p.edge ; alignment = map_alignment f p.alignment }))

(* Add max number of mismatches as argument  *)
let align ?(mub=max_int) g index search_seq =
  let open Nodes in
  let align_against_seq ~search_pos ~node_seq ~node_offset =
    (*let () =
      Printf.printf "align_against_seq search_pos %d \t node_seq %s \t node_offset %d\n"
              search_pos node_seq node_offset
    in*)
    let reached_end_of =
      align_sequences ~s1:search_seq ~o1:search_pos ~s2:node_seq ~o2:node_offset
    in
    match reached_end_of with
    | `Both m
    | `First (m, _)   -> `Finished m      (* end of the search string *)
    | `Second (m, so) -> `GoOn (m, so)
  in
  let over_all_sequence_nodes node m0 ~search_pos =
    let rec descend node m0 ~search_pos =
      (*let () = Printf.printf "descend m0 %d search_pos %d\n" m0 search_pos in *)
      G.fold_succ_e (fun (_, edge, vs) acc ->
          match vs with
          | S _             -> invalid_argf "Another start?"
          | E _             -> acc
          | B _             -> (* This is a weird case since all prev nodes point at a boundary
                                  the acc in this case should be [] *)
                               assert (acc = []);
                               descend vs m0 ~search_pos
          | N (_, node_seq) ->
              match align_against_seq ~search_pos ~node_seq ~node_offset:0 with
              | `Finished mismatches    ->
                  let tm = m0 + mismatches in
                  if tm > mub then acc else
                    { mismatches ; edge; alignment = Leaf tm } :: acc
              | `GoOn (mismatches, search_pos)  ->
                  let tm = m0 + mismatches in
                  if tm > mub then acc else
                    let alignment = Partial (descend vs tm ~search_pos) in
                    { mismatches ; edge; alignment } :: acc)
        g node []
    in
    descend node m0 ~search_pos
  in
  (*let k = Kmer_table.k index in*)
  starting_with index search_seq >>= (fun lst ->
      List.filter_map lst ~f:(fun (p, node_seq, no) ->
          (* The current index table only tells you where a k-mer starts;
             not the full k-mer path. Consequently, there can be other paths
             when aligning:
             match align_against_seq ~search_pos:k ~node_seq ~node_offset:(no + k) with
             is wrong!
          *)
          match align_against_seq ~search_pos:0 ~node_seq ~node_offset:no with
          | `Finished mismatches  ->
              if mismatches > mub then None else Some (Leaf mismatches)
          | `GoOn (m, search_pos) ->
              if m > mub then None else
                Some (Partial (over_all_sequence_nodes (N (p, node_seq)) m ~search_pos)))
      |> fun als -> Ok als)

let name_edges_in_alignment aindex =
  map_alignment (Alleles.Set.to_human_readable aindex)

let rec best_of_paths = function
  | Leaf n      -> n
  | Partial pl  ->
      let bps = List.map ~f:(fun pa -> best_of_paths pa.alignment) pl in
      List.fold_left ~f:min ~init:max_int bps
      (* if pl = [] then there are no alignments that satisfy the mub *)

let to_weights lst =
  let flst = List.map ~f:float_of_int lst in
  let ilst = List.map ~f:(fun x -> 1. /. (1. +. x)) flst in
  let s = List.fold_left ~f:(+.) ~init:0. ilst in
  List.map ~f:(fun x -> x /. s) ilst

let init_alingment_map aindex =
  Alleles.Map.make aindex 0.

(* Weighing Alignments ... inference *)
let alignments_to_weights ?(base_weight=1.) amap al =
  let rec descend w0 = function
  | Leaf _     -> ()
  | Partial pl -> let weights = to_weights (List.map pl ~f:(fun pa -> pa.mismatches)) in
                  List.iter2 pl weights ~f:(fun pa w ->
                    let w_i = w *. w0 in
                    Alleles.Map.update_from pa.edge amap ((+.) w_i);
                    descend w_i pa.alignment)
  in
  descend base_weight al

let most_likely aindex amap =
  Alleles.Map.fold aindex ~init:[] ~f:(fun acc v allele ->
    if v > 0. then (v,allele) :: acc else acc) amap
  |> List.sort ~cmp:(fun ((v1 : float), _) (v2,_) -> compare v2 v1)

  (*
let most_likely { Allele_set.to_allele; _} foo =

  Array.fold_left to_allele ~init:(0,[]) ~f:(fun (i,acc) a ->
      let v = foo.(i) in
      if v > 0. then (i+1, (v,a) :: acc) else (i + 1, acc))
  |> snd
  |> List.sort ~cmp:(fun ((v1 : float), _) (v2,_) -> compare v2 v1)
  *)


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
