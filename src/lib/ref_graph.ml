
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

type start = int * string [@@deriving eq, ord]
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
  type t = BitSet.t
  let hash = Hashtbl.hash
  let compare = BitSet.compare
  let equal = BitSet.equals
  let default = BitSet.empty ()
end


module G = Imperative.Digraph.ConcreteLabeled(Nodes)(Edges)

exception Found of Nodes.t

let next_node_along aset allele g ~from =
  try
    G.fold_succ_e (fun (_, bt, vs) n ->
        if Allele_set.is_set aset bt allele then raise (Found vs) else n)
      g from None
  with Found v ->
    Some v

let fold_along_allele aset ~start allele g ~f ~init =
  let next = next_node_along aset allele g in
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
        let i = Pattern.encode_sub s ~pos:index ~len:k in     (* TODO: use consistent args *)
        let nfull = Kmer_table.update_index f state.full posit i in
        { state with full = nfull }
    | `Part p ->
        let is = Pattern.encode_sub s ~pos:index ~len:p in
        { state with partial = (k - p, is, posit) :: state.partial }
  in
  let fpartial state kseqlst p s =
    let l = String.length s in
    List.fold_left kseqlst ~init:(state, [])
      ~f:(fun (state, acc) (krem, curp, posit) ->
            if krem <= l then
              let pn = Pattern.encode_extend s ~pos:0 ~len:krem curp in
              let nf = Kmer_table.update_index f state.full posit pn in
              { state with full = nf }, acc
            else
              let pn = Pattern.encode_extend s ~pos:0 ~len:l curp in
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

let starting_with index s =
  let k = Kmer_table.k index in
  match String.sub s ~index:0 ~length:k with
  | None -> error "Not long enough: %s vs %d" s k
  | Some ss -> Ok (Kmer_table.lookup index ss)

let mismatches ~s1 ~o1 ~s2 ~o2 =
  let l1 = String.length s1 in
  let l2 = String.length s2 in
  let rec loop i m =
    let i1 = i + o1 in
    let i2 = i + o2 in
    if i1 >= l1 then
      if i2 >= l2 then
        `Both m
      else
        `First m
    else if i2 >= l2 then
      `Second m
    else
      loop (i + 1) (if s1.[i1] = s2.[i2] then m else m + 1)
  in
  loop 0 0

(*
let align g index s =
  starting_with index s >>= (fun lst ->
      List.map lst ~f:(fun (p, s, offset) ->
          let node = N (p, s) in


      (*| [] -> error "Not found"  *)
   *)

(** Output **)

(* TODO:
  - When constructing the dot files, it would be nice if the alleles (edges),
    were in some kind of consistent order. *)

let output_dot ?short fname (aset, g) =
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
      let edge_attributes e = [`Label (Allele_set.to_human_readable aset (G.E.label e))]

    end)
  in
  Dot.output_graph oc g;
  close_out oc

let output ?(pdf=true) ?(open_=true) ~short fname (aset, g) =
  output_dot ~short (fname ^ ".dot") (aset, g);
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
