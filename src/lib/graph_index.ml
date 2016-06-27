(** Graph Indexing

  This methodology uses a Kmer-table of positions where the kmer-ends to
  provide an index.

  Positions point to the last element of a kmer (ie. the "A" in 4-mer "TGCA"),
  for two reasons:
    1. Pointing to the beginning implies that the aligning step will have to
       rescan the graph from the beginning of the kmer, even though we know the
       kmer precisely.
    2. Pointing _past_ the last element opens the door to ambiguity, that position
       might be found in more than one sequence node.

  Two avoid unnecessary creating unnecessary sub-strings and string
  concatenation we perform all of the checks in the 'Kmer_to_int' encoded
  spaces. Specifically, the [encode ~ext] code allows us to convert a k-mer into
  an l-mer where (l > k). *)

open Util
open Graph
open Ref_graph
module Tg = Topological.Make_stable (G)

type ksublength =
  [ `Whole
  | `Part of int
  ]

type 'a kmer_substring =
  { index  : int
  ; length : 'a
  }

let whole ~index = { index; length = `Whole }
let part ~index l = { index; length = `Part l }

(* Fold over whole k-mers in s, but return an assoc of length remaining
   and suffixes that are not matched. *)
let fold_over_kmers_in_string ~k ~f ~init s =
  let l = String.length s in
  let rec loop index acc =
    if index = l then acc else
      let sub_length = l - index in
      if sub_length >= k then
        loop (index + 1) (f acc (whole ~index))
      else
        loop (index + 1) (f acc (part ~index sub_length))
  in
  loop 0 init

type ('a, 'b) kmer_fold_state =
  { full    : 'a
  ; partial : (int * 'b) list
  }

let fold_over_kmers_in_graph ~k ~f ~init ~extend ~close g =
  let open Nodes in
  let rec fill_partial_matches node state =
    match state.partial with
    | [] -> state
    | _l -> fold_partials_on_successor node state
  and fold_partials_on_successor node state =
    G.fold_succ (fun n s -> partial_state_on_successors s n) g node state
  and partial_state_on_successors state = function
    | S _             -> invalid_argf "Start should not be a successor!"
    | E _             -> state
    | B _ as n        -> fold_partials_on_successor n state (* skip node *)
    | (N (p, s) as n) ->
        let l = String.length s in
        let f state (krem, bstate) =
          if krem <= l then
            let nfull = close state.full (p, s) (part ~index:0 krem) bstate in
            { state with full = nfull }
          else
            let np = (krem - l, extend (p, s) (part ~index:0 l) (Some bstate)) :: state.partial in
            { state with partial = np }
        in
        List.fold_left state.partial ~f ~init:{ state with partial = [] }
        |> fill_partial_matches n
  in
  let proc node state =
    match node with
    | S _ | E _ | B _ -> state
    | N (p, s)        ->
        let f state ss =
          match ss with
          | { index; length = `Whole  } ->
              let nfull = f state.full (p, s) (whole ~index) in
              { state with full = nfull }
          | { index; length = `Part l } ->
              { state with partial = (k - l, extend (p, s) (part ~index l) None) :: state.partial }
        in
        fold_over_kmers_in_string s ~k ~f ~init:{ state with partial = [] }
        |> fill_partial_matches node
  in
  Tg.fold proc g { full = init; partial = [] }

let kmer_counts ~k g =
  let init = Kmer_table.make k 0 in
  let f tbl (_al, sequence) { index; length = `Whole } =
    let i = Kmer_to_int.encode sequence ~pos:index ~len:k in
    Kmer_table.update ((+) 1) tbl i;
    tbl
  in
  let extend (_al, sequence) { index; length = `Part len } = function
    | None     -> Kmer_to_int.encode sequence ~pos:index ~len
    | Some ext -> Kmer_to_int.encode sequence ~pos:index ~len ~ext
  in
  let close tbl (_al, sequence) { index; length = `Part len} ext =
    let i = Kmer_to_int.encode sequence ~pos:index ~len ~ext in
    Kmer_table.update ((+) 1) tbl i;
    tbl
  in
  let fs = fold_over_kmers_in_graph ~k g ~f ~close ~extend ~init in
  fs.full

type position =
  { alignment   : alignment_position
  ; sequence    : sequence
  ; km1_offset  : int       (* into sequence *)
  }

(* A Graph index *)
type t = position list Kmer_table.t

let create ~k g =
  let init = Kmer_table.make k [] in
  let f tbl (alignment, sequence) { index; length = `Whole } =
    let i = Kmer_to_int.encode sequence ~pos:index ~len:k in
    let p = { alignment; sequence; km1_offset = index + k - 1 } in
    Kmer_table.update (fun lst -> p :: lst) tbl i;
    tbl
  in
  let extend (_al, sequence) { index; length = `Part len } = function
    | None     -> Kmer_to_int.encode sequence ~pos:index ~len
    | Some ext -> Kmer_to_int.encode sequence ~pos:index ~len ~ext
  in
  let close tbl (alignment, sequence) { index; length = `Part len} ext =
    let i = Kmer_to_int.encode sequence ~pos:index ~len ~ext in
    let p = { alignment; sequence; km1_offset = index + len - 1 } in
    Kmer_table.update (fun lst -> p :: lst) tbl i;
    tbl
  in
  let fs = fold_over_kmers_in_graph ~k g ~f ~close ~extend ~init in
  fs.full

let starting_with index s =
  let k = Kmer_table.k index in
  match String.sub s ~index:0 ~length:k with
  | None    -> error "Not long enough %s for index size %d" s k
  | Some ss -> Ok (Kmer_table.lookup index ss)

let k = Kmer_table.k
