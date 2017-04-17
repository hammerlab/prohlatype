(** Graph Indexing

  This methodology uses a Kmer-table of kmer starts provide an index.

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

let fold_over_kmers_in_graph ~k ~init ~absorb ~extend ~close g =
  let open Nodes in
  let rec fill_partial_matches node state = function
    | [] -> state
    | ls -> fold_partials_on_successor node state ls
  and fold_partials_on_successor node state plst =
    G.fold_succ (fun n s -> partial_state_on_successors s plst n) g.g node state
  and partial_state_on_successors state plst = function
    | S _             -> invalid_argf "Start should not be a successor!"
    | E _             -> state
    | B _ as n        -> fold_partials_on_successor n state plst (* skip node *)
    | (N (p, s) as n) ->
        let l = String.length s in
        let f (state, pacc) (krem, bstate) =
          if krem <= l then
            let nstate = close state (p, s) (part ~index:0 krem) bstate in
            nstate, pacc
          else
            let npacc = (krem - l, extend (p, s) (part ~index:0 l) (Some bstate)) :: pacc in
            state, npacc
        in
        let state, nplst = List.fold_left plst ~f ~init:(state, []) in
        fill_partial_matches n state nplst
  in
  let proc node state =
    match node with
    | S _ | E _ | B _ -> state
    | N (p, s)        ->
        let f (state, pacc) ss =
          match ss with
          | { index; length = `Whole  } ->
              absorb state (p, s) (whole ~index), pacc
          | { index; length = `Part l } ->
              state, (k - l, extend (p, s) (part ~index l) None) :: pacc
        in
        let state, plst = fold_over_kmers_in_string s ~k ~f ~init:(state, []) in
        fill_partial_matches node state plst
  in
  Tg.fold proc g.g init

let fold_over_biological_kmers_in_graph ~k ~init ~absorb ~extend ~close g =
  let open Nodes in
  let rec fill_partial_matches node state = function
    | [] -> state
    | ls -> fold_partials_on_successor node state ls
  and fold_partials_on_successor node state plst =
    G.fold_succ_e (partial_state_on_successors plst) g.g node state
  and partial_state_on_successors plst (_, in_edge, node) state =
    match node with
    | S _      -> invalid_argf "Start should not be a successor!"
    | E _      -> state
    | B _      ->
        (* skip the node BUT filter on the alleles going to this node! *)
        let nplst =
          List.filter_map plst ~f:(fun (k, ke, bstate) ->
            let kei = Alleles.Set.inter in_edge ke in
            if Alleles.Set.cardinal kei > 0 then
              Some (k, kei, bstate)
            else None)
        in
        begin match nplst with
        | [] -> state
        | lst -> fold_partials_on_successor node state lst
        end
    | N (p, s) ->
        let l = String.length s in
        let f (state, pacc) (krem, k_edge, bstate) =
          let kei = Alleles.Set.inter in_edge k_edge in
          if Alleles.Set.cardinal kei > 0 then begin
            if krem <= l then
              close state (p, s) (part ~index:0 krem) bstate, pacc
            else
              state, (krem - l, kei, extend (p, s) (part ~index:0 l) (Some bstate)) :: pacc
          end else
            state, pacc
        in
        let state, nplst = List.fold_left plst ~f ~init:(state, []) in
        fill_partial_matches node state nplst
  in
  (*let all_incoming n =
    G.fold_pred_e (fun (_,e,_) a -> Alleles.Set.union e a) g.g n
      (Alleles.Set.init g.aindex)
  in *)
  let everything = Alleles.Set.(complement (init g.aindex)) in
  let proc node state =
    match node with
    | S _ | E _ | B _ -> state
    | N (p, s)        ->
        (*let in_edges = all_incoming node in *)
        let f (state, pacc) ss =
          match ss with
          | { index; length = `Whole  } ->
              absorb state (p, s) (whole ~index), pacc
          | { index; length = `Part l } ->
              state, (k - l, everything, extend (p, s) (part ~index l) None) :: pacc
        in
        let state, plst = fold_over_kmers_in_string s ~k ~f ~init:(state, []) in
        fill_partial_matches node state plst
  in
  Tg.fold proc g.g init

let kmer_counts ~biological ~k g =
  let init = Kmer_table.make k 0 in
  let absorb tbl (_al, sequence) { index; length = `Whole } =
    let ilst = Kmer_to_int.encode_N_tolerant sequence ~pos:index ~len:k in
    List.iter ilst ~f:(Kmer_table.update ((+) 1) tbl);
    tbl
  in
  let extend (_al, sequence) { index; length = `Part len } = function
    | None      -> Kmer_to_int.encode_N_tolerant sequence ~pos:index ~len
    | Some exts -> Kmer_to_int.encode_N_tolerant sequence ~pos:index ~len ~exts
  in
  let close tbl (_al, sequence) { index; length = `Part len} exts =
    let ilst = Kmer_to_int.encode_N_tolerant sequence ~pos:index ~len ~exts in
    List.iter ilst ~f:(Kmer_table.update succ tbl);
    tbl
  in
  if biological then
    fold_over_biological_kmers_in_graph ~k g ~absorb ~close ~extend ~init
  else
    fold_over_kmers_in_graph ~k g ~absorb ~close ~extend ~init

type position =
  { alignment : alignment_position
  ; sequence  : sequence
  ; offset    : int       (* into sequence *)
  } [@@deriving eq,ord,show]

let specific_position g allele pos =
  let open Nodes in
  Ref_graph.find_node_at g ~allele ~pos >>= function
    | S _      -> error "%d for %s is at Start " pos allele
    | E _      -> error "%d for %s is at End " pos allele
    | B _      -> error "%d for %s is at Boundary " pos allele
    | N (p, s) -> Ok { alignment = p ; sequence = s ; offset = pos - p }

(* A Graph index *)
type t = position list Kmer_table.t

(* The failed method were we store the penultimate position of the Kmer *)
let create_at_penultimate ~k g =
  let init = Kmer_table.make k [] in
  let absorb tbl (alignment, sequence) { index; length = `Whole } =
    let ilst = Kmer_to_int.encode_N_tolerant sequence ~pos:index ~len:k in
    let p = { alignment; sequence; offset = index + k - 1 } in
    List.iter ilst ~f:(Kmer_table.update (fun lst -> p :: lst) tbl);
    tbl
  in
  let extend (al, sequence) { index; length = `Part len } ext_opt =
    match ext_opt with
    | None      -> Kmer_to_int.encode_N_tolerant ~pos:index ~len sequence
    | Some exts -> Kmer_to_int.encode_N_tolerant sequence ~pos:index ~len ~exts
  in
  let close tbl (alignment, sequence) { index; length = `Part len} exts =
    let ilst = Kmer_to_int.encode_N_tolerant sequence ~pos:index ~len ~exts in
    let p = { alignment; sequence; offset = index + len - 1 } in
    List.iter ilst ~f:(Kmer_table.update (fun lst -> p :: lst) tbl);
    tbl
  in
  fold_over_biological_kmers_in_graph ~k g ~absorb ~close ~extend ~init

let create ~k g =
  let init = Kmer_table.make k [] in
  let absorb tbl (alignment, sequence) { index; length = `Whole } =
    let ilst = Kmer_to_int.encode_N_tolerant sequence ~pos:index ~len:k in
    let p = { alignment; sequence; offset = index } in
    List.iter ilst ~f:(Kmer_table.update (fun lst -> p :: lst) tbl);
    tbl
  in
  let extend (alignment, sequence) { index; length = `Part len } ext_opt =
    match ext_opt with
    | None     ->
        let nj = Kmer_to_int.encode_N_tolerant sequence ~pos:index ~len in
        let p  = { alignment; sequence; offset = index } in
        (nj, p)
    | Some (exts, p) ->
        let nj = Kmer_to_int.encode_N_tolerant sequence ~pos:index ~len ~exts in
        nj, p
  in
  let close tbl (_al, sequence) { index; length = `Part len} (exts, p) =
    let ilst = Kmer_to_int.encode_N_tolerant sequence ~pos:index ~len ~exts in
    List.iter ilst ~f:(Kmer_table.update (fun lst -> p :: lst) tbl);
    tbl
  in
  fold_over_biological_kmers_in_graph ~k g ~absorb ~close ~extend ~init

let lookup ?(distance=0) index s =
  Kmer_table.lookup_kmer_neighbors ~d:distance index s >>= fun pos_list_arr ->
    (* Since kmer neighbors are very similar the lookup positions will also
       be repetitive. It is helpful to take a step and trim this list. *)
    Ok (List.concat (Array.to_list pos_list_arr) |> List.dedup)

let k = Kmer_table.k
