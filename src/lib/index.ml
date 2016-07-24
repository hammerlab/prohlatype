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

let fold_over_kmers_in_graph ~k ~f ~init ~extend ~close g =
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
        (* I'm not thrilled with f's lexical hiding (see below) but haven't
           thought of a good name for the external 'f' similar to
           'extend' and 'close' *)
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
              let nstate = f state (p, s) (whole ~index) in
              nstate, pacc
          | { index; length = `Part l } ->
              state, (k - l, extend (p, s) (part ~index l) None) :: pacc
        in
        let state, plst = fold_over_kmers_in_string s ~k ~f ~init:(state, []) in
        fill_partial_matches node state plst
  in
  Tg.fold proc g.g init

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
  fold_over_kmers_in_graph ~k g ~f ~close ~extend ~init

type position =
  { alignment : alignment_position
  ; sequence  : sequence
  ; offset    : int       (* into sequence *)
  }

let specific_position g allele alp =
  let open Nodes in
  Ref_graph.find_node_at g allele alp >>= function
    | S _      -> error "%d for %s is at Start " alp allele
    | E _      -> error "%d for %s is at End " alp allele
    | B _      -> error "%d for %s is at Boundary " alp allele
    | N (p, s) -> Ok { alignment = p ; sequence = s ; offset = alp - p }

(* A Graph index *)
type t = position list Kmer_table.t

let create_at_penultimate ~k g =
  let init = Kmer_table.make k [] in
  let f tbl (alignment, sequence) { index; length = `Whole } =
    let i = Kmer_to_int.encode sequence ~pos:index ~len:k in
    let p = { alignment; sequence; offset = index + k - 1 } in
    (*let () = printf "Adding %d \t %s \t %d\n" alignment (Ref_graph.index_string sequence p.offset) p.offset in *)
    Kmer_table.update (fun lst -> p :: lst) tbl i;
    tbl
  in
  let extend (al, sequence) { index; length = `Part len } ext_opt =
    (*let () = printf "At %d, %s extend at %d len %d" al (Ref_graph.index_string sequence index) index len in *)
    match ext_opt with
    | None     ->
        let nj = Kmer_to_int.encode sequence ~pos:index ~len in
        (*let () = printf " cur: None \t nj %d \n" nj in *)
        nj
    | Some ext ->
        let nj = Kmer_to_int.encode sequence ~pos:index ~len ~ext in
        (*let () = printf " cur: Some %d \t nj %d\n" ext nj in *)
        nj
  in
  let close tbl (alignment, sequence) { index; length = `Part len} ext =
    let i = Kmer_to_int.encode sequence ~pos:index ~len ~ext in
    let p = { alignment; sequence; offset = index + len - 1 } in
    (*let () = printf "Closing %d \t %s \t %d\t ext: %d \n"
      alignment (Ref_graph.index_string sequence p.offset) p.offset ext
    in *)
    Kmer_table.update (fun lst -> p :: lst) tbl i;
    tbl
  in
  fold_over_kmers_in_graph ~k g ~f ~close ~extend ~init

let create ~k g =
  let init = Kmer_table.make k [] in
  let f tbl (alignment, sequence) { index; length = `Whole } =
    let i = Kmer_to_int.encode sequence ~pos:index ~len:k in
    let p = { alignment; sequence; offset = index } in
    (*let () = printf "Adding %d \t %s \t %d\n" alignment (Ref_graph.index_string sequence p.offset) p.offset in *)
    Kmer_table.update (fun lst -> p :: lst) tbl i;
    tbl
  in
  let extend (alignment, sequence) { index; length = `Part len } ext_opt =
    (*let () = printf "At %d, %s extend at %d len %d" alignment
          (Ref_graph.index_string sequence index) index len in *)
    match ext_opt with
    | None     ->
        let nj = Kmer_to_int.encode sequence ~pos:index ~len in
        let p  = { alignment; sequence; offset = index } in
        (*let () = printf " cur: None \t nj %d \n" nj in *)
        (nj, p)
    | Some (ext, p) ->
        let nj = Kmer_to_int.encode sequence ~pos:index ~len ~ext in
        (*let () = printf " cur: Some %d \t nj %d\n" ext nj in *)
        nj, p
  in
  let close tbl (_al, sequence) { index; length = `Part len} (ext, p) =
    let i = Kmer_to_int.encode sequence ~pos:index ~len ~ext in
    (*let () = printf "Closing %d \t %s \t %d\t ext: %d \n"
      alignment (Ref_graph.index_string sequence p.offset) p.offset ext
    in *)
    Kmer_table.update (fun lst -> p :: lst) tbl i;
    tbl
  in
  fold_over_kmers_in_graph ~k g ~f ~close ~extend ~init

let starting_with index s =
  let k = Kmer_table.k index in
  match String.sub s ~index:0 ~length:k with
  | None    -> error "Not long enough %s for index size %d" s k
  | Some ss -> Ok (Kmer_table.lookup_kmer index ss)

let lookup ?(n=0) index s =
  let k = Kmer_table.k index in
  match String.sub s ~index:0 ~length:k with
  | None    -> error "Not long enough %s for index size %d" s k
  | Some ss ->
    if n = 0 then
      Ok (Kmer_table.lookup_kmer index ss)
    else if n = 1 then
      let ns = Kmer_to_int.neighbors ~k (Kmer_to_int.encode ss) in
      let init = Kmer_table.lookup_kmer index ss in
      Array.fold_left ~init ~f:(fun acc n ->
        List.append acc (Kmer_table.lookup index n)) ns
      |> fun lst -> Ok lst
    else
      invalid_argf "not implemented n: %d" n




let k = Kmer_table.k
