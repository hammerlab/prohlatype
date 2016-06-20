
open Util
open Graph
open Ref_graph
module Tg = Topological.Make_stable (G)

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

type position = alignment_position * sequence * int (* offset into sequence *)

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

type t = position list Kmer_table.t

let kmer_list ~k g =
  let f p acc = p :: acc in
  create_kmer_table ~k g f []

let create ~k g =
  (kmer_list ~k g).full

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


