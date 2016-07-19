
open Ref_graph

let adjacents_fold ~f ~init g node =
  G.fold_pred (G.fold_succ_e (fun (_, e, v) a ->
      if v <> node then f a e v else a) g) 
    g node init

let up g f = G.fold_pred f g

module Ons = Set.Make (struct
  type t = int * Nodes.t
  let compare (i1, n1) (i2, n2) = 
    let c1 = compare i1 i2 in
    if c1 = 0 then Nodes.compare n1 n2 else c1
end)

let up_n n f g v init =
  let rec loop s =
    let ((min_height, node) as mn) = Ons.min_elt s in
    if min_height < n then
      let ss = Ons.remove mn s in
      let ns =
        G.fold_pred (fun v s -> Ons.add (min_height + 1, v) s)
          g node ss
      in
      loop ns
    else
      Ons.fold (fun (_, n) a -> f n a) s init
   in
   loop (Ons.singleton (0, v))

module NodesSet = MoreLabels.Set.Make (struct
  type t = Nodes.t
  let compare = Nodes.compare
end)

let adjacents_n max_height ~f ~init g node =
  let add_if_new ~cur n ss =
    if NodesSet.mem n cur then ss else
      NodesSet.add n ss
  in
  let add_and_fold_if_new ~cur (_, e, n) ((ss, acc) as st) =
    if NodesSet.mem n cur then st else
      NodesSet.add n ss, (f e n acc)
  in
  let rec up i acc cur tl =
    if i > max_height then
      acc, cur :: tl
    else
      let parents = NodesSet.fold cur ~f:(G.fold_pred NodesSet.add g) ~init:NodesSet.empty in
      let children, grandkids, nacc = down acc parents cur tl in
      up (i + 1) nacc parents (children :: grandkids)
  and down acc parents cur = function
    | [] ->
        let with_siblings, nacc =
          NodesSet.fold parents ~f:(G.fold_succ_e (add_and_fold_if_new ~cur) g)
            ~init:(cur, acc)
        in
        with_siblings, [], nacc
    | (h :: t) ->
        let children =
          NodesSet.fold parents ~f:(G.fold_succ (add_if_new ~cur) g) ~init:cur
        in
        let grandkids, grand_grand_list, nacc = down acc children h t in
        children, (grandkids :: grand_grand_list), nacc
  in
  let start = NodesSet.singleton node in
  up 0 init start []

let adjacents_until (type a) ?max_height ~f ~init g node =
  let module M = struct exception F of a end in
  let add_if_new ~cur n ss =
    if NodesSet.mem n cur then ss else
      NodesSet.add n ss
  in
  let add_and_fold_if_new ~cur (_, e, n) ((ss, acc) as st) =
    if NodesSet.mem n cur then st else
      match f e n acc with
      | `Stop r     -> raise (M.F r)
      | `Continue r -> NodesSet.add n ss, r
  in
  let rec up i acc cur tl =
    match max_height with
    | Some mh when i > mh -> acc, cur :: tl
    | None | Some _ ->
        let parents = NodesSet.fold cur ~f:(G.fold_pred NodesSet.add g) ~init:NodesSet.empty in
        try 
          let children, grandkids, nacc = down acc parents cur tl in
          up (i + 1) nacc parents (children :: grandkids)
        with M.F r ->
          r, (parents :: cur :: tl)
  and down acc parents cur = function
    | [] ->
        let with_siblings, nacc =
          NodesSet.fold parents ~f:(G.fold_succ_e (add_and_fold_if_new ~cur) g)
            ~init:(cur, acc)
        in
        with_siblings, [], nacc
    | (h :: t) ->
        let children =
          NodesSet.fold parents ~f:(G.fold_succ (add_if_new ~cur) g) ~init:cur
        in
        let grandkids, grand_grand_list, nacc = down acc children h t in
        children, (grandkids :: grand_grand_list), nacc
  in
  let start = NodesSet.singleton node in
  up 0 init start []


(* An edge is a triple of previous vertex, edge label (set) next vertex.
  Since we'll fold over the next vertex, paramterize down to have the same
  signature! 

let down g f (_, _, v) = G.fold_succ_e f g v

let down_n n g f a v =
  let rec loop i k =
    if i = n then k
    else loop (i + 1) (fun f -> k (k f))
  in
  let down_enough = loop 1 (down g) in
  down_enough f a v
*)

