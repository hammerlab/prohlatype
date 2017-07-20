
open Util
open Graph

(* TODO:
  - Hashcons the sequences
  - Turn fold_succ_e from O(n) into something better
*)

type alignment_position = int [@@deriving eq, ord, show]
type start =
  alignment_position * (Alleles.allele [@equal Alleles.equal] [@compare Alleles.compare])
    [@@deriving eq, ord]
type end_ = alignment_position [@@deriving eq, ord]
type sequence = string [@@deriving eq, ord, show]

(* start end pairs *)
type sep = { start : start ; end_ : end_ } [@@deriving eq, ord]

module Nodes = struct

  type t =
    | S of start
    | E of end_
    | B of alignment_position * int       (* Boundary of position and count *)
    | N of alignment_position * sequence  (* Sequences *)
    [@@deriving eq, ord]

  let vertex_name ?(short=true) = function
    | S (n, s)  -> sprintf "S%d-%s" n s
    | E n       -> sprintf "E%d" n
    | B (p, n)  -> sprintf "B%d-%d" n p
    | N (n, s)  -> sprintf "%d%s" n (if short then short_seq s else s)

  let position = function
    | S (p, _) | E p | B (p, _) | N (p, _)  -> p

  let inside pos = function
    | S (p, _) | E p | B (p, _) -> p = pos
    | N (p, s) -> p <= pos && pos < p + String.length s

  let inside_seq p s ~pos =
    p < pos && pos < p + String.length s

  let compare_by_position_first n1 n2 =
    let compare_int (i1 : int) (i2 : int) = Pervasives.compare i1 i2 in
    let r = compare_int (position n1) (position n2) in
    if r = 0 then compare n1 n2 else r

  let hash = Hashtbl.hash

  let is_seq_or_boundary = function
    | N _ | B _ -> true
    | S _ | E _ -> false

end (* Nodes *)

module Edges = struct

  type t = Alleles.set
  let hash = Hashtbl.hash
  let compare = Alleles.set_compare
  let equal = Alleles.set_equals
  let default = Alleles.set_empty

end (* Edges *)

module G = Imperative.Digraph.ConcreteLabeled(Nodes)(Edges)

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

type by_position =
  | NL of Nodes.t list  (* A list of nodes that start at a given position. *)
  | Redirect of int     (* A redirect (index into array - NOT position) of where
                           to find the previous Node list. In
                           [by_position_array] below, that location should
                           contain the nodes that span this position. *)

and by_position_array = by_position array

module NodeSet = struct
  include Set.Make (Nodes)
  let to_string s =
    fold ~f:(fun n a -> Nodes.vertex_name n :: a) ~init:[] s
    |> String.concat ~sep:"; "
    |> sprintf "{%s}"
end (* NodeSet *)

module EdgeNodeSet = struct

  include Set.Make (struct
    (* edges first since these edges point into the respective node. *)
    type t = Edges.t * Nodes.t

    let compare (e1, n1) (e2, n2) =
      let r = Nodes.compare n1 n2 in
      if r = 0 then Edges.compare e1 e2 else r
  end)

  (*let to_string ?max_length ?complement s =
    fold ~f:(fun (e, n) a ->
      (sprintf "(%s -> %s)"
        (Alleles.Set.to_human_readable ?max_length ?complement e)
        (Nodes.vertex_name n)) :: a) ~init:[] s
    |> String.concat ~sep:"; "
    |> sprintf "{%s}"

  let to_table ?max_length ?complement s =
    fold ~f:(fun p l -> p :: l) ~init:[] s
    |> List.sort ~cmp:(fun (_,n1) (_,n2) -> Nodes.compare n1 n2)
    |> List.map ~f:(fun (e, n) ->
        sprintf "%s <- %s"
          (Nodes.vertex_name n)
          (Alleles.Set.to_human_readable ?max_length ?complement e))
    |> String.concat ~sep:"\n"
    |> sprintf "%s" *)

end (* EdgeNodeSet *)

type adjacent_info =
  { edge_node_set : EdgeNodeSet.t
  ; seen_alleles  : Alleles.set
  }

type t =
  { align_date    : string                   (* When the alignment was created by IMGT. *)
  ; reference     : string
  ; g             : G.t                                            (* The actual graph. *)
  ; aindex        : Alleles.index               (* The allele index, for Sets and Maps. *)
  ; aset          : (module Alleles.Set)
  ; amap          : (module Alleles.Map)
  ; bounds        : sep list Alleles.map    (* Map of where the alleles start and stop. *)
  ; posarr        : by_position_array
  ; adjacents_arr : adjacent_info array
  ; offset        : int
  ; merge_map     : (string * string) list  (* Left empty if not a merged|imputed graph *)
  (* TODO: Type the merge_map logic. Strings are too loose and the impute -> merge logic
           isn't uniform. *)
  }

(** Output **)

(* TODO:
  - When constructing the dot files, it would be nice if the alleles (edges),
    were in some kind of consistent order. *)

(** [starts_by_position g] returns an associated list of alignment_position and
    an a set of alleles that start at that position. *)
let starts_by_position { bounds; amap; aset; _ } =
  let module AS = (val aset : Alleles.Set) in
  let module AM = (val amap : Alleles.Map) in
  AM.fold bounds ~init:[] ~f:(fun asc sep_lst allele ->
    List.fold_left sep_lst ~init:asc ~f:(fun asc sep ->
      let pos = fst sep.start in
      try
        let bts, rest = remove_and_assoc pos asc in
        (pos, AS.set bts allele) :: rest
      with Not_found ->
        (pos, AS.singleton allele) :: asc))

let add_allele_edge g aset pv nv allele =
  let module AS = (val aset : Alleles.Set) in
  let new_edge =
    try
      let eset = G.find_edge g pv nv |> G.E.label in
      G.remove_edge g pv nv;
      G.E.create pv (AS.set eset allele) nv
    with Not_found ->
      G.E.create pv (AS.singleton allele) nv
  in
  G.add_edge_e g new_edge

(** Compress the start nodes; join all the start nodes that have the same
    alignment position into one node. *)
let create_compressed_starts t =
  let module AS = (val t.aset : Alleles.Set) in
  let start_asc = starts_by_position t in
  let ng = G.copy t.g in
  let open Nodes in
  List.iter start_asc ~f:(fun (pos, allele_set) ->
    if AS.cardinal allele_set > 1 then begin
      let a_str = AS.to_string ~compress:true allele_set in
      let node = G.V.create (S (pos, a_str)) in
      G.add_vertex ng node;
      AS.iter allele_set ~f:(fun allele ->
        let rm = S (pos, allele) in
        G.iter_succ (fun sv -> add_allele_edge ng t.aset node sv allele) ng rm;
        G.remove_vertex ng rm)
    end);
  { t with g = ng }

let output_dot ?(human_edges=true) ?(compress_edges=true) ?(compress_start=true)
  ?(insert_newlines=true) ?short ?max_length fname t =
  let module AS = (val t.aset : Alleles.Set) in
  let { g; _} = if compress_start then create_compressed_starts t else t in
  let oc = open_out fname in
  let insert_newline = insert_chars ['\n'] in
  let module Dot = Graphviz.Dot (
    struct
      include G
      let graph_attributes _g = []
      let default_vertex_attributes _g = []
      let vertex_name v =
        let s = sprintf "\"%s\"" (Nodes.vertex_name ?short v) in
        if insert_newlines then insert_newline s else s

      let vertex_attributes _v = [`Shape `Box]
      let get_subgraph _v = None

      let default_edge_attributes _t = [`Color 4711]
      let edge_attributes e =
        let compress = compress_edges in
        let s =
          if human_edges then
            AS.to_human_readable ~compress ?max_length (G.E.label e)
          else
            AS.to_string ~compress (G.E.label e)
        in
        if insert_newlines then
          [`Label ( insert_newline s) ]
            else
          [`Label s ]

    end)
  in
  Dot.output_graph oc g;
  close_out oc

let output ?human_edges ?compress_edges ?compress_start ?insert_newlines
  ?max_length ?(pdf=true) ?(open_=true) ~short fname t =
  output_dot ?human_edges ?compress_edges ?compress_start ?max_length ~short (fname ^ ".dot") t;
  if pdf then begin
    let r = Sys.command (sprintf "dot -Tpdf %s.dot -o %s.pdf" fname fname) in
    if r = 0 && open_ then
      Sys.command (sprintf "open %s.pdf" fname)
    else
      r
  end else
    0

let save fname g =
  let oc = open_out fname in
  Marshal.to_channel oc g [];
  close_out oc

let load fname =
  let ic = open_in fname in
  let g : G.t = (Marshal.from_channel ic) in
  close_in ic;
  g

(** Construction *)

let relationship pos v =
  let open Nodes in
  let end_pos p s = p + String.length s in
  match v with
  | S _ | E _                             -> `Ignore
  | B (p, _) when pos < p                 -> `Before p
  | N (p, _) when pos < p                 -> `Before p
  | B (p, _) when pos = p                 -> `Exact
  | N (p, s) when pos = p                 -> `Exact
  | B _                                   -> `After
  | N (p, s) when pos < end_pos p s       -> `In (p, s)
  | N _ (*p, s) when pos >= end_pos p s*) -> `After

(*first_start, last_end, end_to_next_start_assoc *)
let reference_starts_and_ends lst =
  match lst with
  | []                  -> invalid_argf "Reference has no start and ends"
  | {start; end_} :: [] -> start, end_, []
  | {start; end_} :: t  ->
    let rec loop ep acc = function
      | []                  -> invalid_argf "stop before empty"
      | {start; end_} :: [] -> end_, (ep, start) :: acc
      | {start; end_} :: t  -> loop end_ ((ep, start) :: acc) t
    in
    let e, l = loop end_ [] t in
    start, e, l

let add_reference_elems g aset allele ref_elems =
  let module AS = (val aset : Alleles.Set) in
  let open Nodes in
  let bse () = AS.singleton allele in
  let add_start start_pos lst =
    let st = start_pos, allele in
    `Started (st, G.V.create (S st)) :: lst
  in
  let add_end end_pos ~st ~prev lst =
    G.add_edge_e g (G.E.create prev (bse ()) (G.V.create (E end_pos)));
    `Ended (st, end_pos) :: lst
  in
  let add_boundary ~st ~prev ~idx ~pos lst =
    let boundary_node = G.V.create (B (pos, idx)) in
    G.add_edge_e g (G.E.create prev (bse ()) boundary_node);
    `Started (st, boundary_node) :: lst
  in
  let add_seq ~st ~prev start s lst =
    let sequence_node = G.V.create (N (start, s)) in
    G.add_edge_e g (G.E.create prev (bse ()) sequence_node);
    `Started (st, sequence_node) :: lst
  in
  List.fold_left ref_elems ~init:[] ~f:(fun state e ->
    let open MSA_parser in
    match state, e with
    | []            , Start start_pos               -> add_start start_pos state
    | `Ended _ :: _ , Start start_pos               -> add_start start_pos state
    | []            , Boundary _
    | `Ended _ :: _ , Boundary _                    -> state        (* ignore *)
    | []            , al_el
    | `Ended _ :: _ , al_el                         ->
        invalid_argf "Unexpected %s after end for %s"
          (al_el_to_string al_el) allele
    | `Started (st, prev) :: tl, End end_pos          -> add_end end_pos ~st ~prev tl
    | `Started (st, prev) :: tl, Boundary {idx; pos } -> add_boundary ~st ~prev ~idx ~pos tl
    | `Started (st, prev) :: tl, Sequence {start; s } -> add_seq ~st ~prev start s tl
    | `Started (_, _) :: _,      Gap _                -> state       (* ignore gaps *)
    | `Started (_, _) :: tl,     Start sp             ->
        invalid_argf "Unexpected second start at %d for %s" sp allele)
  |> List.map ~f:(function
      | `Started _ -> invalid_argf "Still have a Started in %s ref" allele
      | `Ended (start, end_) -> { start; end_})
  |> List.sort ~cmp:(fun s1 s2 -> compare_start s1.start s2.start)

let test_consecutive_elements allele =
  let open MSA_parser in
  let open Nodes in
  function
  | `Gap close_pos      ->
      begin function
      | (End end_pos) :: tl when end_pos = close_pos        -> `End (close_pos, tl)
      | Sequence {start; s} :: tl when start = close_pos    ->
          let new_close = start + String.length s in
          `Continue (Some (N (start, s)), (`Sequence new_close), tl)
      | Start _ :: _                                        ->
          invalid_argf "For %s another start after a Gap closes %d." allele close_pos
      | []                                                  ->
          invalid_argf "For %s Empty list after Gap close %d." allele close_pos
      | Boundary _ :: _
      | Sequence _ :: _
      | Gap _ :: _
      | End _ :: _                                          -> `Close close_pos
      end
  | `Sequence close_pos ->
      begin function
      | End end_pos :: tl when end_pos = close_pos          -> `End (close_pos, tl)
      | Gap { gstart; length} :: tl when gstart = close_pos ->
          let new_close = gstart + length in
          `Continue (None, (`Gap new_close), tl)
      | Start _ :: _                                        ->
          invalid_argf "For %s another start after a Sequence close %d." allele close_pos
      | []                                                  ->
          invalid_argf "For %s empty list after Sequence close %d." allele close_pos
      | Boundary _ :: _
      | Sequence _ :: _
      | Gap _ :: _
      | End _ :: _                                          -> `Close close_pos
      end

exception Found of Nodes.t

let next_node_along g aset allele ~from =
  let module AS = (val aset : Alleles.Set) in
  try
    G.fold_succ_e (fun (_, bt, vs) n ->
        if AS.is_set bt allele then raise (Found vs) else n)
      g from None
  with Found v ->
    Some v

let add_non_ref g aset reference (first_start, last_end, end_to_next_start_assoc) allele alt_lst =
  let module AS = (val aset : Alleles.Set) in
  let open MSA_parser in
  let open Nodes in
  let first_start_node = S first_start in
  let last_end_node = E last_end in
  let end_to_start_nodes = List.map ~f:(fun (e, s) -> E e, S s) end_to_next_start_assoc in
  let next_reference ~msg from =
    match next_node_along g aset reference ~from with
    | Some n -> n
    | None   -> match List.Assoc.get from end_to_start_nodes with
                | Some n -> n
                | None   -> invalid_arg msg
  in
  let add_allele_edge pv nv = add_allele_edge g aset pv nv allele in
  let do_nothing _ _ = () in
  let advance_until ~visit ~prev ~next pos =
    let rec forward node msg =
      loop node (next_reference ~msg node)
    and loop prev next =
      if next = first_start_node then
        forward next "No next after reference Start!"
      else if next = last_end_node then
        `AfterLast prev
      else
        match relationship pos next with
        | `Ignore    -> forward next (sprintf "Skipping %s" (vertex_name next))
        | `Before ap -> `InGap (prev, next, ap)
        | `Exact     -> `AtNext (prev, next)
        | `In (p, s) -> `InsideNext (prev, next, p, s)
        | `After     -> visit prev next;
                        forward next "Not at End, should have next!"
    in
    loop prev next
  in
  let split_in ~prev ~visit ~next pos =
    let split_and_rejoin p s node =
      let open Nodes in
      let index = pos - p in
      let fs, sn = String.split_at s ~index in
      let pr = G.pred_e g node in
      let su = G.succ_e g node in
      G.remove_vertex g node;
      let v1 = N (p, fs) in
      G.add_vertex g v1;
      let s_inter =
        List.fold_left pr ~init:(AS.init ())
            ~f:(fun bta (p, bt, _) ->
                  G.add_edge_e g (G.E.create p bt v1);
                  AS.union bt bta)
      in
      let v2 = N (pos, sn) in
      G.add_vertex g v2;
      let s_inter = AS.clear s_inter allele in
      G.add_edge_e g (G.E.create v1 s_inter v2);
      List.iter su ~f:(fun (_, e, s) -> G.add_edge_e g (G.E.create v2 e s));
      (v1, v2)
    in
    match advance_until ~prev ~next ~visit pos with
    | `InsideNext (pv, nv, p, s) ->
        let v1, v2 = split_and_rejoin p s nv in
        visit pv v1;
        `AtNext (v1, v2)
    | `AfterLast _ as al -> al
    | `InGap _ as ig     -> ig
    | `AtNext _ as an    -> an
  in
  let rec advance_until_boundary ~visit ~prev ~next pos idx =
    let rec forward node msg =
      loop node (next_reference ~msg node)
    and loop pv nv =
      match nv with
      | S _ | E _ -> forward nv "Skipping start End"
      | B (p, c) when p = pos ->
          if c <> idx then
            invalid_argf "Boundary at %d position diff from reference %d count %d"
              p c idx
          else
            pv, nv
      | B (p, _)
      | N (p, _) when p < pos ->
          visit pv nv;
          forward nv (sprintf "Trying to find B %d %d after %d" pos idx p)
      | B (p, c) (*when p > pos*) ->
          invalid_argf "Next Boundary %d %d after desired boundary %d %d"
            p c pos idx
      | N (p, _) ->
          invalid_argf "Next Sequence position: %d at or after desired boundary pos %d (idx %d) %s"
            p pos idx allele
    in
    loop prev next
  in
  (* How we close back with the reference *)
  let rec rejoin_after_split ~prev ~next split_pos state ~new_node lst =
    match split_in ~prev ~next ~visit:do_nothing split_pos with
    | `AfterLast _          -> solo_loop state new_node lst
    | `InGap (pv, next, ap) -> ref_gap_loop state ~prev:new_node next ap lst
    | `AtNext (_pv, next)   -> add_allele_edge new_node next;
                               main_loop state ~prev:new_node ~next lst
  (* In the beginning we have not 'Start'ed ->
    Loop through the alignment elemends:
      - discarding Boundaries and Gaps
      - on a Start find the position in reference and start correct loop
      - complaining on anything other than a Start *)
  and start_loop previous_starts_and_ends lst =
    let rec find_start_loop = function
      | Boundary _ :: t
      | Gap _ :: t      -> find_start_loop t (* Ignore Gaps & Boundaries before Start *)
      | Start p :: t    -> Some ( p, t)
      | []              ->
        begin
          match previous_starts_and_ends with
          | [] -> invalid_argf "Failed to find start for %s." allele
          | ls -> None
        end
      | s :: _          -> invalid_argf "Encountered %s in %s instead of Start"
                            (al_el_to_string s) allele
    in
    match find_start_loop lst with
    | None -> previous_starts_and_ends (* fin *)
    | Some (start_pos, tl) ->
        let start = start_pos, allele in
        let state = start, previous_starts_and_ends in
        let new_node = G.V.create (S start) in
        (* we can think of a start as as a Gap *)
        close_position_loop state ~prev:first_start_node ~next:first_start_node
          ~allele_node:new_node (`Gap start_pos) tl
  and add_end (start, os) end_ prev tl =
    add_allele_edge prev (G.V.create (E end_));
    let ns = { start; end_ } :: os in
    if tl = [] then
      ns  (* fin *)
    else
      start_loop ns tl
  (* When the only thing that matters is the previous node. *)
  and solo_loop state prev = function
    | []                        -> invalid_argf "No End at allele sequence: %s" allele
    | Start p :: _              -> invalid_argf "Another start %d in %s allele sequence." p allele
    | End end_pos :: tl         -> add_end state end_pos prev tl
    | Boundary { idx; pos} :: t -> let boundary_node = G.V.create (B (pos, idx)) in
                                   add_allele_edge prev boundary_node;
                                   solo_loop state boundary_node t
    | Sequence { start; s} :: t -> let sequence_node = G.V.create (N (start, s)) in
                                   add_allele_edge prev sequence_node;
                                   solo_loop state sequence_node t
    | Gap _ :: t                -> solo_loop state prev t
  (* When traversing a reference gap. We have to keep check allele elements
    position to check when to join back with the next reference node. *)
  and ref_gap_loop state ~prev ref_node ref_pos = function
    | []                                -> invalid_argf "No End at allele sequence: %s" allele
    | Start p :: _                      -> invalid_argf "Another start %d in %s allele sequence." p allele
    | (End end_pos :: tl) as lst        ->
        if end_pos <= ref_pos then
          add_end state end_pos prev tl
        else (* end_pos > ref_pos *)
          let () = add_allele_edge prev ref_node in
          main_loop state ~prev ~next:ref_node lst
    | (Boundary { idx; pos} :: tl) as l ->
        if pos < ref_pos then
          invalid_argf "Allele %s has a boundary %d at %d that is in ref gap ending %d."
            allele idx pos ref_pos
        else if pos = ref_pos then
          if ref_node = B (pos, idx) then
            let () = add_allele_edge prev ref_node in
            main_loop state ~prev ~next:ref_node tl
          else
            invalid_argf "Allele %s has a boundary %d at %d where ref gap ends %d."
              allele idx pos ref_pos
        else
          let () = add_allele_edge prev ref_node in
          main_loop state ~prev ~next:ref_node l
    | (Sequence { start; s} :: tl) as l ->
        if start > ref_pos then begin
          add_allele_edge prev ref_node;
          main_loop state ~prev ~next:ref_node l
        end else
          let new_node = G.V.create (N (start, s)) in
          let () = add_allele_edge prev new_node in
          let close_pos = start + String.length s in
          if close_pos <= ref_pos then
            ref_gap_loop state ~prev:new_node ref_node ref_pos tl
          else (* close_pos > ref_pos *)
            rejoin_after_split ~prev:ref_node ~next:ref_node close_pos state
              ~new_node tl
    | (Gap { gstart; length} :: tl) as l ->
        if gstart > ref_pos then begin
          add_allele_edge prev ref_node;
          main_loop state ~prev ~next:ref_node l
        end else
          let close_pos = gstart + length in
          if close_pos <= ref_pos then
            ref_gap_loop state ~prev ref_node ref_pos tl
          else (* close_pos > ref_pos *)
            rejoin_after_split ~prev:ref_node ~next:ref_node close_pos state
              ~new_node:prev tl
  (* The "run-length"-encoding nature of alignment implies that our edges have
     to link to the next sequence element iff they are consecutive
     (ex. "A..C", "...A..C" ) as opposed to linking back to the reference! *)
  and close_position_loop state ~prev ~next ~allele_node pe lst =
    match test_consecutive_elements allele pe lst with
    | `End (end_pos, tl)                    ->
        add_end state end_pos allele_node  tl
    | `Close pos                            ->
        rejoin_after_split pos state ~prev ~next ~new_node:allele_node lst
    | `Continue (new_node_opt, new_pe, lst) ->
        let allele_node =
          match new_node_opt with
          | None    -> allele_node
          | Some nn -> add_allele_edge allele_node nn; nn
        in
        close_position_loop state ~prev ~next ~allele_node new_pe lst
  (* When not in a reference gap *)
  and main_loop state ~prev ~next = function
    | []              -> invalid_argf "No End at allele sequence: %s" allele
    | Start p :: _    -> invalid_argf "Another start %d in %s allele sequence." p allele
    | End end_pos :: tl  ->
        let prev =
          match split_in ~prev ~next ~visit:add_allele_edge end_pos with
          | `AfterLast p
          | `AtNext (p, _)
          | `InGap (p, _, _) -> p
        in
        add_end state end_pos prev tl
    | Boundary { idx; pos} :: t ->
        let prev, next =
          advance_until_boundary ~prev ~next ~visit:add_allele_edge
            pos idx
        in
        let () = add_allele_edge prev next in
        main_loop state ~prev ~next t
    | Sequence { start; s} :: t ->
        let new_node = G.V.create (N (start, s)) in
        let open_res = split_in ~prev ~next ~visit:add_allele_edge start in begin
        match open_res with
        | `AfterLast prev         ->
            let () = add_allele_edge prev new_node in
            solo_loop state new_node t
        | `InGap (prev, next, _)
        | `AtNext (prev, next)    ->
            let () = add_allele_edge prev new_node in
            let close_pos = start + String.length s in
            close_position_loop ~prev ~next ~allele_node:new_node state (`Sequence close_pos) t
        end
    | Gap {gstart; length} :: t ->
        let open_res = split_in ~prev ~next ~visit:add_allele_edge gstart in begin
        match open_res with
        | `AfterLast prev         ->
            solo_loop state next t
        | `InGap (prev, next, _)
        | `AtNext (prev, next)    ->
            let close_pos = gstart + length in
            close_position_loop ~prev ~next ~allele_node:prev state (`Gap close_pos) t
        end
  in
  start_loop [] alt_lst

(* TODO: use a real heap/pq ds?
  - The one in Graph doesn't allow you to control not putting in duplicates
*)
module NodeQueue = struct

  include Set.Make(struct
      type t = Nodes.t
      let compare = Nodes.compare_by_position_first
    end)

  let add_successors g = G.fold_succ add g

  let at_min_position q =
    let rec loop p q acc =
      if is_empty q then
        q, acc
      else
        let me = min_elt q in
        match acc with
        | [] -> loop (Nodes.position me) (remove me q) [me]
        | _  -> if Nodes.position me = p then
                  loop p (remove me q) (me :: acc)
                else
                  q, acc
    in
    loop min_int q []

end

module FoldAtSamePosition = struct

  let after_start_nodes g amap bounds =
    let module AM = (val amap : Alleles.Map) in
    let open Nodes in
    AM.fold bounds ~init:NodeQueue.empty
      ~f:(fun q sep_lst allele ->
            List.fold_left sep_lst ~init:q
              ~f:(fun q sep ->
                    NodeQueue.add_successors g (S (fst sep.start, allele)) q))

  let step g q =
    let without_old, amp = NodeQueue.at_min_position q in
    let withnew =
      List.fold_left amp ~init:without_old
        ~f:(fun q n ->
              if Nodes.is_seq_or_boundary n then
                NodeQueue.add_successors g n q
              else
                q)
    in
    withnew, amp

  let fold_from g q ~f ~init =
    let rec loop a q =
      if NodeQueue.is_empty q then
        a
      else
        let nq, amp = step g q in
        let na = f a amp in
        loop na nq
    in
    loop init q


  let fold_after_starts g amap bounds ~f ~init =
    fold_from g (after_start_nodes g amap bounds) ~f ~init

  let node_queue_with_start_and_successors_nodes g amap bounds =
    let module AM = (val amap : Alleles.Map) in
    let open Nodes in
    AM.fold bounds ~init:NodeQueue.empty
      ~f:(fun q sep_lst allele ->
            List.fold_left sep_lst ~init:q
              ~f:(fun q sep ->
                let s = S (fst sep.start, allele) in
                let nq = NodeQueue.add s q in
                NodeQueue.add_successors g s nq))

  let f g amap bounds ~f ~init =
    fold_from g (node_queue_with_start_and_successors_nodes g amap bounds) ~f ~init

end (* FoldAtSamePosition *)

let range_pr amap bounds =
  let module AM = (val amap : Alleles.Map) in
  AM.fold_wa bounds ~init:(max_int, min_int)
    ~f:(fun p sep_lst ->
          List.fold_left sep_lst ~init:p ~f:(fun (st, en) sep ->
            (min st (fst sep.start)), (max en sep.end_)))

(** [range g] returns the minimum and maximum alignment position in [g]. *)
let range { amap; bounds; _ } =
  range_pr amap bounds

let create_by_position g amap bounds =
  let st, en = range_pr amap bounds in
  let rec redirects dest num acc =
    if num <= 0 then
      acc
    else
      redirects dest (num - 1) ((Redirect dest) :: acc)
  in
  FoldAtSamePosition.(f g amap bounds ~init:(0, []) ~f:(fun (i,acc) lst ->
    let p = Nodes.position (List.hd_exn lst) in
    let j = p - st in
    let num_redirects = j - (i + 1) in
    (j, NL lst :: (redirects i num_redirects acc))))
  |> snd
  |> List.rev
  |> Array.of_list

module JoinSameSequencePaths = struct

  (* Alignment and sequence pair, set used as queue *)
  module Apsq =
    Set.Make(struct
      type t = alignment_position * sequence
      let compare (a1, s1) (a2, s2) =
        let r = compare_alignment_position a1 a2 in
        if r = 0 then
          compare_sequence s1 s2
        else
          r
    end)

  let peak_min_position q =
    if Apsq.is_empty q then
      None
    else
      Some (fst (Apsq.min_elt q))

  let at_min_position q =
    let rec loop p q acc =
      if Apsq.is_empty q then
        q, acc
      else
        let me = Apsq.min_elt q in
        match acc with
        | [] -> loop (fst me) (Apsq.remove me q) [me]
        | _  -> if fst me = p then
                  loop p (Apsq.remove me q) (me :: acc)
                else
                  q, acc
    in
    loop min_int q []

  let same_debug = ref false

  let rec add_node_successors_only g v q =
    let open Nodes in
    match v with
    | N (a, s) -> Apsq.add (a, s) q
    | E _      -> q
    | S _
    | B _      -> G.fold_succ (add_node_successors_only g) g v q

  let after_starts g amap bounds =
    let module AM = (val amap : Alleles.Map) in
    let open Nodes in
    let add_successors = G.fold_succ (add_node_successors_only g) g in
    AM.fold bounds ~init:Apsq.empty
      ~f:(fun q sep_lst allele ->
            List.fold_left sep_lst ~init:q
              ~f:(fun q sep ->
                    add_successors (S (fst sep.start, allele)) q))

  let do_it g aset amap bounds =
    let module AS = (val aset : Alleles.Set) in
    let open Nodes in
    let add_successors = G.fold_succ (add_node_successors_only g) g in
    let qstart = after_starts g amap bounds in
    let unite_edges pv nv e =
      try
        let ecur = G.find_edge g pv nv in
        let into = G.E.label ecur in
        AS.unite ~into e
      with Not_found ->
        G.add_edge_e g (G.E.create pv e nv)
    in
    let split_and_rejoin ~index q p s =
      let node = N (p, s) in
      let pr = G.pred_e g node in
      let su = G.succ_e g node in
      G.remove_vertex g node;
      let fs, sn = String.split_at s ~index in
      let p1 = p + index in
      let v1 = N (p, fs) in
      if not (G.mem_vertex g v1) then G.add_vertex g v1;
      let s_inter =
        List.fold_left pr ~init:(AS.init ())
          ~f:(fun bta (p, bt, _) ->
                unite_edges p v1 bt;
                AS.union bt bta)
      in
      let v2 = N (p1, sn) in
      if not (G.mem_vertex g v2) then G.add_vertex g v2;
      unite_edges v1 v2 s_inter;
      List.iter su ~f:(fun (_, e, s) -> unite_edges v2 s e);
      Apsq.add (p1, sn) q
    in
    (* assume the list isn't empty *)
    let find_highest_diff ~must_split lst =
      let init = Option.value must_split ~default:max_int in
      let max_compare_length =
        List.fold_left lst ~init ~f:(fun l (_p, s) ->
          min l (String.length s))
      in
      let arr = [| 0; 0; 0; 0; 0 |] in
      let clear_arr () = for i = 0 to 4 do arr.(i) <- 0 done in
      let fill_arr index =
        List.iter lst ~f:(fun (_p, s) ->
          let c = String.get_exn s ~index in
          let j = if c = 'N' then 4 else Kmer_to_int.char_to_int c in
          arr.(j) <- arr.(j) + 1)
      in
      let rec diff_loop index =
        if index >= max_compare_length then
          index
        else begin
          clear_arr ();   (* clear A,C,G,T,N counts in the beginning! *)
          fill_arr index;
          match arr with
          | [| n; 0; 0; 0; 0|] when n >= 1 -> diff_loop (index + 1)
          | [| 0; n; 0; 0; 0|] when n >= 1 -> diff_loop (index + 1)
          | [| 0; 0; n; 0; 0|] when n >= 1 -> diff_loop (index + 1)
          | [| 0; 0; 0; n; 0|] when n >= 1 -> diff_loop (index + 1)
          | [| 0; 0; 0; 0; n|] when n >= 1 -> diff_loop (index + 1)
          | _                            ->
              if !same_debug then
                eprintf "found diff %d: [| %d; %d; %d; %d; %d |]\n"
                  index arr.(0) arr.(1) arr.(2) arr.(3) arr.(4);
              if index > 0 then index else same_loop index
        end
      and same_loop index =
        let nindex = index + 1 in
        if nindex >= max_compare_length then
          index
        else begin
          clear_arr ();
          fill_arr nindex;
          let mx = Array.fold_left ~init:0 ~f:max arr in
          if !same_debug then
            eprintf "found diff %d: [| %d; %d; %d; %d; %d |]\n"
              index arr.(0) arr.(1) arr.(2) arr.(3) arr.(4);
          if mx = 1 then
            same_loop nindex
          else
            index
        end
      in
      diff_loop 0
    in
    let flatten ~next_pos q = function
      | []                -> q
      | (p, s) :: []      ->
          begin
            if !same_debug then
              eprintf "one %d %s next_pos %d\n%!"
                p s (Option.value next_pos ~default:(-1));
            match next_pos with
            | Some np when inside_seq p s ~pos:np ->
                split_and_rejoin ~index:(np - p) q p s
            | None | Some _ ->
                add_successors (N (p, s)) q
          end
      | (p, _) :: _ as ls ->
          let must_split = Option.map ~f:(fun np -> np - p) next_pos in
          let index = find_highest_diff ~must_split ls in
          if !same_debug then begin
            eprintf "index:\t%d must_split:\t%d\n" index (Option.value must_split ~default:(-1));
            List.iter ls ~f:(fun (p, s) -> eprintf "%d: %s\n" p s)
          end;
          if index = 0 then
            List.fold_left ls ~init:q ~f:(fun q (p, s) ->
              if String.length s <= 1 then
                add_successors (N (p, s)) q
              else
                split_and_rejoin ~index:1 q p s)
          else
            List.fold_left ls ~init:q ~f:(fun q (p, s) ->
            if String.length s = index then begin
              if !same_debug then
                eprintf "Not splitting %d %s because length is less than %d\n" p s index;
              (* Don't split to avoid an empty Node!*)
              add_successors (N (p, s)) q
            end else begin
              if !same_debug then
                eprintf "splitting %d %s at %d\n" p (index_string s index) index;
              split_and_rejoin ~index q p s
            end)
    in
    let rec loop q =
      if q = Apsq.empty then
        ()
      else
        let nq, amp = at_min_position q in
        let next_pos = peak_min_position nq in
        if !same_debug then
          eprintf "popping [%s] peaking at %d\n"
            (List.map ~f:(fun (p, s) -> sprintf "(%d,%s)" p s) amp
             |> String.concat ~sep:";")
            (Option.value next_pos ~default:(-1));
        let nq = flatten ~next_pos nq amp in
        loop nq
    in
    loop qstart

end (* JoinSameSequencePaths *)

let find_nodes_at_private offset posarr ~pos =
  (* TODO: At the end-of-the-day this array is still a bit hacky for this
     data-structure. I know ahead of time that for each valid position there
     must be a node! I should encode that in the types. *)
  let rec lookup_at redirect j =
    match posarr.(j) with
    | NL []         -> failwithf "Empty node list at %d" j
    | NL (h :: tl)  -> h, tl
    | Redirect prev ->
        if redirect then
          lookup_at false prev
        else
          failwithf "Double redirect at %d" j
  in
  lookup_at true (pos - offset)

let find_node_at_private ?allele ?ongap ~pos g aset offset posarr =
  let module AS = (val aset : Alleles.Set) in
  let module M = struct exception Found end in
  let open Nodes in
  let all_str, test =
    match allele with
    | None   -> "", (fun n -> Nodes.is_seq_or_boundary n && Nodes.inside pos n)
    | Some a ->
        let is_along_allele node =
          try
            G.fold_pred_e (fun (_, e, _) _b ->
              if AS.is_set e a then
                raise M.Found else false) g node false
          with M.Found ->
            true
        in
        (" " ^ a)
        , (fun n -> Nodes.is_seq_or_boundary n && is_along_allele n && Nodes.inside pos n)
  in
  let h, tl = find_nodes_at_private offset posarr ~pos in
  let lst = h :: tl in
  match List.find lst ~f:test with
  | Some n -> Ok n
  | None   ->
      begin match ongap with
      | None   -> error "%sin a gap at %d" all_str pos
      | Some f -> f lst
      end

let adjacents_debug_ref = ref false

(* Methods for finding adjacent nodes/edges combinations.:
   Nodes with the same (or greater; due to gaps) alignment position.

   Aside from checking the bounds and looking for Boundary positions,
   there is no way of knowing when we have found all adjacents. These
   algorithms expand a set of successively seen nodes at progressively
   farther distance (measured by edges) from the root node. At each
   step they recurse down and explore new adjacents. The two main methods
   {until} and {check_by_levels} provide a fold operation over these
   adjacents (via ~f and ~init) as they are discovered and also returns:
   1. A EdgeNodeSet of actual adjacents. adjacent nodes can have multiple
      edges leading into them.
   2. The final accumulator value.
   3. A stack of encountered nodes (probably not useful). *)
module Adjacents = struct

  (* [add_if_new cur node sibs] Add [node] to [kids] if it isn't in [cur]rent.
     Although [kids] and [cur] may represent nodes at different alignment
     position, to prevent traversing/recursing down the same paths as we
     discover new kids, we keep track of two sets.

     (check) + (add) -> O(log n) + O(log n) = O(2 log n
     as opposed to O(log n) for just adding but then we redo lots (how much?)
     work and can't add the is_empty stop condition in [down].  *)
  let add_if_new ~cur n kids =
    if NodeSet.mem n cur then kids else NodeSet.add n kids

  let add_if_new_and_above ~cur ~pos n kids =
    (* Start Nodes will have an equal node position to pos but will be "above"
       and therefore potentially point at adjacents.*)
    if NodeSet.mem n cur then kids else
      let open Nodes in
      match n with
      | S (np, _) when np = pos -> NodeSet.add n kids
      | _                       -> if Nodes.position n >= pos then kids else
                                     NodeSet.add n kids

  let siblings_and_adjacents pos ~if_new g ~new_nodes ~cur ~adjacents acc =
    let add ((_pn, _el, n) as e) (kids, adjs, acc) =
      if Nodes.position n >= pos then
        let is_new, (nadjs, nacc) = if_new e (adjs, acc) in
        if is_new then
          G.fold_pred (add_if_new_and_above ~cur:kids ~pos) g n kids, nadjs, nacc
        else
          kids, nadjs, nacc
      else
        let nkids = add_if_new ~cur n kids in
        nkids, adjs, acc
    in
    let init = NodeSet.empty, adjacents, acc in
    NodeSet.fold new_nodes ~f:(G.fold_succ_e add g) ~init

  let rec down if_new g pos acc ~adjacents ~new_nodes cur =
    let newer_nodes, nadjacents, nacc =
      siblings_and_adjacents pos ~if_new g ~new_nodes ~cur ~adjacents acc
    in
    let new_cur = NodeSet.union cur new_nodes in
    if NodeSet.is_empty newer_nodes then
      (* no new nodes: don't recurse! *)
      nadjacents, nacc, new_cur
    else
      down if_new g pos nacc ~adjacents:nadjacents ~new_nodes:newer_nodes new_cur

  let look_above g cur =
    NodeSet.fold cur ~f:(G.fold_pred (add_if_new ~cur) g) ~init:NodeSet.empty

  (* A general combination of [up] and [down] that is a parameterized
     search that should tell us when to stop recursing. *)
  let up_and ?max_height ?prev_edge_node_set ~init ~if_new ~down g ~pos node =
    let rec up i adjacents acc ~new_nodes cur =
      match max_height with
      | Some h when i >= h  -> adjacents, acc, NodeSet.union new_nodes cur
      | None | Some _       ->
          match down acc ~adjacents ~new_nodes cur with
          | `Stop t                            -> t
          | `Continue (nadjacents, nacc, ncur) ->
              let newer_nodes = look_above g new_nodes in
              if !adjacents_debug_ref then
                eprintf "Adding new newer_nodes: %s\n from new_nodes: %s\n"
                  (NodeSet.to_string newer_nodes) (NodeSet.to_string new_nodes);
              up (i + 1) nadjacents nacc ~new_nodes:newer_nodes ncur
    in
    let wrap_if_new e a = snd (if_new e a) in
    let initial_edge_set =
      Option.value prev_edge_node_set
        ~default:EdgeNodeSet.empty
    in
    let adj_strt, nacc = G.fold_pred_e wrap_if_new g node (initial_edge_set, init) in
    let new_nodes = G.fold_pred NodeSet.add g node NodeSet.empty in
    let cur = NodeSet.singleton node in
    up 0 adj_strt nacc ~new_nodes cur

  (* A more general method that checks whether to continue after each new
     adjacent. Since we throw an exception to terminate early, the search
     stack isn't correct.

    TODO: Need to benchmark to figure out which method is ultimately faster. *)
  let until (type a) ?max_height ~f ~init g node =
    let module M = struct exception F of a end in
    let if_new (_, e, n) ((adjacents, acc) as s) =
      let en = e, n in
      if EdgeNodeSet.mem en adjacents then false, s else
        match f e n acc with
        | `Stop r     -> raise (M.F r)
        | `Continue r -> true, (EdgeNodeSet.add en adjacents, r)
    in
    let pos = Nodes.position node in
    let downc = down if_new g pos  in
    let down acc ~adjacents ~new_nodes cur =
      try `Continue (downc acc ~adjacents ~new_nodes cur)
      with M.F r -> `Stop (adjacents, r, (NodeSet.union new_nodes cur))
    in
    up_and ?max_height ~init ~if_new ~down g ~pos node

  (* A less general method that requires an extra predicate function to check
     the stopping condition. The advantage is that we check less frequently and
     probably when it matters: when we increase the level of how far back in
     the graph we look for adjacents. *)
  let check_by_levels ?max_height ?prev_edge_node_set ~f ~stop ~init g node =
    let if_new (pn, e, n) ((adjacents, acc) as s) =
      if !adjacents_debug_ref then
        eprintf "if_new check of %s -> %s\n"
          (Nodes.vertex_name pn) (Nodes.vertex_name n);
      let en = e, n in
      if EdgeNodeSet.mem en adjacents then false, s else
        true, (EdgeNodeSet.add en adjacents, (f e n acc))
    in
    let pos = Nodes.position node in
    (* TODO: There is a corner cases where down has redundant calls to if_new.*)
    let downc = down if_new g pos in
    let down acc ~adjacents ~new_nodes cur =
      let (nadjacents, nacc, ptl) as state =
        downc acc ~adjacents ~new_nodes cur
      in
      if stop nacc then
        `Stop (nadjacents, nacc, (NodeSet.union new_nodes cur))
      else
        `Continue state
    in
    let init =
      Option.value_map prev_edge_node_set ~default:init
        ~f:(EdgeNodeSet.fold ~init ~f:(fun (e, n) a -> f e n a))
    in
    up_and ?max_height ?prev_edge_node_set ~init ~if_new ~down g ~pos node

end (* Adjacents *)

let within ?(include_ends=true) pos sep =
  (fst sep.start) <= pos &&
    ((include_ends && pos <= sep.end_) || pos < sep.end_)

let alleles_with_data_private ?include_ends aset amap bounds ~pos =
  let module AM = (val amap : Alleles.Map) in
  let module AS = (val aset : Alleles.Set) in
  AM.fold bounds ~init:(AS.init ()) ~f:(fun es sep_lst allele ->
    List.fold_left sep_lst ~init:es ~f:(fun es sep ->
      if within ?include_ends pos sep then
        AS.set es allele
      else
        es))

(* At or past *)
let adjacents_at_private ?max_edge_debug_length ?(max_height=10000)
  ?prev_edge_node_set ~pos g aset amap offset posarr bounds =
  let module AS = (val aset : Alleles.Set) in
  let max_length = max_edge_debug_length in
  let all_edges = alleles_with_data_private ~include_ends:false aset amap bounds ~pos in
  let prev_edge_node_set =
    Option.map prev_edge_node_set ~f:(EdgeNodeSet.fold ~init:EdgeNodeSet.empty
        ~f:(fun ((edge, node) as en) nes ->
            if (Nodes.position node > pos || Nodes.inside pos node)
            (* Can be optimized with an non-empty intersect method *)
            && AS.(cardinal (inter all_edges edge)) > 0 then
              EdgeNodeSet.add en nes
            else
              nes))
  in
  let rootn, _ = find_nodes_at_private offset posarr ~pos in
  let stop es_acc =
    if es_acc = all_edges then true else
      begin
        if !adjacents_debug_ref then
          eprintf "Still missing\n1:%s\n2:%s\n"
            (AS.to_human_readable ?max_length
              (AS.diff all_edges es_acc))
            (AS.to_human_readable ?max_length
              (AS.diff es_acc all_edges));
          false
      end
  in
  let f edge node edge_set =
    match node with
    | Nodes.E p when p <= pos -> edge_set
    | _ ->
      if !adjacents_debug_ref then
        eprintf "Adding %s <- %s.\n"
          (Nodes.vertex_name node)
          (AS.to_human_readable ?max_length edge);
      AS.union edge edge_set
  in
  let init = AS.init () in
  Adjacents.check_by_levels ?prev_edge_node_set ~max_height ~init ~stop g rootn ~f

let create_adjacents_arr g aset amap offset posarr bounds =
  let st, en = range_pr amap bounds in
  let len = en - st + 1 in
  let pensr = ref None in
  Array.init len ~f:(fun i ->
    let pos = offset + i in
    let (edge_node_set, seen_alleles, _) =
      adjacents_at_private ?prev_edge_node_set:!pensr g aset amap offset posarr bounds ~pos
    in
    pensr := Some edge_node_set;
    {edge_node_set; seen_alleles})

type construction_arg =
  { selectors           : Alleles.Selectors.t list
  ; join_same_sequence  : bool
  }

let default_construction_arg =
  { selectors          = []
  ; join_same_sequence = true
  }

let construction_arg_to_string
  { selectors; join_same_sequence } =
    sprintf "%s_%b"
      (Alleles.Selectors.list_to_string selectors)
      join_same_sequence

let construct_from_parsed ?(merge_map=[]) ?(arg=default_construction_arg) r =
  let open MSA_parser in
  let { selectors ; join_same_sequence; } = arg in
  let { align_date; reference; ref_elems; alt_elems} = r in
  let alt_elems = List.sort ~cmp:(fun (n1, _) (n2, _) -> Alleles.compare n1 n2) alt_elems in
  let alt_alleles = Alleles.Selectors.apply_to_assoc selectors alt_elems in
  let num_alleles = List.length alt_alleles in
  let ref_length = List.length ref_elems in
  let g = G.create ~size:(ref_length * num_alleles) () in
  let aindex = Alleles.index (reference :: List.map ~f:fst alt_alleles) in
  let module Aset = Alleles.MakeSet (struct let index = aindex end) in
  let aset = (module Aset : Alleles.Set) in
  let module Amap = Alleles.MakeMap (struct let index = aindex end) in
  let amap = (module Amap : Alleles.Map) in
  let refs_start_ends = add_reference_elems g aset reference ref_elems in
  let fs_ls_st_assoc = reference_starts_and_ends refs_start_ends in
  let start_and_stop_assoc =
    List.fold_left alt_alleles ~init:[ reference, refs_start_ends]
      ~f:(fun acc (allele, lst) ->
            let start_and_stops = add_non_ref g aset reference fs_ls_st_assoc allele lst in
            (allele, start_and_stops) :: acc)
  in
  let bounds =
    Amap.init (fun allele ->
      match List.Assoc.get allele start_and_stop_assoc with
      | Some alst -> List.rev alst
      | None      -> assert false)
  in
  if join_same_sequence then
    JoinSameSequencePaths.do_it g aset amap bounds; (* mutates g *)
  let offset = fst (range_pr amap bounds) in
  let posarr = create_by_position g amap bounds in
  let adjacents_arr = create_adjacents_arr g aset amap offset posarr bounds in
  { align_date
  ; reference
  ; g
  ; aindex
  ; aset
  ; amap
  ; bounds
  ; offset
  ; posarr
  ; merge_map
  ; adjacents_arr
  }

let construct ?arg input =
  Alleles.Input.construct input >>= fun (mp, merge_map) ->
      Ok (construct_from_parsed ~merge_map ?arg mp)

(* More powerful accessors *)
let all_bounds { bounds; amap; _} allele =
  let module AM = (val amap : Alleles.Map) in
  AM.get bounds allele

let find_bound g allele ~pos =
  all_bounds g allele
  |> List.find ~f:(fun sep -> (fst sep.start) <= pos && pos <= sep.end_)
  |> Option.value_map ~default:(error "%d is not in any bounds" pos)
        ~f:(fun s -> Ok s)

let fold_along_allele g aset allele ~start ~f ~init =
  if not (G.mem_vertex g start) then
    invalid_argf "%s is not a vertex in the graph." (Nodes.vertex_name start)
  else
    let next = next_node_along g aset allele in
    let rec loop from (acc, stop) =
      match stop with
      | `Stop     -> acc
      | `Continue ->
          match next ~from with
          | None    -> acc
          | Some vs -> loop vs (f acc vs)
    in
    loop start (f init start)

(* This methods advantage is that [next_after] works ! *)
let find_position_old ?(next_after=false) t ~allele ~pos =
  let open Nodes in
  let next_after o v = if next_after then Some (v, false) else o in
  find_bound t allele ~pos >>= fun bound ->
    let start = S (fst bound.start, allele) in
    fold_along_allele t.g t.aset allele ~start ~init:None
      ~f:(fun o v ->
            match v with
            | S (p, _) when p <= pos                   -> o, `Continue
            | S (p, _) (*   p > pos *)                 -> next_after o v, `Stop
            | E p                                      -> next_after o v, `Stop
            | B (p, _) when p = pos                    -> Some (v, true), `Stop
            | B (p, _) when p > pos                    -> next_after o v, `Stop
            | B (p, _) (*   p < pos*)                  -> o, `Continue
            | N (p, s) when p = pos                    -> Some (v, true), `Stop
            | N (p, s) when inside_seq p s ~pos        -> Some (v, false), `Stop
            | N (p, s) when p < pos                    -> o, `Continue
            | N (p, s)  (* p + String.length s > pos*) -> next_after o v, `Stop)
    |> Option.value_map ~default:(error "%d in a gap" pos)
        ~f:(fun x -> Ok x)

let find_node_at ?allele ?ongap ~pos { g; aset; offset; posarr; _} =
  find_node_at_private ?allele ?ongap ~pos g aset offset posarr

(* - A list of nodes to start searching from. If there is more than one it
     is a list of the start nodes.
   - A function to modify the accumulated string. For example when we wish
     to start the sequence at a position inside of a node.
   - An initialization of a count of how many characters we've seen. This
     accounts for starting a node that we're going to trim with the function
     of the previous argument.  *)
let parse_start_arg g allele =
  let module AM = (val g.amap : Alleles.Map) in
  let id x = x in
  let open Nodes in function
  | Some (`Node s)        ->
      if G.mem_vertex g.g s then
        Ok ([s], id, 0)
      else
        error "%s vertex not in graph" (vertex_name s)
  | Some (`AtPos pos)     ->
      find_node_at ~allele ~pos g >>= fun v ->
        let pv = position v in
        let index = pos - pv in
        Ok ([v], String.drop ~index, -index)
  | Some (`AtNext pos)    ->
      (* find_node_at is O(1) but it doesn't currently accomodate
         finding the next node if necessary. *)
      find_position_old ~next_after:true ~allele ~pos g >>= fun (v, precise) ->
        let pv = position v in
        let index = pos - pv in
        if precise || index < 0 then
          Ok ([v], id, 0)
        else
          Ok ([v], String.drop ~index, -index)
  | Some (`Pad pos)       ->
      find_position_old ~next_after:true ~allele ~pos g >>= fun (v, precise) ->
        let pv = position v in
        let index = pos - pv in
        if precise then
          Ok ([v], id, 0)
        else if index < 0 then
          Ok ([v], (fun s -> (String.make (-index) '.') ^ s), index)
        else
          Ok ([v], String.drop ~index, -index)
  | None                  ->
      begin match AM.get g.bounds allele with
      | exception Not_found -> error "Allele %s not found in graph!" allele
      | []                  -> error "Allele %s not found in graph!" allele
      | spl                 -> Ok (List.map spl ~f:(fun sep ->
                                 Nodes.S (fst sep.start, allele)), id, 0)
      end

let parse_stop_arg ?(count=0) =
  let stop_if b = if b then `Stop else `Continue in
  let open Nodes in function
  | None             -> (fun _ -> `Continue)
                        , (fun x -> x)
  | Some (`AtPos p)  -> (fun n -> stop_if (position n >= p))
                        , (fun x -> x)
  | Some (`Length n) -> let r = ref count in
                        (function
                          | S _ | E _ | B _ -> `Continue
                          | N (_, s) ->
                              r := !r + String.length s;
                              stop_if (!r >= n))
                        , String.take ~index:n
  | Some (`Pad n)    -> let r = ref count in
                        (function
                          | S _ | E _ | B _ -> `Continue
                          | N (_, s) ->
                              r := !r + String.length s;
                              stop_if (!r >= n))
                        , (fun s ->
                            let l = String.length s in
                            if l < n then
                              s ^ String.make (n - l) '.'
                            else
                              String.take ~index:n s)

(** Accessors. *)
let sequence ?(boundaries=false) ?start ?stop gt allele =
  let open Nodes in
  parse_start_arg gt allele start
    >>= fun (start, pre, count) ->
      let stop, post = parse_stop_arg ~count stop in
      List.fold_left start ~init:[] ~f:(fun acc start ->
        fold_along_allele gt.g gt.aset allele ~start ~init:acc
          ~f:(fun clst node ->
                match node with
                | N (_, s)             -> (s :: clst, stop node)
                | B _  when boundaries -> ("|" :: clst, stop node)
                | S _ | E _ | B _      -> (clst, stop node)))
      |> List.rev
      |> String.concat
      |> pre
      |> post
      |> fun s -> Ok s

let alleles_with_data ?include_ends { aset; amap; bounds; _} ~pos =
  alleles_with_data_private ?include_ends aset amap  bounds ~pos

(* This method is super awkward, since it is mirroring the
  gap searching logic of Adjacents, but not as carefully/exhaustively. *)
let search_through_gap g node ~pos =
  let preds =
    G.fold_pred_e (fun (pn, e, _) l -> (Nodes.position pn, pn, e) :: l) g node []
    |> List.sort ~cmp:compare (* should be precise about comparing by position. *)
  in
  let back_look =
    List.find_map preds ~f:(fun (p, n, e) ->
      if Nodes.inside pos n then Some (n, e) else None)
  in
  match back_look with
  | Some v -> Ok v
  | None   -> (* try looking forward again! *)
      let frds =
        List.fold_left preds ~init:[] ~f:(fun acc (_, n, _) ->
          G.fold_succ_e (fun (_, e, n) l -> (Nodes.position n, n, e) :: l)
            g n [])
        |> List.sort ~cmp:compare
      in
      let forward_look =
        List.find_map frds ~f:(fun (p, n, e) ->
          if Nodes.inside pos n then Some (n, e) else None)
      in
      match forward_look with
      | Some v -> Ok v
      | None   -> error "Couldn't work through gap before %s for %d"
                    (Nodes.vertex_name node) pos

let adjacents_at ?max_edge_debug_length ?max_height ~pos
    { offset; adjacents_arr; _} =
  adjacents_arr.(pos - offset)
