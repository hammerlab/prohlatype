
open Util
open Graph

(* TODO:
  - Hashcons the sequences
  - Turn fold_succ_e from O(n) into something better
*)

type start =
  MSA.position *
  (Alleles.allele [@equal Alleles.equal] [@compare Alleles.compare])
    [@@deriving eq, ord]
type sequence = string [@@deriving eq, ord, show]

(* start end pairs *)
type sep = { start : start ; end_ : MSA.position } [@@deriving eq, ord]

let blts = MSA.boundary_label_to_string

module Nodes = struct

  type t =
    | S of start
    | E of MSA.position
    | B of MSA.position * MSA.boundary_label(* Boundary of position and label *)
    | N of MSA.position * sequence                               (* Sequences *)
    [@@deriving eq, ord]

  let vertex_name ?(short=true) = function
    | S (n, s)  -> sprintf "S%d-%s" n s
    | E n       -> sprintf "E%d" n
    | B (p, l)  -> sprintf "B-%s-%d" (blts l) p
    | N (p, s)  -> sprintf "%d-%s(%d)"
                    p (if short then short_seq s else s)
                    (p + String.length s)

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

  let is_boundary = function
    | B _             -> true
    | N _ | S _ | E _ -> false

  let is_seq = function
    | N _             -> true
    | S _ | E _ | B _ -> false

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

(* We want to be able to have a 'by-position' array of all the nodes of the
   graph that intersect that position. Since some positions might contain the
   same set of nodes we distinguish two cases: *)
type by_position =
  | NL of Nodes.t list     (* A list of nodes that start at a given position. *)
  | Redirect of int        (* A redirect (index into array - NOT position) of
                              where to find the previous Node list. In
                              [by_position_array] below, that location should
                              contain the nodes that span this position. *)

and by_position_array = by_position array

module NodeSet = struct
  include Set.Make (Nodes)
  let to_string s =
    sprintf "{%s}" (string_of_list (elements s) ~sep:"; " ~f:Nodes.vertex_name)

end (* NodeSet *)

module EdgeNodeSet = struct

  include Set.Make (struct
    (* edges first since these edges point into the respective node. *)
    type t = Edges.t * Nodes.t

    let compare (e1, n1) (e2, n2) =
      let r = Nodes.compare n1 n2 in
      if r = 0 then Edges.compare e1 e2 else r
  end)

  (* Expose these once we revert the Alleles.Set functorization
  let to_string ?max_length ?complement s =
    sprintf "{%s}"
      (string_of_list (elements s) ~sep:"; "
        ~f:(fun (e, n) ->
              sprintf "(%s -> %s)"
                (Alleles.Set.to_human_readable ?max_length ?complement e)
                (Nodes.vertex_name n)))


  let to_table ?max_length ?complement s =
    let uns = elements s in
    let lst = List.sort uns ~cmp:(fun (_,n1) (_,n2) -> Nodes.compare n1 n2)  in
    sprintf "%s"
      (string_of_list lst ~sep:"\n"
        ~f:(fun (e, n) ->
              sprintf "%s <- %s"
                (Nodes.vertex_name n)
                (Alleles.Set.to_human_readable ?max_length ?complement e))) *)

end (* EdgeNodeSet *)

type t =
  { align_date    : string                   (* When the alignment was created by IMGT. *)
  ; reference     : string                                (* Name of reference allelle. *)
  ; g             : G.t                                            (* The actual graph. *)
  ; aindex        : Alleles.index               (* The allele index, for Sets and Maps. *)
  ; aset          : (module Alleles.Set)
  ; amap          : (module Alleles.Map)
  ; bounds        : sep list Alleles.map    (* Map of where the alleles start and stop. *)
  ; posarr        : by_position_array                          (* Per position lookups. *)
  ; merge_map     : (string * MSA.Alteration.t list) list
  ; offset        : int
  }

let number_of_alleles t = 
  t.aindex.Alleles.size

(** Some accessors that are necessary for output (and construction) but placed
    here for convenience. *)

(** [starts_by_position g] returns an associated list of positions and
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
  (* UGH. Unfortunately a necessary guard against
     1. Not typing this logic better.
     2. OCamlgraph find_edge stack-overflow *)
  if pv = nv then () else
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

(** Output

     TODO: When constructing the dot files, it would be nice if the alleles (edges), were
           in some kind of consistent order. *)


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
  output_dot ?human_edges ?compress_edges ?compress_start ?max_length ~short
    (fname ^ ".dot") t;
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
  | B _ (*   when pos > p *)              -> `After
  | N (p, s) when pos < end_pos p s       -> `In (p, s)
  | N _ (*p, s) when pos >= end_pos p s*) -> `After

exception GraphConstruction of string

let gc fmt =
  ksprintf (fun s -> raise (GraphConstruction s))
    fmt

(*first_start, last_end, end_to_next_start_assoc *)
let reference_starts_and_ends lst =
  match lst with
  | []                  -> gc "Reference has no start and ends"
  | {start; end_} :: [] -> start, end_, []
  | {start; end_} :: t  ->
    let rec loop ep acc = function
      | []                  -> gc "stop before empty"
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
  let add_boundary ~st ~prev ~label ~pos lst =
    let boundary_node = G.V.create (B (pos, label)) in
    G.add_edge_e g (G.E.create prev (bse ()) boundary_node);
    `Started (st, boundary_node) :: lst
  in
  let add_seq ~st ~prev start s lst =
    let sequence_node = G.V.create (N (start, s)) in
    G.add_edge_e g (G.E.create prev (bse ()) sequence_node);
    `Started (st, sequence_node) :: lst
  in
  List.fold_left ref_elems ~init:[] ~f:(fun state e ->
    let open MSA in
    match state, e with
    | []            , Start start_pos               -> add_start start_pos state
    | `Ended _ :: _ , Start start_pos               -> add_start start_pos state
    | []            , Boundary _
    | `Ended _ :: _ , Boundary _                    -> state        (* ignore *)
    | []            , al_el
    | `Ended _ :: _ , al_el                         ->
        gc "Unexpected %s after end for %s" (al_el_to_string al_el) allele
    | `Started (st, prev) :: tl, End end_pos            -> add_end end_pos ~st ~prev tl
    | `Started (st, prev) :: tl, Boundary {label; pos } -> add_boundary ~st ~prev ~label ~pos tl
    | `Started (st, prev) :: tl, Sequence {start; s }   -> add_seq ~st ~prev start s tl
    | `Started (_, _) :: _,      Gap _                  -> state       (* ignore gaps *)
    | `Started (_, _) :: tl,     Start sp               ->
        gc "Unexpected second start at %d for %s" sp allele)
  |> List.map ~f:(function
      | `Started _ -> gc "Still have a Started in %s ref" allele
      | `Ended (start, end_) -> { start; end_})
  |> List.sort ~cmp:(fun s1 s2 -> compare_start s1.start s2.start)

let test_consecutive_elements allele =
  let open MSA in
  let open Nodes in
  function
  | `Gap close_pos      ->
      begin function
      | (End end_pos) :: tl when end_pos = close_pos        -> `End (close_pos, tl)
      | Sequence {start; s} :: tl when start = close_pos    ->
          let new_close = start + String.length s in
          `Continue (Some (N (start, s)), (`Sequence new_close), tl)
      | Start _ :: _                                        ->
          gc "For %s another start after a Gap closes %d." allele close_pos
      | []                                                  ->
          gc "For %s Empty list after Gap close %d." allele close_pos
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
          gc "For %s another start after a Sequence close %d." allele close_pos
      | []                                                  ->
          gc "For %s empty list after Sequence close %d." allele close_pos
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
  let open MSA in
  let open Nodes in
  let first_start_pos = fst first_start in
  let first_start_node = S first_start in
  let last_end_node = E last_end in
  let end_to_start_nodes = List.map ~f:(fun (e, s) -> E e, S s) end_to_next_start_assoc in
  let next_reference ~msg from =
    match next_node_along g aset reference ~from with
    | Some n -> n
    | None   -> match List.Assoc.get from end_to_start_nodes with
                | Some n -> n
                | None   -> gc "In %s: %s" allele msg
  in
  let add_allele_edge pv nv =
    (*let () = printf "adding allele edge %s -- (%s) -> %s"
      (vertex_name pv) allele (vertex_name nv)
    in *)
    add_allele_edge g aset pv nv allele in
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
                        forward next
                          (sprintf "At pos: %d, next: %s not at End, should have next!"
                              pos (Nodes.vertex_name next))
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
  let rec advance_until_boundary ~visit ~prev ~next pos label =
    let rec forward node msg =
      loop node (next_reference ~msg node)
    and loop pv nv =
      match nv with
      | S _ | E _ -> forward nv "Skipping start End"
      | B (p, c) when p = pos ->
          if c <> label then
            gc "Boundary at %d position diff from reference %s label %s"
              p (blts c) (blts label)
          else
            pv, nv
      | B (p, _)
      | N (p, _) when p < pos ->
          visit pv nv;
          forward nv (sprintf "Trying to find B %d %s after %d"
            pos (blts label) p)
      | B (p, c) (*when p > pos*) ->
          gc "Next Boundary %d %s after desired boundary %d %s"
            p (blts c) pos (blts label)
      | N (p, _) ->
          gc "Next Sequence position: %d at or after desired boundary pos %d (label %s) %s"
            p pos (blts label) allele
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
      - on a Start find the position in reference and start construction loop
      - complaining on anything other than a Start *)
  and start_loop previous_starts_and_ends lst =
    let rec find_start_loop = function
      | Boundary _ :: t
      | Gap _ :: t    -> find_start_loop t (* Ignore Gaps & Boundaries before Start *)
      | Start p :: t  -> Some ( p, t)
      | []            ->
        begin
          match previous_starts_and_ends with
          | [] -> gc "Failed to find start for %s." allele
          | ls -> None
        end
      | s :: _        -> gc "Encountered %s in %s instead of Start"
                            (al_el_to_string s) allele
    in
    match find_start_loop lst with
    | None -> previous_starts_and_ends (* fin *)
    | Some (start_pos, tl) ->
        let start = start_pos, allele in
        let state = start, previous_starts_and_ends in
        let new_node = G.V.create (S start) in
        if start_pos < first_start_pos then
          ref_gap_loop state ~prev:new_node first_start_node first_start_pos tl
        else
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
    | []                          -> gc "No End at allele sequence: %s" allele
    | Start p :: _                -> gc "Another start %d in %s allele sequence." p allele
    | End end_pos :: tl           -> add_end state end_pos prev tl
    | Boundary { label; pos} :: t -> let boundary_node = G.V.create (B (pos, label)) in
                                     add_allele_edge prev boundary_node;
                                     solo_loop state boundary_node t
    | Sequence { start; s} :: t   -> let sequence_node = G.V.create (N (start, s)) in
                                     add_allele_edge prev sequence_node;
                                     solo_loop state sequence_node t
    | Gap _ :: t                  -> solo_loop state prev t
  (* When traversing a reference gap. We have to keep track of allele element's
     position to check when to join back with the next reference node. *)
  and ref_gap_loop state ~prev ref_node ref_pos lst = match lst with
    | []                          -> gc "No End at allele sequence: %s" allele
    | Start p :: _                -> gc "Another start %d in %s allele sequence." p allele
    | End end_pos :: tl           ->
        if end_pos <= ref_pos then
          add_end state end_pos prev tl
        else (* end_pos > ref_pos *)
          let () = add_allele_edge prev ref_node in
          main_loop state ~prev ~next:ref_node lst
    | Boundary { label; pos} :: tl  ->
        if pos < ref_pos then begin
          if ref_pos = first_start_pos then begin
            let bn = G.V.create (B (pos, label)) in
            let () = add_allele_edge prev bn in
            ref_gap_loop state ~prev:bn ref_node ref_pos tl
          end else
            gc "Allele %s has a boundary %s at %d that is in ref gap ending %d."
              allele (blts label) pos ref_pos
        end else if pos = ref_pos then
          if ref_node = B (pos, label) then
            let () = add_allele_edge prev ref_node in
            main_loop state ~prev ~next:ref_node tl
          else
            gc "Allele %s has a boundary %s at %d where ref gap ends %d."
              allele (blts label) pos ref_pos
        else (* pos > ref_pos *)
          let () = add_allele_edge prev ref_node in
          main_loop state ~prev ~next:ref_node lst
    | Sequence { start; s} :: tl  ->
        if start > ref_pos then begin
          add_allele_edge prev ref_node;
          main_loop state ~prev ~next:ref_node lst
        end else
          let new_node = G.V.create (N (start, s)) in
          let () = add_allele_edge prev new_node in
          let close_pos = start + String.length s in
          if close_pos < ref_pos then begin
            ref_gap_loop state ~prev:new_node ref_node ref_pos tl
          end else (* close_pos => ref_pos *)
            rejoin_after_split ~prev:ref_node ~next:ref_node close_pos state
              ~new_node tl
    | Gap { gstart; length} :: tl ->
        if gstart > ref_pos then begin
          add_allele_edge prev ref_node;
          main_loop state ~prev ~next:ref_node lst
        end else
          let close_pos = gstart + length in
          if close_pos < ref_pos then
            ref_gap_loop state ~prev ref_node ref_pos tl
          else (* close_pos >= ref_pos *)
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
    | []              -> gc "No End at allele sequence: %s" allele
    | Start p :: _    -> gc "Another start %d in %s allele sequence." p allele
    | End end_pos :: tl  ->
        let prev =
          match split_in ~prev ~next ~visit:add_allele_edge end_pos with
          | `AfterLast p
          | `AtNext (p, _)
          | `InGap (p, _, _) -> p
        in
        add_end state end_pos prev tl
    | Boundary { label; pos} :: t ->
        let prev, next =
          advance_until_boundary ~prev ~next ~visit:add_allele_edge
            pos label
        in
        let () = add_allele_edge prev next in
        main_loop state ~prev:next ~next t
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

  let at_min_position g q =
    let rec loop p q acc =
      if is_empty q then
        q, acc
      else
        let me = min_elt q in
        if Nodes.position me = p then
          loop p (add_successors g me (remove me q)) (me :: acc)
        else
          q, acc
    and start q =
      let me = min_elt q in
      let nq = remove me q in
      loop (Nodes.position me) (add_successors g me nq) [me]
    in
    start q

end

module FoldAtSamePosition = struct

  (* We know all of the start positions based upon the bounds map; so use that
     to fold upon start nodes to seed a queue. *)
  let after_start_nodes g amap bounds =
    let module AM = (val amap : Alleles.Map) in
    let open Nodes in
    AM.fold bounds ~init:NodeQueue.empty
      ~f:(fun q sep_lst allele ->
            List.fold_left sep_lst ~init:q
              ~f:(fun q sep ->
                    NodeQueue.add_successors g (S (fst sep.start, allele)) q))

  let step = NodeQueue.at_min_position

  let nl_to_string =
    string_of_list ~sep:";" ~f:Nodes.vertex_name

  let nq_to_string nq =
    NodeQueue.fold nq ~init:[] ~f:(fun n acc -> n :: acc)
    |> nl_to_string
    |> sprintf "nq: %s"

  (* Fold from an initialized queue. *)
  let fold_from g ~f ~init q =
    let rec loop a q =
      if NodeQueue.is_empty q then
        a
      else
        let nq, amp = step g q in
        (*if !debug then printf "at: %s\n%!" (nl_to_string amp); *)
        let na = f a amp in
        loop na nq
    in
    loop init q

  let fold_after_starts g amap bounds ~f ~init =
    fold_from g ~f ~init (after_start_nodes g amap bounds)

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
    fold_from g ~f ~init
      (node_queue_with_start_and_successors_nodes g amap bounds)

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

  let debug = ref false

  (* Alignment and sequence pair, set used as queue *)
  module Apsq =
    Set.Make(struct
      type t = MSA.position * sequence
      let compare (a1, s1) (a2, s2) =
        let r = MSA.compare_position a1 a2 in
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

  (* assume the list isn't empty *)
  let find_highest_diff ~must_split lst =
    if !debug then
      printf "find_highest_diff: ms:%s %s\n"
        (match must_split with None -> "None" | Some l -> sprintf "Some %d" l)
        (string_of_list lst ~sep:";" ~f:(fun (p, s) -> sprintf "%d,%s" p s));

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
            if !debug then
              printf "found diff %d: [| %d; %d; %d; %d; %d |]\n"
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
        if !debug then
          printf "found diff %d: [| %d; %d; %d; %d; %d |]\n"
            index arr.(0) arr.(1) arr.(2) arr.(3) arr.(4);
        if mx = 1 then
          same_loop nindex
        else
          index
      end
    in
    diff_loop 0

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
    let flatten ~next_pos q = function
      | []                -> q
      | (p, s) :: []      ->
          begin
            if !debug then
              printf "one %d %s next_pos %d\n%!"
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
          if !debug then begin
            printf "index:\t%d must_split:\t%d\n" index (Option.value must_split ~default:(-1));
            List.iter ls ~f:(fun (p, s) -> printf "%d: %s\n" p s)
          end;
          if index = 0 then
            List.fold_left ls ~init:q ~f:(fun q (p, s) ->
              if String.length s <= 1 then
                add_successors (N (p, s)) q
              else
                split_and_rejoin ~index:1 q p s)
          else
            List.fold_left ls ~init:q ~f:(fun q (p, s) ->
            if String.length s = index + 1 then begin
              if !debug then
                printf "Not splitting %d %s because length is less than %d\n" p s index;
              (* Don't split to avoid an empty Node!*)
              add_successors (N (p, s)) q
            end else begin
              if !debug then
                printf "splitting %d %s at %d\n" p (index_string s index) index;
              split_and_rejoin ~index q p s
            end)
    in
    let rec loop q =
      if q = Apsq.empty then
        ()
      else
        let nq, amp = at_min_position q in
        let next_pos = peak_min_position nq in
        if !debug then
          printf "popping [%s] peaking at %d\n"
            (string_of_list amp ~sep:";" ~f:(fun (p, s) -> sprintf "(%d,%s)" p s))
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
    | None   -> "", (fun n -> Nodes.is_seq n && Nodes.inside pos n)
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
        , (fun n -> Nodes.is_seq n && is_along_allele n && Nodes.inside pos n)
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

(* This used to have more options but they got refactored to other parts of the
   code. We'll preserve the record as a placeholder for other logic that we may
   want, (ex. remove reference). *)
type construction_arg =
  { join_same_sequence : bool
  }

let default_construction_arg = { join_same_sequence = true }

let construction_arg_to_string { join_same_sequence } =
  string_of_bool join_same_sequence

let debug = ref false

let construct_from_parsed ?(arg=default_construction_arg) r =
  let open MSA in
  let { join_same_sequence; } = arg in
  let { Parser.align_date; reference; ref_elems; alt_elems} = r in
  let num_alleles = List.length alt_elems in
  let ref_length = List.length ref_elems in
  let g = G.create ~size:(ref_length * num_alleles) () in
  let aindex = Alleles.index
    (reference :: List.map ~f:(fun a -> a.Parser.allele) alt_elems)
  in
  let module Aset = Alleles.MakeSet (struct let index = aindex end) in
  let aset = (module Aset : Alleles.Set) in
  let module Amap = Alleles.MakeMap (struct let index = aindex end) in
  let amap = (module Amap : Alleles.Map) in
  if !debug then
    printf "adding reference %s to graph.\n%s\n%!"
      reference (al_seq_to_string ~sep:"," ref_elems);
  let refs_start_ends = add_reference_elems g aset reference ref_elems in
  let fs_ls_st_assoc = reference_starts_and_ends refs_start_ends in
  let start_and_stop_assoc =
    List.fold_left alt_elems ~init:[ reference, refs_start_ends]
    ~f:(fun acc { Parser.allele; seq; alters}  ->
            try
              if !debug then
                printf "adding %s (alt: %s) to graph.\n%s\n%!"
                  allele
                  (string_of_list ~sep:";" ~f:Alteration.to_string alters)
                  (al_seq_to_string ~sep:"," seq);
              let start_and_stops =
                add_non_ref g aset reference fs_ls_st_assoc allele seq
              in
              (allele, start_and_stops) :: acc
            with (GraphConstruction s) ->
              eprintf "Graph construction error: %s skipping allele.\n%!" s;
              acc)
  in
  if !debug then printf "calculating bounds.\n%!";
  let bounds =
    Amap.init (fun allele ->
      match List.Assoc.get allele start_and_stop_assoc with
      | Some alst -> List.rev alst
      | None      -> [])
  in
  if !debug then printf "joining same\n%!";
  if join_same_sequence then
    JoinSameSequencePaths.do_it g aset amap bounds; (* mutates g *)
  let offset = fst (range_pr amap bounds) in
  if !debug then printf "creating by positions.\n%!";
  let posarr = create_by_position g amap bounds in
  if !debug then printf "done.\n%!";
  { align_date
  ; reference
  ; g
  ; aindex
  ; aset
  ; amap
  ; bounds
  ; offset
  ; posarr
  ; merge_map = List.map alt_elems ~f:(fun a -> a.Parser.allele, a.Parser.alters)
  }

let construct ?arg input =
  Alleles.Input.construct input >>= fun mp -> Ok (construct_from_parsed ?arg mp)

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

