(** Construction from Mas_parser *)

open Util
open Graph        (* OCamlgraph *)
(*open Ref_graph*)

let inv_argf ?(prefix="") fmt = ksprintf invalid_arg ("%s" ^^ fmt) prefix

(*module SSet = Set.Make (struct
    type t = string
    let compare (x : string) y = compare x y
  end) *)

let relationship pos v =
  let open Ref_graph.Nodes in
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
  let open Ref_graph in
  match lst with
  | []                  -> inv_argf "Reference has no start and ends"
  | {start; end_} :: [] -> start, end_, []
  | {start; end_} :: t  ->
    let rec loop ep acc = function
      | []                  -> inv_argf "stop before empty"
      | {start; end_} :: [] -> end_, (ep, start) :: acc
      | {start; end_} :: t  -> loop end_ ((ep, start) :: acc) t
    in
    let e, l = loop end_ [] t in
    start, e, l

(*let edges_between_to_string g pv nv =
  GE.find_all_edges g pv nv
  |> List.map ~f:(fun (_,l,_) -> l)
  |> String.concat ~sep:"," *)

let add_reference_elems g aset allele ref_elems =
  let open Mas_parser in
  let open Ref_graph in
  let open Nodes in
  let bse = Allele_set.singleton aset allele in
  let add_start start_pos lst =
    let st = start_pos, allele in
    `Started (st, GE.V.create (S st)) :: lst
  in
  let add_end end_pos ~st ~prev lst =
    GE.add_edge_e g (GE.E.create prev bse (GE.V.create (E end_pos)));
    `Ended (st, end_pos) :: lst
  in
  let add_boundary ~st ~prev ~idx ~pos lst =
    let boundary_node = GE.V.create (B (pos, idx)) in
    GE.add_edge_e g (GE.E.create prev bse boundary_node);
    `Started (st, boundary_node) :: lst
  in
  let add_seq ~st ~prev start s lst =
    let sequence_node = GE.V.create (N (start, s)) in
    GE.add_edge_e g (GE.E.create prev bse sequence_node);
    `Started (st, sequence_node) :: lst
  in
  List.fold_left ref_elems ~init:[] ~f:(fun state e ->
      match state, e with
      | []            , Start start_pos               -> add_start start_pos state
      | `Ended _ :: _ , Start start_pos               -> add_start start_pos state
      | []            , Boundary _
      | `Ended _ :: _ , Boundary _                    -> state        (* ignore *)
      | []            , al_el
      | `Ended _ :: _ , al_el                         ->
          inv_argf "Unexpected %s before start for %s"
            (al_el_to_string al_el) allele
      | `Started (st, prev) :: tl, End end_pos          -> add_end end_pos ~st ~prev tl
      | `Started (st, prev) :: tl, Boundary {idx; pos } -> add_boundary ~st ~prev ~idx ~pos tl
      | `Started (st, prev) :: tl, Sequence {start; s } -> add_seq ~st ~prev start s tl
      | `Started (_, _) :: _,      Gap _                -> state       (* ignore gaps *)
      | `Started (_, _) :: tl,     Start sp             ->
          inv_argf "Unexpected second start at %d for %s" sp allele)
  |> List.map ~f:(function
      | `Started _ -> inv_argf "Still have a Started in %s ref" allele
      | `Ended (start, end_) -> { start; end_})
  |> List.sort ~cmp:(fun s1 s2 -> compare_start s1.start s2.start)

let add_non_ref g reference aset (first_start, last_end, end_to_next_start_assoc) allele alt_lst =
  let open Mas_parser in
  let open Ref_graph in
  let open Nodes in
  let first_start_node = S first_start in
  let last_end_node = E last_end in
  let end_to_start_nodes = List.map ~f:(fun (e, s) -> E e, S s) end_to_next_start_assoc in
  let next_reference ~msg from =
    match next_node_along_ge aset reference g ~from with
    | Some n -> n
    | None   -> try List.assoc from end_to_start_nodes
                with Not_found -> invalid_arg msg
  in
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
      let pr = GE.pred_e g node in
      let su = GE.succ_e g node in
      GE.remove_vertex g node;
      let v1 = N (p, fs) in
      GE.add_vertex g v1;
      (* TODO: this set intersect is very clumsy *)
      let pr_edge_s =
        List.fold_left pr ~init:(Allele_set.empty aset)
            ~f:(fun bta (p, bt, _) ->
                  GE.add_edge_e g (GE.E.create p bt v1);
                  Allele_set.union bt bta)
       (* List.fold_left pr ~init:SSet.empty ~f:(fun set (p, l, _) ->
          GE.add_edge_e g (GE.E.create p l v1);
            SSet.add l set)*)
      in
      let v2 = N (pos, sn) in
      GE.add_vertex g v2;
      let su_edge_s =
        List.fold_left su ~init:(Allele_set.empty aset)
            ~f:(fun bta (_, bt, s) ->
                  GE.add_edge_e g (GE.E.create v2 bt s);
                  Allele_set.union bta bt)
        (*List.fold_left su ~init:SSet.empty ~f:(fun set (_, l, s) ->
          GE.add_edge_e g (GE.E.create v2 l s);
            SSet.add l set) *)
      in
      let s_inter = Allele_set.inter pr_edge_s su_edge_s in
      GE.add_edge_e g (GE.E.create v1 s_inter v2);
      (*SSet.iter (fun l -> GE.add_edge_e g (GE.E.create v1 l v2)) s_inter; *)
      (v1, v2)
    in
    match advance_until ~prev ~next ~visit pos with
    | `InsideNext (pv, nv, p, s) ->
        visit pv nv;
        let v1, v2 = split_and_rejoin p s nv in
        `AtNext (v1, v2)
    | `AfterLast _ as al -> al
    | `InGap _ as ig     -> ig
    | `AtNext _ as an    -> an
  in
  let add_allele_edge pv nv =
    try
      Allele_set.set_allele aset (GE.find_edge g pv nv |> GE.E.label) allele
    with Not_found ->
      let bt = Allele_set.singleton aset allele in
      GE.add_edge_e g (GE.E.create pv bt nv)
  in
  let rec advance_until_boundary ~visit ~prev ~next pos idx =
    let rec forward node msg =
      loop node (next_reference ~msg node)
    and loop pv nv =
      match nv with
      | S _ | E _ -> forward nv "Skipping start End"
      | B (p, c) when p = pos ->
          if c <> idx then
            inv_argf "Boundary at %d position diff from reference %d count %d"
              p c idx
          else
            pv, nv
      | B (p, _)
      | N (p, _) when p < pos ->
          visit pv nv;
          forward nv (sprintf "Trying to find B %d %d after %d" pos idx p)
      | B (p, c) (*when p > pos*) ->
          inv_argf "Next Boundary %d %d after desired boundary %d %d"
            p c pos idx
      | N (p, _) ->
          inv_argf "Next Sequence position: %d at or after desired boundary pos %d (idx %d)"
            p pos idx
    in
    loop prev next
  in
  (* How we close back with the reference *)
  let rec rejoin_after_split ~prev ~next split_pos state ~new_node lst =
    match split_in ~prev ~next ~visit:do_nothing split_pos with
    | `AfterLast _            -> solo_loop state new_node lst
    | `InGap (_pv, next, ap)  -> ref_gap_loop state new_node next ap lst
    | `AtNext (_pv, next)     -> add_allele_edge new_node next;
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
          | [] -> inv_argf "Failed to find start for %s." allele
          | ls -> None
        end
      | s :: _          -> inv_argf "Encountered %s in %s instead of Start"
                            (al_el_to_string s) allele
    in
    match find_start_loop lst with
    | None -> previous_starts_and_ends (* fin *)
    | Some (start_pos, tl) ->
        let start = start_pos, allele in
        let state = start, previous_starts_and_ends in
        let new_node = GE.V.create (S start) in
        rejoin_after_split ~prev:first_start_node ~next:first_start_node start_pos state
          ~new_node tl
  and add_end (start, os) end_ prev tl =
    add_allele_edge prev (GE.V.create (E end_));
    let ns = { start; end_ } :: os in
    if tl = [] then
      ns  (* fin *)
    else
      start_loop ns tl
  (* When the only thing that matters is the previous node. *)
  and solo_loop state prev = function
    | []                        -> inv_argf "No End at allele sequence: %s" allele
    | Start p :: _              -> inv_argf "Another start %d in %s allele sequence." p allele
    | End end_pos :: tl         -> add_end state end_pos prev tl
    | Boundary { idx; pos} :: t -> let boundary_node = GE.V.create (B (pos, idx)) in
                                  add_allele_edge prev boundary_node;
                                  solo_loop state boundary_node t
    | Sequence { start; s} :: t -> let sequence_node = GE.V.create (N (start, s)) in
                                  add_allele_edge prev sequence_node;
                                  solo_loop state sequence_node t
    | Gap _ :: t                -> solo_loop state prev t
  (* When traversing a reference gap. We have to keep check allele elements
    position to check when to join back with the next reference node. *)
  and ref_gap_loop state prev ref_node ref_pos = function
    | []                                -> inv_argf "No End at allele sequence: %s" allele
    | Start p :: _                      -> inv_argf "Another start %d in %s allele sequence." p allele
    | (End end_pos :: tl) as lst        ->
        if end_pos <= ref_pos then
          add_end state end_pos prev tl
        else (* end_pos > ref_pos *)
          let () = add_allele_edge prev ref_node in
          main_loop state ~prev ~next:ref_node lst
    | (Boundary { idx; pos} :: tl) as l ->
        if pos < ref_pos then
          inv_argf "Allele %s has a boundary %d at %d that is in ref gap ending %d."
            allele idx pos ref_pos
        else if pos = ref_pos then
          if ref_node = B (pos, idx) then
            let () = add_allele_edge prev ref_node in
            main_loop state ~prev ~next:ref_node tl
          else
            inv_argf "Allele %s has a boundary %d at %d where ref gap ends %d."
              allele idx pos ref_pos
        else
          let () = add_allele_edge prev ref_node in
          main_loop state ~prev ~next:ref_node l
    | Sequence { start; s} :: tl        ->
        let new_node = GE.V.create (N (start, s)) in
        let () = add_allele_edge prev new_node in
        let close_pos = start + String.length s in
        if close_pos <= ref_pos then
          ref_gap_loop state new_node ref_node ref_pos tl
        else (* close_pos > ref_pos *)
          rejoin_after_split ~prev:ref_node ~next:ref_node close_pos state
            ~new_node tl
    | Gap { start; length} :: tl        ->
        let close_pos = start + length in
        if close_pos <= ref_pos then
          ref_gap_loop state prev ref_node ref_pos tl
        else (* close_pos > ref_pos *)
          rejoin_after_split ~prev:ref_node ~next:ref_node close_pos state
            ~new_node:prev tl
  (* When not in a reference gap *)
  and main_loop state ~prev ~next = function
    | []              -> inv_argf "No End at allele sequence: %s" allele
    | Start p :: _    -> inv_argf "Another start %d in %s allele sequence." p allele
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
        let new_node = GE.V.create (N (start, s)) in
        let close_pos = start + String.length s in
        let open_res = split_in ~prev ~next ~visit:add_allele_edge start in begin
        match open_res with
        | `AfterLast prev         ->
            let () = add_allele_edge prev new_node in
            solo_loop state new_node t
        | `InGap (prev, next, _)
        | `AtNext (prev, next)    ->
            let () = add_allele_edge prev new_node in
            rejoin_after_split ~prev ~next close_pos state ~new_node t
        end
    | Gap _ :: t -> main_loop state ~prev ~next t (* skip Gaps *)
  in
  start_loop [] alt_lst

type construct_which_args =
  | NumberOfAlts of int
  | SpecificAlleles of string list

let construct_which_args_to_string = function
  | NumberOfAlts n    -> sprintf "N%d" n
  | SpecificAlleles l -> sprintf "S%s" (String.concat ~sep:"_" l)


let construct_from_parsed ?which r =
  let open Ref_graph in
  let open Mas_parser in
  let { reference; ref_elems; alt_elems} = r in
  let alt_elems = List.sort ~cmp:(fun (n1, _) (n2, _) -> compare n1 n2) alt_elems in
  let alt_alleles =
    match which with
    | None ->
        alt_elems
    | Some (NumberOfAlts num_alt_to_add) ->
        List.take alt_elems num_alt_to_add 
    | Some (SpecificAlleles alst)        ->
        List.map alst ~f:(fun name -> name, List.assoc name alt_elems) 
  in
  let num_alleles = List.length alt_alleles in
  let ref_length = List.length ref_elems in
  let g = GE.create ~size:(ref_length * num_alleles) () in
  let aset = Allele_set.construct (reference :: List.map ~f:fst alt_alleles) in
  let refs_start_ends = add_reference_elems g aset reference ref_elems in
  let fs_ls_st_assoc = reference_starts_and_ends refs_start_ends in
  let () =
    List.iter alt_alleles ~f:(fun (allele_name, lst) ->
      ignore (add_non_ref g reference aset fs_ls_st_assoc allele_name lst))
  in
  (aset, g)

let construct_from_file ?which file =
  construct_from_parsed ?which (Mas_parser.from_file file)

type construct_args =
  { alignment_file  : string       (* this is really a path, basename below *)
  ; which           : construct_which_args option
  }

let construct_args_to_string { alignment_file; which } =
  sprintf "%s_%s"
    (Filename.basename alignment_file)
    (match which with
      | None -> "All"
      | Some w -> construct_which_args_to_string w)

let cache_dir = ".ref_graphs"

let construct_from_file =
  let dir = Filename.concat (Sys.getcwd ()) cache_dir in
  disk_memoize ~dir construct_args_to_string
    (fun { alignment_file; which } ->
       construct_from_file ?which alignment_file)
