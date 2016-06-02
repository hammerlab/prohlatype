(** Construction from Mas_parser *)

open Graph        (* OCamlgraph *)
open Nonstd
open Ref_graph
module String = Sosa.Native_string

let inv_argf ?(prefix="") fmt = ksprintf invalid_arg ("%s" ^^ fmt) prefix

let add_reference_elems g allele ref_elems =
  let open Mas_parser in
  let open Nodes in
  let add_start start_pos lst =
    let vs = G.V.create (S (start_pos, allele)) in
    `Started (vs, vs) :: lst
  in
  let add_end end_pos ~vs ~vp lst =
    let ve = G.V.create (E end_pos) in
    G.add_edge_e g (G.E.create vp allele ve);
    `Ended (vs, ve) :: lst
  in
  let add_boundary ~vs ~vp ~idx ~pos lst =
    let vb = G.V.create (B (pos, idx)) in
    G.add_edge_e g (G.E.create vp allele vb);
    `Started (vs, vb) :: lst
  in
  let add_seq ~vs ~vp ~start s lst =
    let vn = G.V.create (N (start, s)) in
    G.add_edge_e g (G.E.create vp allele vn);
    `Started (vs, vn) :: lst
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
      | `Started (vs, vp) :: tl, End end_pos          -> add_end end_pos ~vs ~vp tl
      | `Started (vs, vp) :: tl, Boundary {idx; pos } -> add_boundary ~vs ~vp ~idx ~pos tl
      | `Started (vs, vp) :: tl, Sequence {start; s } -> add_seq ~vs ~vp ~start s tl
      | `Started (_, _) :: _,    Gap _                -> state       (* ignore gaps *)
      | `Started (vs, vp) :: tl, Start sp             ->
          inv_argf "Unexpected second start at %d for %s" sp allele)
  |> List.map ~f:(function
      | `Started _ -> inv_argf "Still have a Started in %s ref" allele
      | `Ended (vs, vp) -> vs, vp)
  |> List.sort ~cmp:(fun (s1,_) (s2,_) -> compare s1 s2)

module SSet = Set.Make (struct
    type t = string
    let compare = compare
  end)

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

let edges_between_to_string g pv nv =
  G.find_all_edges g pv nv
  |> List.map ~f:(fun (_,l,_) -> l)
  |> String.concat ~sep:","

(*first_start, last_end, end_to_next_start_assoc *)
let reference_starts_and_ends = function
  | []           -> inv_argf "Reference has no start and ends"
  | (s, e) :: [] -> s, e, []
  | (s, e) :: t  ->
    let rec loop ep acc = function
      | []           -> inv_argf "stop before empty"
      | (s, e) :: [] -> e, (ep, s) :: acc
      | (s, e) :: t  -> loop e ((ep, s) :: acc) t
    in
    let e, l = loop e [] t in
    s, e, l

let add_non_ref g reference (first_start, last_end, end_to_next_start_assoc) allele alt_lst =
  let open Mas_parser in
  let open Nodes in
  let next_reference ~msg from =
    match next_node_along reference g ~from with
    | Some n -> n
    | None   -> try List.assoc from end_to_next_start_assoc
                with Not_found -> invalid_arg msg
  in
  let do_nothing _ _ = () in
  let advance_until ~visit ~prev ~next pos =
    let rec forward node msg =
      loop node (next_reference ~msg node)
    and loop prev next =
      if next = first_start then
        forward next "No next after reference Start!"
      else if next = last_end then
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
      (* TODO: this set intersect is very clumsy *)
      let pr_edge_s =
        List.fold_left pr ~init:SSet.empty ~f:(fun set (p, l, _) ->
          G.add_edge_e g (G.E.create p l v1);
            SSet.add l set)
      in
      let v2 = N (pos, sn) in
      G.add_vertex g v2;
      let su_edge_s =
        List.fold_left su ~init:SSet.empty ~f:(fun set (_, l, s) ->
          G.add_edge_e g (G.E.create v2 l s);
            SSet.add l set)
      in
      let s_inter = SSet.inter pr_edge_s su_edge_s in
      SSet.iter (fun l -> G.add_edge_e g (G.E.create v1 l v2)) s_inter;
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
  let add_allele_edge pv nv = G.add_edge_e g (G.E.create pv allele nv) in
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
     Loop through the Alignment elemends:
      - discarding Boundary
      - complaining on anything other than a Start
      - on a Start find the position in reference and start correct loop *)
  and start_loop prev_start_ends lst =
    let rec find_start_loop = function
      | Boundary _ :: t 
      | Gap _ :: t      -> find_start_loop t (* Ignore Gaps & Boundaries before Start *)
      | Start p :: t    -> let vs = G.V.create (S (p, allele)) in
                           Some (vs, p, t)
      | []              ->
        begin
          match prev_start_ends with
          | [] -> inv_argf "Failed to find start for %s." allele
          | ls -> None
        end
      | s :: _          -> inv_argf "Encountered %s in %s instead of Start"
                            (al_el_to_string s) allele
    in
    match find_start_loop lst with
    | None -> prev_start_ends (* fin *)
    | Some (new_node, start_pos, tl) ->
        let state = new_node, prev_start_ends in
        rejoin_after_split ~prev:first_start ~next:first_start start_pos state
          ~new_node tl
  (* When the only thing that matters is the previous node. *)
  and solo_loop ((ps, os) as state) vp = function
    | []                        -> inv_argf "No End at allele sequence: %s" allele
    | Start p :: _              -> inv_argf "Another start %d in %s allele sequence." p allele
    | End pos :: []             -> let ve = G.V.create (E pos) in
                                   add_allele_edge vp ve;
                                   (ve,ps) :: os (* fin *)
    | End pos :: lst            -> let ve = G.V.create (E pos) in
                                   add_allele_edge vp ve;
                                   start_loop ((ps, ve) :: os) lst
    | Boundary { idx; pos} :: t -> let vb = G.V.create (B (pos, idx)) in
                                   add_allele_edge vp vb;
                                   solo_loop state vb t
    | Sequence { start; s} :: t -> let vs = G.V.create (N (start, s)) in
                                   add_allele_edge vp vs;
                                   solo_loop state vs t
    | Gap _ :: t                -> solo_loop state vp t
  (* When traversing a reference gap. We have to keep check allele elements
     position to check when to join back with the next reference node. *)
  and ref_gap_loop ((vs, os) as state) prev ref_node ref_pos = function
    | []                                -> inv_argf "No End at allele sequence: %s" allele
    | Start p :: _                      -> inv_argf "Another start %d in %s allele sequence." p allele
    | (End pos :: tl) as l              ->
        if pos <= ref_pos then
          let ve = G.V.create (E pos) in
          let () = add_allele_edge prev ve in
          let ns = (vs,ve) :: os in
          if tl = [] then
            ns (* fin *)
          else
            start_loop ns tl
        else (* pos > ref_pos *)
          let () = add_allele_edge prev ref_node in
          main_loop state ~prev ~next:ref_node l
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
        let new_node = G.V.create (N (start, s)) in
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
  and main_loop ((vs, os) as state) ~prev ~next = function
    | []              -> inv_argf "No End at allele sequence: %s" allele
    | Start p :: _    -> inv_argf "Another start %d in %s allele sequence." p allele
    | End pos :: tl   ->
        let ve = G.V.create (E pos) in
        let pv =
          match split_in ~prev ~next ~visit:add_allele_edge pos with
          | `AfterLast pv
          | `AtNext (pv, _)
          | `InGap (pv, _, _) -> pv
        in
        add_allele_edge pv ve;
        let ns = (vs,ve) :: os in
        if tl = [] then
          ns (* fin *)
        else
          start_loop ns tl
    | Boundary { idx; pos} :: t ->
        let prev, next =
          advance_until_boundary ~prev ~next ~visit:add_allele_edge
            pos idx
        in
        let () = add_allele_edge prev next in
        main_loop state ~prev ~next t
    | Sequence { start; s} :: t ->
        let new_node = G.V.create (N (start, s)) in
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

let construct ?(num_alt_to_add=max_int) allele_lst r =
  let open Mas_parser in
  let { reference; ref_elems; alt_elems} = r in
  let ref_length = List.length ref_elems in
  let num_alleles = min (List.length alt_elems) num_alt_to_add in
  let g = G.create ~size:(ref_length * num_alleles) () in
  let refs_start_ends = add_reference_elems g reference ref_elems in
  let alt_elems =
    List.sort ~cmp:(fun (n1, _) (n2, _) -> compare n1 n2) alt_elems
  in
  let fs_ls_st_assoc = reference_starts_and_ends refs_start_ends in
  let () =
    if allele_lst = [] then
      List.iteri alt_elems ~f:(fun i (allele_name, lst) ->
        if i < num_alt_to_add then
          ignore (add_non_ref g reference fs_ls_st_assoc allele_name lst))
    else
      List.iter allele_lst ~f:(fun allele_name ->
        let lst = List.assoc allele_name alt_elems in
        ignore (add_non_ref g reference fs_ls_st_assoc allele_name lst))
  in
  g
