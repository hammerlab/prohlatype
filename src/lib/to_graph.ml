
open Graph        (* OCamlgraph *)
open Nonstd
open Ref_graph
module String = Sosa.Native_string

let invalid_argf fmt = ksprintf invalid_arg fmt

(** Construction from Mas_parser *)
let seq_elem_to_vertex se =
  let open Mas_parser in
  match se with
  | Start (_, a)    -> G.V.create (Sequences.S a)
  | End _           -> G.V.create Sequences.E
  | Boundary (n, _) -> G.V.create (Sequences.B n)
  | Nuc (p, s)      -> G.V.create (Sequences.N (p, s))
  | Gap _           -> invalid_argf "Not a valid vertex"

let add_reference allele g pv se =
  let open Mas_parser in
  match se with
  | Start _ -> invalid_arg "Should add Start independently!"
  | End _ | Boundary _ | Nuc _ ->
    let v = seq_elem_to_vertex se in
    G.add_vertex g v;
    G.add_edge_e g (G.E.create pv allele v);
    v
  | Gap _ -> pv (* Ignore gaps, since no nodes to create. *)

let add_reference_elems g lst =
  let open Mas_parser in
  match lst with
  | (Start (_, reference) as s) :: t ->
      let v0 = seq_elem_to_vertex s in
      G.add_vertex g v0;
      List.fold_left ~f:(add_reference reference g) ~init:v0 t |> ignore
  | _ :: _ -> invalid_arg "Reference did not start with Start."
  | []     -> invalid_arg "Empty reference seqeunce."

module GapSet = Set.Make (struct
  type t = int * int
  let compare = compare
end)

module SSet = Set.Make (struct
    type t = string
    let compare = compare
  end)
let add_non_ref reference g ref_gaps_set alt_lst =
  let module Ms = Mas_parser in
  let open Sequences in
  let next_reference from =
    next_node_along reference g ~from
    |> Option.value_exn ~msg:(sprintf "Reached end of reference before %s." reference)
  in
  let end_pos n s = n + String.length s in    (* Another function for GADT *)
  let do_nothing _ _ = None in
  let advance_until ?prev_node ?(s=do_nothing) ?(b=do_nothing) ?(n=do_nothing) ~e
      ~pos start_node =
    let rec loop prev_node ref_node =
      let advance f =
        let next_node = next_reference ref_node in
        match f prev_node ref_node with
        | None   -> loop ref_node next_node
        | Some r -> r
      in
      match ref_node with
      | S _      -> advance s
      | E        -> e prev_node
      | B _      -> advance b
      | N (p, s) ->
          if pos = p then                 (* found it exactly no need to split! *)
            `Exact (prev_node, ref_node)
          else if pos < end_pos p s then  (* in here -> split? *)
            `In (prev_node, ref_node, p, s)
          else                            (* keep going *)
            advance n
    in
    loop (Option.value prev_node ~default:start_node) start_node
  in
  let split_reference_n_at ?prev_node ?s ?b ?n ~e ~pos ~join_split start_node =
    match advance_until ?prev_node ?s ?b ?n ~e ~pos start_node with
    | `Exact (pv, nv)     -> (pv, nv)
    | `In (pv, nv, p, s) ->
        join_split pv nv;
        let open Sequences in
        let index = pos - p in
        let fs, sn = String.split_at s ~index in
        let pr = G.pred_e g nv in
        let su = G.succ_e g nv in
        G.remove_vertex g nv;
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
  let use_boundary_if_it_came_before ~before ~after =
    match before with
    | B _ -> before
    | _ -> after
  in
  let insert_alt_sequence prev_node next_ref_node allele alt_lst =
    let add_allele_edge ?(debug=false) pv nv =
      if debug then
        printf "adding allele edge %s to %s\n" (vertex_name pv) (vertex_name nv);
      G.add_edge_e g (G.E.create pv allele nv);
    in
    let add_allele_edge_and_continue ?debug pv nv =
      add_allele_edge ?debug pv nv;
      None
    in
    let invalid_start action _ _ =
        invalid_argf "Found a start of reference %s in %s." action allele
    in
    let invalid_end action _ =
        invalid_argf "Found end of reference before finishing %s in %s."
          action allele
    in
    let insert_gap ~pv ~nv ~open_pos ~gap_length f action =
      let () =
        printf "inserting gap at %s\t%s\t edges:%s \n"
                (vertex_name pv) (vertex_name nv)
          (G.find_all_edges g pv nv
           |> List.map ~f:(fun (_,l,_) -> l)
           |> String.concat ~sep:",")

      in
      let before_gap_open, gap_start =
        split_reference_n_at ~prev_node:pv ~pos:open_pos nv
          ~s:(invalid_start action)
          ~e:(invalid_end action)
          ~b:add_allele_edge_and_continue
          ~n:(add_allele_edge_and_continue ~debug:true)
          ~join_split:(add_allele_edge)
      in
      let () =
        printf "first part of gap at %s\t%s\t edges:%s \n"
          (vertex_name before_gap_open) (vertex_name gap_start)
          (G.find_all_edges g before_gap_open gap_start
           |> List.map ~f:(fun (_,l,_) -> l)
           |> String.concat ~sep:",")
      in
      let bef_gap_close, gap_close =
        split_reference_n_at
          ~prev_node:before_gap_open gap_start
          ~pos:(open_pos + gap_length)
          ~b:add_allele_edge_and_continue    (* Haven't seen gap's that cross boundaries. *)
          ~n:(fun _ _ -> (* Do nothing since we're NOT following this path *) None)
          ~s:(invalid_start action)
          ~e:(invalid_end action)
          ~join_split:(fun _ _ -> ()) (* do NOT add allele edge *)
      in
      let gap_close = use_boundary_if_it_came_before ~before:bef_gap_close ~after:gap_close in
      f ~before_gap_open ~gap_close
    in
    let insert_seq_node pos seq ~before_gap_open ~gap_close =
      printf "insert_seq_node at %s\t%s \t -- %d %s \n"
        (vertex_name before_gap_open) (vertex_name gap_close)
        pos seq;
      let vm = N (pos, seq) in
      G.add_vertex g vm;
      add_allele_edge before_gap_open vm;
      add_allele_edge vm gap_close;
      vm, gap_close
    in
    let adding_advance ~prev_node ~pos nv =
      advance_until ~prev_node ~pos nv
        ~s:(invalid_start "advancing")
        ~e:(invalid_end "advancing")
        ~b:add_allele_edge_and_continue
        ~n:add_allele_edge_and_continue
    in
    let rec loop ~pv ~nv alt_lst =
      let () =
        printf "currently at %s\t%s \t -- trying to add: %s \n"
                (vertex_name pv) (vertex_name nv)
                  (match alt_lst with [] -> "" | h :: _ -> Ms.sequence_element_to_string h)
      in
      let advance_along_reference ?debug l =
        add_allele_edge ?debug pv nv;
        loop ~pv:nv ~nv:(next_reference nv) l
      in
      match alt_lst with
      | (Ms.End endp as h) :: []  ->
          let v = seq_elem_to_vertex h in
          let before, _after =
            split_reference_n_at ~prev_node:pv ~pos:endp nv
              ~s:(invalid_start "adding end")
              ~e:(function
                  | (N (pos, seq) as n) when pos + String.length seq = endp -> `Exact (n, n)
                  | _ -> invalid_argf "adding end for %s %d" allele endp)
              ~b:add_allele_edge_and_continue
              ~n:(add_allele_edge_and_continue ~debug:true)
              ~join_split:(fun _ _ -> ())
          in
          G.add_vertex g v;
          add_allele_edge before v          (* the _ONLY_ way to terminate ! *)
      | (Ms.End _) :: _t          -> invalid_argf "Alt sequences %s did not end on an end." allele
      | []                        -> invalid_argf "Alt sequences %s had no end node." allele
      | Ms.Start _ :: _           -> invalid_argf "Multiple starts: %s!" allele
      | (Ms.Boundary _ as b) :: t ->
          let v = seq_elem_to_vertex b in
          if nv = v then
            advance_along_reference ~debug:true t
          else                                        (* Not at the right node along reference! *)
            advance_along_reference ~debug:true alt_lst
      | Ms.Gap (gp, gl) :: t      ->
          (* If the gap is in reference, then the appropriate break already
             exists and we just need to replicate the same reference edge. *)
          if GapSet.mem (gp, gl) ref_gaps_set then
            match adding_advance ~prev_node:pv ~pos:gp nv with
            | `Exact (pv, nv)   -> loop ~pv ~nv t
            | `In (_, _, _, _)  -> invalid_argf "How could you not find %d pos of %d for allale %s" gp gl allele
          else
            let action = "trying to create gap" in
            let pv, nv =
              insert_gap ~pv ~nv ~open_pos:gp ~gap_length:gl
                (fun ~before_gap_open ~gap_close ->
                   add_allele_edge before_gap_open gap_close;
                   before_gap_open, gap_close)
                action
            in
            loop ~pv ~nv t
      | Ms.Nuc (spos, seq) :: t   ->
          (* If there is a gap in the reference then we don't need to split
             anything to find the right position. *)
          if GapSet.mem (spos, String.length seq) ref_gaps_set then
            match adding_advance ~prev_node:pv ~pos:spos nv with
            | `Exact (pv, nv)   ->
                let re = G.E.create pv allele nv in (* since we 'advance' the label until nv *)
                G.remove_edge_e g re;
                let () = printf "Adding exact from %s to %s, edges_before: %s\n"
                    (vertex_name pv) (vertex_name nv)
                    (G.pred_e g nv |> List.map ~f:(fun (_, l, _ ) -> l) |> String.concat ~sep:", ")
                in
                let pv, np = insert_seq_node spos seq ~before_gap_open:pv ~gap_close:nv in
                loop ~pv ~nv t
            | `In (_, _, _, _)  -> invalid_argf "How could you not find %d pos" spos
          else
            let action = "trying to insert sequence" in
            let pv, nv =
              insert_gap ~pv ~nv ~open_pos:spos ~gap_length:(String.length seq)
                (insert_seq_node spos seq) action
            in
            loop ~pv ~nv t
    in
    loop ~pv:prev_node ~nv:next_ref_node alt_lst
  in
  match alt_lst with
  | (Ms.Start (start_pos, allele) as h) :: t ->
    let v = seq_elem_to_vertex h in
    G.add_vertex g v;
    let bef_start_pos_n_ref, start_pos =
      split_reference_n_at ~pos:start_pos (S reference)
        ~e:(fun _pv ->
              invalid_argf "Couldn't find a start for %s %d, these aren't aligned!"
                allele start_pos)
        ~join_split:(fun _ _ -> ())
    in
    (* If we start right before a boundary, we want the start to point at
       the boundary, so look ahead at a boundary and stop early. *)
    let start = use_boundary_if_it_came_before ~before:bef_start_pos_n_ref ~after:start_pos in
    G.add_edge_e g (G.E.create v allele start);
    insert_alt_sequence v start allele t
  | _ ->
    invalid_argf "Non reference sequences doesn't start with a start!"

let construct ?(num_alt_to_add=max_int) r =
  let open Mas_parser in
  let { reference; ref_elems; alt_elems} = r in
  let ref_length = List.length ref_elems in
  let num_alleles = min (List.length alt_elems) num_alt_to_add in
  let g = G.create ~size:(ref_length * num_alleles) () in
  add_reference_elems g ref_elems;
  let ref_gaps_set =
    List.fold_left ref_elems ~init:GapSet.empty
      ~f:(fun s v ->
            match v with
            | Gap (p, l) -> GapSet.add (p,l) s
            | Start _
            | End _
            | Boundary _
            | Nuc _ -> s)
  in
  List.iteri alt_elems
    ~f:(fun i p ->
          flush_all ();
          let (_allele_name,p) = p in
          if i < num_alt_to_add then begin
            printf "Adding ---------------- %d ----------\n" i;
            add_non_ref reference g ref_gaps_set p
          end);
  g
