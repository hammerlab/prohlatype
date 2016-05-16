
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

let add_non_ref reference g ref_gaps_set alt_lst =
  let module Ms = Mas_parser in
  let open Sequences in
  let next_reference from =
    let () = printf "next from %s \n" (vertex_name from) in
    next_node_along reference g ~from
    |> Option.value_exn ~msg:(sprintf "Reached end of reference before %s." reference)
  in
  let end_pos n s = n + String.length s in    (* Another function for GADT *)
  let do_nothing _ _ = None in
  let advance_until ?prev_node ?(s=do_nothing) ?(b=do_nothing) ?(n=do_nothing) ~e ~pos start_node =
    let rec loop prev_node ref_node =
      let () = printf "--------------\n%!" in
      let advance f =
        let next_node = next_reference ref_node in
        match f ref_node next_node with
        | None   -> loop ref_node next_node
        | Some r -> r (*p, n *)
      in
      match ref_node with
      | S _      -> advance s
      | E        -> e ()
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
  let split_reference_n_at ?prev_node ?s ?b ?n ~e ~pos start_node =
    match advance_until ?prev_node ?s ?b ?n ~e ~pos start_node with
    | `Exact (pv, nv)     -> (pv, nv)
    | `In (_pv, nv, p, s) ->
          let open Sequences in
          let index = pos - p in
          let fs, sn = String.split_at s ~index in
          let pr = G.pred_e g nv in
          let su = G.succ_e g nv in
          G.remove_vertex g nv;
          let v1 = N (p, fs) in
          G.add_vertex g v1;
          List.iter ~f:(fun (p, l, _) -> G.add_edge_e g (G.E.create p l v1)) pr;
          let v2 = N (pos, sn) in
          G.add_vertex g v2;
          G.add_edge_e g (G.E.create v1 reference v2);
          List.iter ~f:(fun (_, l, s) -> G.add_edge_e g (G.E.create v2 l s)) su;
          (v1, v2)
  in
  let insert_alt_sequence prev_node next_ref_node allele alt_lst =
    let add_allele_edge pv nv =
      G.add_edge_e g (G.E.create pv allele nv);
    in
    let add_allele_edge_and_continue pv nv =
      add_allele_edge pv nv;
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
      let vbefore, vgap_start =
        split_reference_n_at ~prev_node:pv ~pos:open_pos nv
          ~s:(invalid_start action)
          ~e:(invalid_end action)
          ~b:add_allele_edge_and_continue
          ~n:add_allele_edge_and_continue
      in
      let vgap_prev, vgap_end =
        split_reference_n_at ~prev_node:vbefore ~pos:(open_pos + gap_length) vgap_start
          ~b:add_allele_edge_and_continue    (* Haven't seen gap's that cross boundaries. *)
          ~n:(fun _ _ -> (* Do nothing since we're NOT following this path *) None)
          ~s:(invalid_start action)
          ~e:(invalid_end action)
      in
      f vbefore vgap_end;
      vgap_prev, vgap_end
    in
    let insert_seq_node pos seq vbefore vgap_end =
      let vm = N (pos, seq) in
      G.add_vertex g vm;
      add_allele_edge vbefore vm;
      add_allele_edge vm vgap_end
    in
    let adding_advance ~prev_node ~pos nv =
      advance_until ~prev_node ~pos nv
        ~s:(invalid_start "advancing")
        ~e:(invalid_end "advancing")
        ~b:add_allele_edge_and_continue
        ~n:add_allele_edge_and_continue
    in
    let rec loop ~pv ~nv alt_lst =
      let () = printf "currently at %s to %s \n" (vertex_name pv) (vertex_name nv) in
      let advance_along_reference l =
        add_allele_edge pv nv;
        loop ~pv:nv ~nv:(next_reference nv) l
      in
      match alt_lst with
      | (Ms.End _ as h) :: []     ->
          let v = seq_elem_to_vertex h in
          G.add_vertex g v;
          add_allele_edge pv v          (* the _ONLY_ way to terminate ! *)
      | (Ms.End _) :: _t          -> invalid_argf "Alt sequences %s did not end on an end." allele
      | []                        -> invalid_argf "Alt sequences %s had no end node." allele
      | Ms.Start _ :: _           -> invalid_argf "Multiple starts: %s!" allele
      | (Ms.Boundary _ as b) :: t ->
          let v = seq_elem_to_vertex b in
          let () = printf "searching for %s \n" (vertex_name v) in
          if nv = v then
            advance_along_reference t
          else                                        (* Not at the right node along reference! *)
            advance_along_reference alt_lst
      | Ms.Gap (gp, gl) :: t      ->
          (* If the gap is in reference, then the appropriate break already
             exists and we just need to replicate the same reference edge. *)
          if GapSet.mem (gp, gl) ref_gaps_set then
            match adding_advance ~prev_node:pv ~pos:gp nv with
            | `Exact (pv, nv)   -> loop ~pv ~nv t
            | `In (_, _, _, _)  -> invalid_argf "How could you not find %d pos" gp
          else
            let action = "trying to create gap" in
            let pv, nv =
              insert_gap ~pv ~nv ~open_pos:gp ~gap_length:gl add_allele_edge action
            in
            loop ~pv ~nv t
      | Ms.Nuc (spos, seq) :: t   ->
          (* If there is a gap in the reference then we don't need to split
             anything to find the right position. *)
          if GapSet.mem (spos, String.length seq) ref_gaps_set then
            match adding_advance ~prev_node:pv ~pos:spos nv with
            | `Exact (pv, nv)   ->
                insert_seq_node spos seq pv nv;
                loop ~pv ~nv t
            | `In (_, _, _, _)  ->
              invalid_argf "How could you not find %d pos" spos
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
    let _before_start_in_reference, start =
      split_reference_n_at ~pos:start_pos (S reference)
        (* If we start right before a boundary, we want the start to point at
           the boundary, so look ahead at a boundary and stop early. *)
        ~b:(fun bnd nv ->
              match nv with
                      (* Ok, since we ignore the '_before_start_in_reference *)
              | N (p, _) when p = start_pos -> Some (`Exact (bnd, bnd)) 
              | _                           -> None)
        ~e:(fun () ->
              invalid_argf "Couldn't find a start for %s %d, these aren't aligned!"
                allele start_pos)
    in
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
          if i < num_alt_to_add then
            let () = add_non_ref reference g ref_gaps_set p
      in
    ());
  g
