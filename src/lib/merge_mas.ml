
open Util
(*open Common *)

let supported_genes = [ "A"; "B"; "C"]

let list_map_consecutives f lst =
  let rec loop acc = function
    | []
    | _ :: []     -> List.rev acc
    | a :: b :: t -> loop (f a b :: acc) (b :: t)
  in
  loop [] lst

type boundary_marker =
  { index     : int     (* Which segment? *)
  ; position  : int     (* Position of the boundary marker. *)
  ; length    : int     (* Length to next boundary, 1 + following segment. *)
  }

let bm ~index ~position = { index; position; length = 0 }
let before_start_boundary = bm ~index:(-1) ~position:(-1)
let before_start bm = bm.index = -1

let to_boundary { index; position; _ } =
  Mas_parser.Boundary { idx = index; pos = position }

let matches_boundary { index; position; _ } = function
  | Mas_parser.Boundary b when b.idx = index &&  b.pos = position -> true
  | _ -> false

let bm_to_string { index; position; length } =
  sprintf "{index: %d; position: %d; length %d }" index position length

let bm_to_boundary_string { index; position; _ } =
  sprintf "{idx: %d; pos: %d}" index position

let bounded lst =
  let open Mas_parser in
  let il bm i = { bm with length = bm.length + i } in
  let rec loop curb curs acc = function
    | []               -> List.rev ((il curb 1, curs) :: acc)
    | Boundary bp :: t -> let newb = bm ~index:bp.idx ~position:bp.pos in
                          loop newb                           ""            ((il curb 1, curs) :: acc)  t
    | Start s :: t     -> let newb = if before_start curb then { curb with position = s - 1 } else curb in
                          loop newb                           curs          acc                         t
    | End _ :: t       -> loop curb                           curs          acc                         t
    | Gap g :: t       -> loop (il curb g.length)             curs          acc                         t
    | Sequence s :: t  -> loop (il curb (String.length s.s))  (curs ^ s.s)  acc                         t
  in
  loop before_start_boundary "" [] lst

type instruction =
  | FillFromGenetic of
    { genetic  : boundary_marker
    ; offset   : int
    ; sequence : string          (* Operationally we don't need this,
                                    but it is helpful for debugging. *)
    }
  | MergeFromNuclear of
    { nuclear         : boundary_marker
    ; genetic         : boundary_marker
    ; nuclear_offset  : int   (* Two offsets because the first element is the *)
    ; genetic_offset  : int   (* genetic boundary followed by nuclear elements. *)
    ; sequence        : string
    }

let instruction_to_string = function
  | FillFromGenetic { genetic; offset; sequence } ->
      sprintf "FillFromGenetic { genentic: %s; offset: %d; sequence: %s}"
        (bm_to_string genetic) offset (short_seq sequence)
  | MergeFromNuclear { nuclear; nuclear_offset; genetic; genetic_offset; sequence } ->
      sprintf "MergeFromNuclear { nuclear: %s; nuclear_offset: %d; \
                genetic: %s; genetic_offset: %d; sequence: %s }"
        (bm_to_string nuclear) nuclear_offset
        (bm_to_string genetic) genetic_offset
        (short_seq sequence)

let bounded_results_to_string lst =
  String.concat ~sep:";"
    (List.map lst ~f:(fun (bm, seq) ->
      sprintf "(%s, %s)" (bm_to_string bm) (short_seq seq)))

let zip_align ~gen ~nuc =
  let rec loop next_bnd_pos acc g n =
    match n with
    | []        ->
        let _p, rest =
          List.fold_left g ~init:(next_bnd_pos, acc)
            ~f:(fun (next_bnd_pos, acc) (genetic, sequence) ->
                  let offset = next_bnd_pos - genetic.position in
                  let next_end_pos = genetic.position + genetic.length in
                  let instr = FillFromGenetic { genetic; sequence; offset } in
                  (next_end_pos, instr :: acc))
        in
        Ok (List.rev rest)
    | (nuclear, nuc_sequence) :: nt  ->
        begin
          match g with
          | []       -> error "Reached end of genetic sequence before nuclear: %s"
                          (bounded_results_to_string n)
          | (genetic, sequence) :: gt ->
              if sequence = nuc_sequence then
                let genetic_offset = next_bnd_pos - genetic.position in
                let nuclear_offset = genetic_offset + genetic.position - nuclear.position in
                let next_end_pos = next_bnd_pos + nuclear.length in
                let instr = MergeFromNuclear { nuclear; genetic; sequence; genetic_offset; nuclear_offset } in
                loop next_end_pos (instr :: acc) gt nt
              else
                let offset = next_bnd_pos - genetic.position in
                let next_end_pos = next_bnd_pos + genetic.length in
                let instr = FillFromGenetic { genetic; sequence; offset } in
                loop next_end_pos (instr :: acc) gt n
        end
  in
  match gen with
  | []           -> error "Empty genetic sequence!"
  | (bm, _) :: _ -> let next_bnd_pos = bm.position in
                    loop next_bnd_pos [] gen nuc

let prefix_from_f s =
  match String.split s ~on:(`Character '_') with
  | [p; _] -> p
  | _      -> invalid_argf "Missing '_' in %s" s

let align_from_prefix prefix_path =
  let gen_mp = Mas_parser.from_file (prefix_path ^ "_gen.txt") in
  let nuc_mp = Mas_parser.from_file (prefix_path ^ "_nuc.txt") in
  let gen = bounded gen_mp.Mas_parser.ref_elems in
  let nuc = bounded nuc_mp.Mas_parser.ref_elems in
  gen_mp, nuc_mp, zip_align ~gen ~nuc

type state =
  { seen_start  : bool
  ; seen_end    : bool
  ; acc         : string Mas_parser.alignment_element list
  }

let start_state =
  { seen_start  = false
  ; seen_end    = false
  ; acc         = []
  }

let latest_position =
  let open Mas_parser in function
    | Start s     -> Ok s
    | End e       -> Ok e
    | Gap g       -> Ok (g.start + g.length)
    | Sequence s  -> Ok (s.start + String.length s.s)
    | Boundary b  -> error "Can't ask for the latest position of boundary! %d, %d"
                      b.idx b.pos

let latest_position_al_elems = function
  | []     -> error "No latest alignment element."
  | h :: _ -> latest_position h

let sequence_start_check ~bm st lst k =
  match bm with
  | `Nothing                                -> k st lst
  | `DontCheckBut kb                        -> k (kb st) lst
  | `CheckAnd (sbm, kb) ->
      match lst with
      | b :: t when matches_boundary sbm b  -> k (kb b st) t
      | []      -> error "Did not start with Boundary %s but empty."
                        (bm_to_boundary_string sbm)
      | v :: _t -> error "Did not start with Boundary %s but with %s."
                        (bm_to_boundary_string sbm)
                        (Mas_parser.al_el_to_string v)

let list_remove_first p lst =
  let rec loop acc = function
    | []     -> List.rev acc
    | h :: t -> if p h then List.rev acc @ t else loop (h :: acc) t
  in
  loop [] lst

let is_end = function | Mas_parser.End _ -> true | _ -> false

let update offset s =
  let open Mas_parser in function
  | Boundary b -> { s with acc = Boundary { b with pos = b.pos + offset} :: s.acc }
  | Start t    -> if s.seen_start then s else
                    { s with acc = Start (t + offset) :: s.acc
                           ; seen_start = true}
  | End e      -> let without_end = if s.seen_end then list_remove_first is_end s.acc else s.acc in
                  { s with acc = End (e + offset) :: without_end
                         ; seen_end = true }
  | Gap g      -> { s with acc = Gap { g with start = g.start + offset } :: s.acc }
  | Sequence t -> { s with acc = Sequence { t with start = t.start + offset} :: s.acc }

let accumulate_up_to_next_boundary_rev ~bm offset s lst =
  let open Mas_parser in
  let rec loop s = function
    | []
    | (Boundary _ :: _) as l  -> Ok (s, l)
    | (Start _ as e) :: tl
    | (End _ as e) :: tl
    | (Sequence _ as e) :: tl
    | (Gap _ as e) :: tl      -> loop (update offset s e) tl
  in
  sequence_start_check ~bm s lst loop

let drop_until_boundary ~bm lst =
  let open Mas_parser in
  let rec loop () = function
    | []                      -> Ok []
    | (Start _) :: t
    | (End _ ) :: t
    | (Sequence _) :: t
    | (Gap _) :: t            -> loop () t
    | (Boundary _ :: t) as l  -> Ok l
  in
  sequence_start_check ~bm () lst loop

let alignment_elements_to_string lst =
  String.concat ~sep:";\n" (List.map ~f:Mas_parser.al_el_to_string lst)

(* Nuc goes into gen!
  - What about if the underlying merging sequences do not have the sequences
    across the boundaries assumed by the instructions generated against the
    reference? *)
let apply_align_instr ?(fail_on_empty=true)  ~gen ~nuc alg =
  let rec start gen nuc = function
    (* Assume that the instructions start with a fill command from genetic sequences *)
    | FillFromGenetic f :: tl when before_start f.genetic ->
        (*first_position_al_elems gen >>= fun start_pos -> *)
          accumulate_up_to_next_boundary_rev f.offset start_state gen ~bm:`Nothing >>=
            fun (state, gen) -> loop state ~gen ~nuc tl
    | []     -> error "Empty instructions"
    | s :: _ -> error "Did not start with a starting FillFromGenetic but %s"
                  (instruction_to_string s)
  and loop state ~gen ~nuc = function
    | MergeFromNuclear m :: tl ->
        let bm =
          let genetic_boundary = to_boundary m.genetic in
          if before_start m.nuclear then
            `DontCheckBut (fun s -> update m.genetic_offset s genetic_boundary)
          else
            `CheckAnd ( m.nuclear
                      , (fun _b s -> update m.genetic_offset s genetic_boundary))
        in
        accumulate_up_to_next_boundary_rev m.nuclear_offset state nuc ~bm >>=
          fun (nstate, ntl) ->
            (* Modulo new gaps due to different sequences, these should be the same.
               Instead of dropping we should check that they are! *)
            drop_until_boundary ~bm:(`CheckAnd (m.genetic, (fun _ a -> a))) gen >>=
              fun gen -> loop nstate ~gen ~nuc:ntl tl
    | FillFromGenetic f :: tl ->
        accumulate_up_to_next_boundary_rev f.offset state gen
          ~bm:(`CheckAnd (f.genetic, (fun b s -> update f.offset s b))) >>=
          fun (nstate, gen) -> loop nstate ~gen ~nuc tl
    | [] ->
        if fail_on_empty && gen <> [] then
          error "After all instructions genetic not empty: %s."
            (alignment_elements_to_string gen)
        else if fail_on_empty && nuc <> [] then
          error "After all instructions nuclear sequence not empty: %s."
            (alignment_elements_to_string nuc)
        else
          Ok (List.rev state.acc)
  in
  start gen nuc alg

let reference_positions_align ?seq lst =
  let open Mas_parser in
  let emst = Option.value (Option.map seq ~f:(sprintf "For seq %s ")) ~default:"" in
  let nerror n p = error "%spositions are not aligned: %s, pos: %d" emst (al_el_to_string n) p in
  let rec loop p = function
    | []                      -> Ok p
    | (Boundary b as v) :: t  -> if b.pos = p then loop (p + 1) t else nerror v p
    | (Start s as v) :: t     -> if s = p then loop p t else nerror v p
    | (End e as v) :: t       -> if e = p then loop p t else nerror v p
    | (Gap g as v) :: t       -> if g.start = p then loop (p + g.length) t else nerror v p
    | (Sequence s as v) :: t  -> if s.start = p then loop (p + (String.length s.s)) t else nerror v p
  and start = function
    | (Boundary b) :: t  -> loop b.pos t
    | (Start s) :: t     -> loop s t
    | (Gap g) :: t       -> loop g.start t
    | x :: _             -> error "Can't start with %s." (al_el_to_string x)
    | []                 -> error "Empty"
  in
  start lst

let list_fold_ok lst ~f ~init =
  let rec loop acc = function
    | []      -> Ok acc
    | h :: t  -> f acc h >>= fun a -> loop a t
  in
  loop init lst

let init_trie elems =
  let open Nomenclature in
  list_fold_ok elems ~init:Trie.empty ~f:(fun trie s ->
    parse s >>= fun (_gene, allele_resolution) ->
      Ok (Trie.add allele_resolution trie))

module RMap = Map.Make (struct
    type t = Nomenclature.resolution
    let compare = compare
    end)

let init_trie_and_map elems =
  let open Nomenclature in
  list_fold_ok elems ~init:(Trie.empty, RMap.empty)
    ~f:(fun (trie, mp) (s, seq) ->
        parse s >>= fun (_gene, allele_resolution) ->
          Ok ( Trie.add allele_resolution trie
             , RMap.add allele_resolution seq mp))

let new_alts ~gen ~nuc instr trie rmp =
  let open Nomenclature in
  list_fold_ok nuc ~init:[] ~f:(fun acc (alt_name, nuc) ->
    (* Check that the gene is always the same? *)
    parse alt_name >>= fun (gene, allele_resolution) ->
      let closest_allele_res = Trie.nearest allele_resolution trie in
      let gen = RMap.find closest_allele_res rmp in
      let closest_allele_str = resolution_to_string ~gene closest_allele_res in
      (*Printf.printf "%s (nuc) -> %s (gen)" alt_name closest_allele_str; *)
      apply_align_instr ~gen ~nuc instr >>= fun new_alts ->
        Ok ((alt_name, new_alts) :: acc))

(* TODO. We need to make a decision about suffixed (ex. with 'N') alleles
   since they may be the closest. *)

let and_check prefix =
  let open Mas_parser in
  let gen_mp, nuc_mp, instr_ok = align_from_prefix prefix in
  instr_ok >>= fun instr ->
    if gen_mp.reference <> nuc_mp.reference then
      error "References don't match %s vs %s" gen_mp.reference nuc_mp.reference
    else if gen_mp.align_date <> nuc_mp.align_date then
      error "Align dates don't match %s vs %s" gen_mp.align_date nuc_mp.align_date
    else
      apply_align_instr ~gen:gen_mp.ref_elems ~nuc:nuc_mp.ref_elems instr >>= fun new_ref_elems ->
        reference_positions_align ~seq:("reference:" ^ gen_mp.reference) new_ref_elems >>= fun _ref_position_check ->
          let seq_assoc = (gen_mp.reference, gen_mp.ref_elems) :: gen_mp.alt_elems in
          init_trie_and_map seq_assoc >>= fun (gtrie, rmap) ->
                new_alts gen_mp.alt_elems nuc_mp.alt_elems instr gtrie rmap >>= fun alt_lst ->
              Ok  { align_date = gen_mp.align_date
                  ; reference  = gen_mp.reference
                  ; ref_elems  = new_ref_elems
                  ; alt_elems  = alt_lst
                  }
