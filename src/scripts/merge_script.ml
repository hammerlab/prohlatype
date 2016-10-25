
open Util
open Common

type boundary_marker =
  { index     : int
  ; position  : int
  }

let bm ~index ~position = { index; position }
let before_start_boundary = bm ~index:(-1) ~position:(-1)
let before_start bm = bm.index = -1

let to_boundary { index; position } =
  Mas_parser.Boundary { idx = index; pos = position }

let matches_boundary { index; position } = function
  | Mas_parser.Boundary b when b.idx = index &&  b.pos = position -> true
  | _ -> false

let bm_to_string { index; position } =
  sprintf "{index: %d; position: %d}" index position

let bm_to_boundary_string { index; position } =
  sprintf "{idx: %d; pos: %d}" index position

let bounded lst =
  let open Mas_parser in
  let rec loop curb curs acc = function
    | []               -> List.rev ((curb, curs) :: acc)
    | Boundary bp :: t -> let newb = bm ~index:bp.idx ~position:bp.pos in
                          loop newb  ""            ((curb, curs) :: acc) t
    | Start s :: t     -> let newb = if before_start curb then { curb with position = s - 1 } else curb in
                          loop newb  curs          acc                   t
    | End _ :: t
    | Gap _ :: t       -> loop curb  curs          acc                   t
    | Sequence s :: t  -> loop curb  (curs ^ s.s)  acc                   t
  in
  loop before_start_boundary "" [] lst

type instruction =
  | FillFromGenetic of
    { genetic  : boundary_marker
    ; sequence : string          (* Operationally we don't need this,
                                    but it is helpful for debugging. *)
    }
  | MergeFromNuclear of
    { nuclear   : boundary_marker
    ; genetic   : boundary_marker
    ; sequence  : string
    }

let instruction_to_string = function
  | FillFromGenetic { genetic; sequence } ->
      sprintf "FillFromGenetic { genentic: %s; sequence: %s}"
        (bm_to_string genetic) (short_seq sequence)
  | MergeFromNuclear { nuclear; genetic; sequence } ->
      sprintf "MergeFromNuclear { nuclear: %s; genetic: %s; sequence: %s }"
        (bm_to_string nuclear) (bm_to_string genetic) (short_seq sequence)

let bounded_results_to_string lst =
  String.concat ~sep:";"
    (List.map lst ~f:(fun (bm, seq) ->
      sprintf "(%s, %s)" (bm_to_string bm) (short_seq seq)))

let zip_align ~gen ~nuc =
  let rec loop acc g n =
    match n with
    | []        ->
        let reversed = List.rev acc in
        let rest     = List.map g ~f:(fun (genetic, sequence) ->
                          FillFromGenetic { genetic; sequence})
        in
        Ok (reversed @ rest)
    | (nuclear, nuc_sequence) :: nt  ->
        begin
          match g with
          | []       -> error "Reached end of genetic sequence before nuclear: %s"
                          (bounded_results_to_string n)
          | (genetic, sequence) :: gt ->
              if sequence = nuc_sequence then
                let instr = MergeFromNuclear { nuclear; genetic; sequence } in
                loop (instr :: acc) gt nt
              else
                let instr = FillFromGenetic { genetic; sequence } in
                loop (instr :: acc) gt n
        end
  in
  loop [] gen nuc

let prefix_from_f s =
  match String.split s ~on:(`Character '_') with
  | [p; _] -> p
  | _      -> invalid_argf "Missing '_' in %s" s

let setup_test (genf, nucf) =
  let prefix = prefix_from_f genf in
  (*assert ((prefix_from_f nucf) = prefix); *)
  let gen_mp = Mas_parser.from_file (to_alignment_file genf) in
  let nuc_mp = Mas_parser.from_file (to_alignment_file nucf) in
  let gen_bp = bounded gen_mp.Mas_parser.ref_elems in
  let nuc_bp = bounded nuc_mp.Mas_parser.ref_elems in
  prefix, gen_bp, nuc_bp

  (*
let test () =
  List.map paired_files ~f:(fun p ->
    let prefix, gen_bp, nuc_bp = setup_test p in
    prefix, zip_align ~gen:gen_bp ~nuc:nuc_bp)
    *)

let ok_prefix =
  [ "A"; "B"; "C"; "DMA"; "DMB"; "DOA"; "DOB"; "DPA1"; "DPB2"; "DQA1"
  ; "DRA" ; "F"; "G"; "HFE"; "H"; "J"; "K"; "L"; "MICA"; "TAP1"; "V"; "Y"
  ]

let err_prefix =
  [ "DRB", (* Comparing DRB1_gen vs DRB_nuc. 2nd exon. *)
  "CACGTTTCTTGTGGCAGCTTAAGTTTGAATGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTTGCTGGAAAGATGCATCTATAACCAAGAGGAGTCCGTGCGCTTCGACAGCGACGTGGGGGAGTACCGGGCGGTGACGGAGCTGGGGCGGCCTGATGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAGGCGGGCCGCGGTGGACACCTACTGCAGACACAACTACGGGGTTGGTGAGAGCTTCACAGTGCAGCGGCGAG
 4                                                                                                                                                                                                                             |                   |  ||
   CACGTTTCTTGTGGCAGCTTAAGTTTGAATGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTTGCTGGAAAGATGCATCTATAACCAAGAGGAGTCCGTGCGCTTCGACAGCGACGTGGGGGAGTACCGGGCGGTGACGGAGCTGGGGCGGCCTGATGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAGGCGGGCCGCGGTGGACACCTATTGCAGACACAACTACGGGGCTGTGGAGAGCTTCACAGTGCAGCGGCGAG"
  ; "DPB1", (* Mismatch in the 2nd exon *)
    "AGAATTACGTGTACCAGGGACGGCAGGAATGCTACGCGTTTAATGGGACACAGCGCTTCCTGGAGAGATACATCTACAACCGGGAGGAGTACGCGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGCTGCGGAGTACTGGAACAGCCAGAAGGACATCCTGGAGGAGAAGCGGGCAGTGCCGGACAGGGTATGCAGACACAACTACGAGCTGGACGAGGCCGTGACCCTGCAGCGCCGAG
  14         | | |                                                                             |  |                                                        |  |                                     |                    | |                      |  | |  |
     AGAATTACCTTTTCCAGGGACGGCAGGAATGCTACGCGTTTAATGGGACACAGCGCTTCCTGGAGAGATACATCTACAACCGGGAGGAGTTCGTGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGATGAGGAGTACTGGAACAGCCAGAAGGACATCCTGGAGGAGGAGCGGGCAGTGCCGGACAGGATGTGCAGACACAACTACGAGCTGGGCGGGCCCATGACCCTGCAGCGCCGAG"
  ; "DQB1", (* Mismatch in the 1st exon, didn't check the rest.*)
    "ATGTCTTGGAAGAAGTCTTTGCGGATCCCCGGAGACCTTCGGGTAGCAACTGTCACCTTGATGCTGGCGATCCTGAGCTCCTCACTGGCTGAGGGCAGAGACTCTCCCG
  10            |   |                  |        |         |            |    |      |  |  |
     ATGTCTTGGAAAAAGGCTTTGCGGATCCCCGGAGGCCTTCGGGCAGCAACTGTGACCTTGATGCTGTCGATGCTGAGCACCCCAGTGGCTGAGGGCAGAGACTCTCCCG"
  ; "E", (* 6th, last, exon. *)
    "GGAGCGACAGTGCCCAGGGGTCTGAGTCTCACAGCTTGTAA
   0
     GGAGCGACAGTGCCCAGGGGTCTGAGTCTCACAGCTTGTAAAG"
  ; "MICB", (* 2nd exon different. *)
   "AGCCCCACAGTCTTCGTTACAACCTCATGGTGCTGTCCCAGGATGAATCTGTGCAGTCAGGGTTTCTCGCTGAGGGACATCTGGATGGTCAGCCCTTCCTGCGCTATGACAGGCAGAAACGCAGGGCAAAGCCCCAGGGACAGTGGGCAGAAGATGTCCTGGGAGCTAAGACCTGGGACACAGAGACCGAGGACTTGACAGAGAATGGGCAAGACCTCAGGAGGACCCTGACTCATATCAAGGACCAGAAAGGAG
  2                                              |                                                                                                                         |
    AGCCCCACAGTCTTCGTTACAACCTCATGGTGCTGTCCCAGGATGGATCTGTGCAGTCAGGGTTTCTCGCTGAGGGACATCTGGATGGTCAGCCCTTCCTGCGCTATGACAGGCAGAAACGCAGGGCAAAGCCCCAGGGACAGTGGGCAGAAGATGTCCTGGGAGCTGAGACCTGGGACACAGAGACCGAGGACTTGACAGAGAATGGGCAAGACCTCAGGAGGACCCTGACTCATATCAAGGACCAGAAAGGAG"
  ; "TAP2", (* 1st exon, didn't check the rest: the C vs T 5th from end. *)
"ATGCGGCTCCCTGACCTGAGACCCTGGACCTCCCTGCTGCTGGTGGACGCGGCTTTACTGTGGCTGCTTCAGGGCCCTCTGGGGACTTTGCTTCCTCAAGGGCTGCCAGGACTATGGCTGGAGGGGACCCTGCGGCTGGGAGGGCTGTGGGGGCTGCTAAAGCTAAGAGGGCTGCTGGGATTTGTGGGGACACTGCTGCTCCCGCTCTGTCTGGCCACCCCCCTGACTGTCTCCCTGAGAGCCCTGGTCGCGGGGGCCTCACGTGCTCCCCCAGCCAGAGTCGCTTCAGCCCCTTGGAGCTGGCTGCTGGTGGGGTACGGGGCTGCGGGGCTCAGCTGGTCACTGTGGGCTGTTCTGAGCCCTCCTGGAGCCCAGGAGAAGGAGCAGGACCAGGTGAACAACAAAGTCTTGATGTGGAGGCTGCTGAAGCTCTCCAGGCCGGACCTGCCTCTCCTCGTTGCCGCCTTCTTCTTCCTTGTCCTTGCTGTCTTGG
1                                                                                                                            |
ATGCGGCTCCCTGACCTGAGACCCTGGACCTCCCTGCTGCTGGTGGACGCGGCTTTACTGTGGCTGCTTCAGGGCCCTCTGGGGACTTTGCTTCCTCAAGGGCTGCCAGGACTATGGCTGGAGGGGACCCTGCGGCTGGGAGGGCTGTGGGGGCTGCTAAAGCTAAGAGGGCTGCTGGGATTTGTGGGGACACTGCTGCTCCCGCTCTGTCTGGCCACCCCCCTGACTGTCTCCCTGAGAGCCCTGGTCGCGGGGGCCTCACGTGCTCCCCCAGCCAGAGTCGCTTCAGCCCCTTGGAGCTGGCTGCTGGTGGGGTACGGGGCTGCGGGGCTCAGCTGGTCACTGTGGGCTGTTCTGAGCCCTCCTGGAGCCCAGGAGAAGGAGCAGGACCAGGTGAACAACAAAGTCTTGATGTGGAGGCTGCTGAAGCTCTCCAGGCCGGACCTGCCTCTCCTCGTTGCCGCCTTCTTCTTCCTTGTCCTTGCTGTTTTGG"
  ; "P", (* No nuc *) ""
  ]

let align_from_prefix p =
  let gen_mp = Mas_parser.from_file (to_alignment_file (p ^ "_gen")) in
  let nuc_mp = Mas_parser.from_file (to_alignment_file (p ^ "_nuc")) in
  let gen = bounded gen_mp.Mas_parser.ref_elems in
  let nuc = bounded nuc_mp.Mas_parser.ref_elems in
  gen_mp, nuc_mp, zip_align ~gen ~nuc

let sequence_start_check ~bm acc lst k =
  match bm with
  | `Nothing                               -> k acc lst
  | `DontCheckBut kb                       -> k (kb acc) lst
  | `CheckAnd (sbm, kb) ->
      match lst with
      | b :: t when matches_boundary sbm b -> k (kb b acc) t
      | []      -> invalid_argf "Did not start with Boundary %s but empty."
                        (bm_to_boundary_string sbm)
      | v :: _t -> invalid_argf "Did not start with Boundary %s but with %s."
                        (bm_to_boundary_string sbm)
                        (Mas_parser.al_el_to_string v)

type state =
  { seen_start  : bool
  ; seen_end    : bool
  ; next_pos    : Mas_parser.position
  ; acc         : string Mas_parser.alignment_element list
  }

let start_state pos =
  { seen_start  = false
  ; seen_end    = false
  ; next_pos    = pos
  ; acc         = []
  }

let update s =
  let open Mas_parser in function
  | Boundary b -> { s with acc = Boundary { b with pos = s.next_pos} :: s.acc
                         ; next_pos = s.next_pos + 1 }
  | Start _    -> if s.seen_start then s else
                    { s with acc = Start s.next_pos :: s.acc
                           ; seen_start = true}
  | End e      -> if s.seen_end then s else
                    { s with acc = End s.next_pos :: s.acc
                           ; seen_end = true }
  | Gap g      -> { s with acc = Gap { g with start = s.next_pos } :: s.acc
                         ; next_pos = s.next_pos + g.length }
  | Sequence t -> { s with acc = Sequence { t with start = s.next_pos} :: s.acc
                         ; next_pos = s.next_pos + (String.length t.s) }

let accumulate_up_to_next_boundary_rev ~bm s lst =
  let open Mas_parser in
  let rec loop s = function
    | []
    | (Boundary _ :: _) as l  -> s, l
    | (Start _ as e) :: tl
    | (End _ as e) :: tl
    | (Sequence _ as e) :: tl
    | (Gap _ as e) :: tl      -> loop (update s e) tl
  in
  sequence_start_check ~bm s lst loop

let drop_until_boundary ~bm lst =
  let open Mas_parser in
  let rec loop () = function
    | []                      -> []
    | (Start _) :: t
    | (End _ ) :: t
    | (Sequence _) :: t
    | (Gap _) :: t            -> loop () t
    | (Boundary _ :: t) as l  -> l
  in
  sequence_start_check ~bm () lst loop

let first_position =
  let open Mas_parser in function
    | Start s     -> Ok s
    | Gap g       -> Ok g.start
    | Boundary b  -> Ok b.pos
    | v           -> error "Can't start with %s" (al_el_to_string v)

let alignment_elements_to_string lst =
  String.concat ~sep:";" (List.map ~f:Mas_parser.al_el_to_string lst)

(* Nuc goes into gen!
  - What about if the underlying merging sequences do not have the sequences
    across the boundaries assumed by the instructions generated against the
    reference? *)
let apply_align_instr ~gen ~nuc start_pos alg =
  let rec start gen nuc = function
    (* Assume that the instructions start with a fill command from genetic sequences *)
    | FillFromGenetic f :: tl when before_start f.genetic ->
        let state, gen =
          accumulate_up_to_next_boundary_rev (start_state start_pos) gen
            ~bm:`Nothing
        in
        loop state ~gen ~nuc tl
    | []     -> error "Empty instructions"
    | s :: _ -> error "Did not start with a starting FillFromGenetic but %s"
                  (instruction_to_string s)
  and loop state ~gen ~nuc = function
    | MergeFromNuclear m :: tl ->
        let bm =
          let genetic_boundary = to_boundary m.genetic in
          if before_start m.nuclear then
            `DontCheckBut (fun s -> update s genetic_boundary)
          else
            `CheckAnd (m.nuclear, (fun _b s -> update s genetic_boundary))
        in
        let nstate, ntl = accumulate_up_to_next_boundary_rev state nuc ~bm in
        (* Modulo new gaps due to different sequences, these should be the same!
           Instead of dropping we should check that they are! *)
        let gen = drop_until_boundary ~bm:(`CheckAnd (m.genetic, (fun _ a -> a))) gen in
        loop nstate ~gen ~nuc:ntl tl
    | FillFromGenetic gs :: tl ->
        let nstate, gen =
            accumulate_up_to_next_boundary_rev state gen
              ~bm:(`CheckAnd (gs.genetic, (fun b s -> update s b)))
        in
        loop nstate ~gen ~nuc tl
    | [] ->
        if gen <> [] then
          error "After all instructions genetic not empty: %s."
            (alignment_elements_to_string gen)
        else if nuc <> [] then
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

let new_alts ~gen ~nuc start_pos instr trie rmp =
  let open Nomenclature in
  list_fold_ok nuc ~init:[] ~f:(fun acc (alt_name, nuc) ->
    (* Check that the gene is always the same? *)
    parse alt_name >>= fun (gene, allele_resolution) ->
      let closest_allele_res = Trie.nearest allele_resolution trie in
      let gen = RMap.find closest_allele_res rmp in
      let closest_allele_str = resolution_to_string ~gene closest_allele_res in
      apply_align_instr ~gen ~nuc start_pos instr >>= fun new_alts ->
        let seq = sprintf "%s -> %s" alt_name closest_allele_str in
        reference_positions_align ~seq new_alts >>= fun _alt_seq_pos_check ->
          Ok ((alt_name, new_alts) :: acc))

(* TODO. We need to make a decision about suffixed (ex. with 'N') alleles
   since they may be the closest. *)

let merge_and_check prefix =
  let open Mas_parser in
  let gen_mp, nuc_mp, instr_ok = align_from_prefix prefix in
  instr_ok >>= fun instr ->
    if gen_mp.reference <> nuc_mp.reference then
      error "References don't match %s vs %s" gen_mp.reference nuc_mp.reference
    else if gen_mp.align_date <> nuc_mp.align_date then
      error "Align dates don't match %s vs %s" gen_mp.align_date nuc_mp.align_date
    else
      match List.hd gen_mp.ref_elems with
      | None             -> error "empty reference elems!"
      | Some hd_ref_elem ->
        first_position hd_ref_elem >>= fun start_pos ->
          apply_align_instr ~gen:gen_mp.ref_elems ~nuc:nuc_mp.ref_elems start_pos instr >>= fun new_ref_elems ->
            reference_positions_align ~seq:"gen_mp.reference" new_ref_elems >>= fun _ref_position_check ->
              let seq_assoc = (gen_mp.reference, gen_mp.ref_elems) :: gen_mp.alt_elems in
              init_trie_and_map seq_assoc >>= fun (gtrie, rmap) ->
                new_alts gen_mp.alt_elems nuc_mp.alt_elems start_pos instr gtrie rmap >>= fun alt_lst ->
                  Ok  { align_date = gen_mp.align_date
                      ; reference  = gen_mp.reference
                      ; ref_elems  = new_ref_elems
                      ; alt_elems  = alt_lst
                      }

let merge_test () =
  ok_prefix
  |> list_fold_ok ~init:() ~f:(fun () p ->
      Printf.printf "testing: %s\n" p;
      merge_and_check p >>= fun _ -> Ok ())
