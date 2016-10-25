
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
    ; offset    : int
    }

(* Nuc goes into gen! *)
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
          | []       -> Error n
          | (genetic, sequence) :: gt ->
              if sequence = nuc_sequence then
                let offset = genetic.position - nuclear.position in
                let instr = MergeFromNuclear { nuclear; genetic; sequence; offset } in
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
    | Start s     -> s
    | Gap g       -> g.start
    | Boundary b  -> b.pos
    | v           -> invalid_argf "Can't start with %s" (al_el_to_string v)

let apply_align_instr ~gen ~nuc start_pos alg =
  let rec start gen nuc = function
    (* Assume that the instructions start with a fill command from genetic sequences *)
    | FillFromGenetic f :: tl when before_start f.genetic ->
        let state, gen =
          accumulate_up_to_next_boundary_rev (start_state start_pos) gen
            ~bm:`Nothing
        in
        loop state ~gen ~nuc tl
    | _ -> assert false
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
        (*
        if gen <> [] then
          invalid_argf "Gen not empty: %s."
            (String.concat ~sep:";" (List.map ~f:Mas_parser.al_el_to_string gen));
        if nuc <> [] then
          invalid_argf "Nuc not empty: %s."
            (String.concat ~sep:";" (List.map ~f:Mas_parser.al_el_to_string nuc)); *)
        List.rev state.acc, gen, nuc
  in
  start gen nuc alg

let reference_positions_align lst =
  let open Mas_parser in
  let rec loop p = function
    | []                      -> Ok p
    | (Boundary b as v) :: t  -> if b.pos = p then loop (p + 1) t else Error (al_el_to_string v, p)
    | (Start s as v) :: t     -> if s = p then loop p t else Error (al_el_to_string v, p)
    | (End e as v) :: t       -> if e = p then loop p t else Error (al_el_to_string v, p)
    | (Gap g as v) :: t       -> if g.start = p then loop (p + g.length) t else Error (al_el_to_string v, p)
    | (Sequence s as v) :: t  -> if s.start = p then loop (p + (String.length s.s)) t else Error (al_el_to_string v, p)
  and start = function
    | (Boundary b) :: t  -> loop b.pos t
    | (Start s) :: t     -> loop s t
    | (Gap g) :: t       -> loop g.start t
    | x :: _ -> invalid_argf "Can't start with %s" (al_el_to_string x)
    | []     -> invalid_argf "Empty"
  in
  start lst

