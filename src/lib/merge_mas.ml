(* Merging of alignment element sequences.

Problem: There are many (10x) more nucleic allele sequences than genetic ones.
But the genetic sequences allow us to match many more reads. We would like to
"impute" sequences for (and hence a graph containing) all of the alleles for
which we have nucleic sequence data. To accomplish this, we use the current
nucleic sequences as scaffolding around which we append information from the
genetic ones.

The algorithm proceeds in several steps. At a high level for a given gene
(ex. HLA-A) we:

  1. Use the reference allele to construct a set of instructions for how to
     zip the exons into the introns.
  2. Afterwards we use those instructions to, merge the reference.
  3. Merge all alleles where we have both nucleic and genetic data.
  4. Finally for the other alleles we use the same instructions to merge the
     them with the closest allele for which we have genetic data.

At a lower level:


TODO: This nuclear vs genetic distinction is stupid it may be redundant
but sometimes we mean exon vs intron we should label the parts accordingly.
*)

open Util
open Mas_parser

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
  ; length    : int     (* Length to next boundary
                        , or 1 + the total length of the following segment. *)
  }

let bm ~index ~position = { index; position; length = 0 }
let before_start_boundary = bm ~index:(-1) ~position:(-1)
let before_start bm = bm.index = -1

let to_boundary ~offset { index; position; _ } =
  Boundary { idx = index; pos = position + offset }

let matches_boundary { index; position; _ } = function
  | Boundary b when b.idx = index &&  b.pos = position -> true
  | _ -> false

let bm_to_string { index; position; length } =
  sprintf "{index: %d; position: %d; length %d }" index position length

let bm_to_boundary_string { index; position; _ } =
  sprintf "{idx: %d; pos: %d}" index position

let bounded lst =
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

type 'a region =
  { bm        : boundary_marker
  ; offset    : int
  ; sequence  : 'a
  }

let region_to_string sequence_to_string { bm; offset; sequence } =
  sprintf "{ bm: %s; offset: %d; sequence: %s}"
    (bm_to_string bm) offset (sequence_to_string sequence)

(* When we want to have nuclear and genetic labels...
   we're going to recurse this record!!!!! *)
type ('n, 'g) ngd =
  { nuclear : 'n
  ; genetic : 'g
  }

let ngd_to_string ns gs { nuclear; genetic }=
  sprintf "{ nuclear: %s; genetic: %s }"
    (ns nuclear) (gs nuclear)

type ('genetic, 'nuclear) instruction =
  | FillFromGenetic of 'genetic region
  | MergeFromNuclear of ('nuclear region, 'genetic region) ngd

let instruction_to_string genetic_seq_to_string nuclear_seq_to_string = function
  | FillFromGenetic r ->
      sprintf "FillFromGenetic %s"
        (region_to_string genetic_seq_to_string r)
  | MergeFromNuclear ngd ->
      sprintf "MergeFromNuclear %s"
        (ngd_to_string (region_to_string nuclear_seq_to_string)
            (region_to_string genetic_seq_to_string) ngd)

let bounded_results_to_string lst =
  String.concat ~sep:";"
    (List.map lst ~f:(fun (bm, seq) ->
      sprintf "(%s, %s)" (bm_to_string bm) (short_seq seq)))

(* Create a framework of instructions based off of, probably the alignment
   in the reference sequence. *)
let zip_align ~gen ~nuc =
  let rec loop next_boundary_pos acc g n =
    match n with
    | []        ->
        let _p, rest =
          List.fold_left g ~init:(next_boundary_pos, acc)
            ~f:(fun (next_boundary_pos, acc) (genetic, sequence) ->
                  let offset = next_boundary_pos - genetic.position in
                  let next_end_pos = genetic.position + genetic.length in
                  let instr = FillFromGenetic { bm = genetic; sequence; offset } in
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
                let genetic  =
                  { bm       = genetic
                  ; offset   = next_boundary_pos - genetic.position
                  ; sequence = sequence
                  }
                in
                let nuclear  =
                  { bm       = nuclear
                  ; offset   = next_boundary_pos (*+ genetic.position*) - nuclear.position
                  ; sequence = nuc_sequence
                  }
                in
                let instr = MergeFromNuclear { genetic ; nuclear } in
                loop (next_boundary_pos + nuclear.bm.length) (instr :: acc) gt nt
              else
                let instr = FillFromGenetic
                  { bm = genetic
                  ; sequence
                  ; offset = next_boundary_pos - genetic.position
                  }
                in
                loop (next_boundary_pos + genetic.length) (instr :: acc) gt n
        end
  in
  match gen with
  | []           -> error "Empty genetic sequence!"
  | (bm, _) :: _ -> let next_boundary_pos = bm.position in
                    loop next_boundary_pos [] gen nuc

let prefix_from_f s =
  match String.split s ~on:(`Character '_') with
  | [p; _] -> p
  | _      -> invalid_argf "Missing '_' in %s" s

let align_from_prefix prefix_path =
  let gen_mp = from_file (prefix_path ^ "_gen.txt") in
  let nuc_mp = from_file (prefix_path ^ "_nuc.txt") in
  let gen = bounded gen_mp.ref_elems in
  let nuc = bounded nuc_mp.ref_elems in
  zip_align ~gen ~nuc >>= fun i -> Ok (gen_mp, nuc_mp, i)

let position = function
  | Start s     -> s
  | End e       -> e
  | Gap g       -> g.start
  | Sequence s  -> s.start
  | Boundary b  -> b.pos

let latest_position = function
  | Start s     -> s
  | End e       -> e
  | Gap g       -> (g.start + g.length)
  | Sequence s  -> (s.start + String.length s.s)
  | Boundary b  -> b.pos

let latest_position_al_elems = function
  | []     -> error "No latest alignment element."
  | h :: _ -> Ok (latest_position h)

let is_sequence = function | Sequence _ -> true | _ -> false
let is_gap      = function | Gap _      -> true | _ -> false
let is_end      = function | End _      -> true | _ -> false
let is_start    = function | Start _    -> true | _ -> false

let shift_al_el offset = function
  | Boundary b -> Boundary { b with pos = b.pos + offset}
  | Start s    -> Start (s + offset)
  | End e      -> End (e + offset)
  | Gap g      -> Gap { g with start = g.start + offset }
  | Sequence s -> Sequence { s with start = s.start + offset}

(*alignment_elements_to_string *)
let al_els_to_string lst =
  if lst = [] then "[]" else
    String.concat ~sep:";\n" (List.map ~f:al_el_to_string lst)

let not_boundary_with_pos pos = function
  | Boundary b when b.pos = pos -> false
  | _ -> true

let list_split_and_map lst ~p ~f =
  let rec loop acc = function
    | h :: t when p h -> loop (f h :: acc) t
    | lst             -> (List.rev acc, lst)
  in
  loop [] lst

let split_at_boundary_and_offset ~offset ~pos =
  list_split_and_map ~p:(not_boundary_with_pos pos) ~f:(shift_al_el offset)

let split_at_end =
  list_split_and_map ~p:(fun e -> not (is_end e)) ~f:id

let split_at_start =
  list_split_and_map ~p:(fun e -> not (is_start e)) ~f:id

(* Droping the boundary from the returned result simplifies the special casing.
   We now know that we (might) have to add (not for the first one!) the
   boundaries from the instructions.
   But we don't have to keep checking for them during the traversals. *)
let split_at_boundary_offset_before_and_drop_it ~offset ~pos lst =
  match split_at_boundary_and_offset ~offset ~pos lst with
  | before, (_bnd :: after) -> before, after
  | before, []              -> before, []

let replace_sequence r s = { r with sequence = s }

(* Special splicing adjustments. *)
(* TODO: Figure out how to impute these from genetic. *)
let a_01_11N r =
  let open Mas_parser in
  if r.bm.index = 1 (*&& r.bm.position = 396 *) then begin
    printf "Adjusting A*01:11N's third exon!\n";
    { r with sequence =
        List.filter r.sequence
          ~f:(function | Gap { start = 1084 } -> false | _ -> true)
        @ [ Sequence { start = 1085; s = "T" }
          ; Gap { start = 1102; length = 17 }
          ; Gap { start = 1122; length = 14 }
          ]
    }
  end else
    r

let drop_splice_failure_from_exon name ~index ~end_pos r =
  if r.bm.index = index then begin
    printf "Adjusting %s's %d exon" name (index + 1);
    { r with sequence =
        List.filter r.sequence ~f:(fun e -> position e < end_pos)
    }
  end else
    r

let known_splice_adjustments = function
  | "A*01:11N"            -> a_01_11N
  (* r.bm.position = 745 *)
  | "A*29:01:01:02N"
  | "A*26:01:01:03N"
  | "A*03:01:01:02N" as n -> drop_splice_failure_from_exon n ~index:2 ~end_pos:2048
  (* r.bm.position = 0 *)
  | "B*15:01:01:02N" as n -> drop_splice_failure_from_exon n ~index:(-1) ~end_pos:75
  (* r.bm.position = 1176 *)
  | "C*04:09N"       as n -> drop_splice_failure_from_exon n ~index:6 ~end_pos:2999
  | _                     -> id

let map_instr_to_alignments ?(fail_on_empty=true) nallele ~gen ~nuc alg =
  let splice_adj = known_splice_adjustments nallele in
  let rec start = function
    | FillFromGenetic f :: tl when before_start f.bm ->
                fill_from_g [] f tl ~gen ~nuc
    | []     -> error "Mapping instructions for %s, empty!" nallele
    | s :: _ -> error "Mapping instructions for %s, did not start with a \
                      \"starting\" FillFromGenetic but %s"
                  nallele (instruction_to_string id id s)
  and loop acc ~gen ~nuc = function
    | MergeFromNuclear {genetic; nuclear} :: tl -> merge_from_n acc genetic nuclear tl ~gen ~nuc
    | FillFromGenetic f :: tl                   -> fill_from_g acc f tl ~gen ~nuc
    | []  ->
        if fail_on_empty && gen <> [] then
          error "Mapping instructions for %s, after all instructions \
            genetic not empty: %s."
            nallele (al_els_to_string gen)
        else if fail_on_empty && nuc <> [] then
          error "Mapping instructions for %s, after all instructions \
            nuclear sequence not empty: %s."
            nallele (al_els_to_string nuc)
        else
          Ok (List.rev acc)                           (* Fin *)
  and fill_from_g acc f tl ~gen ~nuc =
    let pos     = f.bm.position + f.bm.length in
    let before, after = split_at_boundary_offset_before_and_drop_it
      ~pos ~offset:f.offset gen
    in
    let nacc = FillFromGenetic (replace_sequence f before) :: acc in
    loop nacc ~gen:after ~nuc tl
  and merge_from_n acc genetic nuclear tl ~gen ~nuc =
    let gen_pos = genetic.bm.position + genetic.bm.length in
    let gen_els, gen_after = split_at_boundary_offset_before_and_drop_it
      ~pos:gen_pos ~offset:genetic.offset gen
    in
    let nuc_pos = nuclear.bm.position + nuclear.bm.length in
    let nuc_els, nuc_after = split_at_boundary_offset_before_and_drop_it
      ~pos:nuc_pos ~offset:nuclear.offset nuc
    in
    let nacc = MergeFromNuclear
      { genetic = replace_sequence genetic gen_els
      ; nuclear = splice_adj (replace_sequence nuclear nuc_els)
      } :: acc
    in
    loop nacc ~gen:gen_after ~nuc:nuc_after tl
  in
  start alg

(* When we have a nuc:a1, gen:a2 (ie. a1 <> a2):
   take a2's instructions and pair the nuc element! *)
let merge_different_nuc_into_alignment ?(fail_on_empty=true) nallele nuc alg =
  let rec start = function
    | (FillFromGenetic f as fi) :: tl when before_start f.bm ->
        loop [fi] nuc tl
    | []     -> error "Mapping instructions for %s, empty!" nallele
    | s :: _ -> error "Mapping instructions for %s, did not start with a \
                      \"starting\" FillFromGenetic but %s"
                  nallele
                  (instruction_to_string al_els_to_string al_els_to_string s)
  and loop acc nuc = function
    | MergeFromNuclear m :: tl        ->
        let nuclear = m.nuclear in
        let nuc_pos = nuclear.bm.position + nuclear.bm.length in
        let nuc_els, nuc_after = split_at_boundary_offset_before_and_drop_it
          ~pos:nuc_pos ~offset:nuclear.offset nuc
        in
        let nacc = MergeFromNuclear
          { m with nuclear =
            { nuclear with sequence =
              (* Mind blown! *)
              { genetic = nuclear.sequence
              ; nuclear = nuc_els
              }}} :: acc
        in
        loop nacc nuc_after tl
    | (FillFromGenetic _ as f) :: tl  -> loop (f :: acc) nuc tl
    | []                              ->
        if fail_on_empty && nuc <> [] then
          error "In merge different nuc for %s after all instructions nuclear \
            sequence not empty: %s."
            nallele (al_els_to_string nuc)
        else
          Ok (List.rev acc)                           (* Fin *)
  in
  start alg

let without_starts_or_ends =
  List.filter ~f:(function
    | Start _
    | End _       -> false
    | Boundary _
    | Gap _
    | Sequence _  -> true)

let align_reference instr =
  List.map instr ~f:(function
    | FillFromGenetic f ->
        if before_start f.bm then
          f.sequence
        else
          to_boundary ~offset:f.offset f.bm :: f.sequence
    | MergeFromNuclear { genetic; nuclear } ->
        to_boundary ~offset:genetic.offset genetic.bm ::
          (* This will work fine for A,B,C but not for DRB where there are
             restarts in the reference! *)
          without_starts_or_ends nuclear.sequence)
  |> List.concat

let empty_or_only_gaps = List.for_all ~f:is_gap

let handle_split_end ?(verbose=false) name seq next_boundary_pos
  ~no_split ~split_end_at_end ~split_end_middle =
  match split_at_end seq with
  | before, []         ->
      if verbose then printf "In %s split needed because no end.\n" name;
      no_split before
  | before, e :: after ->
    let ep = position e in
    if ep = next_boundary_pos then
      if after <> [] then
        invalid_argf "In %s malformed sequence: %s after end but before next boundary."
          name (al_els_to_string after)
      else begin
        if verbose then
          printf "In %s dropping end at %d before next_boundary_pos %d.\n"
            name ep next_boundary_pos;
        split_end_at_end before
      end
    else if empty_or_only_gaps after then begin
      if verbose then
        printf "In %s found only gaps for an early end at %d before next_boundary_pos %d.\n"
          name ep next_boundary_pos;
      split_end_middle ep   (* Return just the position so that the caller
                               handles what state is returned. *)
    end else
      invalid_argf "In %s Did not find just gaps after end: %s."
        name (al_els_to_string after)

type align_state =
  { acc         : string alignment_element list list
  ; in_nucleic  : bool      (* are inside the nucleic acid sequence?
                                         we have data for it! *)
  ; need_start  : bool      (* need to add a start. *)
  }

let init_align_state =
  { acc        = []
  ; in_nucleic = false
  ; need_start = false
  }

let append s ?in_nucleic ?need_start t =
  { acc        = t :: s.acc
  ; in_nucleic = Option.value in_nucleic ~default:s.in_nucleic
  ; need_start = Option.value need_start ~default:s.need_start
  }

let align_same ?(verbose=false) name instr =
  let check_for_ends genetic_end start nuc_seq next_boundary_pos s =
    if genetic_end then begin (* Trust that the ends are aligned. *)
      if verbose then
        printf "In %s not dropping end!" name;
      append s (start nuc_seq) ~in_nucleic:false ~need_start:true
    end else
      handle_split_end ~verbose name nuc_seq next_boundary_pos
        ~no_split:        (fun b -> append s (start b) ~in_nucleic:true)
        ~split_end_at_end:(fun b -> append s (start b) ~in_nucleic:false)
        (* Keep the actual end (and possible gap's) to signal missing *)
        ~split_end_middle:(fun _p -> append s (start nuc_seq)
                                      ~in_nucleic:false ~need_start:true)
  in
  let latest_position_or_invalid = function
    | []      -> invalid_argf "In %s empty acc, merging right away?" name
    | [] :: _ -> invalid_argf "In %s empty element on acc" name
    | ls :: _ -> List.last ls
                 |> Option.value_exn ~msg:"We checked for non empty list!"
                 |> latest_position
  in
  let { acc; _} =
    List.fold_left instr ~init:init_align_state ~f:(fun state instr ->
      match instr with
      | FillFromGenetic f ->
          if before_start f.bm then
            let need_start = not (List.exists ~f:is_start f.sequence) in
            append state f.sequence ~need_start
          else
            let sbm = to_boundary ~offset:f.offset f.bm in
            let region_has_sequence = List.exists f.sequence ~f:is_sequence in
            if state.need_start && region_has_sequence then
              append state (Start f.bm.position :: sbm :: f.sequence) ~need_start:false
            else
              append state (sbm :: f.sequence)
      | MergeFromNuclear { genetic; nuclear } ->
          let next_boundary_pos = nuclear.bm.position + nuclear.offset + nuclear.bm.length in
          let start_bndr = to_boundary ~offset:genetic.offset genetic.bm in
          let genetic_end = List.exists genetic.sequence ~f:is_end in
          if state.in_nucleic then  (* Within nuclear data *)
            check_for_ends genetic_end (fun l -> start_bndr :: l) nuclear.sequence
              next_boundary_pos state
          else
            begin match split_at_start nuclear.sequence with
            | before, [] ->
                (* Did NOT find a start when we haven't started we have missing data! *)
                if not (empty_or_only_gaps before) then
                  invalid_argf "In %s did not find just gaps before start: %s"
                    name (al_els_to_string before)
                else begin
                  let lp = latest_position_or_invalid state.acc in
                  if verbose then
                    printf "In %s haven't started insert an end at %d \n" name lp;
                  append state [End lp] ~in_nucleic:false ~need_start:true
                end
            | before, s :: after ->
                if not (empty_or_only_gaps before) then
                  invalid_argf "In %s did not find just gaps before start: %s"
                    name (al_els_to_string before)
                else
                  if position s = nuclear.bm.position + nuclear.offset + 1 then
                    check_for_ends genetic_end (fun l -> start_bndr :: l) after
                      next_boundary_pos state
                  else begin
                    let lp = latest_position_or_invalid state.acc in
                    if verbose then
                      printf "In %s about to start insert an end at %d \n" name lp;
                    check_for_ends genetic_end (fun l -> End lp :: s :: start_bndr :: l)
                      after next_boundary_pos state
                  end
            end)
  in
  List.concat (List.rev acc)

let align_different ?(verbose=false) name instr =
  let check_for_ends start ~genx ~nucx next_boundary_pos s =
    handle_split_end ~verbose (name ^ "_nuc") nucx next_boundary_pos
      ~no_split:        (fun b -> append s (start b) ~in_nucleic:true)
      ~split_end_at_end:(fun b -> append s (start b) ~in_nucleic:false)
      ~split_end_middle:(fun ep ->
        (* Have to consider possible ends in the genetic sequence! *)
        let nuc_before = List.filter nucx ~f:(fun a -> position a < ep) in
        let gafter_end = List.filter genx ~f:(fun a -> position a >= ep) in
        let everything_ends gen =
          append s (start (nuc_before @ gen)) ~in_nucleic:false ~need_start:false
        in
        handle_split_end ~verbose (name ^ "_gen") gafter_end next_boundary_pos
          ~no_split:(everything_ends)
          ~split_end_at_end:(everything_ends)
          (* Keep the actual end of genetic seq (and possible gap's) to signal missing *)
          ~split_end_middle:(fun _gp -> everything_ends gafter_end))
  in
  let check_for_ends_gen genetic_end start ~genx next_boundary_pos s =
    if genetic_end then begin (* Trust that the ends are aligned. *)
      if verbose then
        printf "In %s not dropping end because it aligns with gen-gen.!" name;
      append s (start genx) ~in_nucleic:false ~need_start:true
    end else
      handle_split_end ~verbose (name ^ "_gen_gen") genx next_boundary_pos
        ~no_split:        (fun b -> append s (start b))
        ~split_end_at_end:(fun b -> append s (start b))
        (* Keep the actual end (and possible gap's) to signal missing *)
        ~split_end_middle:(fun _p -> append s (start genx) ~need_start:true)
  in
  let { acc; _ } =
    List.fold_left instr ~init:init_align_state ~f:(fun state instr ->
      match instr with
      | FillFromGenetic f ->
          if before_start f.bm then
            let need_start = not (List.exists ~f:is_start f.sequence) in
            append state f.sequence ~need_start
          else
            let sbm = to_boundary ~offset:f.offset f.bm in
            let region_has_sequence = List.exists f.sequence ~f:is_sequence in
            if state.need_start && region_has_sequence then
              append state (Start f.bm.position :: sbm :: f.sequence)
                ~need_start:false
            else
              append state (sbm :: f.sequence)
      | MergeFromNuclear { genetic; nuclear } ->
          let next_boundary_pos =
            nuclear.bm.position + nuclear.offset + nuclear.bm.length
          in
          let start_bndr = to_boundary ~offset:genetic.offset genetic.bm in
          let genx = nuclear.sequence.genetic in
          let nucx = nuclear.sequence.nuclear in
          if state.in_nucleic then  (* Within nuclear data *)
            check_for_ends (fun l -> start_bndr :: l) ~genx ~nucx
              next_boundary_pos state
          else
            let genx = List.filter ~f:(fun e -> not (is_start e)) genx in
            begin match split_at_start nucx with
            | before, [] ->
                (* Did NOT find a start when we haven't started we have missing data! *)
                if not (empty_or_only_gaps before) then
                  invalid_argf "In %s did not find just gaps before start: %s"
                    name (al_els_to_string before)
                else begin
                  if verbose then
                    printf "In %s haven't started inserting genetic sequence \
                      instead!\n" name;
                  let genetic_end = List.exists genetic.sequence ~f:is_end in
                  check_for_ends_gen genetic_end (fun l -> start_bndr :: l) ~genx
                    next_boundary_pos state
                end
            | before, s :: after ->
                if not (empty_or_only_gaps before) then
                  invalid_argf "In %s did not find just gaps before start: %s"
                    name (al_els_to_string before)
                else
                  if position s = nuclear.bm.position + nuclear.offset + 1 then
                    check_for_ends (fun l -> start_bndr :: l) ~genx ~nucx:after
                      next_boundary_pos state
                  else begin
                    (* TODO: what if before_gen is empty? insert Start/End ?*)
                    let before_gen =
                      List.filter genx ~f:(fun a -> position a < position s)
                    in
                    if verbose then
                      printf "In %s before start at %d merging genetic exon data %s.\n" name
                        (position s) (al_els_to_string before_gen);
                    check_for_ends (fun l -> start_bndr :: (before_gen @ l))
                      ~genx ~nucx:after next_boundary_pos state
                  end
            end)
  in
  List.concat (List.rev acc)

let reference_positions_align ?seq lst =
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

let same_alts instr =
  list_fold_ok ~init:([], []) ~f:(fun (instr_acc, alt_acc) (allele, gen, nuc) ->
    map_instr_to_alignments allele ~gen ~nuc instr >>= fun ainstr ->
    let all_elems = align_same ~verbose:true ("same: " ^ allele) ainstr in
      Ok ( (allele, ainstr) :: instr_acc
         , (allele, all_elems) :: alt_acc))

let init_trie elems =
  let open Nomenclature in
  list_fold_ok elems ~init:Trie.empty ~f:(fun trie s ->
    parse s >>= fun (_gene, (allele_resolution, suffix_opt)) ->
      Ok (Trie.add allele_resolution suffix_opt trie))

module RMap = Map.Make (struct
    type t = Nomenclature.resolution * Nomenclature.suffix option
    let compare = compare
    end)

let init_trie_and_map elems =
  let open Nomenclature in
  list_fold_ok elems ~init:(Trie.empty, RMap.empty)
    ~f:(fun (trie, mp) (s, seq) ->
        parse s >>= fun (_gene, (allele_resolution, suffix_opt)) ->
          Ok ( Trie.add allele_resolution suffix_opt trie
             , RMap.add (allele_resolution, suffix_opt) seq mp))

let diff_alts ?(verbose=false) ~trie ~rmap ~na =
  let open Nomenclature in
  list_fold_ok na ~init:([],[]) ~f:(fun (alt_acc, map_acc) (nallele, nuc) ->
    (* Check that the gene is always the same? *)
    parse nallele >>= fun (gene, (allele_resolution, _so)) ->
      let closest_allele_res = Trie.nearest allele_resolution trie in
      let alg = RMap.find closest_allele_res rmap in
      let closest_allele_str = resolution_and_suffix_opt_to_string ~gene closest_allele_res in
      let name = sprintf "%s (nuc) -> %s (gen)" nallele closest_allele_str in
      if verbose then printf "Merging %s.\n" name;
      merge_different_nuc_into_alignment nallele nuc alg >>= fun instr ->
        let new_alts = align_different ~verbose name instr in
        Ok ( (nallele, new_alts) :: alt_acc
           , (nallele, closest_allele_str) :: map_acc))

(* TODO: Figure out a way around this special case. Maybe when we do our
   own alignment? *)

let c0409N_exon_extension =
  "GACAGCTGCCTGTGTGGGACTGAGATGCAGGATTTCTTCACACCTCTCCTTTGTGACTTCAAGAGCCTCTGGCATCTCTTTCTGCAAAGGCATCTGA"

(* TODO. We need to make a decision about suffixed (ex. with 'N') alleles
   since they may be the closest. *)

let same_and_diff ~gen_assoc ~nuc_assoc =
  let gen_sort = List.sort ~cmp:compare gen_assoc in
  let nuc_sort = List.sort ~cmp:compare nuc_assoc in
  let rec loop s dn dg n = function
    | [] -> s, dg, n @ dn
    | ((ga,gl) :: gt) as glst ->
        match n with
        | [] ->
            s, (glst @ dg), dn
        | ((na,nl) :: nt) as nlst ->
            let r = compare na ga in
            if r = 0 then
              loop ((ga, gl, nl) :: s) dn dg nt gt
            else if r < 0 then
              loop s ((na,nl) :: dn) dg nt glst
            else (* r > 0 *)
              loop s dn ((ga,gl) :: dg) nlst gt
  in
  loop [] [] [] nuc_sort gen_sort

let reference_as_diff lst =
  let no_sequences = List.filter ~f:(fun a -> not (is_sequence a)) in
  List.map lst ~f:(function
    | FillFromGenetic f ->
        FillFromGenetic { f with sequence = no_sequences f.sequence }
    | MergeFromNuclear { nuclear; genetic } ->
        MergeFromNuclear
          { nuclear = { nuclear with sequence = no_sequences nuclear.sequence }
          ; genetic = { genetic with sequence = no_sequences genetic.sequence }
          })

let and_check prefix =
  align_from_prefix prefix >>= fun (gen_mp, nuc_mp, instr) ->
    if gen_mp.reference <> nuc_mp.reference then
      error "References don't match %s vs %s" gen_mp.reference nuc_mp.reference
    else if gen_mp.align_date <> nuc_mp.align_date then
      error "Align dates don't match %s vs %s" gen_mp.align_date nuc_mp.align_date
    else
      let nuc = nuc_mp.ref_elems in
      let gen = gen_mp.ref_elems in
      map_instr_to_alignments nuc_mp.reference ~gen ~nuc instr >>= fun ref_instr ->
        let new_ref_elems = align_reference ref_instr in
        reference_positions_align ~seq:("reference:" ^ gen_mp.reference)
          new_ref_elems >>= fun _ref_position_check ->
            let gen_assoc = gen_mp.alt_elems in
            let nuc_assoc = nuc_mp.alt_elems in
            let same, just_gen, just_nuc = same_and_diff ~gen_assoc ~nuc_assoc in
            let () =
              if just_gen <> [] then
                printf "Found these alleles with only genetic data: %s\n"
                  (String.concat ~sep:"; " (List.map just_gen ~f:fst))
            in
            (* Add the same, alleles with both genetic and nucleic data.*)
            same_alts instr same >>= fun (alt_inst, alt_als) ->
              let rdiff = reference_as_diff ref_instr in
              let instr_assoc = (gen_mp.reference, rdiff) :: alt_inst in
              (* Create a trie and map for lookups *)
              init_trie_and_map instr_assoc >>= fun (trie, rmap) ->
                (* Add the alleles with just nucleic data. *)
                diff_alts ~verbose:true ~trie ~rmap ~na:just_nuc >>=
                  fun (diff_alt_lst, diff_map_lst) ->
                    let map_lst = List.map same ~f:(fun (a, _, _) -> (a,a)) in
                    printf "Finished merging!\n%!";
                    Ok ({ align_date = gen_mp.align_date
                        ; reference  = gen_mp.reference
                        ; ref_elems  = new_ref_elems
                        ; alt_elems  = alt_als @ diff_alt_lst
                        } ,
                        map_lst @ diff_map_lst)
