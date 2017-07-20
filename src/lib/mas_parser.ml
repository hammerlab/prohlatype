(* How this works:

This is a line based format. For our purposes, there are three types of lines
as determined by {type line}, where `SeqData` has most of the actual alignment
data. We divide the parsing of sequences into two types: reference and
alternatives. We take special care of keeping track of gaps in the reference so
that positions, when parsing the alternatives, are annotated with regard to
the reference positions.

*)

open Util

(* This refers to the alignment position.  *)
type position = int

(* We'll parse N as a nucleotide, but why does IMGT report this! *)
let is_nucleotide allele  = function
  | 'A' | 'C' | 'G' | 'T' -> true
  | 'N' -> eprintf "IMGT reports that %s has an 'N'as a base!" allele; true
  | _   -> false

let is_amino_acid =
  function
  | 'A'| 'C'| 'D'| 'E'| 'F'| 'G'| 'H'| 'I'| 'K'| 'L'
  | 'M'| 'N'| 'P'| 'Q'| 'R'| 'S'| 'T'| 'V'| 'W'| 'Y' -> true
  | _ -> false

type 'e parse_struct =
  { allele    : string
  (* For now, it makes it easier to have some kind of sane delimiter to align
     and spot check these alignments. The "|" are the 'boundaries'. This keeps
     track of the most recently encountered boundary marker, starting with 0. *)
  ; boundary  : int
  (* Where we are in the sequence, this starts from the first number specified
     in the file and increments as we read characters. *)
  ; position  : position
  ; sequence  : 'e list
  ; in_data   : bool
  }

let init_ps allele position =
  { allele
  ; position
  ; boundary = 0
  ; sequence = []
  ; in_data  = false
  }

let where ps =
  sprintf "allele: %s, position: %d, sequence length: %d"
    ps.allele ps.position (List.length ps.sequence)

type 'sr sequence = { start : position; s : 'sr }
and gap = { gstart : position; length : int }
and boundary = { idx : int; pos : position }
and 'sr alignment_element =
  | Start of position
  | End of position
  | Boundary of boundary
  | Sequence of 'sr sequence
  | Gap of gap

let start_position = function
  | Start s               -> s
  | End e                 -> e
  | Boundary { pos; _ }   -> pos
  | Sequence { start; _ } -> start
  | Gap { gstart; _ }     -> gstart

let end_position sr_length = function
  | Start s     -> s
  | End e       -> e
  | Gap g       -> (g.gstart + g.length)
  | Sequence s  -> (s.start + sr_length s.s)
  | Boundary b  -> b.pos + 1

let al_el_to_string = function
  | Start p                 -> sprintf "Start at %d" p
  | End p                   -> sprintf "End at %d" p
  | Boundary { idx; pos }   -> sprintf "Boundary %d at %d" idx pos
  | Sequence { start; s }   -> sprintf "Sequence %s at %d" s start
  | Gap { gstart; length }  -> sprintf "Gap of %d from %d" length gstart

let is_sequence = function | Sequence _ -> true | _ -> false
let is_gap      = function | Gap _      -> true | _ -> false
let is_end      = function | End _      -> true | _ -> false
let is_start    = function | Start _    -> true | _ -> false

type data_switch =
  | Data
  | NotData
  | NoEffect

let next ps ds f =
  let npos = ps.position + 1 in
  let insert_start_before_boundary = function
    | (Boundary {pos; _} as b) :: tl when pos + 1 = npos ->
             b :: (Start pos) :: tl
    | lst -> (Start npos) :: lst
  in
  let in_data, new_seq =
    match ds, ps.in_data with
    | Data,     true    -> true,       ps.sequence
    | Data,     false   -> true,       insert_start_before_boundary ps.sequence
    | NotData,  true    -> false,      (End npos) :: ps.sequence
    | NotData,  false   -> false,      ps.sequence
    | NoEffect, _       -> ps.in_data, ps.sequence
  in
  { ps with in_data; position = npos; sequence = f npos new_seq }

let insert_boundary ps =
  let nps =
    next ps NoEffect (fun position sequence ->
        Boundary { idx = ps.boundary; pos = position} :: sequence)
  in
  { nps with boundary = ps.boundary + 1 }

let insert_gap ps =
  next ps NoEffect (fun position -> function
    | Gap { gstart; length } :: t when gstart + length = position
                            -> Gap { gstart; length = length + 1 } :: t
    | []
    | End _ :: _
    | Start _ :: _
    | Boundary _ :: _
    | Gap _ :: _
    | Sequence _ :: _ as l  -> Gap { gstart = position; length = 1 } :: l)

let insert_nuc_error fmt =
  invalid_argf ~prefix:"Trying to insert sequence element " fmt

let insert_same ~fail_on_same ps =
  if fail_on_same then
    invalid_argf "Encountered unexpected '-' same char for : %s" (where ps)
  else
    next ps Data (fun _position sequence -> sequence)

let insert_nuc c ps =
  next ps Data (fun position -> function
    | Sequence {start; s} :: t when start + (List.length s) = position
                            -> Sequence { start; s = c :: s} :: t
    | []
    | End _ :: _            -> invalid_argf "Adding a Nuc %c %d %s after End or Empty!"
                                  c position ps.allele
    | Start _ :: _
    | Boundary _ :: _
    | Gap _ :: _
    | Sequence _ :: _ as l  -> Sequence { start = position; s = c :: [] } :: l )

let insert_unknown ps =
  next ps NotData (fun _position sequence -> sequence)

let update ~allele ~dna ~fail_on_same ps s =
  let is_nuc = if dna then is_nucleotide allele else is_amino_acid in
  let rec to_ref_seq_elems_char ps = function
    | []                    -> ps
    | '|' :: t              -> to_ref_seq_elems_char (insert_boundary ps) t
    | '*' :: t              -> to_ref_seq_elems_char (insert_unknown ps) t
    | 'X' :: t when not dna -> to_ref_seq_elems_char (insert_unknown ps) t (* add End *)
    | '.' :: t              -> to_ref_seq_elems_char (insert_gap ps) t
    | '-' :: t              -> to_ref_seq_elems_char (insert_same ~fail_on_same ps) t
    | c :: t when is_nuc c  -> to_ref_seq_elems_char (insert_nuc c ps) t
    | x :: _                -> invalid_argf "Unrecognized char '%c' in %s" x (where ps)
  in
  to_ref_seq_elems_char ps (String.to_character_list s)

let gaps_to_string gps =
  String.concat ~sep:";" (List.map ~f:(fun (p,l) -> sprintf "(%d,%d)" p l) gps)

type 'sr parse_result =
  { dna       : bool        (* DNA or Amino Acid sequence -> diff characters *)
  ; start_pos : position
  ; ref       : string      (* Name of reference. *)
  (* As we parse the alternative tracks, we have to Keep track of the gaps that
     we encounter in the reference, so that all positions are with respect to
     the reference. *)
  ; ref_ps    : 'sr parse_struct
  ; alg_htbl  : (string, 'sr parse_struct) Hashtbl.t
  }

let empty_result ref_allele dna position =
  { dna       = dna
  ; start_pos = position - 1
  ; ref       = ref_allele
  ; ref_ps    = init_ps ref_allele (position - 1)
  ; alg_htbl  = Hashtbl.create 100
  }

let reverse_seq lst =
  let to_string l = String.of_character_list (List.rev l) in
  List.rev lst
  |> List.map ~f:(function
    | Start _
    | End _
    | Boundary _
    | Gap _ as e            -> e
    | Sequence { start; s } -> Sequence { start; s = to_string s })

let normalized_seq ps =
  let rec has_end = function
    | End _ :: _      -> true
    | Boundary _ :: t -> has_end t
    | Gap _ :: t      -> has_end t
    | _               -> false
  in
  if has_end ps.sequence then
    reverse_seq ps.sequence
  else
    reverse_seq (End (ps.position + 1) :: ps.sequence)

type line =
  | Position of bool * int  (* nucleotide or amino acid sequence  *)
  | Dash
  | SeqData of string * string list

(* Assume that it has been trimmed. *)
let parse_data line =
  String.split line ~on:(`Character ' ')
  |> List.filter ~f:((<>) String.empty)
  |> function
      | s :: _ when String.get s 0 = Some '|' -> Dash
      | "AA" :: "codon" :: _         -> Dash (* not really but not modeling this at the moment. *)
      | "gDNA" :: pos :: _           -> Position (true, int_of_string pos)
      | "cDNA" :: pos :: _           -> Position (true, int_of_string pos)
      | "Prot" :: pos :: _           -> Position (false, int_of_string pos)
      | []                           -> invalid_arg "Empty data line!"
      | s :: lst                     -> SeqData (s, lst)

type parse_state =
  | Header
  | Empty
  | Data of line

type result =
  { align_date  : string
  ; reference   : string
  ; ref_elems   : string alignment_element list
  ; alt_elems   : (string * string alignment_element list) list
  }

let report = ref false

let parse_align_date ic =
  try
    let rec loop () =
      let line = input_line ic in
      match String.split ~on:(`String ":") line with
      | ["Sequences Aligned"; ad ] -> Some ad
      | _ -> loop ()
    in
    loop ()
  with End_of_file ->
    None

let from_in_channel ic =
  (*let previous_reference_position = ref min_int in*)
  let latest_reference_position = ref min_int in
  let update x = function
    (* Sometimes, the files position counting seems to disagree with this
       internal count, usually because of multiple boundaries. Not certain
       how to get to the bottom, but my manual string counts lead me to
       believe that there isn't a bug in the parsing code. One possibility is
       that there is no '0' the position in the files; their indices are
       [-1, 1, 2].

       So we don't check for: x.ref_ps.position = p as well. *)
    | Position (dna, p) -> assert (x.dna = dna); x
    | Dash              -> x (* ignore dashes *)
    | SeqData (allele, s) ->
      let fold = update ~allele ~dna:x.dna in
      if x.ref = allele then begin
        (*let prev_pos = x.ref_ps.position in *)
        let nref_ps = List.fold_left ~f:(fold ~fail_on_same:true) ~init:x.ref_ps s in
        latest_reference_position := nref_ps.position;
        { x with ref       = allele
               ; ref_ps    = nref_ps
        }
      end else begin
        let cur_ps =
          try Hashtbl.find x.alg_htbl allele
          with Not_found -> init_ps allele x.start_pos
        in
        let new_ps = List.fold_left ~f:(fold ~fail_on_same:false) ~init:cur_ps s in
        (* Can't make this into an assertion because of sequences such as
            C*04:09N that have sequences extending after the end of the
            reference. *)
        if !report
           (*&& new_ps.in_data <> After*)
           && new_ps.position <> !latest_reference_position then
          printf
            "position mismatch %d vs %d for %s.\n"
            !latest_reference_position new_ps.position new_ps.allele;
        Hashtbl.replace x.alg_htbl allele new_ps;
        x
      end
  in
  let rec loop state acc =
    match input_line ic |> String.strip ~on:`Both with
    | exception End_of_file -> acc
    | line ->
      match state with
      | Header when String.is_empty line -> loop Empty acc
      | Header                           -> loop Header acc
      | Empty  when String.is_empty line -> loop Empty acc
      | Empty                            -> if String.is_prefix line ~prefix:"Please" then
                                              acc
                                            else
                                              let d = parse_data line in
                                              loop (Data d) (update acc d)
      | Data _ when String.is_empty line -> loop Empty acc
      | Data _ ->                           let d = parse_data line in
                                            loop (Data d) (update acc d)
  in
  let rec loop_header state =
    match input_line ic |> String.strip ~on:`Both with
    | exception End_of_file -> invalid_arg "Didn't get to the data!"
    | line ->
      match state with
      | Header when String.is_empty line -> loop_header Empty
      | Header                           -> loop_header Header
      | Empty  when String.is_empty line -> loop_header Empty
      | Empty                            ->
          begin
            let d = parse_data line in
            match d with
            | Position _ -> loop_header (Data d)
            | _          -> invalid_arg "First data not position."
          end
      | Data _ when String.is_empty line -> loop_header state
      | Data (Position (dna, p)) ->
          begin
            match parse_data line with
            | SeqData (allele, _) as d -> let res = empty_result allele dna p in
                                          loop (Data d) (update res d)
            | _                        -> loop_header state
          end
      | Data _ -> loop_header state
  in
  match parse_align_date ic with
  | None ->
      close_in ic;
      invalid_argf "Couldn't extract sequence align date."
  | Some align_date ->
    let reversed = loop_header Header in
    let ref_elems = normalized_seq reversed.ref_ps in
    let alt_elems =
      Hashtbl.fold ~init:[] ~f:(fun ~key:all ~data:ps acc ->
          if ps.sequence = [] then begin
            printf "Dropping empty sequence: %s\n" ps.allele;
            acc
          end else
            (all, normalized_seq ps) :: acc)
        reversed.alg_htbl
    in
    { align_date
    ; reference = reversed.ref
    ; ref_elems
    ; alt_elems
    }

let from_file f =
  let ic = open_in f in
  try
    let r = from_in_channel ic in
    close_in ic;
    r
  with e ->
    close_in ic;
    raise e

let sequence_length l =
  List.fold_left l ~init:(None, 0) ~f:(fun (spo, sum) m ->
      match spo, m with
      | None,   Start s -> (Some s), sum
      | Some p, End e   -> None, (e - p + sum)
      | _,      _       -> spo, sum)
    |> snd

module Boundaries = struct

  type marker =
    { index       : int   (* Which segment? *)
    ; position    : int   (* Position of the boundary marker. *)
    ; length      : int   (* Length to next boundary, or 1 plus the total length
                             of the following segment.*)
    ; seq_length  : int   (* Length of sequences contained in this boundary.*)
    }

  let marker ~index ~position = { index; position; length = 0 ; seq_length = 0}
  let before_start_boundary = marker ~index:(-1) ~position:(-1)
  let before_start m = m.index = -1

  let to_boundary ~offset { index; position; _ } =
    Boundary { idx = index; pos = position + offset }

  let matches_boundary { index; position; _ } = function
    | Boundary b when b.idx = index &&  b.pos = position -> true
    | _ -> false

  let marker_to_string { index; position; length; seq_length } =
    sprintf "{index: %d; position: %d; length %d; seq_length: %d }"
      index position length seq_length

  let marker_to_boundary_string { index; position; _ } =
    sprintf "{idx: %d; pos: %d}" index position

  let fold ~boundary ~start ~end_ ~gap ~seq lst =
    let il s m i =
      if s then { m with length = m.length + i; seq_length = m.seq_length + i }
           else { m with length = m.length + i}
    in
    let rec loop curb curs acc = function
      | []                ->
          (il false curb 1, curs) :: acc
      | Boundary bp :: t  ->
          let newb = marker ~index:bp.idx ~position:bp.pos in
          let nacc = ((il false curb 1, curs) :: acc) in
          loop newb (boundary (Some (curs, bp, newb)))  nacc  t
      | Start s :: t      ->
          let newb = if before_start curb then { curb with position = s - 1 } else curb in
          loop newb (start curs s)                      acc   t
      | End e :: t        ->
          loop curb (end_ curs e)                       acc   t
      | Gap g :: t        ->
          let newb = il false curb g.length in
          loop newb (gap curs g)                        acc   t
      | Sequence s :: t   ->
          let newb = il true curb (String.length s.s) in
          loop newb (seq curs s)                        acc   t
    in
    let init = boundary None in
    loop before_start_boundary init [] lst

  let bounded lst =
    let boundary _ = "" in
    let start s _ = s in
    let end_ s _ = s in
    let gap s _ = s in
    let seq s ss = s ^ ss.s in
    List.rev (fold lst ~boundary ~start ~end_ ~gap ~seq)

end (* Boundaries *)

(* Label these arguments to differentiate that we want a "length/index" into
   contained element, not a position. *)
let split_sequence seq ~pos =
  let index = pos - seq.start in
  let before, after = String.split_at seq.s ~index in
  { start = seq.start;          s = before },
  { start = seq.start + index;  s = after }

let split_gap gap ~pos =
  let length = pos - gap.gstart in
  { gstart = gap.gstart;          length = length },
  { gstart = gap.gstart + length; length = gap.length - length }

(* Applying the alignment elements. *)
module type Data_projection = sig

  type t

  val to_string : t -> string

  (* Initialize the projection. *)
  val of_seq : in_ref:bool -> bool -> string sequence -> t

  val of_gap : in_ref:bool -> bool -> gap -> t

  val of_boundary : in_ref:bool -> bool -> boundary -> t

  (* Mutate the projection. *)
  val add_seq : position -> t -> string sequence -> t

  val add_gap : position -> t -> gap -> t

  (* Pass a "length" into the curreng segment/projection. *)
  val split_due_to_start_stop : int -> bool -> t -> t option * t

end

module MakeZip (R : Data_projection) = struct

  type state =
    { start_pos : int
    ; end_pos   : int
    ; project   : R.t
    ; started   : bool
    }

  let state_to_string s =
    sprintf "{ start_pos = %d; \
             ; end_pos = %d; \
             ; project = %s; \
             ; started = %b }"
    s.start_pos s.end_pos (R.to_string s.project) s.started

  let add_gap cur g =
    { cur with project = R.add_gap cur.start_pos cur.project g }

  let add_seq cur s =
    { cur with project = R.add_seq cur.start_pos cur.project s }

  let split st pos new_start =
    let index = pos - st.start_pos in
    let before_opt, after = R.split_due_to_start_stop index new_start st.project in
    Option.map before_opt ~f:(fun before -> { st with project = before ; end_pos = pos }),
    { start_pos = pos
    ; end_pos   = st.end_pos
    ; project   = after
    ; started   = new_start
    }

  let apply ~reference ~allele =
    let append acc cur =
      if !report then printf "Appending %s \n" (R.to_string cur.project);
      cur.project :: acc
    in
    let rec advance_allele cur acc a =
      match a with
      | []      -> cur.started, (append acc cur), []
      | h :: t  ->
          let sa = start_position h in
          if !report then
            printf "at %s %d when cur is %s.\n"
              (al_el_to_string h) sa (state_to_string cur);
          if sa < cur.start_pos then
            advance_allele cur acc t
          else if sa < cur.end_pos then
            mutate_segment cur acc t h
          else (* sa >= cur.end_pos *)
            cur.started, append acc cur, a
    and mutate_segment cur acc at = function
      | Start pos       -> let b_opt, a = split cur pos true in
                           let nacc = Option.value_map b_opt ~default:acc ~f:(append acc) in
                           advance_allele a nacc at
      | End pos         -> let b_opt, a = split cur pos false in
                           let nacc = Option.value_map b_opt ~default:acc ~f:(append acc) in
                           advance_allele a nacc at
      | Boundary b as e -> if cur.start_pos <> b.pos &&
                              cur.end_pos <> cur.start_pos + 1 then
                              invalid_argf "Found %s in sequence at %d!"
                                (al_el_to_string e) (cur.start_pos)
                           else (* ignore *)
                             advance_allele cur acc at
      | Sequence ss     ->
          let slen = String.length ss.s in
          if ss.start + slen > cur.end_pos then
            let nc, na = split_sequence ss ~pos:cur.end_pos in
            advance_allele (add_seq cur nc) acc (Sequence na :: at)
          else
            let ncur = add_seq cur ss in
            advance_allele ncur acc at
      | Gap g           ->
          if g.gstart + g.length > cur.end_pos then
            let gc, ga = split_gap g ~pos:cur.end_pos in
            advance_allele (add_gap cur gc) acc (Gap ga :: at)
          else
            let ncur = add_gap cur g in
            advance_allele ncur acc at
    in
    let reference_pass =
      Boundaries.fold reference
        ~start:(fun s p -> s)          (* Nothing to do for reference Start's *)
        ~end_:(fun s p -> s)             (* Nothing to do for reference End's *)
        ~boundary:(function
          | None            -> (false, [], allele)
          | Some ((started, acc, allele_lst), b, _m) ->
                let new_state =
                  { start_pos = b.pos
                  ; end_pos = b.pos + 1
                  ; project = R.of_boundary ~in_ref:true started b
                  ; started
                  }
                in
                advance_allele new_state [] allele_lst)
        ~gap:(fun (started, acc, allele_lst) gap ->
                  let new_state =
                    { start_pos = gap.gstart
                    ; end_pos = gap.gstart + gap.length
                    ; project = R.of_gap ~in_ref:true started gap
                    ; started
                    }
                  in
                  advance_allele new_state acc allele_lst)
        ~seq:(fun (started, acc, allele_lst) seq ->
                  let new_state =
                    { start_pos = seq.start
                    ; end_pos = seq.start + String.length seq.s
                    ; project = R.of_seq ~in_ref:true started seq
                    ; started
                    }
                  in
                  advance_allele new_state acc allele_lst)
  in
  let just_bm_and_proj = List.rev_map ~f:(fun (bm, (_, p, _)) -> (bm, p)) in
  match reference_pass with
  | [] -> invalid_argf "Didn't even have a start boundary?"
  | (bm, (final_st, acc, alleles_lst)) :: tl ->
      let final_st, nacc =
        List.fold_left alleles_lst ~init:(final_st, acc) ~f:(fun (st, a) e ->
          match e with
          | End _      -> false, a
          | Start _    -> true, a
          | Boundary b -> st, (R.of_boundary ~in_ref:false st b) :: a
          | Gap g      -> st, (R.of_gap ~in_ref:false st g) :: a
          | Sequence s -> st, (R.of_seq ~in_ref:false st s) :: a)
      in
      just_bm_and_proj ((bm, (final_st, nacc, [])) :: tl)

end (* MakeZip *)

module B = Sosa.Native_bytes

let default_gap_char = '.'

module AlleleSequences = MakeZip (struct

  type t = bool * B.t

  let to_string (s,b) =
    sprintf "%b,%s" s (B.to_native_string b)

  (* Ignore the in_reference parameter, because we want to create the buffer
     for the allele when it is outside the reference. *)
  let of_seq ~in_ref started s =
    match B.of_native_string s.s with
    | `Error (`wrong_char_at d) ->
        invalid_argf "wrong char at %d in %s" d s.s
    | `Ok b ->
        started, b

  let of_gap ~in_ref started g =
    started, B.make g.length default_gap_char

  (* Boundary markers are added by the allele_sequence function. *)
  let of_boundary ~in_ref started b =
    started, B.empty

  let add_seq start_pos (started, buffer) {start; s} =
    let off = start - start_pos in
    for index = 0 to String.length s - 1 do
      B.mutate_exn buffer ~index:(index + off) (String.get_exn s ~index)
    done;
    (started, buffer)

  let add_gap start_pos (stared, buffer) {gstart; length} =
    let off = gstart - start_pos in
    for i = 0 to length - 1 do
      B.mutate_exn buffer ~index:(i + off) default_gap_char
    done;
    (stared, buffer)

  let split_due_to_start_stop index new_start (st, buffer) =
    if index = 0 then
      None,         (new_start, buffer)
    else
      let b, a = B.split_at buffer ~index in
      Some (st, b), (new_start, a)

end)

let allele_sequences ~reference ~allele =
  let ss = AlleleSequences.apply ~reference ~allele in
  List.map ss ~f:(fun (bm, plst) ->
    bm,
    B.(to_native_string (concat
      (List.rev_filter_map plst ~f:(function
        | (false, _) -> None
        | (true, t)  -> Some (filter t ~f:(fun c -> c <> default_gap_char)))))))

let allele_sequence ?boundary_char ~reference ~allele () =
  let sep = Option.value_map ~default:B.empty ~f:B.of_character boundary_char in
  allele_sequences ~reference ~allele
  |> List.map ~f:snd
  |> B.concat ~sep

let reference_sequence_from_ref_alignment_elements ?boundary_char l =
  begin match boundary_char with
  | None    ->
    List.filter_map l ~f:(function | Sequence s -> Some s.s | _ -> None)
  | Some bc ->
    let bs = String.of_character bc in
    List.filter_map l ~f:(function | Sequence s -> Some s.s | Boundary _ -> Some bs | _ -> None)
  end
  |> String.concat ~sep:""

let reference_sequence ?boundary_char mp =
  reference_sequence_from_ref_alignment_elements ?boundary_char mp.ref_elems

let split_by_boundaries_rev lst =
  let rec loop acc gacc = function
    | [] -> (List.rev acc) :: gacc
    | Boundary _ :: t ->
        loop [] ((List.rev acc) :: gacc) t
    | e :: t ->
        loop (e :: acc) gacc t
  in
  loop [] [] lst

module DistanceProjection = struct

  type r =
    { started     : bool
    ; start_pos   : position
    ; length      : int
    ; mismatches  : int
    }

  let r_to_string r =
    sprintf "{%b; %d; %d; %d}"
      r.started r.start_pos r.length r.mismatches

  let split_r length new_start r =
    { r with length },
    { started     = new_start
    ; start_pos   = r.start_pos + length
    ; length      = r.length - length
    ; mismatches  = 0
    }

  type t =
    | G of r
    | S of r
    | B

  let to_string = function
    | G r -> sprintf "G %s" (r_to_string r)
    | S r -> sprintf "S %s" (r_to_string r)
    | B   -> "B"

  let incr m = function
    | S r -> S { r with mismatches = r.mismatches + m }
    | G r -> G { r with mismatches = r.mismatches + m }
    | B   -> invalid_argf "Asked to increase Boundary by %d" m

  let of_seq ~in_ref started seq =
    let start_pos = seq.start in
    let length = String.length seq.s in
    if in_ref then
      S { started; start_pos; length; mismatches = 0 }
    else
      S { started; start_pos; length; mismatches = length }

  let of_gap ~in_ref started (gap : gap) =
    let start_pos = gap.gstart in
    let length = gap.length in
    if in_ref then
      G { started; start_pos; length; mismatches = 0 }
    else
      G { started; start_pos; length; mismatches = gap.length }

  let of_boundary ~in_ref started boundary =
    B

  (* zipping logic makes sure position is inside s *)
  let add_seq _ t {s; _} =
    incr (String.length s) t

  (* seqs in seqs are diffs! *)
  let add_gap _ t ({length; _} : gap) =
    match t with
    | S r -> S { r with mismatches = r.mismatches + length } (* filling a gap *)
    | G _ -> t                                  (* gap in gap, no difference. *)
    | B   -> invalid_argf "Asked to add_gap of %d to Boundary." length

  let split_due_to_start_stop index new_start = function
    | S r when index = 0 -> None, S { r with started = new_start }
    | S r                -> let b, a = split_r index new_start r in
                            Some (S b) , S a
    | G r when index = 0 -> None , G { r with started = new_start }
    | G r                -> let b, a = split_r index new_start r in
                            Some (G b) , G a
    | B                 -> if index = 0 then None, B else
                            invalid_argf "Asked to increase Boundary by %d." index

end

module AlleleToRefDistance = MakeZip (struct
  include DistanceProjection
end)

type allele_segment_relationship =
  | Missing
  | Partial of int    (* sequence length *)
  | Full of int       (* sequence length might be > reference length *)

type 'a segment =
  { reference_seq_length  : int
  ; mismatches            : int
  ; allele_relationships  : 'a
  }

let allele_distances ~reference ~allele =
  let open DistanceProjection in
  let update_state mismatches seq_length = function
    | `Missing,        true  -> `Start (mismatches, seq_length)
    | `Missing,        false -> `NoStart
    | `NoStart,        true  -> `Partial (mismatches, seq_length)
    | `NoStart,        false -> `NoStart
    | `Start (m, l),   true  -> `Start (m + mismatches, l + seq_length)
    | `Start (m, l),   false -> `Partial (m, l)
    | `Partial (m, l), true  -> `Partial (m + mismatches, l + seq_length)
    | `Partial (m, l), false -> `Partial (m, l)
  in
  AlleleToRefDistance.apply ~reference ~allele
  |> List.map ~f:(fun (bm, lst) ->
      (* `Missing : the default state when we have observed 0 projections
         `Start   : we have observed at least 1 started projection and no not-started
         `NoStart : we have observed at least 1 not-started projection and no started
         `Partial : we have observed both. *)
      let final_state = List.fold_left lst ~init:`Missing ~f:(fun state p ->
        match p with
        | B   -> state
        | G r -> update_state r.mismatches 0 (state, r.started)
        | S r -> update_state r.mismatches r.length (state, r.started))
      in
      match final_state with
      | `Missing
      | `NoStart        -> { reference_seq_length = bm.Boundaries.seq_length
                           ; mismatches           = 0
                           ; allele_relationships = Missing
                           }
      | `Start (m, l)   -> { reference_seq_length = bm.Boundaries.seq_length
                           ; mismatches           = m
                           ; allele_relationships = Full l
                           }
      | `Partial (m, l) -> { reference_seq_length = bm.Boundaries.seq_length
                           ; mismatches           = m
                           ; allele_relationships = Partial l
                           })

let split_long_al_els pos = function
  | Start _
  | End _
  | Boundary _ as e -> invalid_argf "Asked to split %s at %d" (al_el_to_string e) pos
  | Sequence seq    -> let b, a = split_sequence seq ~pos in
                       Sequence b, Sequence a
  | Gap gap         -> let b, a = split_gap gap ~pos in
                       Gap b, Gap a

type 'a at_same_pos =
  | Fin
  | Fst of 'a alignment_element * 'a alignment_element list * 'a alignment_element list
  | Snd of 'a alignment_element * 'a alignment_element list * 'a alignment_element list
  | Both of 'a alignment_element * 'a alignment_element list * 'a alignment_element * 'a alignment_element list

let same_pos_and_length_step s1 s2 =
  match s1,   s2 with
  | [],       []       -> Fin
  | h :: t,   []       -> Fst (h, t, [])
  | [],       h :: t   -> Snd (h, [], t)
  | h1 :: t1, h2 :: t2 ->
      let sp1 = start_position h1 in
      let sp2 = start_position h2 in
      let ep1 = end_position String.length h1 in
      let ep2 = end_position String.length h2 in
      if sp1 < sp2 then
        if ep1 <= sp2 then        (* h1 is completely first! *)
          Fst  (h1, t1, s2)
        else                      (* intersect *)
          let b, a = split_long_al_els sp2 h1 in
          Fst  (b, a :: t1, s2)
      else if sp1 = sp2 then
        if ep1 < ep2 then
          if is_start h1 || is_end h1 then      (* prevent splitting because of 0-length starts.*)
            Fst (h1, t1, s2)
          else
            let b, a = split_long_al_els ep1 h2 in
            Both (h1, t1, b, a :: t2)
        else if ep1 = ep2 then
          Both (h1, t1, h2, t2)
        else (* ep1 > ep2 *)
          if is_start h2 || is_end h2 then
            Snd (h2, s1, t2)
          else
            let b, a = split_long_al_els ep2 h1 in
            Both (b, a :: t1, h2, t2)
      else (* sp1 > sp2 *)
        if ep2 <= sp1 then        (* h2 is completely first! *)
          Snd  (h2, s1, t2)
        else
          let b, a = split_long_al_els sp1 h2 in
          Snd (b, s1, a :: t2)

let same_position_fold ~init ~f s1 s2 =
  let rec loop acc r =
    match r with
    | Fin                 -> acc
    | Fst (e, s1, s2)     -> loop (f acc (`Fst e))       (same_pos_and_length_step s1 s2)
    | Snd (e, s1, s2)     -> loop (f acc (`Snd e))       (same_pos_and_length_step s1 s2)
    | Both (g, s1, h, s2) -> loop (f acc (`Both (g, h))) (same_pos_and_length_step s1 s2)
  in
  loop init (same_pos_and_length_step s1 s2)

let same_position_print =
  same_position_fold ~init:() ~f:(fun _ -> function
    | `Fst e       -> printf "Fst %s \n" (al_el_to_string e)
    | `Snd e       -> printf "Snd %s \n" (al_el_to_string e)
    | `Both (g, h) -> printf "Both %s %s \n" (al_el_to_string g) (al_el_to_string h))

(* an end position so compare for equality. *)
let same_pos_and_length_step_upto end_pos s1 s2 =
  let relationship = function
    | Start s     -> if s >= end_pos then `After else `Inside false
    | End e       -> if e >= end_pos then `After else `Inside false
    | Boundary b  -> if b.pos >= end_pos then `After else `Inside false
    | Gap g       -> if g.gstart >= end_pos then `After else
                      `Inside (g.gstart < end_pos && end_pos < g.gstart + g.length)
    | Sequence s  -> if s.start >= end_pos then `After else
                      `Inside (s.start < end_pos && end_pos < s.start + (String.length s.s))
  in
  let r = same_pos_and_length_step s1 s2 in
  match r with
  | Fin                   -> Some r
  | Fst (e, s1, s2)       -> begin match relationship e with
                             | `After        -> None
                             | `Inside false -> Some r
                             | `Inside true  -> let b, a = split_long_al_els end_pos e in
                                                Some (Fst (b, a :: s1, s2))
                             end
  | Snd (e, s1, s2)       -> begin match relationship e with
                             | `After        -> None
                             | `Inside false -> Some r
                             | `Inside true  -> let b, a = split_long_al_els end_pos e in
                                                Some (Snd (b, s1, a :: s2))
                             end
  | Both (e1, s1, e2, s2) -> begin match relationship e1 with
                             | `After        -> None
                             | `Inside false -> Some r
                             | `Inside true  -> let b1, a1 = split_long_al_els end_pos e1 in
                                                let b2, a2 = split_long_al_els end_pos e2 in
                                                Some (Both (b1, a1 :: s1, b2, a2 :: s2))
                             end

module Zip2 = struct

  type state =
    { started1    : bool
    ; started2    : bool
    ; is_seq1     : bool          (* to keep track of an allele's seq length. *)
    ; is_seq2     : bool          (* to keep track of an allele's seq length. *)
    ; is_boundary : bool
    ; start_pos   : position
    ; end_pos     : position
    ; mismatches  : int           (* to each other *)
    }

  let state_of_boundary ~started1 ~started2 b =
    { started1 ; started2 ; mismatches = 0
    ; is_seq1 = false; is_seq2 = false
    ; is_boundary = true
    ; start_pos = b.pos ; end_pos = b.pos + 1
    }

  let state_of_gap ?(in_ref=true) ~started1 ~started2 gap =
    let mismatches = if in_ref then 0 else gap.length in
    { started1 ; started2 ; mismatches
    ; is_seq1 = false; is_seq2 = false
    ; is_boundary = false
    ; start_pos = gap.gstart
    ; end_pos = gap.gstart + gap.length
    }

  let state_of_sequence ?(in_ref=true) ~started1 ~started2 seq =
    let n = String.length seq.s in
    let mismatches = if in_ref then 0 else n in
    { started1 ; started2 ; mismatches
    ; is_seq1 = true; is_seq2 = true
    ; is_boundary = false
    ; start_pos = seq.start
    ; end_pos = seq.start + String.length seq.s
    }

  let state_to_string s =
    sprintf "{%b; %b; %b; %b; %d; %d; %d}"
      s.started1 s.started2 s.is_seq1 s.is_seq2 s.start_pos s.end_pos s.mismatches

  let split ?new_start1 ?new_start2 cur pos =
    let new_start1 = Option.value ~default:cur.started1 new_start1 in
    let new_start2 = Option.value ~default:cur.started2 new_start2 in
    if pos = cur.start_pos then
      None, ( { cur with started1 = new_start1 ; started2 = new_start2 })
    else
      Some ( { cur with end_pos = pos }),
      ( { cur with started1 = new_start1
                 ; started2 = new_start2
                 ; start_pos = pos
                 ; mismatches = 0})

  let add_mismatches ?is_seq1 ?is_seq2 cur n =
    let is_seq1 = Option.value ~default:cur.is_seq1 is_seq1 in
    let is_seq2 = Option.value ~default:cur.is_seq2 is_seq2 in
    { cur with mismatches = cur.mismatches + n ; is_seq1 ; is_seq2 }

  let add_seq ?is_seq1 ?is_seq2 cur seq =
    add_mismatches ?is_seq1 ?is_seq2 cur (String.length seq.s)

  let add_seq_wrap_one cur seq = function
    | `Fst  -> add_mismatches ~is_seq1:true cur (String.length seq.s)
    | `Snd  -> add_mismatches ~is_seq2:true cur (String.length seq.s)

  let set_seq ~is_seq1 ~is_seq2 cur =
    { cur with is_seq1; is_seq2 }

  let add_gap ?is_seq1 ?is_seq2 cur gap =
    add_mismatches ?is_seq1 ?is_seq2 cur gap.length

  let add_gap_wrap_one cur gap = function
    | `Fst  -> add_mismatches ~is_seq1:false cur gap.length
    | `Snd  -> add_mismatches ~is_seq2:false cur gap.length

  let zip2 ~reference ~allele1 ~allele2 =
    let append acc state =
      if !report then printf "Appending %s \n" (state_to_string state);
      state :: acc
    in
    let rec split_wrap cur pos ?new_start1 ?new_start2 acc s1 s2 =
      let b_opt, a = split cur pos ?new_start1 ?new_start2 in
      let nacc = Option.value_map b_opt ~default:acc ~f:(append acc) in
      advance_allele a nacc s1 s2

    and split_wrap_only_one cur pos acc s1 s2 new_start = function
      | `Fst -> split_wrap cur pos ~new_start1:new_start acc s1 s2
      | `Snd -> split_wrap cur pos ~new_start2:new_start acc s1 s2

    and advance_allele cur acc a1 a2 =
      match same_pos_and_length_step_upto cur.end_pos a1 a2 with
      | None                  -> cur.started1, cur.started2, append acc cur, a1, a2
      | Some st               ->
          begin match st with
          | Fin                   -> cur.started1, cur.started2, append acc cur, a1, a2
          | Fst (e, s1, s2)       -> one_side cur acc s1 s2 `Fst e
          | Snd (e, s1, s2)       -> one_side cur acc s1 s2 `Snd e
          | Both (e1, s1, e2, s2) -> two_side cur acc s1 s2 e1 e2
          end

    and one_side cur acc s1 s2 which e =
      if !report then
        printf "one side at %s when cur is %s.\n"
          (al_el_to_string e) (state_to_string cur);
      if start_position e < cur.start_pos then
        advance_allele cur acc s1 s2 (* skip element *)
      (* The same_pos_and_length_step_upto logic prevents us form seeing
         elements after the current position. *)
      else
        match e with
        | Start pos       -> split_wrap_only_one cur pos acc s1 s2 true which
        | End pos         -> split_wrap_only_one cur pos acc s1 s2 false which
        | Boundary b as e ->
            if cur.start_pos <> b.pos && cur.end_pos <> cur.start_pos + 1 then
              invalid_argf "Found %s in sequence at %d!"
                (al_el_to_string e) (cur.start_pos)
            else (* ignore *)
              advance_allele cur acc s1 s2
        | Sequence seq    -> advance_allele (add_seq_wrap_one cur seq which) acc s1 s2
        | Gap gap         -> advance_allele (add_gap_wrap_one cur gap which) acc s1 s2

    and two_side cur acc s1 s2 e1 e2 =
      if !report then
        printf "two side at %s %s when cur is %s.\n"
          (al_el_to_string e1) (al_el_to_string e2) (state_to_string cur);
      if start_position e1 < cur.start_pos then
        advance_allele cur acc s1 s2
      else
        match e1, e2 with
        | Start s, Start _ -> split_wrap cur s acc s1 s2 ~new_start1:true ~new_start2:true
        | Start s, End _   -> split_wrap cur s acc s1 s2 ~new_start1:true ~new_start2:false
        | End e,   Start _ -> split_wrap cur e acc s1 s2 ~new_start1:false ~new_start2:true
        | End e,   End _   -> split_wrap cur e acc s1 s2 ~new_start1:false ~new_start2:false
        | Start _, _
        | End _,   _
        | _,       Start _
        | _,       End _   -> invalid_argf "%s paired with non-zero-length al-el %s"
                                  (al_el_to_string e1) (al_el_to_string e2)

        | Boundary b, Boundary _ -> (* make sure we're at a boundary in the ref *)
            if cur.start_pos <> b.pos && cur.end_pos <> cur.start_pos + 1 then
              invalid_argf "Boundaries don't align %s %s in sequence at %d!"
                (al_el_to_string e1) (al_el_to_string e2) (cur.start_pos)
            else (* ignore *)
              advance_allele cur acc s1 s2
        | Boundary _, _
        | _,          Boundary _ -> invalid_argf "Boundaries not aligned %s %s in sequence at %d!"
                                      (al_el_to_string e1) (al_el_to_string e2) (cur.start_pos)

        (* Same *)
        | Gap _,       Gap _                        -> advance_allele (set_seq cur ~is_seq1:false ~is_seq2:false) acc s1 s2
        | Sequence q1, Sequence q2 when q1.s = q2.s -> advance_allele (set_seq cur ~is_seq1:true ~is_seq2:true) acc s1 s2

        (* Diff *)
        | Sequence s,  Sequence _  (*q1.s <> q2.s*) -> advance_allele (add_seq cur s ~is_seq1:true ~is_seq2:true) acc s1 s2
        | Sequence _,  Gap g                        -> advance_allele (add_gap cur g ~is_seq1:true ~is_seq2:false) acc s1 s2
        | Gap g,       Sequence _                   -> advance_allele (add_gap cur g ~is_seq1:false ~is_seq2:true) acc s1 s2
    in
    let reference_pass =
      Boundaries.fold reference
        ~start:(fun state _ -> state)          (* Nothing to do for reference Start's *)
        ~end_:(fun state _ -> state)             (* Nothing to do for reference End's *)
        ~boundary:(function
          | None                                             ->
              (false, false, [], allele1, allele2)
          | Some ((started1, started2, acc, a1, a2), b, _bm) ->
              (* reset acc to [] *)
              let new_state = state_of_boundary ~started1 ~started2 b in
              advance_allele new_state [] a1 a2)
        ~gap:(fun (started1, started2, acc, a1, a2) gap ->
              let new_state = state_of_gap ~started1 ~started2 gap in
              advance_allele new_state acc a1 a2)
        ~seq:(fun (started1, started2, acc, a1, a2) seq ->
              let new_state = state_of_sequence ~started1 ~started2 seq in
              advance_allele new_state acc a1 a2)
  in
  let just_bm_and_state = List.rev_map ~f:(fun (bm, (_, _, p, _, _)) -> (bm, p)) in
  match reference_pass with
  | [] -> invalid_argf "Didn't even have a start boundary for zip2?"
  | (bm, (fs1, fs2, acc, a1, a2)) :: tl ->
      let final1, final2, nacc =
        same_position_fold a1 a2 ~init:(fs1, fs2, acc)
          ~f:(fun (started1, started2, acc) ns ->
              match ns with
              | `Fst e ->
                begin match e with
                | Start _    -> (true, started2, acc)
                | End _      -> (false, started2, acc)
                | Boundary b -> let ns = state_of_boundary ~started1 ~started2 b in
                                (started1, started2, ns :: acc)
                | Gap g      -> let ns = state_of_gap ~started1 ~started2 ~in_ref:false g in
                                (started1, started2, ns :: acc)
                | Sequence s -> let ns = state_of_sequence ~started1 ~started2 ~in_ref:false s in
                                (started1, started2, ns :: acc)
                end
              | `Snd e ->
                begin match e with
                | Start _    -> (started1, true, acc)
                | End _      -> (started1, false, acc)
                | Boundary b -> let ns = state_of_boundary ~started1 ~started2 b in
                                (started1, started2, ns :: acc)
                | Gap g      -> let ns = state_of_gap ~started1 ~started2 ~in_ref:false g in
                                (started1, started2, ns :: acc)
                | Sequence s -> let ns = state_of_sequence ~started1 ~started2 ~in_ref:false s in
                                (started1, started2, ns :: acc)
                end
              | `Both (e1, e2) ->
                begin match e1, e2 with
                | Start _, Start _ -> (true,  true,  acc)
                | Start _, End _   -> (true,  false, acc)
                | End _,   Start _ -> (false, true,  acc)
                | End _,   End _   -> (false, false, acc)
                | Start _, _
                | End _,   _
                | _,       Start _
                | _,       End _   ->
                    invalid_argf "%s paired with non-zero-length al-el %s, past reference"
                      (al_el_to_string e1) (al_el_to_string e2)

                | Boundary b, Boundary _  ->
                      let ns = state_of_boundary ~started1 ~started2 b in
                      (started1, started2, ns :: acc)
                | Boundary _, _
                | _,          Boundary _  ->
                    invalid_argf "Boundaries not aligned %s %s in sequence, past reference!"
                      (al_el_to_string e1) (al_el_to_string e2)

                (* Same *)
                | Gap _,       Gap _                        -> (started1, started2, acc)
                | Sequence q1, Sequence q2 when q1.s = q2.s -> (started1, started2, acc)

                (* Diff *)
                | Sequence s,  Sequence _ ->
                    let ns = state_of_sequence ~started1 ~started2 ~in_ref:false s in
                    (started1, started2, ns :: acc)
                | Sequence _,  Gap g      ->
                    let ns = state_of_gap ~started1 ~started2 ~in_ref:false g in
                    (started1, started2, ns :: acc)
                | Gap g,       Sequence _ ->
                    let ns = state_of_gap ~started1 ~started2 ~in_ref:false g in
                    (started1, started2, ns :: acc)
                end)
      in
      just_bm_and_state ((bm, (final1, final2, nacc, [], [])) :: tl)

end (* Zip2 *)

let allele_distances_between ~reference ~allele1 ~allele2 =
  let open Zip2 in
  let update_state_of_one { is_boundary; end_pos; start_pos; _} is_seq p =
    if is_boundary then (* ignore boundaries. *)
      fst p
    else
      match p with
      | `Missing,  true   -> `Start (if is_seq then end_pos - start_pos else 0)
      | `Missing,  false  -> `NoStart
      | `Start sl, true   -> `Start (if is_seq then sl + end_pos - start_pos else sl)
      | `Start sl, false  -> `Partial sl (* check that end = start_pos or len 1? *)
      | `NoStart, true    -> `Partial (if is_seq then end_pos - start_pos else 0)
      | `NoStart, false   -> `NoStart
      | `Partial p, true  -> `Partial (if is_seq then p + end_pos - start_pos else p)
      | `Partial p, false -> `Partial p
  in
  Zip2.zip2 ~reference ~allele1 ~allele2
  |> List.map ~f:(fun (bm, slst) ->
      let (fs1, fs2, mismatches) =
        List.fold_left slst ~init:(`Missing, `Missing, 0)
          ~f:(fun (st1, st2, m) zs ->
                let nst1 = update_state_of_one zs zs.is_seq1 (st1, zs.started1) in
                let nst2 = update_state_of_one zs zs.is_seq2 (st2, zs.started2) in
                nst1, nst2, m + zs.mismatches)
      in
      { reference_seq_length  = bm.Boundaries.seq_length
      ; mismatches
      ; allele_relationships =
        match (fs1, fs2) with
        | (`Missing, `Missing)
        | (`Missing, `NoStart)
        | (`NoStart, `Missing)
        | (`NoStart, `NoStart)      -> Missing, Missing

        | (`Missing, `Start l)      -> Missing, Full l
        | (`Missing, `Partial l)    -> Missing, Partial l
        | (`NoStart, `Start l)      -> Missing, Full l
        | (`NoStart, `Partial l)    -> Missing, Partial l

        | (`Start l, `Missing)      -> Full l, Missing
        | (`Start l, `NoStart)      -> Full l, Missing
        | (`Start l, `Start s)      -> Full l, Full s
        | (`Start l, `Partial p)    -> Full l, Partial p

        | (`Partial l, `Missing)    -> Partial l, Missing
        | (`Partial l, `NoStart)    -> Partial l, Missing
        | (`Partial l, `Start s)    -> Partial l, Full s
        | (`Partial l, `Partial p)  -> Partial l, Partial p
      })
