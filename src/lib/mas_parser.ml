(* How this works:

This is a line based format. For our purposes, there are three types of lines
as determined by {type line}, where `SeqData` has most of the actual alignment
data. We divide the parsing of sequences into two types: reference and
alternatives. We take special care of keeping track of gaps in the reference so
that positions, when parsing the alternatives, are annotated with regard to
the reference positions.

We return
*)

open Util

(* This refers to the alignment position.  *)
type position = int

(* We'll parse N as a nucleotide. *)
let is_nucleotide = function 'A' | 'C' | 'G' | 'T' | 'N' -> true | _ -> false
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
  | Start p               -> p
  | End p                 -> p
  | Boundary { pos; _ }   -> pos
  | Sequence { start; _ } -> start
  | Gap { gstart; _ }     -> gstart

let end_position sr_length = function
  | Start s     -> s
  | End e       -> e
  | Gap g       -> (g.gstart + g.length)
  | Sequence s  -> (s.start + sr_length s.s)
  | Boundary b  -> b.pos

let al_el_to_string = function
  | Start p                 -> sprintf "Start at %d" p
  | End p                   -> sprintf "End at %d" p
  | Boundary { idx; pos }   -> sprintf "Boundary %d at %d" idx pos
  | Sequence { start; s }   -> sprintf "Sequence %s at %d" s start
  | Gap { gstart; length }  -> sprintf "Gap of %d from %d" length gstart

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

let update ~dna ~fail_on_same ps s =
  let is_nuc = if dna then is_nucleotide else is_amino_acid in
  let rec to_ref_seq_elems_char ps = function
    | []                    -> ps
    | '|' :: t              -> to_ref_seq_elems_char (insert_boundary ps) t
    | '*' :: t              -> to_ref_seq_elems_char (insert_unknown ps) t
    | 'X' :: t when not dna -> to_ref_seq_elems_char (insert_unknown ps) t (* add End *)
    | '.' :: t              -> to_ref_seq_elems_char (insert_gap ps) t
    | '-' :: t              -> to_ref_seq_elems_char (insert_same ~fail_on_same ps) t
    | c :: t when is_nuc c  -> to_ref_seq_elems_char (insert_nuc c ps) t
    | x :: _                -> invalid_argf "Unrecognized char %c in %s" x (where ps)
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
      | s :: _ when s.[0] = Some '|' -> Dash
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
      if x.ref = allele then begin
        (*let prev_pos = x.ref_ps.position in *)
        let nref_ps =
          List.fold_left ~f:(update ~dna:x.dna ~fail_on_same:true) ~init:x.ref_ps s
        in
        latest_reference_position := nref_ps.position;
        { x with ref       = allele
               ; ref_ps    = nref_ps
        }
      end else begin
        let cur_ps =
          try Hashtbl.find x.alg_htbl allele
          with Not_found -> init_ps allele x.start_pos
        in
        let new_ps =
          List.fold_left ~f:(update ~dna:x.dna ~fail_on_same:false) ~init:cur_ps s
        in
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
      Hashtbl.fold (fun all ps acc ->
          if ps.sequence = [] then begin
            printf "Dropping empty sequence: %s\n" ps.allele;
            acc
          end else
            (all, normalized_seq ps) :: acc)
        reversed.alg_htbl []
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
          List.rev ((il false curb 1, curs) :: acc)
      | Boundary bp :: t  ->
          let newb = marker ~index:bp.idx ~position:bp.pos in
          let nacc = ((il false curb 1, curs) :: acc) in
          loop newb                                (boundary newb)       nacc  t
      | Start s :: t      ->
          let newb = if before_start curb then { curb with position = s - 1 } else curb in
          loop newb                                (start curs s) acc   t
      | End e :: t        ->
          loop curb                                (end_ curs e) acc   t
      | Gap g :: t        ->
          loop (il false curb g.length)            (gap curs g) acc   t
      | Sequence s :: t   ->
          loop (il true curb (String.length s.s))  (seq curs s) acc   t
    in
    loop before_start_boundary (boundary before_start_boundary) [] lst

  let bounded lst =
    let boundary _ = "" in
    let start s _ = s in
    let end_ s _ = s in
    let gap s _ = s in
    let seq s ss = s ^ ss.s in
    fold lst ~boundary ~start ~end_ ~gap ~seq

end (* Boundaries *)

(* Applying the alignment elements. *)
module type Data_projection = sig

  type t

  val to_string : t -> string

  (* Initialize the projection. *)
  val of_seq : in_ref:bool -> bool -> string sequence -> t

  val of_gap : in_ref:bool -> bool -> gap -> t

  (* Mutate the projection. *)
  val add_seq : position -> t -> string sequence -> t

  val add_gap : position -> t -> gap -> t

  (* Use [int] instead of [position] because we're passing where in
     the current segment the allele starts/stops. *)
  val add_start : t -> int -> t

  val add_end : t -> int -> t

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

  let left st pos =
    { st with project = R.add_end st.project (pos - st.start_pos)
            ; started = false
            ; end_pos = pos
    }

  let right st pos =
    let index = pos - st.start_pos in
    { st with project = R.add_start st.project index
            ; started = true
            ; start_pos = st.start_pos + index
    }

  let new_cur started ~n ~b ~s = function
    | Start _
    | End _       -> n ()
    | Boundary bb -> b bb
    | Gap g       -> s { start_pos = g.gstart
                       ; end_pos = g.gstart + g.length
                       ; project = R.of_gap ~in_ref:true started g
                       ; started
                       }
    | Sequence ss -> let len = String.length ss.s in
                     s { start_pos = ss.start
                       ; end_pos = ss.start + len
                       ; project = R.of_seq ~in_ref:true started ss
                       ; started
                       }

  let apply ~reference ~allele =
    let rec append acc project =
      if !report then printf "Appending %s \n" (R.to_string project);
      (`P project) :: acc
    and new_ref_segment started acc r a =
      match r with
      | []     ->
          let _final_started, nacc =
            List.fold_left a ~init:(started, acc) ~f:(fun (st, a) e ->
              match e with
              | End _       -> false, a
              | Start _     -> true, a
              | Boundary bb -> st, `B bb :: a
              | Gap g       -> st, append a (R.of_gap ~in_ref:false st g)
              | Sequence s  -> st, append a (R.of_seq ~in_ref:false st s))
          in
          if !report then printf "Ending after no more reference: %d!\n" (List.length nacc);
          List.rev nacc (* Fin. *)
      | h :: t ->
          new_cur started h
            ~n:(fun () -> new_ref_segment started acc t a)
            ~b:(fun bb -> new_ref_segment started (`B bb :: acc) t a)
            ~s:(fun cur -> loop cur acc t a)
    and allele_ended cur acc r =
      if !report then printf "Ending after no more allele! %d\n" (List.length acc);
      let init = append acc cur.project in
      let acc =
        List.fold_left r ~init ~f:(fun a e ->
          match e with
          | End _      -> a
          | Start _    -> a
          | Boundary b -> `B b :: a
          | Gap g      -> append a (R.of_gap ~in_ref:true cur.started g)
          | Sequence s -> append a (R.of_seq ~in_ref:true cur.started s))
      in
      List.rev acc (* Fin *)
    and mutate_segment cur acc r at = function
      | Start p                         -> loop (right cur p) acc r at
      | End p                           -> let nst = left cur p in
                                           (* What if we restart in the same segment? *)
                                           new_ref_segment false (append acc nst.project) r at
      | Boundary _ as e                 -> invalid_argf "Found %s in sequence at %d"
                                            (al_el_to_string e) (cur.start_pos)
      | Sequence ss                     ->
          let slen = String.length ss.s in
          if ss.start + slen > cur.end_pos then
            let index = cur.end_pos - ss.start in
            let na =
              Sequence { start = ss.start + index
                       ; s = String.drop ss.s ~index
                       }
            in
            let ncur = add_seq cur { ss with s = String.take ss.s ~index } in
            loop ncur acc r (na :: at)
          else
            let ncur = add_seq cur ss in
            loop ncur acc r at
      | Gap g                           ->
          if g.gstart + g.length > cur.end_pos then
            let length = cur.end_pos - g.gstart in
            let na =
              Gap { gstart = g.gstart + length
                  ; length = g.length - length
                  }
            in
            let ncur = add_gap cur { g with length } in
            loop ncur acc r (na :: at)
          else
            let ncur = add_gap cur g in
            loop ncur acc r at
    and at_end_of_segment cur acc r at = function
      | Boundary _ -> loop cur acc r at   (* Add bs based on reference. *)
      | Start _    -> new_ref_segment true (append acc cur.project) r at
      | End _      -> new_ref_segment false (append acc cur.project) r at
      | g_or_s     -> new_ref_segment cur.started (append acc cur.project) r (g_or_s :: at)
    and loop cur acc r a =
      match a with
      | []        -> allele_ended cur acc r
      | ah :: at  ->
          let sa = start_position ah in
          if !report then
            printf "at %s %d when cur is %s.\n"
              (al_el_to_string ah) sa (state_to_string cur);
          if sa < cur.start_pos then
            loop cur acc r at
          else if sa > cur.end_pos then
            new_ref_segment cur.started (append acc cur.project) r a
          else if sa = cur.end_pos then
            at_end_of_segment cur acc r at ah
          else
            mutate_segment cur acc r at ah
    in
    new_ref_segment false [] reference allele

end (* MakeZip *)

module B = Sosa.Native_bytes

let gap_char = '.'

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
    started, B.make g.length gap_char

  let add_seq start_pos (started, buffer) {start; s} =
    let off = start - start_pos in
    for index = 0 to String.length s - 1 do
      B.mutate_exn buffer ~index:(index + off) (String.get_exn s ~index)
    done;
    (started, buffer)

  let add_gap start_pos (stared, buffer) {gstart; length} =
    let off = gstart - start_pos in
    for i = 0 to length - 1 do
      B.mutate_exn buffer ~index:(i + off) gap_char
    done;
    (stared, buffer)

  let add_start (_, buffer) index =
    true, B.drop buffer ~index

  (* return true because we're keeping the 'in allele' part
    TODO: the start / end logic is once again awkward. It seems like add_end
    should return a pair of in/out sequence segments. *)
  let add_end (_, buffer) index =
    true, B.take buffer ~index

end)

let allele_sequence ?boundary_char ~reference ~allele () =
  let ss = AlleleSequences.apply ~reference ~allele in
  let nacc =
    match boundary_char with
    | None   ->
        List.filter_map ss ~f:(function
          | `P (true, t)  -> Some (B.filter ~f:(fun c -> c <> gap_char) t)
          | `P (false, _) -> None
          | `B _          -> None)
    | Some c ->
        let bs = String.of_character c in
        List.filter_map ss ~f:(function
          | `P (true, t)  -> Some (B.filter ~f:(fun c -> c <> gap_char) t)
          | `P (false, _) -> None
          | `B _          -> Some bs)
  in
  B.(to_native_string (concat nacc))

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

  type t =
    | G of { started : bool; gap : gap; mismatches : int }
    | S of { started : bool; seq : string sequence; mismatches : int }

  let length = function
    | G r -> r.gap.length
    | S r -> String.length r.seq.s

  let mismatches = function
    | G r -> r.mismatches
    | S r -> r.mismatches

end

module AlleleToRefDistance = MakeZip (struct

  open DistanceProjection
  (* what is in this part and # of mismatches *)
  type t = DistanceProjection.t

  let to_string = function
    | G {started; mismatches; _} -> sprintf "G %b %d" started mismatches
    | S {started; mismatches; _} -> sprintf "S %b %d" started mismatches

  let incr m = function
    | S r -> S { r with mismatches = r.mismatches + m }
    | G r -> G { r with mismatches = r.mismatches + m }

  let of_seq ~in_ref started seq =
    if in_ref then S { started; seq; mismatches = 0 }
    else S { started; seq; mismatches = String.length seq.s }

  let of_gap ~in_ref started gap =
    if in_ref then G { started; gap; mismatches = 0 }
    else G { started; gap; mismatches = gap.length }

  (* zipping logic makes sure position is inside s *)
  let add_seq _ t {s; _} =
    incr (String.length s) t

  (* seqs in seqs are diffs! *)
  let add_gap _ t {length; _} =
    match t with
    | S r -> S { r with mismatches = r.mismatches + length } (* filling a gap *)
    | G _ -> t                                  (* gap in gap, no difference. *)

  let add_start t msm =
    incr msm t

  let add_end t len =
    incr (DistanceProjection.length t - len) t

end)

let allele_distances ~reference ~allele =
  let rec loop sum acc = function
    | []        -> List.rev (sum :: acc)
    | `B _ :: t -> loop 0 (sum :: acc) t
    | `P d :: t -> loop (sum + DistanceProjection.mismatches d) acc t
  in
  loop 0 [] (AlleleToRefDistance.apply ~reference ~allele)

