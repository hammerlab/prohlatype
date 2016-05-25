(* How this works:

This is a line based format. For our purposes, there are three types of lines
as determined by {type line}, where `SeqData` has most of the actual alignment
data. We divide the parsing of sequences into two types: reference and
alternatives. We take special care of keeping track of gaps in the reference so
that positions, when parsing the alternatives, are annotated with regard to
the reference positions.

We return
*)

open Printf
module String = Sosa.Native_string

let invalid_argf ?(prefix="") fmt = ksprintf invalid_arg ("%s " ^^ fmt) prefix

(* This refers to the alignment position.  *)
type position = int

let is_nucleotide = function 'A' | 'C' | 'G' | 'T' -> true | _ -> false
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
  }

let init_ps allele position =
  { allele
  ; position
  ; boundary = 0
  ; sequence = []
  }

let where ps =
  sprintf "allele: %s, position: %d, sequence length: %d"
    ps.allele ps.position (List.length ps.sequence)

type 'sr alignment_element =
  | Boundary of { idx : int; pos : position }
  | Sequence of { start : position; s : 'sr }
  | Gap of { start : position; length : int }
  | Unknown of { start : position; length : int }

let to_string srts = function
  | Boundary { idx; pos }     -> sprintf "Boundary %d at %d" idx pos
  | Sequence { start; s }     -> sprintf "Sequence %s at %d" (srts s) start 
  | Gap { start; length }     -> sprintf "Gap of %d from %d" length start
  | Unknown {start; length }  -> sprintf "Unknown of %d from %d" length start

let next ps = { ps with position = ps.position + 1 }
let with_next ps f = f (next ps)

let insert_boundary ps =
  with_next ps (fun ps ->
    { ps with sequence = Boundary { idx = ps.boundary; pos = ps.position} :: ps.sequence
            ; boundary = ps.boundary + 1
    })

let insert_gap ps =
  with_next ps (fun ps ->
    match ps.sequence with
    | Gap { start; length } :: t when start + length = ps.position
                            -> { ps with sequence = Gap { start; length = length + 1 } :: t }
    | []
    | Boundary _ :: _
    | Gap _ :: _
    | Unknown _ :: _
    | Sequence _ :: _ as l  -> { ps with sequence = Gap { start = ps.position; length = 1 } :: l })

let insert_unknown ps =
  with_next ps (fun ps ->
    match ps.sequence with
    | Unknown { start; length } :: t when start + length = ps.position
                            -> { ps with sequence = Unknown { start; length = length + 1 } :: t }
    | []
    | Boundary _ :: _
    | Gap _ :: _
    | Unknown _ :: _
    | Sequence _ :: _ as l  -> { ps with sequence = Unknown { start = ps.position; length = 1 } :: l })

let insert_nuc_error fmt =
  invalid_argf ~prefix:"Trying to insert sequence element" fmt

let insert_same ~fail_on_same ps =
  if fail_on_same then
    invalid_argf "Encountered unexpected '-' same char for : %s" (where ps)
  else
    next ps

let insert_nuc c ps =
  with_next ps (fun ps ->
    match ps.sequence with
    | Sequence {start; s} :: t when start + (List.length s) = ps.position
                            -> { ps with sequence = Sequence { start; s = c :: s} :: t }
    | []
    | Boundary _ :: _
    | Gap _ :: _
    | Unknown _ :: _
    | Sequence _ :: _ as l  -> { ps with sequence = Sequence { start = ps.position; s = c :: []} :: l })

let update ~dna ~fail_on_same ps s =
  let is_nuc = if dna then is_nucleotide else is_amino_acid in
  let rec to_ref_seq_elems_char ps = function
    | []                    -> ps
    | '|' :: t              -> to_ref_seq_elems_char (insert_boundary ps) t
    | '*' :: t              -> to_ref_seq_elems_char (insert_unknown ps) t
    | 'X' :: t when not dna -> ps (* don't add anything else *)
    | '.' :: t              -> to_ref_seq_elems_char (insert_gap ps) t
    | '-' :: t              -> to_ref_seq_elems_char (insert_same ~fail_on_same ps) t
    | c :: t when is_nuc c  -> to_ref_seq_elems_char (insert_nuc c ps) t
    | x :: _                -> invalid_argf "Unrecognized char %c in %s" x (where ps)
  in
  to_ref_seq_elems_char ps (String.to_character_list s)


(*
type 'sr sequence_element =
  | Start of int * string
  | End of int
  | Boundary of (int * int)     (* idx and pos *)
  | Nuc of (int * 'sr * int)
  | Gap of (int * int)          (* pos and length *)
let sequence_element_to_string_g ~sr_to_string = function
  | Start (p, allele) -> sprintf "Start %s at %d" allele p
  | End p             -> sprintf "End %d" p
  | Boundary (i, p)   -> sprintf "Boundary %d at %d" i p
  | Nuc (p, s, o)     -> sprintf "Nucleotide %s at %d offset end by %d" (sr_to_string s) p o
  | Gap (p, l)        -> sprintf "Gap of length %d starting at %d" l p

let sequence_element_to_string =
  sequence_element_to_string_g ~sr_to_string:(fun x -> x)
   *)


(*let position = function
  | Start (p, _)
  | End p
  | Boundary (_, p)
  | Nuc (p, _, _)
  | Gap (p, _) -> p *)
(*
let increment_position ps = { ps with position = ps.position + 1 }

let new_status ps datum =
  match ps.in_data, ps.sequence, datum with
  | Before, [], true  -> { ps with in_data = In ; sequence = Start (ps.position, ps.allele)}
  | Before, ls, true  -> invalid_argf "Trying to start when we already have data: %s %d" ps.allele (List.length ls)
  | Before,  _, false -> ps
  | In,      _, false -> { ps with in_data = After ; sequence = End ps.position :: ps.sequence }
  | In,      _, true  -> ps
  | After,   _, true  -> invalid_argf "Adding a datum to a closed sequence: %s %d" ps.allele ps.position
  | After,   _, false -> ps

(*
let update_ps ?incr ~datum insert_seq ps =
  let ps, closed = new_status ps datum in
  if closed then ps else
  match incr with
  | None -> { ps with sequence = insert_seq ps.sequence }
  | Some `Pos ->
    begin
      match ps.gaps with
      | [] ->
          { ps with sequence = insert_seq ps.sequence; position = ps.position + 1 }
      | (p0, gl) :: tl ->
        if ps.position < p0 then
          { ps with sequence = insert_seq ps.sequence; position = ps.position + 1 }
        else
          { ps with sequence = insert_seq ps.sequence
                  ; position = ps.position        (* do NOT increase position! *)
                  ; gaps     = if gl = 1 then tl else (p0, gl - 1) :: tl }
    end
  | Some `Bnd ->
    { ps with sequence = insert_seq ps.sequence; boundary = ps.boundary + 1 }
       *)

(*let update_gaps_in_ps new_gaps ps =
  let adj_gaps =
    (* Adjust gaps that might be spanning the current break boundary. *)
    match ps.sequence with
    | Gap (pos, len) :: _ ->
      List.map (fun ((p, l) as g) -> if p = pos then (p, l - len) else g) new_gaps
    | _ -> new_gaps
  in
  { ps with gaps = adj_gaps } *)

(*** Inserting Sequences ***)
let to_forward_pos p len = function
  | `Ref                            -> p + len
  | `Alt ((gp, _) :: _) when gp = p -> p
  | `Alt _                          -> p + len

let insert_nuc pos ~seq c = function
  | End _ :: _              -> invalid_argf "Trying to insert nucleotide %c at %d past end" c pos
  | []
  | Start _ :: _
  | Boundary _ :: _
  | Gap _ :: _ as l         -> Nuc (pos, c :: [], 0) :: l
  | (Nuc (p, ss, o) :: t) as l ->
      let len = List.length ss in
      let forward_pos = to_forward_pos p len seq in
      if forward_pos = pos then
        Nuc (p, c :: ss, o) :: t
      (* there is a gap in the reference where we have a alternate sequence! *)
      else if forward_pos > pos then
        Nuc (p, c :: ss, o + (forward_pos - pos)) :: t
      else
        Nuc (pos, c :: [], 0) :: l

let gaps_to_string gps =
  String.concat ~sep:";" (List.map (fun (p,l) -> sprintf "(%d,%d)" p l) gps)

let output_debug action state ps seq =
  printf "%s %s, p: %d, b %d, %s, forward_pos %d head %s\n"
    action state ps.position ps.boundary
    (gaps_to_string ps.gaps)
    (match ps.sequence with
       | Nuc (p, ss, _) :: _ -> let len = List.length ss in to_forward_pos p len seq
       | _                   -> -10000)
    (sequence_element_to_string_g ~sr_to_string:(String.of_character_list) (List.hd ps.sequence))

let insert_nuc_s ~seq c ps =
  (*let () = if ps.allele = "DRB1*01:01:01" then output_debug "inserting nuc" "before" ps seq in *)
  let r = update_ps ~incr:`Pos (insert_nuc ~seq ps.position c) ps in
  (*let () = if r.allele = "DRB1*01:01:01" then output_debug "inserting nuc" "after" r seq in *)
  r

(*** Inserting Gaps ***)
let to_forward_pos2 p len = function
  | `Ref                            -> p
  | `Alt ((gp, _) :: _) when gp = p -> p
  | `Alt _                          -> p + len

let debug_gap_ref = ref false

let insert_gap in_reference_gap pos = function
  | End _ :: _               -> invalid_argf "Trying to insert gap at %d past end" pos
  | []
  | Start _ :: _
  | (Boundary _ :: _)
  | (Nuc _ :: _) as l        -> Gap (pos, 1) :: l
  | (Gap (p, len) :: t) as l ->
      if in_reference_gap then
        if p = pos then
          Gap (p, len + 1) :: t
        else
          Gap (pos
      else
        Gap (p, len + 1) :: t

(* let forward_pos = to_forward_pos2 p len seq in
      if forward_pos = pos then begin
        if !debug_gap_ref then Printf.printf "forward pos %d equals pos %d p: %d len %d \n" forward_pos pos p len;
        Gap (p, len + 1) :: t
      end else if p + len > pos then  begin (* Our gap spans the reference gap *)
        if !debug_gap_ref then Printf.printf "p: %d + len %d greater than pos %d \n" p len pos;
        l
      end else begin
        if !debug_gap_ref then Printf.printf "tail case\n";
        Gap (pos, 1) :: l
      end *)

let insert_gap_s ?incr ~seq ps =
  let () = if ps.position > 2345 && ps.position < 2355 && ps.allele = "DRB1*01:01:01" then begin output_debug "inserting gap" "before" ps seq; debug_gap_ref := true end in
  let r = update_ps ?incr (insert_gap ~seq ps.position) ps in
  let () = if ps.position > 2345 && ps.position < 2355 && r.allele = "DRB1*01:01:01" then begin output_debug "inserting gap" "after" r seq; debug_gap_ref := false end in
  r

let insert_missing ps = update_ps ~datum:false ~incr:`Pos (fun l -> l) ps

let insert_boundary ps =
  update_ps ~incr:`Bnd
    (fun s -> Boundary (ps.boundary + 1, ps.position) :: s ) ps

let to_ref_seq_elems dna ps s =
  (*let () = output_debug "ref seq" s ps in*)
  let is_nuc = if dna then is_nucleotide else is_amino_acid in
  (* do not increase the position based on gaps in reference sequence.*)
  let rec to_ref_seq_elems_char ps = function
    | []                          -> ps
    | b :: t when is_boundary b   -> to_ref_seq_elems_char (insert_boundary ps) t
    | m :: t when is_missing m    -> to_ref_seq_elems_char (insert_missing ps) t
    | g :: t when is_gap g        -> to_ref_seq_elems_char (insert_gap_s ~seq:`Ref ps) t
    | c :: t when is_nuc c        -> to_ref_seq_elems_char (insert_nuc_s ~seq:`Ref c ps) t
    | x :: _                      -> invalid_argf "unrecognized char in reference ps: %c" x
  in
  to_ref_seq_elems_char ps (String.to_character_list s)

(* Sequence has not been reversed! *)
let find_reference_gaps until seq =
  let rec loop acc = function
    | Gap (p, l) :: t when p >= until   -> loop ((p,l) :: acc) t
    | h :: t when (position h) >= until -> loop acc t
    | _ -> acc
  in
  loop [] seq

let update_seq dna ps s =
  (*let () =
    if ps.allele = "A*68:02:02" then
      output_debug "update_seq" s ps
  in *)
  let is_nuc = if dna then is_nucleotide else is_amino_acid in
  let rec update_seq_char ps = function
    | []                          -> ps
    | b :: t when is_boundary b   -> update_seq_char (insert_boundary ps) t
    | m :: t when is_missing m    -> update_seq_char (insert_missing ps) t
    | g :: t when is_gap g        -> update_seq_char (insert_gap_s ~seq:(`Alt ps.gaps) ~incr:`Pos ps) t
    | c :: t when is_nuc c        -> update_seq_char (insert_nuc_s ~seq:(`Alt ps.gaps) c ps) t
    | s :: t when is_same s       -> update_seq_char (update_ps ~incr:`Pos (fun l -> l) ps) t
    | 'X' :: _ when (not dna)     -> ps (* signals an 'end' for AA's *)
    | x :: _                      -> invalid_argf "unrecognized char in seq: %c" x
  in
  update_seq_char ps (String.to_character_list s)
*)

type 'sr parse_result =
  { dna       : bool        (* DNA or Amino Acid sequence -> diff characters *)
  ; start_pos : int
  ; ref       : string      (* Name of reference. *)
  (* As we parse the alternative tracks, we have to Keep track of the gaps that
     we encounter in the reference, so that all positions are with respect to
     the reference. *)
  ; ref_ps    : 'sr parse_struct
  ; alg_htbl  : (string, 'sr parse_struct) Hashtbl.t
  }

let empty_result ref_allele dna position =
  { dna       = dna
  ; start_pos = position
  ; ref       = ref_allele
  ; ref_ps    = init_ps ref_allele position
  ; alg_htbl  = Hashtbl.create 100
  }

let reverse_seq lst =
  let to_string l = String.of_character_list (List.rev l) in 
  List.rev lst |> List.map (function
    (*| Start _ 
    | End _  *)
    | Boundary _
    | Unknown _
    | Gap _ as e            -> e
    | Sequence { start; s } -> Sequence { start; s = to_string s })

(*let normalize_seq ps =
  match ps.sequence with
  | End _ :: _ -> reverse_seq ps.sequence
  | _ ->          reverse_seq (End ps.position :: ps.sequence) *)

type line =
  | Position of bool * int  (* nucleotide or amino acid sequence  *)
  | Dash
  | SeqData of string * string list

(* Assume that it has been trimmed. *)
let parse_data line =
  String.split line ~on:(`Character ' ')
  |> List.filter ((<>) String.empty)
  |> function
      | "|" :: _                 -> Dash
      | "AA" :: "codon" :: _    -> Dash (* not really but not modeling this at the moment. *)
      | "gDNA" :: pos :: _      -> Position (true, int_of_string pos)
      | "cDNA" :: pos :: _      -> Position (true, int_of_string pos)
      | "Prot" :: pos :: _      -> Position (false, int_of_string pos)
      | []                      -> invalid_arg "Empty data line!"
      | s :: lst                -> SeqData (s, lst)

type parse_state =
  | Header
  | Empty
  | Data of line

type result =
  { reference : string
  ; ref_elems : string alignment_element list
  ; alt_elems : (string * string alignment_element list) list
  }

let report = ref false

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
        let nref_ps = List.fold_left (update ~dna:x.dna ~fail_on_same:true) x.ref_ps s in
        latest_reference_position := nref_ps.position;
        { x with ref       = allele
               ; ref_ps    = nref_ps
        }
      end else begin
        let cur_ps =
          try Hashtbl.find x.alg_htbl allele
          with Not_found -> init_ps allele x.start_pos
        in
        let new_ps = List.fold_left (update ~dna:x.dna ~fail_on_same:false) cur_ps s in
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
  let reversed = loop_header Header in
  let ref_elems = reverse_seq reversed.ref_ps.sequence in
  let alt_elems =
    Hashtbl.fold (fun all ps acc -> (all, reverse_seq ps.sequence) :: acc)
      reversed.alg_htbl []
  in
  { reference = reversed.ref ; ref_elems ; alt_elems }

let from_file f =
  let ic = open_in f in
  try
    let r = from_in_channel ic in
    close_in ic;
    r
  with e ->
    close_in ic;
    raise e
