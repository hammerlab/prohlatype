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

let invalid_argf fmt = ksprintf invalid_arg fmt

(* This 'could' be 2 separate types (for reference and non). *)
type 'sr sequence_element =
  | Start of int * string
  | End of int
  | Boundary of (int * int)   (* idx and pos *)
  | Nuc of (int * 'sr)        (* pos and what *)
  | Gap of (int * int)        (* pos and length *)

let sequence_element_to_string_g ~sr_to_string = function
  | Start (p, allele) -> sprintf "Start %s at %d" allele p
  | End p             -> sprintf "End %d" p
  | Boundary (i, p)   -> sprintf "Boundary %d at %d" i p
  | Nuc (p, s)        -> sprintf "Nucleotide %s at %d" (sr_to_string s) p
  | Gap (p, l)        -> sprintf "Gap %d starting at %d" l p

let sequence_element_to_string =
  sequence_element_to_string_g ~sr_to_string:(fun x -> x)

let position = function
  | Start (p, _)
  | End p
  | Boundary (_, p)
  | Nuc (p, _)
  | Gap (p, _) -> p

type data_pos =
  | Before
  | In
  | After

type 'sr parse_struct =
  { allele    : string
  (* For now, it makes it easier to have some kind of sane delimiter to align
     and spot check these alignments. The "|" are the 'boundaries'. This keeps
     track of the most recently encountered boundary marker, starting with 0. *)
  ; boundary  : int
  (* Where we are in the sequence, this starts from the first number specified
     in the file and increments as we read characters. This value must be with
     regard to the _reference_ position. *)
  ; position  : int
  ; sequence  : 'sr sequence_element list
  (* pos, size : also used to indicate how much of the gap left to traverse *)
  ; gaps      : (int * int) list
  (* Sequences use '*' to indicate missing information at the beginning or end
     of the alignment. This flag (starts as false), determines if we're past
     the parts that are missing. *)
  ; in_data   : data_pos
  }

let init_ps allele position =
  { allele
  ; position
  ; boundary = 0
  ; sequence = []
  ; gaps     = []
  ; in_data  = Before
  }

let new_status ps datum =
  match ps.in_data, datum with
  | Before, true  -> { ps with in_data = In
                             ; sequence = Start (ps.position, ps.allele) :: ps.sequence }
                     , false
  | Before, false -> ps, false
  | In, false     -> { ps with in_data = After
                             ; sequence = End ps.position :: ps.sequence }
                     , false
  | In, true      -> ps, false
  | After, true
  | After, false  -> ps, true

let update_ps ?incr ?(datum=true) insert_seq ps =
  let ps, closed = new_status ps datum in
  if closed then ps else
  match incr with
  | None ->
    { ps with sequence = insert_seq ps.sequence }
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

let update_gaps_in_ps new_gaps ps =
   { ps with gaps = List.filter (fun (p, _) -> p >= ps.position) (new_gaps @ ps.gaps) }

(* Define special characters *)
let boundary_char = '|'
let gap_char  = '.'
let same_char = '-'
let missing_char = '*'

let is_boundary c = c = boundary_char
let is_gap c      = c = gap_char
let is_same c     = c = same_char
let is_missing c  = c = missing_char
let is_nucleotide = function 'A' | 'C' | 'G' | 'T' -> true | _ -> false
let is_amino_acid =
  function
  | 'A'| 'C'| 'D'| 'E'| 'F'| 'G'| 'H'| 'I'| 'K'| 'L'
  | 'M'| 'N'| 'P'| 'Q'| 'R'| 'S'| 'T'| 'V'| 'W'| 'Y' -> true
  | _ -> false

let to_forward_pos p len = function
  | `Ref                            -> p + len
  | `Alt ((gp, _) :: _) when gp = p -> p
  | `Alt _                          -> p + len

let to_forward_pos2 p len = function
  | `Ref                            -> p
  | `Alt ((gp, _) :: _) when gp = p -> p
  | `Alt _                          -> p + len

let insert_nuc pos ~seq c = function
  | End _ :: _              -> invalid_argf "Trying to insert nucleotide %c at %d past end" c pos
  | []
  | Start _ :: _
  | Boundary _ :: _
  | Gap _ :: _ as l         -> Nuc (pos, c :: []) :: l
  | (Nuc (p, ss) :: t) as l ->
      let forward_pos = to_forward_pos p (List.length ss) seq in
      if forward_pos = pos then
        Nuc (p, c :: ss) :: t
      else
        Nuc (pos, c :: []) :: l

let insert_nuc_s ~seq c ps =
  update_ps ~incr:`Pos (insert_nuc ~seq ps.position c) ps

let insert_gap ~seq pos = function
  | End _ :: _               -> invalid_argf "Trying to insert gap at %d past end" pos
  | []
  | Start _ :: _
  | (Boundary _ :: _)
  | (Nuc _ :: _) as l        -> Gap (pos, 1) :: l
  | (Gap (p, len) :: t) as l ->
    (*if p = pos then  Gaps do not extend the position. *)
      let forward_pos = to_forward_pos2 p len seq in
      if forward_pos = pos then
        Gap (p, len + 1) :: t
      else
        Gap (pos, 1) :: l

let insert_gap_s_f ?incr ~seq ps = update_ps ?incr (insert_gap ~seq ps.position) ps

let output_debug action state ps =
  let p0, gl = match ps.gaps with | (x,y) :: _ -> (x, y) | [] -> (-1, -1) in
  printf "%s %s, p: %d, b %d, %s, (%d,%d) \n"
    action state ps.position ps.boundary
    (try sequence_element_to_string_g ~sr_to_string:String.of_character_list (List.hd ps.sequence)
     with Failure _ -> "empty hd")
        p0 gl

let insert_gap_s ?incr ~seq ps =
  (*let () = if ps.allele = "A*68:02:02" then output_debug "inserting gap" "before" ps in *)
  let r = insert_gap_s_f ?incr ~seq ps in
  (*let () = if r.allele = "A*68:02:02" then output_debug "inserting gap" "after" r in *)
  r

let insert_missing ps = update_ps ~datum:false ~incr:`Pos (fun l -> l) ps

let insert_boundary_s ps =
  update_ps ~incr:`Bnd
    (fun s -> Boundary (ps.boundary + 1, ps.position) :: s ) ps

let to_ref_seq_elems dna ps s =
  (*let () = output_debug "ref seq" s ps in*)
  let is_nuc = if dna then is_nucleotide else is_amino_acid in
  (* do not increase the position based on gaps in reference sequence.*)
  let rec to_ref_seq_elems_char ps = function
    | []                          -> ps
    | b :: t when is_boundary b   -> to_ref_seq_elems_char (insert_boundary_s ps) t
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
    | b :: t when is_boundary b   -> update_seq_char (insert_boundary_s ps) t
    | m :: t when is_missing m    -> update_seq_char (insert_missing ps) t
    | g :: t when is_gap g        -> update_seq_char (insert_gap_s ~seq:(`Alt ps.gaps) ~incr:`Pos ps) t
    | c :: t when is_nuc c        -> update_seq_char (insert_nuc_s ~seq:(`Alt ps.gaps) c ps) t
    | s :: t when is_same s       -> update_seq_char (update_ps ~incr:`Pos (fun l -> l) ps) t
    | 'X' :: _ when (not dna)     -> ps (* signals an 'end' for AA's *)
    | x :: _                      -> invalid_argf "unrecognized char in seq: %c" x
  in
  update_seq_char ps (String.to_character_list s)

type 'sr parse_result =
  { dna       : bool        (* DNA or Amino Acid sequence -> diff characters *)
  ; start_pos : int
  ; ref       : string      (* Name of reference. *)
  (* As we parse the alternative tracks, we have to Keep track of the gaps that
     we encounter in the reference, so that all positions are with respect to
     the reference. *)
  ; ref_gaps  : (int * int) list  (* Recently encountered reference gaps. *)
  ; ref_ps    : 'sr parse_struct
  ; alg_htbl  : (string, 'sr parse_struct) Hashtbl.t
  }

let empty_result ref_allele dna position =
  { dna       = dna
  ; start_pos = position
  ; ref       = ref_allele
  ; ref_gaps  = []
  ; ref_ps    = init_ps ref_allele position
  ; alg_htbl  = Hashtbl.create 100
  }

let reverse_seq lst =
  List.rev lst
  |> List.map (function
      | Nuc (n, sl)   -> Nuc (n, List.rev sl |> String.of_character_list)
      | Boundary b    -> Boundary b
      | Gap g         -> Gap g
      | Start (p, s)  -> Start (p, s)
      | End e         -> End e)

let normalize_seq ps =
  match ps.sequence with
  | End _ :: _ -> reverse_seq ps.sequence
  | _ ->          reverse_seq (End ps.position :: ps.sequence)

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
  ; ref_elems : string sequence_element list
  ; alt_elems : (string * string sequence_element list) list
  }

let gaps_to_string gps =
  String.concat ~sep:";" (List.map (fun (p,l) -> sprintf "(%d,%d)" p l) gps)

let report = ref false

let from_in_channel ic =
  let previous_reference_position = ref min_int in
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
        let prev_pos = x.ref_ps.position in
        previous_reference_position := !latest_reference_position;
        let nref_ps = List.fold_left (to_ref_seq_elems x.dna) x.ref_ps s in
        latest_reference_position := nref_ps.position;
        let new_gaps = find_reference_gaps prev_pos nref_ps.sequence in
        let new_gaps =
          (* are we crossing a gap across parse boundaries? *)
          match new_gaps with
          | (p, new_length) :: tl when p = prev_pos ->
              begin try
                let _, prev_length = List.find (fun (p, _) -> p = prev_pos) x.ref_gaps in
                (p, new_length - prev_length) :: tl
                with Not_found -> new_gaps
              end
          | _ -> new_gaps
        in
        { x with ref       = allele
               ; ref_ps    = nref_ps
               ; ref_gaps  = new_gaps
        }
      end else begin
        let cur_ps =
          try Hashtbl.find x.alg_htbl allele
          with Not_found -> init_ps allele x.start_pos
        in
        let new_ps =
          List.fold_left (update_seq x.dna) (update_gaps_in_ps x.ref_gaps cur_ps) s
        in
        (* Can't make this into an assertion because of sequences such as
            C*04:09N that have sequences extending after the end of the
            reference. *)
        if !report
           && new_ps.in_data <> After
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
  let ref_elems = normalize_seq reversed.ref_ps in
  let alt_elems =
    Hashtbl.fold (fun all ps acc -> (all, normalize_seq ps) :: acc)
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

