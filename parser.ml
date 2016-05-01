(*
HLA-A Genomic Sequence Alignments
IPD-IMGT/HLA Release: 3.24.0
Sequences Aligned: 2016 April 15
Steven GE Marsh, Anthony Nolan Research Institute.
Please see http://hla.alleles.org/terms.html for terms of use.

 gDNA              -300
                   |
 A*01:01:01:01     CAGGAGCAGA GGGGTCAGGG CGAAGTCCCA GGGCCCCAGG CGTGGCTCTC AGGGTCTCAG GCCCCGAAGG CGGTGTATGG ATTGGGGAGT CCCAGCCTTG
 A*01:01:01:02N    ********** ********** ********** ********** ********** ********** ********** *********- ---------- ----------
 A*01:01:01:03     ********** ********** ********** ********** ********** ********** ********** ********** ********** **********
   *)

open Printf
module String = Sosa.Native_string

let invalid_argf fmt = ksprintf invalid_arg fmt

type data =
  | Position of int
  | Dash
  | SeqData of string * string list

(* Assume that it has been trimmed. *)
let parse_data line =
  String.split line ~on:(`Character ' ')
  |> List.filter ((<>) String.empty)
  |> function
      | ["|"]                   -> Dash
      | "AA" :: "codon" :: _    -> Dash (* not really but not modeling this at the moment. *)
      | ["gDNA"; pos]           -> Position (int_of_string pos)
      | ["cDNA"; pos]           -> Position (int_of_string pos)
      | []                      -> invalid_arg "Empty data line!"
      | s :: lst (*when lst <> []*) -> SeqData (s, lst)
      (*| _                       -> invalid_arg ("Unexpected data line:" ^  line) *)

(* This 'could' be 2 separate types (for reference and non). *)
type 'sr seq_elems =
  | Start of int * string
  | End of int
  | Boundary of (int * int)   (* idx and pos *)
  | Nuc of (int * 'sr)        (* pos and what *)
  | Gap of (int * int)        (* pos and length *)

let position = function
  | Start (p, _)
  | End p
  | Boundary (_, p)
  | Nuc (p, _)
  | Gap (p, _) -> p

type 'sr parse_struct =
  { allele    : string
  (* For now, it makes it easier to have some kind of sane delimiter to align
     and spot check these alignments. The "|" are the 'boundaries' *)
  ; boundary  : int
  (* Where we are in the sequence, this starts from the first number specified
     in the file and increments as we read characters. *)
  ; position  : int
  ; sequence  : 'sr seq_elems list
  ; gaps      : (int * int) list    (* pos and size *)
  ; in_data   : bool
  }

let init_ps allele position =
  { allele
  ; position
  ; boundary = 0
  ; sequence = []
  ; gaps     = []
  ; in_data  = false
  }

let new_status ps datum =
  match ps.in_data, datum with
  | false, true  -> { ps with in_data = true
                            ; sequence = Start (ps.position, ps.allele) :: ps.sequence }
  | false, false -> ps
  | true, false  -> { ps with in_data = false
                            ; sequence = End ps.position :: ps.sequence }
  | true, true   -> ps

let update_ps ?incr ?(datum=true) insert_seq ps =
  let ps = new_status ps datum in
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

(* reference mode *)
let insert_nuc pos c = function
  | End _ :: _ -> invalid_argf "Trying to insert nucleotide %c at %d past end" c pos
  | []
  | Start _ :: _
  | Boundary _ :: _
  | Gap _ :: _ as l   -> Nuc (pos, c :: []) :: l
  | Nuc (p, ss) :: tl -> Nuc (p, c :: ss) :: tl

let insert_nuc_s c ps =
  update_ps ~incr:`Pos (insert_nuc ps.position c) ps

let insert_gap n = function
  | End _ :: _ -> invalid_argf "Trying to insert gap at %d past end" n
  | []
  | Start _ :: _
  | (Boundary _ :: _)
  | (Nuc _ :: _) as l -> Gap (n, 1) :: l
  | Gap (p, l) :: tl  -> Gap (p, l + 1) :: tl

let insert_gap_s ?incr ps = update_ps ?incr (insert_gap ps.position) ps

let insert_missing ps = update_ps ~datum:false ~incr:`Pos (fun l -> l) ps

let insert_boundary_s ps =
  update_ps ~incr:`Bnd
    (fun s -> Boundary (ps.boundary + 1, ps.position) :: s ) ps

let to_ref_seq_elems ps s =
  (* do not increase the position based on gaps in reference sequence.*)
  let rec to_ref_seq_elems_char ps = function
    | []                          -> ps
    | b :: t when is_boundary b   -> to_ref_seq_elems_char (insert_boundary_s ps) t
    | m :: t when is_missing m    -> to_ref_seq_elems_char (insert_missing ps) t
    | g :: t when is_gap g        -> to_ref_seq_elems_char (insert_gap_s ps) t
    | c :: t when is_nucleotide c -> to_ref_seq_elems_char (insert_nuc_s c ps) t
    | x :: _                      -> invalid_arg (sprintf "unrecognized char in reference ps: %c" x)
  in
  to_ref_seq_elems_char ps (String.to_character_list s)

(* Sequence has not been reversed! *)
let find_reference_gaps until seq =
  let rec loop acc = function
    | Gap (p, l) :: t when p > until   -> loop ((p,l) :: acc) t
    | h :: t when (position h) > until -> loop acc t
    | _ -> acc
  in
  loop [] seq

let update_seq ps s =
  let rec update_seq_char ps = function
    | []                          -> ps
    | b :: t when is_boundary b   -> update_seq_char (insert_boundary_s ps) t
    | m :: t when is_missing m    -> update_seq_char (insert_missing ps) t
    | g :: t when is_gap g        -> update_seq_char (insert_gap_s ~incr:`Pos ps) t
    | c :: t when is_nucleotide c -> update_seq_char (insert_nuc_s c ps) t
    | '-' :: t                    -> update_seq_char (update_ps ~incr:`Pos (fun l -> l) ps) t
    | x :: _                      -> invalid_arg (sprintf "unrecognized char in seq: %c" x)
  in
  update_seq_char ps (String.to_character_list s)

type 'sr result =
  { start_pos : int
  ; reference : string
  ; ref_gaps  : (int * int) list
  ; gap_until : int                     (* This is such ugly hackery for poor reasons. *)
  ; ref_ps    : 'sr parse_struct
  ; alg_htbl  : (string, 'sr parse_struct) Hashtbl.t
  }

let empty_result ref_allele position =
  { start_pos = position
  ; reference = ref_allele
  ; ref_gaps  = []
  ; gap_until = min_int
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

type parse_state =
  | Header
  | Empty
  | Data of data

let parse_ic ic =
  let update x = function
    | Position p     -> { x with start_pos = p }
    | Dash           -> x (* ignore dashes *)
    | SeqData (allele, s) ->
      if x.reference = allele then
        let nref_ps = List.fold_left to_ref_seq_elems x.ref_ps s in
        (* Since the reference sequence is updated first, we can for the
           most recent gaps. *)
        let new_gaps = find_reference_gaps x.gap_until nref_ps.sequence in
        let new_until = match new_gaps with | [] -> x.gap_until | (h,_) :: _ -> h in
        printf "new gaps: until %d \n" x.gap_until;
        List.iter (fun (p,l) -> printf "%d - %d\n" p l) new_gaps;
        { x with reference = allele
               ; ref_ps    = nref_ps
               ; ref_gaps  = new_gaps
               ; gap_until = new_until
        }
      else
        let cur_ps =
          try Hashtbl.find x.alg_htbl allele
          with Not_found -> init_ps allele x.start_pos
        in
        let cur_ps' = List.fold_left update_seq { cur_ps with gaps = x.ref_gaps } s in
        Hashtbl.replace x.alg_htbl allele cur_ps';
        x
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
            | Position p -> loop_header (Data d)
            | _          -> invalid_arg "First data not position."
          end
      | Data _ when String.is_empty line -> loop_header state 
      | Data (Position p) -> 
          begin
            match parse_data line with
            | SeqData (allele, _) as d -> let res = empty_result allele p in
                                          loop (Data d) (update res d)
            | _                        -> loop_header state
          end
      | Data _ -> loop_header state
  in
  let reversed = loop_header Header in
  let ref_seq = normalize_seq reversed.ref_ps in
  let all_seq =
    Hashtbl.fold (fun all ps acc -> (all, normalize_seq ps) :: acc)
      reversed.alg_htbl []
  in
  reversed.reference, ref_seq, all_seq

let parse_f f =
  let ic = open_in f in
  try
    let r = parse_ic ic in
    close_in ic;
    r
  with e ->
    close_in ic;
    raise e


(*
let do_it () =
  let f = "/Users/leonidrozenberg/Documents/code/foreign/IMGTHLA/alignments/A_gen.txt" in
  let ic = open_in f in
  let r = parse ic in
  close_in ic;
  r
*)
