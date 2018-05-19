(* Multiple Sequences Alignment.

  Parse them and other operations.
*)

open Util
open Biology

(* This refers to the alignment position.  *)
type position = int [@@deriving eq, ord, show, yojson]

type 'sr sequence = { start : position; s : 'sr }
and gap = { gstart : position; length : int }
and boundary =
  { label : Gene_region.t
  ; pos : position
  }
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
  | Boundary b  -> b.pos

let al_el_to_string = function
  | Start p                 -> sprintf "Start at %d" p
  | End p                   -> sprintf "End at %d" p
  | Boundary { pos; label } -> sprintf "Boundary %s at %d"
                                (Gene_region.to_string label) pos
  | Sequence { start; s }   -> sprintf "Sequence %s at %d" s start
  | Gap { gstart; length }  -> sprintf "Gap of %d from %d" length gstart

let al_seq_to_string ?(enclose=true) ?(sep=";") lst =
  if lst = [] then begin
    if enclose then "[]" else ""
  end else
    let s = string_of_list ~sep ~f:al_el_to_string lst in
    if enclose then
      sprintf "[%s]" s
    else
      sprintf "%s" s

let is_sequence = function | Sequence _ -> true | _ -> false
let is_gap      = function | Gap _      -> true | _ -> false
let is_end      = function | End _      -> true | _ -> false
let is_start    = function | Start _    -> true | _ -> false

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

type 'a alignment_sequence = 'a alignment_element list


(** Information about how we alter an allele's sequence.

    We define Alterations here because they're fundamental to some of the later
    algorithms so it is helpful to include the type now. But we dont' add any
    alterations to the Parsed allele information in this module. *)
module Alteration = struct

  type per_segment =
    { full  : bool
    ; type_ : Gene_region.t
    ; start : position
    ; end_  : position
    } [@@deriving eq, ord, show, yojson]

  let per_segment_to_string p =
    sprintf "[%s:%s:%d,%d)"
      (Gene_region.to_short p.type_)
      (if p.full then "full" else "partial")
      p.start p.end_

  let per_segment_list_to_string l =
    string_of_list ~sep:";" ~f:per_segment_to_string l

  type t =
    { allele    : string
    ; why       : string
    ; distance  : float
    ; positions : per_segment list
    } [@@deriving show, yojson]

  let to_string { allele; why; distance; positions } =
    sprintf "%s %s (d: %f) positions: [%s]"
      allele why distance (per_segment_list_to_string positions)

end (* Alteration. *)

(* How this works:

  This is a line based format. For our purposes, there are three types of lines
  as determined by {type line}, where `SeqData` has most of the actual alignment
  data. We divide the parsing of sequences into two types: reference and
  alternatives. We take special care of keeping track of gaps in the reference
  so that positions, when parsing the alternatives alleles, are annotated with
  regard to this global-alignment position. *)
module Parser = struct

  exception Error of string

  let perrorf fmt =
    ksprintf (fun s -> raise (Error s)) fmt

  type 'e parse_struct =
    { allele_n        : string
    (* For now, it makes it easier to have some kind of sane delimiter to align
       and spot check these alignments. The "|" are the 'boundaries'. This keeps
       track of the most recently encountered boundary marker, starting with 0.
     *)
    ; sequence_type   : Sequence.type_
    ; boundary_label  : Gene_region.t
    (* Where we are in the sequence, this starts from the first number specified
       in the file and increments as we read characters. *)
    ; position        : position
    ; sequence        : 'e list
    ; in_data         : bool
    }

  (* This method is here, as opposed to Biology.Gene_region, because it is
   * specific to the way that the MSA file is parsed and what it represents.
   *)
  let first_gene_region st pos =
    match st with
    | Sequence.GDNA ->
        if pos < 1 then
          Gene_region.UTR5
        else
          Gene_region.Exon 1
    | Sequence.CDNA
    | Sequence.Protein -> Gene_region. Exon 1

  let init_ps sequence_type allele_n position =
    let boundary_label = first_gene_region sequence_type position in
    { allele_n
    ; position
    ; sequence_type
    ; boundary_label
    ; sequence = [ Boundary { label = boundary_label; pos = position } ]
    ; in_data  = false
    }

  let where ps =
    sprintf "allele_n: %s, position: %d, sequence length: %d"
      ps.allele_n ps.position (List.length ps.sequence)

  (* What type of character have we seen? *)
  type type_of_character =
    | SequenceCharacter                                   (* new data: 'ACGT' *)
    | MetaCharacter                                      (* gap/boundary '|.' *)
    | UnknownChar                                   (* End of knowledge: '*'. *)

  (* Method to make writing 'inserts' that add an element and advance parser
     struct, incrementing position, easier. *)
  let next ps ds ~f =
    (* The ONLY way to insert a Start *)
    let insert_start_before_boundary = function    (* List is stored reversed *)
      | Boundary b :: tl when b.pos = ps.position -> Boundary b :: (Start b.pos) :: tl
      | lst                                       -> (Start ps.position) :: lst
    in
    let now_in_data, new_seq =
      match ds, ps.in_data with
      | SequenceCharacter,  true  -> true,       ps.sequence
      | SequenceCharacter,  false -> true,       insert_start_before_boundary ps.sequence
      | UnknownChar,        true  -> false,      (End ps.position) :: ps.sequence
      | UnknownChar,        false -> false,      ps.sequence
      | MetaCharacter,      _     -> ps.in_data, ps.sequence
    in
    { ps with in_data = now_in_data
            ; position = ps.position + 1
            ; sequence = f ps.position new_seq
    }

  (* Specific alignment element inserts *)

  let insert_boundary ps =
    let next_label = Gene_region.next ps.sequence_type ps.boundary_label in
    let nps =
      next ps MetaCharacter ~f:(fun position sequence ->
          Boundary { label = next_label; pos = position} :: sequence)
    in
    { nps with boundary_label = next_label
             (* Do not advance position due to boundaries.*)
             ; position = ps.position }

  let insert_gap ps =
    next ps MetaCharacter ~f:(fun position l -> match l with
      | Gap { gstart; length } :: t when gstart + length = position
                        -> Gap { gstart; length = length + 1 } :: t
      | []
      | End _ :: _
      | Start _ :: _
      | Boundary _ :: _
      | Gap _ :: _
      | Sequence _ :: _ -> Gap { gstart = position; length = 1 } :: l)

  (* We want to fail_on_same when parsing the reference,
     but not for alternate alleles *)
  let insert_same ~fail_on_same ps =
    if fail_on_same then
      perrorf "Encountered unexpected '-' same char for : %s" (where ps)
    else
      next ps SequenceCharacter ~f:(fun _position sequence -> sequence)

  let insert_nuc c ps =
    next ps SequenceCharacter ~f:(fun position l -> match l with
      | Sequence {start; s} :: t when start + (List.length s) = position
                        -> Sequence { start; s = c :: s} :: t
      | []              -> perrorf "Adding char %c %d %s at Empty state!"
                              c position ps.allele_n
      | End _ :: _      -> perrorf "Adding char %c %d %s after End!"
                              c position ps.allele_n
      | Start _ :: _
      | Boundary _ :: _
      | Gap _ :: _
      | Sequence _ :: _ -> Sequence { start = position; s = c :: [] } :: l )

  let insert_unknown ps =
    next ps UnknownChar ~f:(fun _position sequence -> sequence)

  (* Advance the parser struct across a string of characters. *)
  let update ~allele ~st ~fail_on_same ps s =
    let is_vc = Sequence.is_valid_character ~s:allele st in
    let rec loop ps = function
      | []                    -> ps
      | '|' :: t              -> loop (insert_boundary ps) t
      | '*' :: t              -> loop (insert_unknown ps) t
      | 'X' :: t when st = Protein -> loop (insert_unknown ps) t   (* add End *)
      | '.' :: t              -> loop (insert_gap ps) t
      | '-' :: t              -> loop (insert_same ~fail_on_same ps) t
      | c :: t when is_vc c   -> loop (insert_nuc c ps) t
      | x :: _                -> perrorf "Unrecognized char '%c' in %s"
                                   x (where ps)
    in
    loop ps (String.to_character_list s)

  type 'sr parse_result =
    { start_pos : position
    ; ref       : string                                (* Name of reference. *)
    (* As we parse the alternative tracks, we have to keep track of the gaps
       that we encounter in the reference, so that all positions are with
       respect to the reference. *)
    ; ref_ps    : 'sr parse_struct
    ; alt_htbl  : (string, 'sr parse_struct) Hashtbl.t         (* Alternates. *)
    }

  let init_parse_result sequence_type ref_allele position =
    { start_pos = position
    ; ref       = ref_allele
    ; ref_ps    = init_ps sequence_type ref_allele position
    ; alt_htbl  = Hashtbl.create 100
    }

  let reverse_seq ~boundary_swap lst =
    let to_string l = String.of_character_list (List.rev l) in
    List.rev_map lst ~f:(function
      | Start _
      | End _
      | Gap _ as e            -> e
      | Boundary b            -> Boundary (boundary_swap b)
      | Sequence { start; s } -> Sequence { start; s = to_string s })

  (* We can only have Boundaries, Gap's or End's at the end, if the sequence
     hasn't terminated with an UnknownChar then we need to explicitly add an
     End in this normalization step. *)
  let normalized_seq ~boundary_swap ps =
    let empty_seq =
      [ Boundary { label = first_gene_region ps.sequence_type ps.position
                 ; pos = ps.position } ]
    in
    if ps.sequence = empty_seq then
      []
    else
      let rec has_end = function
        | End _ :: _      -> true
        | Boundary _ :: t -> has_end t
        | Gap _ :: t      -> has_end t
        | _               -> false
      in
      if has_end ps.sequence then
        reverse_seq ~boundary_swap ps.sequence
      else
        reverse_seq ~boundary_swap (End ps.position :: ps.sequence)

  type line =
    | Position of Sequence.type_ * int             (* nucleotide or amino acid sequence *)
    | Dash
    | SeqData of string * string list

  (* Assume that it has been trimmed. *)
  let parse_data line =
    let open Sequence in
    String.split line ~on:(`Character ' ')
    |> List.filter ~f:((<>) String.empty)
    |> function
        | s :: _ when String.get s 0 = Some '|' -> Some Dash
        | "AA" :: "codon" :: _  -> Some Dash   (* not modeling this at the moment. *)
        | "gDNA" :: pos :: _    -> Some (Position (GDNA, int_of_string pos))
        | "cDNA" :: pos :: _    -> Some (Position (CDNA, int_of_string pos))
        | "Prot" :: pos :: _    -> Some (Position (Protein, int_of_string pos))
        | []                    -> None
        | s :: lst              -> Some (SeqData (s, lst))

  type parse_state =
    | Header
    | Empty
    | Data of line

  type alt =
    { allele : string
    ; seq    : string alignment_element list
    ; alters : Alteration.t list
    }

  let sort_alts_by_nomenclature lst =
    List.map lst ~f:(fun a -> Nomenclature.parse_to_resolution_exn a.allele, a)
    |> List.sort ~cmp:(fun (n1, _) (n2, _) -> Nomenclature.compare_by_resolution n1 n2)
    |> List.map ~f:snd

  type sequence_alignment =
    { release     : string
    ; align_date  : string
    ; locus       : Nomenclature.locus
    ; reference   : string
    ; ref_elems   : string alignment_element list
    ; alt_elems   : alt list
    }

  let report = ref false

  let find_header_lines ic =
    let parse_before_3_32 first_line =
      let colon_split = String.split ~on:(`Character ':') in
      try
        let rec find_release line =
          match colon_split line with
          | [ "IPD-IMGT/HLA Release"; release ] -> find_start release (input_line ic)
          | _                                   -> find_release (input_line ic)
        and find_start release line =
          match colon_split line with
          | ["Sequences Aligned"; ad ] -> Some (release, ad)
          | _                          -> find_start release (input_line ic)
        in
        find_release first_line
      with End_of_file ->
        None
    in
    let parse_after_3_32 first_line =
      let colon_split_after_number_sign s =
        String.split ~on:(`Character ':') (String.drop s ~index:1)
      in
      try
        let rec find_start line =
          match colon_split_after_number_sign line with
          | [ " date"; ad ] -> find_release ad (input_line ic)
          | _               -> find_start (input_line ic)
        and find_release start line =
          match colon_split_after_number_sign line with
          | [ " version"; release ] -> Some (release, start)
          | _                       -> find_release start (input_line ic)
        in
        find_start first_line
      with End_of_file ->
        None
    in
    let first_line = input_line ic in
    match first_line.[0] with
    | None      -> None (* Empty first line! *)
    | Some '#'  -> parse_after_3_32 first_line
    | Some _    -> parse_before_3_32 first_line

  let from_in_channel sequence_type locus ic =
    let latest_reference_position = ref min_int in
    let update x = function
      (* Sometimes, the files position counting seems to disagree with this
        internal count, usually because of multiple boundaries. Not certain
        how to get to the bottom, but my manual string counts lead me to
        believe that there isn't a bug in the parsing code.

        One possibility is that there is no '0' the position in the files;
        their indices are [-1, 1, 2].

        Furthermore, x.ref_ps.position now refers to the alignment position
        as opposed to the reference, ie gaps advance the position.

        So we don't check for: x.ref_ps.position = p as well. *)
      | Position (st, p)    -> assert (sequence_type = st); x
      | Dash                -> x                             (* ignore dashes *)
      | SeqData (allele, s) ->
          let fold_sequence ~fail_on_same init =
            List.fold_left ~init ~f:(update ~allele ~st:sequence_type ~fail_on_same)
          in
          if allele = x.ref then begin
            let nref_ps = fold_sequence ~fail_on_same:true x.ref_ps s in
            latest_reference_position := nref_ps.position;
            { x with ref   = allele ; ref_ps = nref_ps }
          end else begin
            let cur_ps =
              try Hashtbl.find x.alt_htbl allele
              with Not_found -> init_ps sequence_type allele x.start_pos
            in
            let new_ps = fold_sequence ~fail_on_same:false cur_ps s in
            (* Can't make this into an assertion because of sequences such as
              C*04:09N that have sequences extending after the end of the
              reference. *)
            if !report && new_ps.position <> !latest_reference_position then
              eprintf "position mismatch %d vs %d for %s.\n"
                !latest_reference_position new_ps.position new_ps.allele_n;
            Hashtbl.replace x.alt_htbl allele new_ps;
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
        | Empty ->
            begin match parse_data line with
            | None   -> loop Empty acc
            | Some (SeqData ("Please", _)) -> acc   (* Fin *)
            | Some d -> loop (Data d) (update acc d)
            end
        | Data  _ ->
            begin match parse_data line with
            | None   -> loop Empty acc
            | Some d -> loop (Data d) (update acc d)
            end
    in
    let rec loop_header state =
      match input_line ic |> String.strip ~on:`Both with
      | exception End_of_file -> perrorf "Didn't get to the data!"
      | line ->
        match state with
        | Header when String.is_empty line -> loop_header Empty
        | Header                           -> loop_header Header
        | Empty  when String.is_empty line -> loop_header Empty
        | Empty                            ->
            begin
              match parse_data line with
              | Some ((Position _) as d)  -> loop_header (Data d)
              | _                         -> perrorf "First data not position."
            end
        | Data _ when String.is_empty line -> loop_header state
        | Data (Position (st, p)) ->
            if st <> sequence_type then
              perrorf "Different sequence type %s from expected: %s"
                (Sequence.show_type_ st) (Sequence.show_type_ sequence_type)
            else begin
              match parse_data line with
              | Some (SeqData (allele, _) as d) ->
                  let res = init_parse_result sequence_type allele p in
                  loop (Data d) (update res d)
              | _                        ->
                  loop_header state
            end
        | Data _ -> loop_header state
    in
    match find_header_lines ic with
    | None ->
        close_in ic;
        perrorf "Couldn't extract sequence align date."
    | Some (release, align_date) ->
        let reversed = loop_header Header in
        let boundary_swap =
          let open Sequence in
          let open Gene_region in
          match sequence_type with
          | CDNA
          | Protein -> id
          | GDNA    ->
              let last_reference_boundary = reversed.ref_ps.boundary_label in
              begin match last_reference_boundary with
              | Intron _n -> fun b ->
                                if b.label = last_reference_boundary then
                                  { b with label = UTR3}
                                else
                                  b
              | rb        -> perrorf "Strangely: parsing gDNA encoded file \
                              didn't end in an Intron to replace with UTR."
              end
        in
        let ref_elems = normalized_seq ~boundary_swap reversed.ref_ps in
        let alt_elems =
          Hashtbl.fold ~init:[] ~f:(fun ~key:all ~data:ps acc ->
              if ps.sequence = [] then begin
                printf "Dropping empty sequence: %s\n" ps.allele_n;
                acc
              end else
                { allele = ps.allele_n
                ; seq    = normalized_seq ~boundary_swap ps
                ; alters = []
                } :: acc)
            reversed.alt_htbl
        in
        { release
        ; align_date
        ; locus
        ; reference = reversed.ref
        ; ref_elems
        ; alt_elems
        }

  let from_file ?sequence_type f =
    let extension_less = Filename.chop_extension (Filename.basename f) in
    let locus, last_after_underscore =
      match String.split ~on:(`Character '_') extension_less with
      | locus_s :: lau :: [] ->
          begin match Nomenclature.parse_locus locus_s with
          | Ok l    -> (l, lau)
          | Error e ->
              (* DRB_nuc hack! *)
              if locus_s = "DRB" && lau = "nuc" then
                (Nomenclature.DRB1, "nuc")
              else
                perrorf "%s" e
          end
      | _ -> perrorf "Filename: %s doesn't match the \"locus_[prot|nuc|gene].txt format."
              f
    in
    let sequence_type =
      match sequence_type with
      | Some bs -> bs
      | None    ->
          match last_after_underscore with
          | "prot" -> Sequence.Protein
          | "nuc"  -> Sequence.CDNA
          | "gen"  -> Sequence.GDNA
          | _           ->
              perrorf "Unrecognized alignment file name: %s, unable to infer \
                            Boundary scheme." f
    in
    let ic = open_in f in
    try
      let r = from_in_channel sequence_type locus ic in
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

  let in_order_invariant l =
    let rec loop = function
      | [] -> Ok ()
      | h :: [] -> Ok ()
      | a :: b :: t ->
          if end_position String.length a > start_position b then
            Error (a, b)
          else
            loop (b :: t)
    in
    loop l

end (* Parser *)

module Boundaries = struct

  type marker =
    { label       : Gene_region.t
    ; position    : position
    ; length      : int
    ; seq_length  : int
    }
    [@@deriving show]

  let marker ~label ~position =
    { label
    ; position
    ; length = 0
    ; seq_length = 0
    }

  let marker_to_string = show_marker

  let of_boundary { label; pos } =
    marker ~label ~position:pos

  let to_boundary ?(offset=0) { label; position; } =
    Boundary { label = label; pos = position + offset }

  let matches_boundary { label; position; _ } = function
    | Boundary b when b.label = label &&  b.pos = position -> true
    | _ -> false

  let all_boundaries_before_start_or_end =
    let rec loop = function
      | Start s :: Boundary bp :: tl when s = bp.pos ->
        Boundary bp :: Start s :: loop tl
      | End e :: Boundary bp :: tl when e = bp.pos ->
        Boundary bp :: End e :: loop tl
      | [] -> []
      | h :: t -> h :: loop t
    in
    loop

  let first_boundary_before_start = function
    | Start s :: Boundary bp :: tl when s = bp.pos ->
        Boundary bp :: Start s :: tl
    | lst -> lst

  (* Does NOT reverse the accumulator, that's up to the caller. There are
     cases where the last segment can require special treatment such as
     in a splice variant. *)
  let fold ~boundary ~start ~end_ ~gap ~seq lst ~init =
    let advance_marker sequence m i =
      if sequence then
        { m with length = m.length + i; seq_length = m.seq_length + i }
      else
        { m with length = m.length + i}
    in
    let rec loop cm cur acc = function
      | []                -> Ok ((cm, cur) :: acc)
      | Boundary bp :: t  -> boundary cur bp >>= fun ob ->
                                loop (of_boundary bp) ob ((cm, cur) :: acc) t
      | Start s :: t      -> start cur s >>= fun os -> loop cm os acc t
      | End e :: t        -> end_ cur e >>= fun oe -> loop cm oe  acc t
      | Gap g :: t        -> let nb = advance_marker false cm g.length in
                             gap cur g >>= fun og -> loop nb og acc t
      | Sequence s :: t   -> let nb = advance_marker true cm (String.length s.s) in
                             seq cur s >>= fun os -> loop nb os  acc t
    and init_loop = function
      | Boundary bp :: t  -> boundary init bp >>= fun ob -> loop (of_boundary bp) ob [] t
      | []                -> error "empty list"
      | e :: _            ->
          error "Alignment sequence didn't start with Boundary, but %s."
            (al_el_to_string e)
    in
    init_loop lst

  let grouped lst =
    let boundary _prev _ = Ok [] in       (* Start new acc on every boundary. *)
    let start l s = Ok (Start s :: l) in
    let end_ l e = Ok (End e :: l) in
    let gap l g = Ok (Gap g :: l) in
    let seq l s = Ok (Sequence s :: l) in
    fold ~boundary ~start ~end_ ~gap ~seq ~init:[]
      (first_boundary_before_start lst) >>= fun lst ->
        Ok (List.rev_map lst ~f:(fun (b, l) -> (b, List.rev l)))

  let ungrouped lst =
    List.map lst ~f:(fun (b, l) -> to_boundary b :: l)
    |> List.concat

end (* Boundaries *)

(* Applying the alignment elements. We want to create methods that take the
   reference and alternate allele alignment_sequence's and act upon aligned
   information. We think of this problem as one where we "initialize" the
   projection based upon the reference. And then we change it based upon
   an allele's sequence/gap information. *)
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

  (* Pass a "length" into the current segment/projection. *)
  val split_due_to_start_stop : int -> bool -> t -> t option * t

end

module MakeZip (R : Data_projection) = struct

  let debug = ref false

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
             ; started = %b \
             }"
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
      if !debug then printf "Appending %s \n" (R.to_string cur.project);
      cur.project :: acc
    in
    let rec advance_allele cur acc a =
      match a with
      | []      -> Ok (cur.started, (append acc cur), [])
      | h :: t  ->
          let sa = start_position h in
          if !debug then
            printf "at %s %d when cur is %s.\n"
              (al_el_to_string h) sa (state_to_string cur);
          if sa < cur.start_pos then
            advance_allele cur acc t
          else if sa < cur.end_pos then
            mutate_segment cur acc t h
          else (* sa >= cur.end_pos *)
            Ok (cur.started, append acc cur, a)
    and mutate_segment cur acc at = function
      | Start pos       -> let b_opt, a = split cur pos true in
                           let nacc = Option.value_map b_opt ~default:acc ~f:(append acc) in
                           advance_allele a nacc at
      | End pos         -> let b_opt, a = split cur pos false in
                           let nacc = Option.value_map b_opt ~default:acc ~f:(append acc) in
                           advance_allele a nacc at
      | Boundary b as e -> if cur.start_pos <> b.pos &&
                              cur.end_pos <> cur.start_pos + 1 then
                              error "Found %s in sequence at %d!"
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
    let reference_pass_r =
      let in_ref = true in
      Boundaries.(fold (first_boundary_before_start reference)
        ~start:(fun s p -> Ok s)          (* Nothing to do for reference Start's *)
        ~end_:(fun s p -> Ok s)             (* Nothing to do for reference End's *)
        ~boundary:(fun (started, acc, allele_lst) b ->
            let new_state =
              { start_pos = b.pos
              ; end_pos = b.pos
              ; project = R.of_boundary ~in_ref started b
              ; started
              }
            in
            advance_allele new_state [] allele_lst)
        ~gap:(fun (started, acc, allele_lst) gap ->
            let new_state =
              { start_pos = gap.gstart
              ; end_pos = gap.gstart + gap.length
              ; project = R.of_gap ~in_ref started gap
              ; started
              }
            in
            advance_allele new_state acc allele_lst)
        ~seq:(fun (started, acc, allele_lst) seq ->
            let new_state =
              { start_pos = seq.start
              ; end_pos = seq.start + String.length seq.s
              ; project = R.of_seq ~in_ref started seq
              ; started
              }
            in
            advance_allele new_state acc allele_lst)
        ~init:(false, [], (first_boundary_before_start allele)))
  in
  let just_bm_and_proj = List.rev_map ~f:(fun (bm, (_, p, _)) -> (bm, p)) in
  reference_pass_r >>= function
    | [] -> error "Didn't even have a start boundary?"
    (* It is possible to have elements of the alternate allele _after_ the
      reference as occurs with some splice variants. *)
    | (bm, (final_st, acc, alleles_lst)) :: tl ->
        let final_st, nacc =
          let in_ref = false in
          List.fold_left alleles_lst ~init:(final_st, acc)
            ~f:(fun (st, a) e ->
                  match e with
                  | End _      -> false, a
                  | Start _    -> true, a
                  | Boundary b -> st, (R.of_boundary ~in_ref st b) :: a
                  | Gap g      -> st, (R.of_gap ~in_ref st g) :: a
                  | Sequence s -> st, (R.of_seq ~in_ref st s) :: a)
        in
        Ok (just_bm_and_proj ((bm, (final_st, nacc, [])) :: tl))

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
        local_errorf "wrong char at %d in %s" d s.s
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

end)  (* AlleleSequences *)

let allele_sequences ~reference ~allele =
  try
    AlleleSequences.apply ~reference ~allele
    >>| List.map ~f:(fun (bm, plst) ->
      bm,
      B.(to_native_string (concat
        (List.rev_filter_map plst ~f:(function
          | (false, _) -> None
          | (true, t)  -> Some (filter t ~f:(fun c -> c <> default_gap_char)))))))
  with Local_error e ->
    Error e

let allele_sequence ?boundary_char ~reference ~allele () =
  let sep = Option.value_map ~default:String.empty ~f:String.of_character boundary_char in
  allele_sequences ~reference ~allele
  >>| List.map ~f:snd
  >>| String.concat ~sep

let reference_sequence_from_ref_alignment_elements ?boundary_char l =
  begin match boundary_char with
  | None    ->
    List.filter_map l ~f:(function
        | Sequence s -> Some s.s
        | _ -> None)
  | Some bc ->
    let bs = String.of_character bc in
    List.filter_map l ~f:(function
        | Sequence s -> Some s.s
        | Boundary _ -> Some bs
        | _ -> None)
  end
  |> String.concat ~sep:""

let reference_sequence ?boundary_char mp =
  reference_sequence_from_ref_alignment_elements
    ?boundary_char mp.Parser.ref_elems

(* Compute a hamming distance between sequences. *)
module Hamming_distance_projection = struct

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
    | G of r             (* Gap *)
    | S of r        (* Sequence *)
    | B             (* Boundary *)

  let to_string = function
    | G r -> sprintf "G %s" (r_to_string r)
    | S r -> sprintf "S %s" (r_to_string r)
    | B   -> "B"

  let incr m = function
    | S r -> S { r with mismatches = r.mismatches + m }
    | G r -> G { r with mismatches = r.mismatches + m }
    | B   -> local_errorf "Asked to increase Boundary by %d" m

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
    | B   -> local_errorf "Asked to add_gap of %d to Boundary." length

  let split_due_to_start_stop index new_start = function
    | S r when index = 0 -> None, S { r with started = new_start }
    | S r                -> let b, a = split_r index new_start r in
                            Some (S b) , S a
    | G r when index = 0 -> None , G { r with started = new_start }
    | G r                -> let b, a = split_r index new_start r in
                            Some (G b) , G a
    | B                 -> if index = 0 then None, B else
                            local_errorf "Asked to increase Boundary by %d." index

end

module Allele_reference_hamming_distance = MakeZip (Hamming_distance_projection)

module Segments = struct

  type relationship =
    | Missing
    | Partial of int    (* sequence length *)
    | Full of int       (* sequence length might be > reference length *)

  type 'a t =
    { seq_length    : int
    ; mismatches    : int
    ; relationship  : 'a
    }

  let distances ~reference ~allele =
    let open Hamming_distance_projection in
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
    try
      Allele_reference_hamming_distance.apply ~reference ~allele
      >>| List.map  ~f:(fun (bm, lst) ->
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
          | `NoStart        -> { seq_length   = bm.Boundaries.seq_length
                              ; mismatches   = 0
                              ; relationship = Missing
                              }
          | `Start (m, l)   -> { seq_length   = bm.Boundaries.seq_length
                              ; mismatches   = m
                              ; relationship = Full l
                              }
          | `Partial (m, l) -> { seq_length   = bm.Boundaries.seq_length
                              ; mismatches   = m
                              ; relationship = Partial l
                              })
    with (Local_error e) ->
      Error e

  let split_long_al_els pos = function
    | Start _
    | End _
    | Boundary _ as e -> local_errorf "Asked to split %s at %d" (al_el_to_string e) pos
    | Sequence seq    -> let b, a = split_sequence seq ~pos in
                         Sequence b, Sequence a
    | Gap gap         -> let b, a = split_gap gap ~pos in
                         Gap b, Gap a

  type 'a at_same_pos =
    | Fin
    | Fst of 'a alignment_element * 'a alignment_element list * 'a alignment_element list
    | Snd of 'a alignment_element * 'a alignment_element list * 'a alignment_element list
    | Both of 'a alignment_element
            * 'a alignment_element list
            * 'a alignment_element
            * 'a alignment_element list

  (* Advance (zip) two alignment_sequence's such that we return Both if two
     elements start at the same position and have the same length. *)
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
        (* TODO: Refactor this code to use the Interval logic. *)
        if sp1 < sp2 then
          if ep1 <= sp2 then        (* h1 is completely first! *)
            Fst (h1, t1, s2)
          else                      (* intersect *)
            let b, a = split_long_al_els sp2 h1 in
            Fst (b, a :: t1, s2)
        else if sp1 = sp2 then
          if ep1 < ep2 then
            match h1 with
            | Start _ | End _ | Boundary _ -> Fst (h1, t1, s2)
            | Gap _ | Sequence _           -> let b, a = split_long_al_els ep1 h2 in
                                              Both (h1, t1, b, a :: t2)
          else if ep1 = ep2 then
            Both (h1, t1, h2, t2)
          else (* ep1 > ep2 *)
            match h2 with
            | Start _ | End _ | Boundary _ -> Snd (h2, s1, t2)
            | Gap _ | Sequence _           -> let b, a = split_long_al_els ep2 h1 in
                                              Both (b, a :: t1, h2, t2)
        else (* sp1 > sp2 *)
          if ep2 <= sp1 then        (* h2 is completely first! *)
            Snd (h2, s1, t2)
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
    same_position_fold ~init:() ~f:(fun () s ->
      match s with
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
    | Fin                     -> r
    | Fst (e, s1, s2)       ->
        begin match relationship e with
        | `After              -> Fin
        | `Inside false       -> r
        | `Inside true        -> let b, a = split_long_al_els end_pos e in
                                 Fst (b, a :: s1, s2)
        end
    | Snd (e, s1, s2)       ->
        begin match relationship e with
        | `After              -> Fin
        | `Inside false       -> r
        | `Inside true        -> let b, a = split_long_al_els end_pos e in
                                 Snd (b, s1, a :: s2)
        end
    | Both (e1, s1, e2, s2) ->
        begin match relationship e1 with
        | `After              -> Fin
        | `Inside false       -> r
        | `Inside true        -> let b1, a1 = split_long_al_els end_pos e1 in
                                 let b2, a2 = split_long_al_els end_pos e2 in
                                 Both (b1, a1 :: s1, b2, a2 :: s2)
        end

  (* We 'zip' 2 allele alignment_sequence's together against the reference. *)
  module Zip2 = struct

    let debug = ref false

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
      { started1
      ; started2
      ; mismatches = 0
      ; is_seq1 = false
      ; is_seq2 = false
      ; is_boundary = true
      ; start_pos = b.pos
      ; end_pos = b.pos
      }

    let state_of_gap ?(in_ref=true) ~started1 ~started2 gap =
      let mismatches = if in_ref then 0 else gap.length in
      { started1
      ; started2
      ; mismatches
      ; is_seq1 = false
      ; is_seq2 = false
      ; is_boundary = false
      ; start_pos = gap.gstart
      ; end_pos = gap.gstart + gap.length
      }

    let state_of_sequence ?(in_ref=true) ~started1 ~started2 seq =
      let n = String.length seq.s in
      let mismatches = if in_ref then 0 else n in
      { started1
      ; started2
      ; mismatches
      ; is_seq1 = true
      ; is_seq2 = true
      ; is_boundary = false
      ; start_pos = seq.start
      ; end_pos = seq.start + String.length seq.s
      }

    let state_to_string s =
      sprintf "{%b; %b; %b; %b; %d; %d; %d}"
        s.started1 s.started2 s.is_seq1 s.is_seq2 s.start_pos s.end_pos
        s.mismatches

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
        if !debug then printf "Appending %s \n" (state_to_string state);
        state :: acc
      in
      (* Advance the 2 allele alignment sequences until we have an element,
         from one (if so which one) or both at the next position.*)
      let rec adv_alleles cur acc a1 a2 =
        match same_pos_and_length_step_upto cur.end_pos a1 a2 with
        | Fin                   -> Ok (cur.started1, cur.started2, append acc cur, a1, a2)
        | Fst (e, s1, s2)       -> one_side cur acc s1 s2 `Fst e
        | Snd (e, s1, s2)       -> one_side cur acc s1 s2 `Snd e
        | Both (e1, s1, e2, s2) -> two_side cur acc s1 s2 e1 e2

      and split_wrap cur pos ?new_start1 ?new_start2 acc s1 s2 =
        let b_opt, a = split cur pos ?new_start1 ?new_start2 in
        let nacc = Option.value_map b_opt ~default:acc ~f:(append acc) in
        adv_alleles a nacc s1 s2

      and split_wrap_only_one cur pos acc s1 s2 new_start = function
        | `Fst -> split_wrap cur pos ~new_start1:new_start acc s1 s2
        | `Snd -> split_wrap cur pos ~new_start2:new_start acc s1 s2

      and one_side cur acc s1 s2 which e =
        if !debug then
          printf "one side at %s when cur is %s.\n"
            (al_el_to_string e) (state_to_string cur);
        if start_position e < cur.start_pos then
          adv_alleles cur acc s1 s2 (* skip element *)
        (* The same_pos_and_length_step_upto logic prevents us form seeing
          elements after the current position. *)
        else
          match e with
          | Start pos       -> split_wrap_only_one cur pos acc s1 s2 true which
          | End pos         -> split_wrap_only_one cur pos acc s1 s2 false which
          | Boundary b as e ->
              if cur.start_pos <> b.pos && cur.end_pos <> cur.start_pos + 1 then
                error "Found %s in sequence at %d!"
                  (al_el_to_string e) (cur.start_pos)
              else (* ignore *)
                adv_alleles cur acc s1 s2
          | Sequence seq    -> adv_alleles (add_seq_wrap_one cur seq which) acc s1 s2
          | Gap gap         -> adv_alleles (add_gap_wrap_one cur gap which) acc s1 s2

      and two_side cur acc s1 s2 e1 e2 =
        if !debug then
          printf "two side at %s %s when cur is %s.\n"
            (al_el_to_string e1) (al_el_to_string e2) (state_to_string cur);
        if start_position e1 < cur.start_pos then
          adv_alleles cur acc s1 s2
        else
          match e1,     e2 with
          | Start s,    Start _    -> split_wrap cur s acc s1 s2 ~new_start1:true ~new_start2:true
          | Start s,    End _      -> split_wrap cur s acc s1 s2 ~new_start1:true ~new_start2:false
          | Start s,    Boundary _ -> split_wrap cur s acc s1 (e2 :: s2) ~new_start1:true

          | End e,      Start _    -> split_wrap cur e acc s1 s2 ~new_start1:false ~new_start2:true
          | Boundary _, Start s    -> split_wrap cur s acc (e1 :: s1) s2 ~new_start2:true
          | End e,      End _      -> split_wrap cur e acc s1 s2 ~new_start1:false ~new_start2:false

          | Start _,    _
          | End _,      _
          | _,          Start _
          | _,          End _      -> error "%s paired with non-zero-length al-el %s"
                                        (al_el_to_string e1) (al_el_to_string e2)

          | Boundary b, Boundary _ -> (* make sure we're at a boundary in the ref *)
              if cur.start_pos <> b.pos && cur.end_pos <> cur.start_pos then
                error "Boundaries don't align %s %s in sequence at %d!"
                  (al_el_to_string e1) (al_el_to_string e2) (cur.start_pos)
              else (* ignore *)
                adv_alleles cur acc s1 s2
          | Boundary _, _
          | _,          Boundary _ -> error "Boundaries not aligned %s %s in sequence at %d!"
                                        (al_el_to_string e1) (al_el_to_string e2) (cur.start_pos)

          (* Same *)
          | Gap _,       Gap _                        -> adv_alleles (set_seq cur ~is_seq1:false ~is_seq2:false) acc s1 s2
          | Sequence q1, Sequence q2 when q1.s = q2.s -> adv_alleles (set_seq cur ~is_seq1:true ~is_seq2:true) acc s1 s2

          (* Diff *)
          | Sequence s,  Sequence _  (*q1.s <> q2.s*) -> adv_alleles (add_seq cur s ~is_seq1:true ~is_seq2:true) acc s1 s2
          | Sequence _,  Gap g                        -> adv_alleles (add_gap cur g ~is_seq1:true ~is_seq2:false) acc s1 s2
          | Gap g,       Sequence _                   -> adv_alleles (add_gap cur g ~is_seq1:false ~is_seq2:true) acc s1 s2
      in
      let reference_pass_r =
        Boundaries.(fold (first_boundary_before_start reference)
          ~init:( false                                                     (* started1 *)
                , false                                                     (* started2 *)
                , []                                                     (* accumulator *)
                , (first_boundary_before_start allele1)
                , (first_boundary_before_start allele2))
          ~start:(fun state _ -> Ok state)       (* Nothing to do for reference Start's *)
          ~end_:(fun state _ -> Ok state)          (* Nothing to do for reference End's *)
          ~boundary:(fun (started1, started2, acc, a1, a2) b ->
                let new_state = state_of_boundary ~started1 ~started2 b in
                adv_alleles new_state [] a1 a2)                      (* reset acc to [] *)
          ~gap:(fun (started1, started2, acc, a1, a2) gap ->
                let new_state = state_of_gap ~started1 ~started2 gap in
                adv_alleles new_state acc a1 a2)
          ~seq:(fun (started1, started2, acc, a1, a2) seq ->
                let new_state = state_of_sequence ~started1 ~started2 seq in
                adv_alleles new_state acc a1 a2))
    in
    let just_bm_and_state = List.rev_map ~f:(fun (bm, (_, _, p, _, _)) -> (bm, p)) in
    reference_pass_r >>= function
    | [] -> error "Didn't even have a start boundary for zip2?"
    | (final_bm, (fs1, fs2, acc, a1, a2)) :: tl ->
        if !debug then
          printf "After reference pass still have:\nAllele1: %s\nAllele2: %s\n"
            (al_seq_to_string a1)
            (al_seq_to_string a2);
        try
          let final1, final2, nacc =
            same_position_fold a1 a2 ~init:(fs1, fs2, acc)
              ~f:(fun (started1, started2, acc) ns ->
                  match ns with
                  | `Fst (Start _)    -> (true, started2, acc)
                  | `Fst (End _)      -> (false, started2, acc)
                  | `Snd (Start _)    -> (started1, true, acc)
                  | `Snd (End _)      -> (started1, false, acc)
                  | `Fst (Boundary b)
                  | `Snd (Boundary b) -> let ns = state_of_boundary ~started1 ~started2 b in
                                        (started1, started2, ns :: acc)
                  | `Fst (Gap g)
                  | `Snd (Gap g)      -> let ns = state_of_gap ~started1 ~started2 ~in_ref:false g in
                                        (started1, started2, ns :: acc)
                  | `Fst (Sequence s)
                  | `Snd (Sequence s) -> let ns = state_of_sequence ~started1 ~started2 ~in_ref:false s in
                                        (started1, started2, ns :: acc)

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
                        local_errorf "%s paired with non-zero-length al-el %s, past reference"
                          (al_el_to_string e1) (al_el_to_string e2)

                    | Boundary b, Boundary _  ->
                          let ns = state_of_boundary ~started1 ~started2 b in
                          (started1, started2, ns :: acc)
                    | Boundary _, _
                    | _,          Boundary _  ->
                        local_errorf "Boundaries not aligned %s %s in sequence, past reference!"
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
          Ok (just_bm_and_state ((final_bm, (final1, final2, nacc, [], [])) :: tl))
        with (Local_error e) ->
          Error e

  end (* Zip2 *)

  let distances_between ~reference ~allele1 ~allele2 =
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
    try
      Zip2.zip2 ~reference ~allele1 ~allele2
      >>| List.map ~f:(fun (bm, slst) ->
          let (fs1, fs2, mismatches) =
            List.fold_left slst ~init:(`Missing, `Missing, 0)
              ~f:(fun (st1, st2, m) zs ->
                    let nst1 = update_state_of_one zs zs.is_seq1 (st1, zs.started1) in
                    let nst2 = update_state_of_one zs zs.is_seq2 (st2, zs.started2) in
                    nst1, nst2, m + zs.mismatches)
          in
          { seq_length  = bm.Boundaries.seq_length
          ; mismatches
          ; relationship =
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
    with (Local_error e) ->
      Error e

end (* Segments *)
