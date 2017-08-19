(* Altering (Merging/Imputing) of cDNA and gDNA alignment sequences. *)

open Util
open MSA

(* The set of genes for which I've tested these algorithms. *)
let supported_genes = [ "A"; "B"; "C"]

(* Filter out the string sequence elements. *)
let no_sequences = List.filter ~f:(fun a -> not (is_sequence a))

(* Specify string containers. *)
let end_position = end_position String.length

(* start end pairs, TODO: Unify this with the sep in Ref_graph *)
type sep2 =
  { start : position
  ; end_  : position
  } [@@deriving eq, ord]

let sep2_to_string s = sprintf "[%d,%d)" s.start s.end_

type info =
  | FullSequence                         (* Nothing added same as original. *)
  | Added of
    { alternate_allele : string
    ; distance         : float
    ; positions        : sep2 list
    }

let info_to_string = function
  | FullSequence -> ""
  | Added { alternate_allele; distance; positions} ->
      sprintf "Added %s (%f) [%s]" alternate_allele distance
        (String.concat ~sep:";" (List.map positions ~f:sep2_to_string))


let full_sequence name l =
  match l with
  | Start s :: _ ->
      begin match List.last l with
      | Some (End e) -> Ok (s, e)
      | _  -> error "%s doesn't end at the end!: %s" name (al_seq_to_string ~sep:"," l)
      end
  | _      -> error "%s doesn't start right away!" name


(* Logic that describes how we split the sequences. *)
module Split : sig

  val al_el_at
      : pos:position
      -> excluding:bool
      -> acc:string alignment_sequence
      -> string alignment_sequence
      -> string alignment_sequence * string alignment_sequence

  val at_start
      : ?acc:'a alignment_sequence
      -> 'a alignment_sequence
      -> 'a alignment_sequence * 'a alignment_sequence

  val at_end
      : ?acc:'a alignment_sequence
      -> 'a alignment_sequence
      -> 'a alignment_sequence * 'a alignment_sequence

end (* Split *) = struct

(* Fold through a list (though I mean an alignment_sequence) searching for a
    position where [f] can tell us to stop, and return the elements seen
    before and the remaining.

    The elements 'before' are pushed onto `acc`.

    [f] describes the where the queried element belongs in relationship to the
    split position.

  `Later        -> Haven't found the right position/element, keep going.
  `Split (b,a)  -> stop and push b(efore) onto the accumulator and a(fter)
                    back onto the list.
  `After        -> element is after the position but stop. *)
  let rec at ?(acc=[]) ~f = function
    | []     -> acc, []
    | h :: t ->
        begin
          match f h with
          | `Split (b, a) -> (b :: acc), (a :: t)
          | `After        -> acc,        (h :: t)
          | `Later        -> at ~f ~acc:(h :: acc) t
        end

  let al_el_at ~pos ~excluding ~acc =
    let single p =
      if pos < p then
        `After
      else if excluding && pos = p then
        `After
      else (* pos >~ p *)
        `Later
    in
    let split_at_pos = function
      | Start p
      | End p        -> single p
      | Boundary b   -> single b.pos
      | Sequence seq -> if pos >= seq.start + String.length seq.s then
                          `Later
                        else if pos <= seq.start then
                          `After
                        else let b, a = split_sequence seq ~pos in
                          `Split (Sequence b, Sequence a)
      | Gap gap      -> if pos >= gap.gstart + gap.length then
                          `Later
                        else if pos <= gap.gstart then
                          `After
                        else let b, a = split_gap gap ~pos in
                          `Split (Gap b, Gap a)
    in
    at ~acc ~f:split_at_pos

  let at_start ?acc l =
    at ?acc ~f:(function | Start _ -> `After | _ -> `Later) l

  let at_end ?acc l=
    at ?acc ~f:(function | End _ -> `After | _ -> `Later) l

end (* Split *)

module Impute : sig

  val debug : bool ref

  (* TODO: Hide this function *)
  val impute : bigger:string alignment_sequence
             -> smaller:string alignment_sequence
             -> (string alignment_sequence * sep2 list, string) result

  val do_it : Distances.logic
            -> Parser.result
            -> (Parser.result * (StringMap.key * info) list, string) result

end (* Impute *) = struct

  let debug = ref false

  (* Keep track of the positions that we're imputing into the smaller sequence. *)
  let add_new_positions ~start ~end_ pos_acc =
    if start = end_ then pos_acc else { start; end_} :: pos_acc

  let reverse_and_remove_boundary_duplicates lst =
    let rec loop a l = match l with
      | Boundary b1 :: ((Boundary b2 :: _) as tl) when b1 = b2 -> loop a tl
      | h :: t  -> loop (h :: a) t
      | []      -> a
    in
    loop [] lst

  let impute ~bigger ~smaller =
    full_sequence "impute bigger sequence" bigger >>= fun (bigger_start, bigger_end) ->
      (* smaller is just past it's start position. *)
      let rec merge_until_end merged pos_acc smaller_start_pos bigger ~smaller =
        if smaller_start_pos >= bigger_end then
          error "smaller start pos %d is at the bigger end %d"
            smaller_start_pos bigger_end
        else
          (* Take bigger elements up to the start position, excluding boundaries. *)
          let start_with, rest_b = Split.al_el_at ~pos:smaller_start_pos ~excluding:false ~acc:merged bigger in
          let npos_acc = add_new_positions ~start:bigger_start ~end_:smaller_start_pos pos_acc in
          (* Take smaller elements up to, excluding, next end. *)
          let without_end, rest_s = Split.at_end ~acc:start_with smaller in
          let smaller_end_pos = end_position (List.hd_exn rest_s) in
          (* Drop bigger elements up to the smaller end. *)
          let _ignore_me, b_after_s_end = Split.al_el_at ~pos:smaller_end_pos ~excluding:true ~acc:[] rest_b in
          (* Lastly, move smaller upto a potential next start. *)
          let _before_another_start, smaller_at_start_or_empty = Split.at_start  rest_s in
          start_merging smaller_end_pos without_end npos_acc
            ~smaller:smaller_at_start_or_empty
            ~bigger:b_after_s_end
          (* Recursive to handle the End -> Start (Unkown-Gap in sequence) cases. *)
      and start_merging last_smaller_pos acc pos_acc ~smaller ~bigger =
        match smaller with
        | Start p :: t -> merge_until_end acc pos_acc p bigger ~smaller:t
        | []           -> if !debug then
                            printf "ending impute with: last_smaller_pos: %d bigger: %s\n%!"
                              last_smaller_pos (al_seq_to_string ~sep:";" bigger);
                          let npos_acc =
                            match bigger with
                            | []     -> pos_acc
                            | h :: t -> let end_ = Option.value_map (List.last t)
                                                    ~f:end_position
                                                    ~default:(end_position h)
                                        in
                                        add_new_positions ~start:last_smaller_pos ~end_ pos_acc
                          in
                          let final_seq =
                            List.rev_append bigger acc
                            |> reverse_and_remove_boundary_duplicates
                          in
                          Ok (final_seq, (List.rev npos_acc))                    (* Fin *)
        | ae :: _      -> invalid_argf "Not at Start or Empty: %s"
                            (al_el_to_string ae)
      in
      let _before_smaller_start, at_smaller_start = Split.at_start smaller in
      match at_smaller_start with
      | Start p :: t -> merge_until_end [] [] p bigger ~smaller:t
      | []           -> invalid_arg "Empty smaller sequence!"
      | ae :: _      -> invalid_argf "Not at Start: %s" (al_el_to_string ae)

  (* TODO: This logic is much simpler if the Distance logic is reference:
     move an allele's Start/End's to the references's. We don't need to
     go the round trip of saying that the distance is 0. *)
  let do_it logic mp =
    let open Parser in
    let module Sm = StringMap in
    let reflength = sequence_length mp.ref_elems in
    let to_fill, full = List.partition_map mp.alt_elems ~f:(fun (a, s) ->
      let l = sequence_length s in
      if l < reflength then `Fst (l, a, s) else `Snd (a, s))
    in
    (* Strip the sequences from the reference since for alternate alleles
       they represent a difference to the reference; therefore without the
       sequences the reference has no difference to itself. *)
    let init = Sm.singleton mp.reference (no_sequences mp.ref_elems) in
    let candidates =
      List.fold_left full ~init ~f:(fun m (key, data) ->
        Sm.add ~key ~data m)
    in
    let merge_assoc = List.map full ~f:(fun (a, _) -> (a, FullSequence)) in
    let to_distances =
      Distances.one ~reference:mp.reference ~reference_sequence:mp.ref_elems
        logic
    in
    List.sort to_fill ~cmp:(fun (l1,_,_) (l2,_,_) -> compare l2 l1)
    |> list_fold_ok ~init:(candidates, merge_assoc)
      ~f:(fun (cm, ma) (_length, a_name, a_seq) ->
            to_distances ~candidates ~allele:(a_name, a_seq) >>= function
              | []            ->
                  error "Did not return _any_ distances to %s allele" a_name
              | (key, d) :: _ ->
                  let bigger = Sm.find key candidates in
                  impute ~smaller:a_seq ~bigger >>= fun (merged, positions) ->
                    let merge_info =
                      Added { alternate_allele = key
                            ; distance = d
                            ; positions
                            } in
                    let nma = (a_name, merge_info) :: ma in
                    Ok (Sm.add ~key:a_name ~data:merged cm, nma))
      >>= fun (merged_alleles, merge_assoc) ->
        let nmp =
          let nalt_elems = Sm.bindings (Sm.remove mp.reference merged_alleles) in
          { mp with alt_elems = nalt_elems}
        in
        Ok (nmp, merge_assoc)

end (* Impute *)

(*
Problem: There are many (10x) more cDNA derived (occasionally referred to as
nucleic following IMGT's naming these alignment files Gene_nuc.txt) allele
sequences than gDNA (occasionally referred to as genetic also following IMGT's
convention of naming these alignment files Gene_gen.txt) ones. But since the
gDNA sequences are longer they allow us to match more reads. We would like
to "impute" sequences for (and hence a graph or PHMM containing) all of the
alleles for which we have some cNDA sequence data. To accomplish this, we use
the current cDNA sequences as scaffolding around which we append information
from the gDNA ones.

The algorithm proceeds in several steps. At a high level for a given gene
(ex. HLA-A) we:

  1. Use the reference allele to construct a set of instructions for how to
     zip the exons into the introns.
  2. Afterwards we use those instructions to, merge the reference.
  3. Merge all alleles where we have both cDNA and gDNA.
  4. Finally for the other alleles we use the same instructions to merge them
     with the closest allele for which we have gDNA data.

At a lower level:

  1. Use the reference allele to construct a set of instructions for how to
     zip the exons into the introns. This will calculate a list of
     instructions that will alternate between an intron or UTR and an exon:
     [align_mp_into_instructions].

     A lot of downstream code is special cased to start with an UTR5
     (aka FillFromGenetic). This step also calculates all of the necessary
     offsetting logic so that alignment elements that are contained in these
     instructions can be moved to the appropriate position for a consistent
     merged alignment.

  2. Merging the reference is fairly straightforward, but one might object as
     to WHY we would even need to merge them! Specifically shouldn't we have
     the same data for the reference? We have the same _sequence_ information
     for the reference (and we check this, see
     [reference_sequence_from_ref_alignment_elements] in [zip_align]), but the
     position of the gaps, with respect to the global alignment is different.
     When we're in an exon we want to use the cDNA alignment but switch to the
     gDNA alignment in an intron/UTR.

  3. Before merging the other alleles we separate them into 3 cases using
      [same_and_diff_warn_on_just_gen]:
      - Alleles where we have only gDNA data (4).
      - Alleles where we have both cDNA and gDNA data (6).
      - Alleles where we have only cDNA data (6).

  4. It would be surprising if we had alleles where we have only gDNA data,
     (since one could presumably infer the nucleic components), so an error
     message is written, but no action is taken.

  5. When merging the alleles where we have both gDNA and cDNA data,
     [same_merge] we take the instructions formed from the reference and
     replace sequence elements with those from the gDNA and cDNA data
     ([map_instr_to_alignments]) and then we "execute" (alternately insert
     them into the list) them.

  6. For alleles where we only have cDNA data we first compute the distances
     to the other alleles (Distances.compute) and follow a similar merging
     logic [diff_merge]. Insert them into the instructions
     ([map_instr_to_alignments]), check to for Start's/End's inside of the
     Exon part and replace them with repositioned gDNA sequence
     ([insert_missing_data_in_merge_nuclear]) and finally execute them.
*)

module Merge : sig

  val debug : bool ref

  type instructions

  val align_mp_into_instructions
      : gen_mp:Parser.result
      -> nuc_mp:Parser.result
      -> (instructions, string) result

  val map_instr_to_alignments
      : string
      -> gen:string alignment_sequence
      -> nuc:string alignment_sequence
      -> instructions
      -> (instructions, string) result

  val insert_missing_data_in_merge_nuclear
      : string
      -> instructions
      -> (instructions, string) result

  val exec_instructions
      : instructions
      -> (string alignment_sequence, string) result

  val to_distance_arguments
      : Parser.result
      -> string list
      -> Distances.arg

  val do_it
      : gen_mp:Parser.result
      -> nuc_mp:Parser.result
      -> Distances.logic
      -> (Parser.result * 'a list, string) result

end (* Merge *) = struct

  let debug = ref false

  type 'a region =
    { bm      : Boundaries.marker
    ; nbp     : position
    (* New boundary position. A convenience for debugging, not necessary. *)
    ; offset  : int
    (* How much to shift elements in this region to conform with the merged
       global alignment. *)
    ; al_seq  : 'a
    (* In general this field is meant to hold ['a alignment_sequence] but I'll
       keep it fully abstract as the basic instructions don't need actual
       sequences and this will be unit. *)
    }

  let region_to_string al_seq_to_string { bm; nbp; offset; al_seq } =
    sprintf "{ bm: %s; nbp: %d; offset: %d; al_seq: %s}"
      (Boundaries.marker_to_string bm) nbp offset (al_seq_to_string al_seq)

  let region_to_final_al_seq r =
    (Boundaries.to_boundary ~offset:r.offset r.bm) :: r.al_seq

  (* When we want to have nuclear and genetic labels. *)
  type 'a ngd =
    { nuclear : 'a                                       (* Coming from cDNA. *)
    ; genetic : 'a                                       (* Coming from gDNA. *)
    }

  let ngd_to_string ts { nuclear; genetic } =
    sprintf "{ nuclear: %s; genetic: %s }" (ts nuclear) (ts genetic)

  type 'a instruction =
    | FillFromGenetic of 'a region
    (* Take the elements from the genetic region and ignore the nuclear region.
       This is the case for intron's or UTR's; where we don't have any
       alignment sequence data to construct a sequence. We still have to
       adjust to the correct position based off of the region's offset.

       All sequences that can be used as a source _must_ have their full
       sequence imputed beforehand; otherwise [map_instr_to_alignments] will
       fail. *)
    | MergeFromNuclear of { reference : 'a region ngd
                          ; to_merge  : 'a region ngd
                          }
    (* Generally, take the data from the nuclear sequence. This is the case
       when we're in an exon region. But what happens if we are in an Exon
       but the original cDNA sequence doesn't have data (ex exon 6 or 7
       of a class I allele where we generally only have exon's 2 & 3); in this
       case use the genetic sequence data that we're imputing from.

       For these regions we also need to remember the scaffolding of the
       reference to determine unique Start/End positions that may occur in
       the merging nuclear sequence. *)

  let instruction_to_string seq_to_string = function
    | FillFromGenetic r ->
        sprintf "FillFromGenetic %s"
          (region_to_string seq_to_string r)
    | MergeFromNuclear { reference; to_merge} ->
        sprintf "MergeFromNuclear {reference: %s; to_merge: %s}"
          (ngd_to_string (region_to_string seq_to_string) reference)
          (ngd_to_string (region_to_string seq_to_string) to_merge)

  let instruction_to_string_seqless =
    instruction_to_string (fun _ -> "")

  let bounded_results_to_string lst =
    String.concat ~sep:";"
      (List.map lst ~f:(fun (bm, seq) ->
        sprintf "(%s, %s)" (Boundaries.marker_to_string bm) (short_seq seq)))

  let shift_al_el offset = function
    | Boundary b -> Boundary { b with pos = b.pos + offset}
    | Start s    -> Start (s + offset)
    | End e      -> End (e + offset)
    | Gap g      -> Gap { g with gstart = g.gstart + offset }
    | Sequence s -> Sequence { s with start = s.start + offset}

  let shift_all_al_els offset =
    List.map ~f:(shift_al_el offset)

  let new_region al_seq new_boundary_pos bm =
    let offset = new_boundary_pos - bm.Boundaries.position in
    { bm = bm
    ; nbp = new_boundary_pos
    ; offset
    ; al_seq = shift_all_al_els offset al_seq
    }

  type s_seq = bytes MSA.alignment_sequence

  type instructions = s_seq instruction list

  (* Create a framework of instructions based off of, the alignment in the
     reference sequence. *)
  let zip_align ~gen ~nuc =
    let ref_seq = reference_sequence_from_ref_alignment_elements in
    let open Boundaries in
    (* Convert the reference alignment_sequence's to Boundary delimited strings,
       this way we can make sure that the sequences are the same. *)
    let gblst = grouped (all_boundaries_before_start_or_end gen) in
    let nblst = grouped (all_boundaries_before_start_or_end nuc) in
    let rec need_exon new_boundary_pos acc gl nl = match gl, nl with
      | ({ label = Exon gx} as gb, gseq) :: gtl
      , ({ label = Exon nx} as nb, nseq) :: ntl ->
        let gss = ref_seq gseq in
        let nss = ref_seq nseq in
        if gx <> nx then
          error "Genetic Exon %d doesn't line up with Nuclear exon %d" gx nx
        else if gss <> nss then
          error "Genetic sequence disagrees with Nucleic at Exon %d: %s" gx
            (manual_comp_display ~width:160 ~labels:("Genetic", "Nuclear")
              gss nss)
        else
          let genetic = new_region gseq new_boundary_pos gb in
          let nuclear = new_region nseq new_boundary_pos nb in
          let ngd = {genetic; nuclear} in
          let instr = MergeFromNuclear { reference = ngd; to_merge = ngd } in
          need_intron (new_boundary_pos + nb.length) (instr :: acc)
            gtl ntl
      | (gb, _) :: _, (nb, _) :: _ ->
          error "Genetic boundary %s doesn't line up with Nuclear boundary %s"
            (marker_to_string gb) (marker_to_string nb)
      | [], _   -> error "Genetic segments ended before nuclear."
      |  _, []  -> error "Nuclear segments ended before genetic."
    and need_intron new_boundary_pos acc gl nl = match gl, nl with
      | ({ label = Intron gx} as gb, gseq) :: gtl
      , ({ label = Exon nx}, _) :: _ ->
        if gx + 1 <> nx then
          error "Genetic Intron %d doesn't line up with Nuclear exon %d" gx nx
        else
          let instr = FillFromGenetic (new_region gseq new_boundary_pos gb) in
          need_exon (new_boundary_pos + gb.length) (instr :: acc) gtl nl
      | ({ label = UTR3 } as gb, gseq) :: [], [] ->
          let instr = FillFromGenetic (new_region gseq new_boundary_pos gb) in
          Ok (List.rev (instr :: acc))                                (* Fin. *)

      | (gb, _) :: _, (nb, _) :: _ ->
          error "Genetic boundary %s doesn't line up with Nuclear boundary %s"
            (marker_to_string gb) (marker_to_string nb)
      | [], _   -> error "Genetic segments ended before nuclear."
      |  _, []  -> error "Nuclear segments ended before genetic."
    in
    match gblst, nblst with
    | ({ label = UTR5 } as gb1, gseq) :: gtl
    , ({ label = Exon 1}, _) :: _ ->
        need_exon (gb1.length + gb1.position)
          [ FillFromGenetic (new_region gseq gb1.position gb1)] gtl nblst
    | (gb1, _) :: _
    , (nb1, _) :: _ -> error "First genetic marker isn't UTR5: %s or First nucleic marker %s isn't Exon 1"
                        (marker_to_string gb1) (marker_to_string nb1)
    | [], _         -> error "Empty genetic alignment_sequence"
    | _, []         -> error "Empty nuclear alignment_sequence"

  let strip_nuclear_start_ends_region { nuclear; genetic } =
    let nnuclear =
      { nuclear with al_seq =
                 List.filter nuclear.al_seq ~f:(fun a ->
                    (is_gap a) || (is_sequence a))
      }
    in
    { genetic; nuclear = nnuclear }

  let strip_nuclear_start_ends allele_name instr =
    List.map instr ~f:(fun e ->
      match e with
      | FillFromGenetic _ -> e
      | MergeFromNuclear { reference; to_merge } ->
          let reference = strip_nuclear_start_ends_region reference in
          let to_merge = strip_nuclear_start_ends_region to_merge in
          MergeFromNuclear {reference; to_merge })

  (* Public *)
  let align_mp_into_instructions ~gen_mp ~nuc_mp =
    let open Parser in
    zip_align ~gen:gen_mp.ref_elems ~nuc:nuc_mp.ref_elems
    >>= fun instr -> Ok (strip_nuclear_start_ends gen_mp.reference instr)

  let replace_al_seq r s =
    { r with al_seq = shift_all_al_els r.offset s }

  (* Public *)
  (* Take sequences (gen, nuc):
      1. insert them into their relevant instructions
      2. shift them as necessary.
  *)
  let map_instr_to_alignments allele_name ~gen ~nuc instr =
    let open Boundaries in
    (* Convert the reference alignment_sequence's to Boundary delimited strings,
       this way we can make sure that the sequences are the same. *)
    full_sequence "Merge.map_instr_to_alignments genetic sequence" gen >>=
      fun (_gen_start_pos, _gen_end_pos) ->
      let gblst = grouped (all_boundaries_before_start_or_end gen) in
      let nblst = grouped (all_boundaries_before_start_or_end nuc) in
      let rec fill_from_g acc fr itl ~nuc ~gen = match gen with
        | (gbm, gseq) :: gtl when gbm.label = fr.bm.label ->
            let nacc = FillFromGenetic (replace_al_seq fr gseq) :: acc in
            begin match itl with
            | []  -> Ok (List.rev nacc)                                (* Fin *)
            | _   -> need_merge nacc itl ~gen:gtl ~nuc
            end
        | []            -> error "Empty genetic sequence"
        | (gbm, _) :: _ -> error "Genetic boundary label %s doesn't match \
                            instruction's: %s"
                            (boundary_label_to_string gbm.label)
                            (boundary_label_to_string fr.bm.label)

      and fill_from_m reference acc itl ~genetic ~nuclear ~gen ~nuc = match gen, nuc with
        | (gbm, gseq) :: gtl, (nbm, nseq) :: ntl
            (* the genetic and nuclear labels already match *)
            when gbm.label = genetic.bm.label && nbm.label = genetic.bm.label ->
            let nacc =
              let genetic = replace_al_seq genetic gseq in
              let nuclear = replace_al_seq nuclear nseq in
              let to_merge = { genetic; nuclear }  in
              MergeFromNuclear {reference; to_merge }:: acc
            in
            need_fill nacc itl ~gen:gtl ~nuc:ntl
        | [],  _        -> error "%s has premature empty genetic sequence"
                            allele_name
        |  _, []        -> error "%s has premature empty nuclear sequence"
                            allele_name
        | (gbm, _) :: _
        , (nbm, _) :: _ -> error "Genetic %s or Nuclear %s boundary label \
                            doesn't match instruction's: %s"
                            (boundary_label_to_string gbm.label)
                            (boundary_label_to_string nbm.label)
                            (boundary_label_to_string genetic.bm.label)

      and need_fill acc ilst ~gen ~nuc = match ilst with
        | FillFromGenetic f :: tl
                -> fill_from_g acc f tl ~gen ~nuc
        | []     -> error "Premature end to instructions for %s" allele_name
        | s :: _ -> error "Mapping instructions for %s, did not start with a \
                          \"starting\" FillFromGenetic but %s"
                      allele_name (instruction_to_string_seqless s)
      and need_merge acc ilst ~gen ~nuc = match ilst with
        | MergeFromNuclear { reference; to_merge}  :: tl ->
            let { genetic; nuclear } = to_merge in
            fill_from_m reference acc tl ~genetic ~nuclear ~gen ~nuc
        | i :: _ -> error "Didn't find a merge from nuclear instruction: %s"
                      (instruction_to_string_seqless i)
        | []     -> error "Premature end to instructions for %s" allele_name

      in
      need_fill [] instr ~gen:gblst ~nuc:nblst

  (* When we try to MergeFromNuclear we may have a case where the nuclear (cDNA
     sourced) sequence has a start (end) in the middle of an exon. In this case
     we need to merge the elements before (after) this position from the genetic
     (gDNA) sequence. But we need to project this nuclear position into the
     genetic alignment.

     nuclear position -> genetic position

     Since each instruction has already grouped the elements by Boundary and
     the genetic sequence is full, this nuclear position can be in one of 4
     states:

     1. Inside a sequence.

     2. Inside a gap. Though this would be highly suspect: the sequence would
        start (end) inside an area where we know that there are no bases:
        CAGT....ACGT
        ----..******        ---> Error.

        CAGT....ACGT
        ******..----        ---> Error.

     3. Before a gap:
        CAGT.....ACGT
        ****.....----       ---> Error.

        For an end:
        CAGT.....ACGT
        ----*********       ---> Ok.

     4. Start after a gap:
        CAGT.....ACGT
        *********....       ---> Ok.

        End after a gap:
        CAGT.....ACGT
        ----.....****       ---> Error.

    For the cases that are Ok, when we need to map from one positions to the
    other a sensible algorithm is:
    1. Figure out the reference sequence (we've already checked that they're
       the same: see [ref_seq] above) length at that position: Fold through
       the nuclear sequence elements adding up their lengths until we reach
       this point.
    2. Traverse the same sequence length in the nuclear elements, adding up
       our total position.

    nuclear position ->
    reference sequence length (from nuclear) ->
    genetic position (from reference sequence length in genetic).
  *)

  (* 1. Is the passed position OK (according to above conditions)?
     2. What is the sequence length to reach the position? *)
  let reference_sequence_position an se pos nuclear =
    let rec loop sp = function
      | []             -> error "In %s ran out of elements looking for %d" an pos
      | h :: tl ->
          begin match h with
          | Boundary b -> invalid_argf "%s should be removed!"
                            (al_el_to_string h)
          | Start _
          | End _      -> invalid_argf "%s at position search: %d"
                            (al_el_to_string h) pos
          | Gap g      ->
              let gend = g.gstart + g.length in
              begin match se with
              | `Start  ->
                if g.gstart = pos then
                  error "In %s start %d at start of gap %d" an pos g.gstart
                else if pos < gend then
                  error "In %s start %d inside alignment gap (%d,%d]" an pos g.gstart gend
                else if pos = gend then
                  Ok sp                             (* no sequence across gap. *)
                else (* pos > gend *)
                  loop sp tl
              | `End    ->
                if g.gstart = pos then
                  Ok sp
                else if pos <= gend then
                  error "In %s end %d inside or at end alignment gap (%d,%d]"
                    an pos g.gstart gend
                else (* pos > gend *)
                  loop sp tl
              end
          | Sequence s -> let n = String.length s.s in
                          if s.start = pos then
                            Ok sp
                          else if pos <= s.start + n then
                            Ok (sp + pos - s.start)
                          else (* pos > s.start + n *)
                            loop (sp + n) tl
          end
    in
    loop 0 nuclear.al_seq

  let nuclear_position_in_reference an ~length genetic =
    let rec loop sp = function
      | []             -> error "In %s not enough elements!" an
      | h :: tl ->
          begin match h with
          | Boundary b -> error "In %s %s should be removed!"
                            an (al_el_to_string h)
          | Start _
          | End _      -> error "In %s %s in reference!" an (al_el_to_string h)
          | Gap _      -> loop sp tl
          | Sequence s -> let n = String.length s.s in
                          let d = sp + n - length in
                          if d >= 0 then
                            Ok (s.start + d + genetic.offset)
                          else (* d < 0  -> sp + n < length *)
                            loop (sp + n) tl
          end
    in
    loop 0 genetic.al_seq

  let nuclear_position_to_genetic an se pos region =
    let { nuclear; genetic } = region in
    reference_sequence_position an se pos nuclear >>=
      fun length -> nuclear_position_in_reference an ~length genetic

  let update_hd has_data b l =
    List.fold_left l ~init:has_data ~f:(fun h e ->
      match e with
      | Start _     -> true
      | End _       -> false
      | Boundary _
      | Gap _
      | Sequence _  -> h)

  let drop_before p =
    let rec loop l = match l with
      | []      -> []
      | h :: t  -> if start_position h >= p then
                      l
                   else
                     loop t
    in
    loop

  let take_upto p =
    let rec loop l = match l with
      | []      -> []
      | h :: t  -> if start_position h > p then
                      []
                   else
                     h :: loop t
    in
    loop

  let between a b s =
    take_upto b (drop_before a s)

  let start_check an ~last_end ~reference ~genetic ~nuclear =
    match Split.at_start nuclear.al_seq with
    | [], (Start s :: tl) when s = last_end ->  (* Start immediately. *)
        let nnuclear = { nuclear with al_seq = tl } in
        Ok (Some nnuclear)
    | b, (Start s :: tl) ->
        if s = last_end then
          error "In %s found %s in between last end position and Start at %d"
            an (al_seq_to_string ~sep:"," b) s
        else
          (* b can contain gaps, but they're 'unknown' so lets drop them for
             what the imputed sequence contains. *)
          nuclear_position_to_genetic an `Start s reference >>= fun gstart ->
            nuclear_position_to_genetic an `End last_end reference >>= fun gend ->
              let nal_seq = between gend gstart genetic.al_seq in
              let shifted = shift_all_al_els nuclear.offset nal_seq in
              let nnuclear = { nuclear with al_seq = shifted @ tl } in
              Ok (Some nnuclear)
    | _, (e :: _) ->
        error "In %s Asked to split at Start not: %s" an (al_el_to_string e)
    | _n, [] -> Ok None

  let end_check allele_name ~nuclear =
    match Split.at_end nuclear.al_seq with
    | b, (End e :: tl) ->
        Ok (Some (b, e, {nuclear with al_seq = tl}))
    | _, (e :: _) ->
        error "In %s Asked to split at End not: %s"
          allele_name (al_el_to_string e)
    | ns, []  -> Ok None

  let check_for_start_ends allele_name has_data ~reference ~to_merge =
    let {nuclear; genetic } = to_merge in
    let rec check ?last_end has_data nuclear =
      if has_data then
        end_check allele_name ~nuclear >>= function
          | None -> Ok (true, nuclear)       (* Didn't find End *)
          | Some (before, last_end, after) ->
              if !debug then
                printf "an end! b: %s le: %d a: %s\n"
                  (al_seq_to_string ~sep:"," before)
                  last_end
                  (al_seq_to_string ~sep:"," after.al_seq);
              check ~last_end false after
                >>= fun (has_data, anuclear) ->
                  Ok (has_data
                     , { anuclear with
                          al_seq = (List.rev before) @ anuclear.al_seq })
      else (* not has_data *)
        let last_end = Option.value last_end ~default:nuclear.nbp in
        start_check allele_name ~last_end ~reference ~genetic ~nuclear >>= function
          | Some nuclear -> check true nuclear
          | None -> Ok (false, nuclear)                  (* Didn't find Start *)
    in
    check has_data nuclear

  (* Public *)
  (* Fold/map through the instructions and in MergeFromNuclear instructions
     where the nuclear sequence hasn't started insert from the genetic exon.

    Genetic, our source sequences _must_ have their full sequence imputed
    before being used. This is assured by [map_instr_to_alignments] this
    way we know that they have no start/end's inside the genetic instructions.  *)
  let insert_missing_data_in_merge_nuclear allele_name instr =
    let open Boundaries in
    let rec n has_data acc l = match l with
      | [] ->
          error "Mapping instructions for %s, empty!" allele_name
      | (FillFromGenetic f as h) :: tl when f.bm.label = UTR3 ->
          Ok (List.rev (h :: acc))                                     (* Fin *)
      | (FillFromGenetic _ as h) :: tl ->
          m has_data (h :: acc) tl
      | e :: _ ->
          error "Expected a Fill but received: %s"
            (instruction_to_string_seqless e)
    and m has_data acc l = match l with
      | MergeFromNuclear { reference; to_merge } :: tl ->
          check_for_start_ends allele_name has_data ~reference ~to_merge >>=
            fun (has_data, nuclear) ->
              let to_merge = {to_merge with nuclear } in
              n has_data (MergeFromNuclear { reference; to_merge} :: acc) tl
      | e :: _ ->
          error "Expected a Merge but received: %s"
            (instruction_to_string_seqless e)
      | [] ->
          error "Mapping instructions for %s, empty!" allele_name
    in
    n false [] instr

  let start_before_first_boundary = function
    | Boundary b :: Start s :: tl when b.pos = s ->
        Ok (Start s :: Boundary b :: tl)
    | lst ->
        error "Final instructions don't have Boundary Start, but: %s"
          (al_seq_to_string ~sep:";" lst)

  let exec_instructions l =
    List.map l ~f:(function
        | FillFromGenetic g               -> region_to_final_al_seq g
        | MergeFromNuclear {to_merge; _ } -> region_to_final_al_seq to_merge.nuclear)
    |> List.concat
    |> start_before_first_boundary

  (* Split the allele association lists into
     1. Elements in both.
     2. Elements just in "gen_assoc"
     3. Elements just in "nuc_assoc" *)
  let same_and_diff ~gen_assoc ~nuc_assoc =
    let gen_sort = List.sort ~cmp:compare gen_assoc in
    let nuc_sort = List.sort ~cmp:compare nuc_assoc in
    let rec loop s dg dn nlst glst = match glst with
      | []                -> s, dg, nlst @ dn
      | g :: gt ->
          let ga, gl = g in
          begin match nlst with
          | []            -> s, (glst @ dg), dn
          | n :: nt ->
              let na, nl = n in
              let r = compare na ga in
              if r = 0 then                                           (* same *)
                loop ((ga,gl,nl) :: s) dg dn nt gt
              else if r > 0 then            (* gen before nuc -> just gen ?!? *)
                loop s (g :: dg) dn nlst gt
              else (* r < 0 *)                  (* nuc before gen -> just nuc *)
                loop s dg (n :: dn) nt glst
          end
    in
    loop [] [] [] nuc_sort gen_sort

  let same_and_diff_warn_on_just_gen ~gen_mp ~nuc_mp =
    let same, just_gen, just_nuc =
        same_and_diff ~gen_assoc:gen_mp.Parser.alt_elems
          ~nuc_assoc:nuc_mp.Parser.alt_elems
    in
    if just_gen <> [] then
      eprintf "Found these alleles with only genetic data: %s\n"
        (String.concat ~sep:"; " (List.map just_gen ~f:fst));
    same, just_nuc

  let same_merge allele_name ~gen ~nuc instr =
    if !debug then printf "same_merge of %s\n%!" allele_name;
    map_instr_to_alignments allele_name ~gen ~nuc instr >>= fun i ->
      exec_instructions (strip_nuclear_start_ends allele_name i) >>= fun l ->
        Ok (allele_name, l)

  let same_alts instr same =
    list_fold_ok same ~init:[]
      ~f:(fun acc (allele_name, gen, nuc) ->
            same_merge allele_name ~gen ~nuc instr >>= fun p ->
              if !debug then
                match Parser.in_order_invariant (snd p) with
                | Ok ()       -> Ok (p :: acc)
                | Error (a,b) -> error "same_alts: %s out of order %s before %s"
                                  allele_name
                                  (al_el_to_string a) (al_el_to_string b)
              else
                Ok (p :: acc))

  let diff_merge allele_name ~gen ~nuc instr =
    map_instr_to_alignments allele_name ~gen ~nuc instr >>=
      insert_missing_data_in_merge_nuclear allele_name >>= fun i ->
        exec_instructions i >>= fun l ->
          if !debug then
            match Parser.in_order_invariant l with
            | Ok ()       -> Ok (allele_name, l)
            | Error (a,b) -> error "diff_merge %s out of order %s before %s"
                               allele_name (al_el_to_string a)
                               (al_el_to_string b)
          else
            Ok (allele_name, l)

  let to_distance_arguments nuc_mp gen_alt_alleles =
    let gs = string_set_of_list gen_alt_alleles in
    let nm = string_map_of_assoc nuc_mp.Parser.alt_elems in
    let candidates_without_ref, targets =
      StringMap.partition nm ~f:(fun key _al_seq -> StringSet.mem key gs)
    in
    let reference = nuc_mp.Parser.reference in
    let reference_sequence = nuc_mp.Parser.ref_elems in
    let candidates =
      StringMap.add candidates_without_ref ~key:reference
        ~data:(no_sequences reference_sequence)
    in
    { Distances.reference; reference_sequence; targets; candidates }

  let diff_alts dl instr nuc_mp gen_mp =
    let open Distances in
    let gen_alt_alleles = List.map ~f:fst gen_mp.Parser.alt_elems in
    let darg = to_distance_arguments nuc_mp gen_alt_alleles in
    let candidates_g =
      let gma = string_map_of_assoc gen_mp.Parser.alt_elems in
      StringMap.mapi darg.candidates ~f:(fun key _nal_seq ->
        if key = darg.reference then
          no_sequences gen_mp.Parser.ref_elems
        else
          StringMap.find key gma)
    in
    compute darg dl >>= fun dmap ->
      StringMap.bindings dmap
      |> list_fold_ok ~init:[] ~f:(fun acc (nuc_allele_name, dlst) ->
            let (gen_allele_name, _distance) = List.hd_exn dlst in
            let gen = StringMap.find gen_allele_name candidates_g in
            let nuc = StringMap.find nuc_allele_name darg.targets in
            if !debug then
              printf "diff_merge of %s with %s's introns.\n%!"
                nuc_allele_name gen_allele_name;
            diff_merge nuc_allele_name ~gen ~nuc instr >>=
              fun p -> Ok (p :: acc))

  let do_it ~gen_mp ~nuc_mp dl =
    let open Parser in
    align_mp_into_instructions ~gen_mp ~nuc_mp >>= fun instr ->
    exec_instructions instr >>= fun new_ref_elems ->
      if !debug then
        printf "new refererence: %s\n%!" (al_seq_to_string ~sep:"," new_ref_elems);
      Impute.do_it dl gen_mp >>= fun (igen_mp, _igen_assoc) ->
        let same, _just_nuc = same_and_diff_warn_on_just_gen ~gen_mp:igen_mp ~nuc_mp in
        same_alts instr same >>= fun same_alts ->
        diff_alts dl instr nuc_mp igen_mp >>= fun diff_alts ->
          Ok ({ align_date = gen_mp.align_date
              ; reference  = gen_mp.reference
              ; ref_elems  = new_ref_elems
              ; alt_elems  = same_alts @ diff_alts
              },
              [])

end (* Merge *)
