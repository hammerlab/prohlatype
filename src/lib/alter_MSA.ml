(* Altering (Merging/Imputing) of cDNA and gDNA alignment sequences. *)

open Util
open MSA
open Biology

(* The set of loci for which I've tested these algorithms.  *)
let supported_loci_lst =
  Nomenclature.[ A; B; C
               ; E ; F ; G ; HFE ; H ; J ; K ; L ; MICA ; MICB ; P ; TAP1 ; TAP2 ; T ; V ; W ; Y
               ]


let supported_loci l =
    List.mem ~set:supported_loci_lst l

(* Filter out the string sequence elements. *)
let no_sequences = List.filter ~f:(fun a -> not (is_sequence a))

(* Specify string containers. *)
let end_position = end_position String.length

let full_sequence name l =
  match l with
  | Start s :: _ ->
      begin match List.last l with
      | Some (End e) -> Ok (s, e)
      | _  -> error "%s doesn't end at the end!: %s"
                name (al_seq_to_string l)
      end
  | _      -> error "%s doesn't start right away!: %s"
                name (al_seq_to_string l)

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

  val one : bigger:string alignment_sequence
          -> smaller:string alignment_sequence
          -> (string alignment_sequence * Alteration.per_segment list
             , string) result

  val do_it : Distances.logic
            -> Parser.sequence_alignment
            -> (Parser.sequence_alignment, string) result

end (* Impute *) = struct

  let debug = ref false

  type blst = boundary list

  let lookup blst final_end ~start ~end_ =
    let rec loop acc start = function
      | []             -> error "Failed to find a boundary position for %d" start
      | b :: []        -> if end_ > final_end then
                            error "Asking to create a per_segment after %d final_end %d"
                              end_ final_end
                          else
                            Ok ({ Alteration.full = start = b.pos && end_ = final_end
                                ; type_ = b.label
                                ; start ; end_ } :: acc)
      | b1 :: b2 :: tl -> if b1.pos <= start && start < b2.pos then begin
                            if end_ <= b2.pos then
                              Ok ({ Alteration.full = start = b1.pos && end_ = b2.pos
                                  ; type_ = b1.label
                                  ; start ; end_ } :: acc)
                            else (* end_ > b2.pos *)
                              let nacc =
                                { Alteration.full = start = b1.pos
                                ; type_ = b1.label
                                ; start
                                ; end_ = b2.pos } :: acc
                              in
                              loop nacc b2.pos (b2 :: tl)
                          end else
                            loop acc start (b2 :: tl)
    in
    loop [] start blst

  (* Keep track of the positions that we're imputing into the smaller sequence. *)
  let add_per_segment blst final_end ~start ~end_ ps_acc =
    let open Alteration in
    let first = (List.hd_exn blst).pos in
    if start = end_ then begin
      let not_start_nor_end = start <> first && start <> final_end in
      if not_start_nor_end && !debug then
        printf "Tried to add empty per_segment at %d\n" start;
      Ok ps_acc
    end else
      lookup blst final_end ~start ~end_ >>= fun nacc ->
        Ok (nacc @ ps_acc)

  let reverse_and_remove_boundary_duplicates lst =
    let rec loop a l = match l with
      | Boundary b1 :: ((Boundary b2 :: _) as tl) when b1 = b2 -> loop a tl
      | h :: t  -> loop (h :: a) t
      | []      -> a
    in
    loop [] lst

  let trim_same acc = function
    | []  -> acc
    | h :: _  ->
        match acc with
        | ha :: tl when ha = h -> tl
        | _  -> acc

  let one ~bigger ~smaller =
    let blst =
      List.filter_map bigger ~f:(function
          | Boundary b                           -> Some b
          | Start _ | End _ | Gap _ | Sequence _ -> None)
    in
    full_sequence "impute bigger sequence" bigger >>= fun (bigger_start, bigger_end) ->
      let aps = add_per_segment blst bigger_end in                             (* Curry *)
      (* smaller is just past it's start position. *)
      let rec merge_until_end merged ps_acc smaller_start_pos bigger ~smaller =
        if smaller_start_pos >= bigger_end then
          error "smaller start pos %d is at the bigger end %d"
            smaller_start_pos bigger_end
        else
          (* Take bigger elements up to the start position, excluding boundaries. *)
          let start_with, rest_b =
            Split.al_el_at ~pos:smaller_start_pos ~excluding:false ~acc:merged bigger in
          aps ~start:bigger_start ~end_:smaller_start_pos ps_acc >>= fun nps_acc ->
            (* Take smaller elements up to, excluding, next end. *)
            (*let smaller_without_same = trim_same start_with smaller in*)
            let without_end, rest_s = Split.at_end ~acc:start_with smaller in
            let smaller_end_pos = end_position (List.hd_exn rest_s) in
            (* Drop bigger elements up to the smaller end. *)
            let _ignore_me, b_after_s_end =
              Split.al_el_at ~pos:smaller_end_pos ~excluding:true ~acc:[] rest_b in
            (* Lastly, move smaller upto a potential next start. *)
            let _before_another_start, smaller_at_start_or_empty = Split.at_start rest_s in
            start_merging smaller_end_pos without_end nps_acc
              ~smaller:smaller_at_start_or_empty
              ~bigger:b_after_s_end
            (* Recursive to handle the End -> Start (Unkown-Gap in sequence) cases. *)
      and start_merging last_smaller_pos acc ps_acc ~smaller ~bigger =
        match smaller with
        | Start p :: t -> merge_until_end acc ps_acc p bigger ~smaller:t
        | []           -> if !debug then
                            printf "ending impute with: last_smaller_pos: %d bigger: %s\n%!"
                              last_smaller_pos (al_seq_to_string bigger);
                          let final_seq =
                            List.rev_append bigger acc
                            |> reverse_and_remove_boundary_duplicates
                          in
                          begin
                            match bigger with
                            | []     ->
                                Ok (final_seq, List.rev ps_acc)
                            | h :: t ->
                                let end_ = Option.value_map (List.last t)
                                            ~f:end_position
                                              ~default:(end_position h)
                                in
                                aps ~start:last_smaller_pos ~end_ ps_acc >>= fun nps_acc ->
                                  Ok (final_seq, (List.rev nps_acc))               (* Fin *)
                          end
        | ae :: _      -> error "Not at Start or Empty: %s"
                            (al_el_to_string ae)
      in
      let _before_smaller_start, at_smaller_start = Split.at_start smaller in
      match at_smaller_start with
      | Start p :: t -> merge_until_end [] [] p bigger ~smaller:t
      | []           -> error "Empty smaller sequence!"
      | ae :: _      -> error "Not at Start: %s" (al_el_to_string ae)

  (* TODO: This logic is much simpler if the Distance logic is reference:
     move an allele's Start/End's to the references's. We don't need to
     go the round trip of saying that the distance is 0. *)
  let do_it logic mp =
    let open Parser in
    let module Sm = StringMap in
    let reflength = sequence_length mp.ref_elems in
    let to_fill, full = List.partition_map mp.alt_elems ~f:(fun alt ->
      let l = sequence_length alt.seq in
      if l < reflength then `Fst (l, alt) else `Snd alt)
    in
    (* Strip the sequences from the reference since for alternate alleles
       they represent a difference to the reference; therefore without the
       sequences the reference has no difference to itself. *)
    let init = Sm.singleton mp.reference (no_sequences mp.ref_elems) in
    let candidates =
      List.fold_left full ~init ~f:(fun m alt ->
        Sm.add ~key:alt.allele ~data:alt.seq m)
    in
    let to_distances =
      Distances.one ~reference:mp.reference ~reference_sequence:mp.ref_elems
        logic
    in
    List.sort to_fill ~cmp:(fun (l1,_) (l2,_) -> compare l2 l1)
    |> list_fold_ok ~init:full ~f:(fun alt_list (length, alt) ->
        to_distances ~candidates ~allele:(alt.allele, alt.seq) >>= function
        | []            ->
            error "Did not return _any_ distances to %s allele" alt.allele
        | (key, d) :: _ ->
            if !debug then
              printf "Imputing %s with %s d: %f because %d < %d\n%!"
                alt.allele key d length reflength;
            let bigger = Sm.find key candidates in
            one ~smaller:alt.seq ~bigger >>= fun (merged, positions) ->
              let a =
                { Alteration.allele = key
                ; why = "imputing missing gDNA"
                ; distance = d
                ; positions
                }
              in
              Ok ({ alt with alters = [a]; seq = merged } :: alt_list))
      >>= fun nalt_elems ->
        Ok ({ mp with alt_elems = sort_alts_by_nomenclature nalt_elems})

end (* Impute *)

(*
Problem: There are many (10x) more cDNA derived (occasionally referred to as
nucleic following IMGT's naming these alignment files Gene_nuc.txt) allele
sequences than gDNA (occasionally referred to as genetic also following IMGT's
convention of naming these alignment files Gene_gen.txt) ones. But since the
gDNA sequences are longer they allow us to match more reads. We would like
to "impute" sequences for (and hence a graph or PHMM containing) all of the
alleles for which we have some cDNA sequence data. To accomplish this, we use
the current cDNA sequences as scaffolding around which we append information
from the gDNA ones.

The algorithm proceeds in several steps. At a high level for a given gene
(ex. HLA-A) we:

  1. Use the reference allele to construct a set of instructions for how to
     zip the exons into the introns.
  2. Afterwards we use those instructions to merge the reference.
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

  3. Before merging the other alleles we separate and check 3 cases using
      [same_warn_on_just_gen]:
      - Alleles where we have only gDNA data (4).
      - Alleles where we have both cDNA and gDNA data (6).
      - Alleles where we have only cDNA data (6).

  4. It would be surprising if we had alleles where we have only gDNA data,
     (since one could presumably infer the nucleic components), so an error
     message is written, but no action is taken, and those alleles are ignored.

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
      : gen_mp:Parser.sequence_alignment
      -> nuc_mp:Parser.sequence_alignment
      -> (instructions, string) result

  val map_instr_to_alignments
      : bigger_gDNA:string alignment_sequence
      -> bigger_cDNA:string alignment_sequence
      -> smaller_cDNA:string alignment_sequence
      -> string
      -> instructions
      -> (instructions, string) result

  val insert_missing_data_in_merge_nuclear
      : string
      -> instructions
      -> (instructions, string) result

  val exec_instructions
      : instructions
      -> ((string alignment_sequence * Alteration.per_segment list), string) result

  val to_distance_arguments
      : Parser.sequence_alignment
      -> string list
      -> Distances.arg

  val do_it
      : gen_mp:Parser.sequence_alignment
      -> nuc_mp:Parser.sequence_alignment
      -> Distances.logic
      -> (Parser.sequence_alignment, string) result

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

  let region_to_per_segment r =
    { Alteration.start = r.nbp
    ; end_ = r.nbp + r.bm.Boundaries.length
    ; type_ = r.bm.Boundaries.label
    ; full = true
    }

  let region_to_string al_seq_to_string { bm; nbp; offset; al_seq } =
    sprintf "{ bm: %s; nbp: %d; offset: %d; al_seq: %s}"
      (Boundaries.marker_to_string bm) nbp offset (al_seq_to_string al_seq)

  let region_to_final_al_seq r =
    (Boundaries.to_boundary ~offset:r.offset r.bm) :: r.al_seq

  type 'a instruction =
    | FillFromGenetic of 'a region
    (* Take the elements from the genetic region and ignore the nuclear region.
       This is the case for intron's or UTR's; where we don't have any
       alignment sequence data to construct a sequence. We still have to
       adjust to the correct position based off of the region's offset.

       All sequences that can be used as a source _must_ have their full
       sequence imputed beforehand; otherwise [map_instr_to_alignments] will
       fail. *)
    | MergeFromNuclear of { nuclear : 'a region          (* Coming from cDNA. *)
                          ; genetic : 'a region          (* Coming from gDNA. *)
                          ; pss     : Alteration.per_segment list
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
    | MergeFromNuclear { nuclear; genetic; pss} ->
        sprintf "MergeFromNuclear {nuclear: %s; genetic: %s; pss: %s}"
          (region_to_string seq_to_string nuclear)
          (region_to_string seq_to_string genetic)
          (Alteration.per_segment_list_to_string pss)

  let instruction_to_string_seqless =
    instruction_to_string (fun _ -> "")

  let bounded_results_to_string =
    string_of_list ~sep:";"
      ~f:(fun (bm, seq) -> sprintf "(%s, %s)"
            (Boundaries.marker_to_string bm) (short_seq seq))

  let shift_al_el offset = function
    | Boundary b -> Boundary { b with pos = b.pos + offset}
    | Start s    -> Start (s + offset)
    | End e      -> End (e + offset)
    | Gap g      -> Gap { g with gstart = g.gstart + offset }
    | Sequence s -> Sequence { s with start = s.start + offset}

  let shift_all_al_els offset =
    List.map ~f:(shift_al_el offset)

  let new_region ?offset al_seq new_boundary_pos bm =
    let offset =
      Option.value offset ~default:(new_boundary_pos - bm.Boundaries.position)
    in
    { bm = bm
    ; nbp = new_boundary_pos
    ; offset
    ; al_seq = shift_all_al_els offset al_seq
    }

  type s_seq = string MSA.alignment_sequence

  type instructions = s_seq instruction list

  (* Create a framework of instructions based off of, the alignment in the
     reference sequence. *)
  let zip_align ~gen ~nuc =
    let ref_seq = reference_sequence_from_ref_alignment_elements in
    let open Boundaries in
    (* Convert the reference alignment_sequence's to Boundary delimited strings,
       this way we can make sure that the sequences are the same. *)
    grouped (all_boundaries_before_start_or_end gen) >>= begin fun gblst ->
      grouped (all_boundaries_before_start_or_end nuc) >>= begin fun nblst ->
        let rec need_exon new_boundary_pos acc gl nl = match gl, nl with
          | ({ label = Exon gx} as gb, gseq) :: gtl
          , ({ label = Exon nx} as nb, nseq) :: ntl ->
            let gss = ref_seq gseq in
            let nss = ref_seq nseq in
            if gx <> nx then
              error "Genetic Exon %d doesn't line up with Nuclear exon %d" gx nx
            else if gss <> nss then
              error "Genetic sequence disagrees with Nucleic at Exon %d:\n %s" gx
                (manual_comp_display ~width:160 ~labels:("Genetic ", "Nuclear ")
                  gss nss)
            else
              let nuclear = new_region nseq new_boundary_pos nb in
              (* When constructing other sequences's we're going to take the nuclear
                (from cDNA) sequence of the bigger allele (to represent the true
                "genetic" region. This is a pretty bad mis-use of names that needs
                to be fixed with better naming. But for now it is important to
                remember that we want the nuclear offset. *)
              let genetic = new_region ~offset:nuclear.offset gseq new_boundary_pos gb in
              let instr =
                MergeFromNuclear { genetic; nuclear; pss = [] }
              in
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
      end
    end

  let strip_start_ends_seq =
    List.filter ~f:(fun a -> (is_gap a) || (is_sequence a))

  let strip_start_ends_region r =
    { r with al_seq = strip_start_ends_seq r.al_seq }

  let strip_nuclear_start_ends allele_name instr =
    List.map instr ~f:(fun e ->
      match e with
      | FillFromGenetic _ -> e
      | MergeFromNuclear { genetic; nuclear; pss } ->
          let nnuclear = strip_start_ends_region nuclear in
          MergeFromNuclear { genetic; nuclear = nnuclear ; pss })

  (* Public *)
  let align_mp_into_instructions ~gen_mp ~nuc_mp =
    let open Parser in
    zip_align ~gen:gen_mp.ref_elems ~nuc:nuc_mp.ref_elems
    >>= fun instr -> Ok (strip_nuclear_start_ends gen_mp.reference instr)

  let replace_al_seq r s =
    { r with al_seq = shift_all_al_els r.offset s }

  let full_sequence_no_op msg s =
    full_sequence msg s >>= fun (_, _) -> Ok ()

  (* Public *)
  (* Take sequences
      bigger_gDNA:  gDNA sequence of Bigger/full sequence:
                     we'll use these Introns/UTR's.
      bigger_cDNA:  cDNA sequence of Bigger/full sequence:
                    we'll use this sequences' Exons if there are
                    Start/End's inside of:
      smaller_cDNA: cDNA sequence of Smaller/incomplete sequence.

    Then
      1. insert them into their relevant instructions
      2. shift them as necessary.
  *)
  let map_instr_to_alignments ~bigger_gDNA ~bigger_cDNA ~smaller_cDNA
    allele_name instr =
    let open Boundaries in
    (* Convert the reference alignment_sequence's to Boundary delimited strings,
       this way we can make sure that the sequences are the same. *)
    let to_msg = sprintf "Merge.map_instr_to_alignments of %s %s sequence" allele_name in
    full_sequence_no_op (to_msg "bigger gDNA") bigger_gDNA >>= fun () ->
    full_sequence_no_op (to_msg "bigger cDNA") bigger_cDNA >>= fun () ->
    grouped (all_boundaries_before_start_or_end bigger_gDNA) >>= fun bigger_gDNA ->
    grouped (all_boundaries_before_start_or_end bigger_cDNA) >>= fun bigger_cDNA ->
    grouped (all_boundaries_before_start_or_end smaller_cDNA) >>= (fun smaller_cDNA ->
      let rec fill_from_g acc fr itl ~bigger_gDNA ~bigger_cDNA ~smaller_cDNA =
        match bigger_gDNA with
        | (gbm, gseq) :: gtl when gbm.label = fr.bm.label ->
            let nacc = FillFromGenetic (replace_al_seq fr gseq) :: acc in
            begin match itl with
            | []  -> Ok (List.rev nacc)                                (* Fin *)
            | _   -> need_merge nacc itl ~bigger_gDNA:gtl ~bigger_cDNA ~smaller_cDNA
            end
        | []            -> error "Empty genetic sequence"
        | (gbm, _) :: _ -> error "Genetic boundary label %s doesn't match \
                            instruction's: %s"
                            (Gene_region.to_string gbm.label)
                            (Gene_region.to_string fr.bm.label)

      and fill_from_m acc itl ~genetic ~nuclear ~bigger_gDNA ~bigger_cDNA ~smaller_cDNA =
        match bigger_gDNA, bigger_cDNA, smaller_cDNA with
        | (bgbm, _bgseq) :: bgtl
        , (bcbm, bcseq) :: bctl
        , (sbm, sseq) :: stl
            (* the genetic and nuclear labels already match *)
            when bgbm.label = genetic.bm.label
              && bcbm.label = genetic.bm.label
              && sbm.label = genetic.bm.label ->
            let nacc =
              let genetic = replace_al_seq genetic (strip_start_ends_seq bcseq) in
              let nuclear = replace_al_seq nuclear sseq in
              MergeFromNuclear { genetic; nuclear; pss = [] }:: acc
            in
            need_fill nacc itl ~bigger_gDNA:bgtl ~bigger_cDNA:bctl ~smaller_cDNA:stl
        | [],  _,  _    -> error "%s has premature empty bigger gDNA sequence"
                            allele_name
        |  _, [],  _    -> error "%s has premature empty bigger cDNA sequence"
                            allele_name
        |  _,  _, []    -> error "%s has premature empty smaller cDNA sequence"
                            allele_name
        | (bgbm, _) :: _
        , (bcbm, _) :: _
        , (sbm, _) :: _ -> error "Bigger gDNA %s or Bigger cDNA %s or Smaller \
                            cDNA %s boundary labels don't match instruction's: %s"
                            (Gene_region.to_string bgbm.label)
                            (Gene_region.to_string bcbm.label)
                            (Gene_region.to_string sbm.label)
                            (Gene_region.to_string genetic.bm.label)

      and need_fill acc ilst ~bigger_gDNA ~bigger_cDNA ~smaller_cDNA = match ilst with
        | FillFromGenetic f :: tl
                -> fill_from_g acc f tl ~bigger_gDNA ~bigger_cDNA ~smaller_cDNA
        | []     -> error "Premature end to instructions for %s" allele_name
        | s :: _ -> error "Mapping instructions for %s, did not start with a \
                          \"starting\" FillFromGenetic but %s"
                      allele_name (instruction_to_string_seqless s)
      and need_merge acc ilst ~bigger_gDNA ~bigger_cDNA ~smaller_cDNA = match ilst with
        | MergeFromNuclear { genetic; nuclear; _ } :: tl ->
            fill_from_m acc tl ~genetic ~nuclear ~bigger_gDNA ~bigger_cDNA ~smaller_cDNA
        | i :: _ -> error "Didn't find a merge from nuclear instruction: %s"
                      (instruction_to_string_seqless i)
        | []     -> error "Premature end to instructions for %s" allele_name
      in
      need_fill [] instr ~bigger_gDNA ~bigger_cDNA ~smaller_cDNA)

   let drop_before pos l =
     let _forget, rest = Split.al_el_at l ~pos ~excluding:false ~acc:[] in
     rest

  let take_upto pos l =
    let reverse_me, _forget = Split.al_el_at l ~pos ~excluding:true ~acc:[] in
    List.rev reverse_me

  let between a b s =
    take_upto b (drop_before a s)

  (* start_check will only trim the Start not a Sequence immediately
     after a Start that might be duplicated. *)
  let check_same_and_append ~before ~after =
    match (List.rev before) with
    | []        -> after
    | lb :: tb  ->
        match after with
        | ha :: _ when ha = lb -> (List.rev tb) @ after
        | _                    -> before @ after

  let check_same_and_append_rev ~before ~after =
    match before with
    | []        -> after
    | hb :: tb  ->
        match after with
        | ha :: _ when ha = hb -> (List.rev tb) @ after
        | _                    -> List.rev before @ after

  let start_check an ~last_end ~genetic ~nuclear =
    match Split.at_start nuclear.al_seq with
    | [], (Start s :: tl) when s = last_end ->  (* Start immediately. *)
        Ok (true, None, { nuclear with al_seq = tl })
    | b, (Start s :: tl) ->
        if s = last_end then
          error "In %s found %s in between last end position and Start at %d"
            an (al_seq_to_string b) s
        else
          (* b can contain gaps, but they're 'unknown' so lets drop them for
             what the imputed sequence contains. *)
          let pss =
            { Alteration.start = last_end
            ; end_  = s
            ; full  = false
            ; type_ = nuclear.bm.Boundaries.label
            }
          in
          let gal_seq = between last_end s genetic.al_seq in
          let nal_seq = check_same_and_append ~before:gal_seq ~after:tl in
          Ok (true, Some pss, { nuclear with al_seq = nal_seq })
    | _, (e :: _) ->
        error "In %s Asked to split at Start not: %s" an (al_el_to_string e)
    | _n, [] -> (* No start borrow everything from genetic. *)
        let gal_seq = drop_before last_end genetic.al_seq in
        let pss =
          { Alteration.start = last_end
          ; end_  = nuclear.nbp + nuclear.bm.Boundaries.length
          ; full  = true
          ; type_ = nuclear.bm.Boundaries.label
          }
        in
        Ok (false, Some pss, { nuclear with al_seq = gal_seq })

  let end_check allele_name ~nuclear =
    match Split.at_end nuclear.al_seq with
    | b, (End e :: tl) ->
        Ok (Some (b, e, { nuclear with al_seq = tl }))
    | _, (e :: _) ->
        error "In %s Asked to split at End not: %s"
          allele_name (al_el_to_string e)
    | ns, []  -> Ok None

  let check_for_start_ends allele_name has_data ~genetic ~nuclear =
    let rec check ?last_end has_data pss_acc nuclear =
      if has_data then
        end_check allele_name ~nuclear >>= function
          | None -> Ok (true, List.rev pss_acc, nuclear)   (* Didn't find End *)
          | Some (before, last_end, after) ->
              if !debug then
                printf "an end in an exon! b: %s le: %d a: %s\n"
                  (al_seq_to_string before)
                  last_end
                  (al_seq_to_string after.al_seq);
              check ~last_end false [] after
                >>= fun (has_data, pss_acc, anuclear) ->
                  let nal_seq = check_same_and_append_rev ~before ~after:anuclear.al_seq in
                  Ok ( has_data
                     , List.rev pss_acc
                     , { anuclear with al_seq = nal_seq })
      else (* not has_data *)
        (* We'll have non-Exon data, so grab only from Boundary position. *)
        let last_end = Option.value last_end ~default:nuclear.nbp  in
        start_check allele_name ~last_end ~genetic ~nuclear >>=
          fun (found_start, pss, nuclear) ->
            let pss_s, npss_acc =
              match pss with
              | None  -> "None", pss_acc
              | Some ps -> (sprintf "Some %s" (Alteration.per_segment_to_string ps))
                           , ps ::pss_acc
            in
            if found_start then begin
              if !debug then
                printf "a start in an exon! le: %d, pss opt: %s, nuclear seq:%s\n"
                  last_end pss_s (al_seq_to_string nuclear.al_seq);
              check true npss_acc nuclear
            end else begin
              if !debug then
                printf "Didn't find start in %s; last_end: %d took everything from genetic:%s? \t %s \n"
                  (Gene_region.to_string nuclear.bm.Boundaries.label)
                  last_end
                  (al_seq_to_string nuclear.al_seq)
                  (al_seq_to_string genetic.al_seq);
              Ok (false, List.rev npss_acc, nuclear)
            end
    in
    check has_data [] nuclear

  (* Public *)
  (* Fold/map through the instructions and in MergeFromNuclear instructions
     where the nuclear sequence hasn't started insert from the genetic exon.

    Genetic, our source sequences _must_ have their full sequence imputed
    before being used. This is assured by [map_instr_to_alignments] this
    way we know that they have no start/end's inside the genetic instructions.  *)
  let insert_missing_data_in_merge_nuclear allele_name instr =
    let open Boundaries in
    let rec n nuc_has_data acc l = match l with
      | [] ->
          error "Mapping instructions for %s, empty!" allele_name
      | (FillFromGenetic f as h) :: tl when f.bm.label = UTR3 ->
          Ok (List.rev (h :: acc))                                     (* Fin *)
      | (FillFromGenetic _ as h) :: tl ->
          m nuc_has_data (h :: acc) tl
      | e :: _ ->
          error "Expected a Fill but received: %s"
            (instruction_to_string_seqless e)
    and m nuc_has_data acc l = match l with
      | MergeFromNuclear { genetic; nuclear } :: tl ->
          check_for_start_ends allele_name nuc_has_data ~genetic ~nuclear >>=
            fun (nuc_has_data, pss, nnuclear) ->
              let nacc = MergeFromNuclear { genetic; nuclear = nnuclear; pss} :: acc in
              n nuc_has_data nacc tl
      | e :: _ ->
          error "Expected a Merge but received: %s"
            (instruction_to_string_seqless e)
      | [] ->
          error "Mapping instructions for %s, empty!" allele_name
    in
    n false [] instr

  let start_before_first_boundary l = match l with
    | Boundary b :: Start s :: tl when b.pos = s ->
        Ok (Start s :: Boundary b :: tl)
    | Start s :: Boundary b :: _ when b.pos = s ->
        Ok l
    | lst ->
        error "First instructions does not have Boundary Start, but: %s"
          (al_seq_to_string lst)

  let exec_instructions l =
    List.map l ~f:(function
        | FillFromGenetic g                -> region_to_final_al_seq g
                                              , [region_to_per_segment g]
        | MergeFromNuclear {nuclear; pss } -> region_to_final_al_seq nuclear
                                              , pss )
    |> List.split
    |> fun (instr, pss) ->
        start_before_first_boundary (List.concat instr) >>= fun l ->
          Ok (l, (List.concat pss))

  let impute_nuc_if_necessary nmp galt n_seq =
    let open Parser in
    let reflength = sequence_length nmp.ref_elems in
    let impute bigger_name distance bigger =
      if !debug then
        printf "imputing nuc of %s with %s bigger seq: %s\n"
          galt.Parser.allele bigger_name (al_seq_to_string bigger);
      Impute.one ~smaller:n_seq ~bigger >>= fun (seq, positions) ->
        let alter =
          { Alteration.allele = bigger_name
          ; why = "imputing missing components"
          ; distance
          ; positions
          }
        in
        start_before_first_boundary seq >>= fun sseq ->
          Ok (sseq, [alter])
    in
    match galt.Parser.alters with
    | []                                     ->
        if !debug then
          printf "Not imputing nuc of %s: %s\n"
            galt.Parser.allele (al_seq_to_string n_seq);
        Ok (n_seq, [])
    | { Alteration.allele; distance; _} :: _ ->
        let nlength = sequence_length n_seq in
        if nlength = reflength then   (* We had to impute gDNA, but not cDNA! *)
          Ok (n_seq, [])
        else
          if allele = nmp.Parser.reference then
            impute allele distance (no_sequences nmp.Parser.ref_elems)
          else
            match List.find nmp.Parser.alt_elems
                    ~f:(fun alt -> alt.Parser.allele = allele)
            with
            | None     -> error "Didn't find %s in nmp" allele
            | Some alt -> impute allele distance alt.Parser.seq

  (* Split the allele association lists into
     1. Elements in both.
     2. Elements just in "gen"
     3. Elements just in "nuc"


  But remember that the sequence in "gen" (from gDNA) are imputed
  (missing UTR's and maybe end Exons'): since they're full we have to make
  sure that their "nuc" (from cDNA) sequences have the same imputation.
  *)
  let same_and_diff ~gen ~nuc =
    let open Parser in
    let rec loop s dg dn ~nlst ~glst = match glst with
      | []                -> Ok (s, dg, nlst @ dn)
      | g :: gt ->
          begin match nlst with
          | []            -> Ok (s, (glst @ dg), dn)
          | n :: nt ->
              let r = compare n.allele g.allele in
              if r = 0 then                                           (* same *)
                impute_nuc_if_necessary nuc g n.seq >>= fun nseq ->
                  loop ((g.allele, g.seq, nseq) :: s) dg dn ~nlst:nt ~glst:gt
              else if r > 0 then            (* gen before nuc -> just gen ?!? *)
                loop s (g :: dg) dn ~nlst ~glst:gt
              else (* r < 0 *)                  (* nuc before gen -> just nuc *)
                loop s dg (n :: dn) ~nlst:nt ~glst
          end
    in
    let cmp alt1 alt2 = compare alt1.allele alt2.allele in
    loop [] [] [] (List.sort ~cmp nuc.alt_elems) (List.sort ~cmp gen.alt_elems)

  let same_warn_on_just_gen ~igen_mp ~nuc_mp =
    same_and_diff ~gen:igen_mp ~nuc:nuc_mp >>= fun (same, just_gen, _just_nuc) ->
      if just_gen <> [] then
        eprintf "Found these alleles with only genetic data: %s\n"
          (string_of_list ~sep:"; " ~f:(fun a -> a.Parser.allele) just_gen);
      Ok same

  let same_merge allele_name instr ~gDNA ~cDNA =
    if !debug then printf "same_merge of %s%!\n" allele_name;
    map_instr_to_alignments allele_name instr
      (* For the same sequence the cDNA difference doesn't matter! *)
      ~bigger_gDNA:gDNA ~bigger_cDNA:cDNA ~smaller_cDNA:cDNA >>= fun i ->
        exec_instructions (strip_nuclear_start_ends allele_name i) >>= fun l ->
        Ok (allele_name, l)

  let check_invariant where allele seq alters acc =
    let new_alt = { Parser.allele = allele; seq; alters } in
    if !debug then
      match Parser.in_order_invariant seq with
      | Error (a,b) -> error "%s: %s out of order %s before %s"
                        where
                        allele (al_el_to_string a) (al_el_to_string b)
      | Ok ()       -> Ok (new_alt :: acc)
    else
      Ok ( new_alt :: acc)

  let same_alts instr same =
    list_fold_ok same ~init:[]
      ~f:(fun acc (allele_name, gDNA, (cDNA, cDNA_alter_lst)) ->
            same_merge allele_name instr ~gDNA ~cDNA >>=
              fun (allele_name, (al_seq, pss_lst)) ->
              let new_alter =
                match pss_lst with
                | []  -> cDNA_alter_lst
                | _   ->
                    let exon_insert = List.exists pss_lst
                      ~f:(function { Alteration.type_ = Exon _ } -> true | _ -> false)
                    in
                    if exon_insert then
                      printf "when same merging %s: exon insert! %s\n"
                        allele_name (Alteration.per_segment_list_to_string pss_lst);
                    { Alteration.allele = allele_name
                    ; why = "reforming alignment"
                    ; distance = 0.
                    ; positions = pss_lst
                    } :: cDNA_alter_lst
              in
              check_invariant "same_alts" allele_name al_seq new_alter acc)

  let diff_merge gen_allele_name distance nuc_allele_name instr
    ~bigger_gDNA ~bigger_cDNA ~smaller_cDNA
    bigger_c_alters
    acc =
    map_instr_to_alignments nuc_allele_name instr ~bigger_gDNA ~bigger_cDNA ~smaller_cDNA >>=
      fun i ->
        if !debug then printf "mapped instructions: %s\n"
          (string_of_list ~sep:";" ~f:(instruction_to_string al_seq_to_string) i);
      insert_missing_data_in_merge_nuclear nuc_allele_name i >>= fun i ->
        exec_instructions i >>= fun (al_seq, pss_lst) ->
          if pss_lst = [] then
            eprintf "empty pss list when merging %s with %s (nuc)!\n"
            gen_allele_name nuc_allele_name;
          let new_alter =
            { Alteration.allele = gen_allele_name
            ; why = "imputing missing non-exonic regions"
            ; distance
            ; positions = pss_lst
            } :: bigger_c_alters
          in
          if !debug then
            printf "seq: %s\nalters: %s\n" (al_seq_to_string al_seq)
              (string_of_list ~sep:";" ~f:Alteration.to_string new_alter);
          check_invariant "diff_merge" nuc_allele_name al_seq new_alter acc

  let to_distance_arguments nuc_mp gen_alt_alleles =
    let gs = string_set_of_list gen_alt_alleles in
    let nm =
      nuc_mp.Parser.alt_elems
      |> List.map ~f:(fun alt -> (alt.Parser.allele, alt.Parser.seq))
      |> string_map_of_assoc
    in
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

  let diff_alts dl instr nuc_mp same ref_g ref_n =
    let open Distances in
    let gen_alt_alleles = List.map ~f:(fun (a, _, _) -> a) same in
    let darg = to_distance_arguments nuc_mp gen_alt_alleles in
    let candidates_g =
      let gma =
        List.map ~f:(fun (a, gs, nsp) -> (a, (gs, nsp))) same
        |> string_map_of_assoc
      in
      StringMap.mapi darg.candidates ~f:(fun key _nal_seq ->
        if key = darg.reference then
          (no_sequences ref_g
          , (no_sequences ref_n, []))
        else
          StringMap.find key gma)
    in
    compute darg dl >>= fun dmap ->
      StringMap.bindings dmap
      |> list_fold_ok ~init:[] ~f:(fun acc (nuc_allele_name, dlst) ->
          let (gen_allele_name, distance) = List.hd_exn dlst in
          let bigger_gDNA, (bigger_cDNA, bigger_c_alters) =
            StringMap.find gen_allele_name candidates_g
          in
          let smaller_cDNA = StringMap.find nuc_allele_name darg.targets in
          if !debug then
            printf "diff_merge of %s with %s's introns.\n%!"
              nuc_allele_name gen_allele_name;
          diff_merge gen_allele_name distance nuc_allele_name instr
            ~bigger_gDNA ~bigger_cDNA ~smaller_cDNA
            bigger_c_alters
            acc)

  let do_it ~gen_mp ~nuc_mp dl =
    let open Parser in
    align_mp_into_instructions ~gen_mp ~nuc_mp >>= fun instr ->
    exec_instructions instr >>= fun (new_ref_elems, ref_pss_lst) ->
      if !debug then
        printf "new refererence: %s\nreference pss lst: %s\n%!"
          (al_seq_to_string new_ref_elems)
          (Alteration.per_segment_list_to_string ref_pss_lst);
      Impute.do_it dl gen_mp >>= fun igen_mp ->
        same_warn_on_just_gen ~igen_mp ~nuc_mp >>= fun same ->
        same_alts instr same >>= fun same_alt_elems ->
          let ref_g = igen_mp.ref_elems in
          let ref_n = nuc_mp.ref_elems in
          if !debug then printf "finished same_alts\n";
          diff_alts dl instr nuc_mp same ref_g ref_n >>= fun diff_alt_elems ->
            let nalt_elems = same_alt_elems @ diff_alt_elems in
            Ok ({ release    = gen_mp.release
                ; align_date = gen_mp.align_date
                ; locus      = gen_mp.locus
                ; reference  = gen_mp.reference
                ; ref_elems  = new_ref_elems
                ; alt_elems  = sort_alts_by_nomenclature nalt_elems
                })

end (* Merge *)
