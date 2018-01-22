
open Util
open ParPHMM
open ParPHMM_drivers

let f_of_yojson =
  (Output.of_yojson Multiple_loci.final_read_info_of_yojson)

let of_json_file f =
  Yojson.Safe.stream_from_file f
  |> Stream.next (* only one element per file. *)
  |> f_of_yojson
  |> unwrap

type read_info =
  | Pr of Alleles_and_positions.t Multiple_loci.paired
  | Soi of Alleles_and_positions.t Multiple_loci.single_or_incremental
  [@@deriving show]

let read_metric of_alp_list =
  let open Pass_result in
  let of_aalp_pr = function
    | Filtered _    -> invalid_argf "read was filtered ?!?"
    | Completed alp -> of_alp_list alp
  in
  let take_regular r c = Alleles_and_positions.descending_cmp r c <= 0 in
  let mlo fp =
    Orientation.most_likely_between ~take_regular fp
    |> Pass_result.map ~f:snd
  in
  function
  | Soi (Multiple_loci.SingleRead or_) ->
      Sp.Single (of_aalp_pr (mlo or_))
  | Soi (Multiple_loci.PairedDependent pd) ->
      let p1 = of_alp_list pd.Multiple_loci.first in
      let p2 = of_aalp_pr pd.Multiple_loci.second in
      Sp.Paired (min p1 p2, max p1 p2)
  | Pr (Multiple_loci.FirstFiltered ff) ->
      Sp.Single (of_aalp_pr (Multiple_loci.(mlo ff.ff_second)))
  | Pr (Multiple_loci.FirstOrientedSecond fos) ->
      let p1 = of_alp_list (fos.Multiple_loci.first) in
      let p2 = of_aalp_pr (fos.Multiple_loci.second) in
      Sp.Paired (min p1 p2, max p1 p2) 

let read_position =
  let of_alp alp = (List.hd_exn alp).position in
  read_metric of_alp

let read_llhd =
  let of_alp alp = (List.hd_exn alp).llhd in
  read_metric of_alp

let compare_sp_snd rp1 rp2 =
  match rp1, rp2 with
  | Sp.Single p1,      Sp.Single p2
  | Sp.Single p1,      Sp.Paired (_, p2)
  | Sp.Paired (_, p1), Sp.Single p2
  | Sp.Paired (_, p1), Sp.Paired (_, p2) -> compare p1 p2

let aggregate_read_positions ?(readsize=100) =
  List.fold_left ~init:[] ~f:(fun acc (_, read_info) ->
    match read_position read_info with
    | Sp.Single end_         -> (end_ - readsize, end_) :: acc
    | Sp.Paired (end1, end2) -> (end1 - readsize, end1) :: 
                                (end2 - readsize, end2) :: acc)


let group_by_boundary_positions bp_lst rlst =
  let rec loop acc bm bp bp_lst rlst =
    match bp_lst with
    | []                ->
        List.rev ((bm, bp, rlst) :: acc)
    | (bm2, bp2) :: tl  ->
        let before, after =
          List.split_while rlst
            ~f:(function
                | (_, Util.Sp.Single p)
                | (_, Util.Sp.Paired (_, p)) -> p < bp2)
        in
        let nacc = (bm, bp, before) :: acc in
        loop nacc bm2 bp2 tl after
  in
  match bp_lst with
  | []                -> []
  | (bm, bp) :: bp_tl -> loop [] bm bp bp_tl rlst
  

let reads_by_loci po =
  let open Multiple_loci in
  List.fold_left po.Output.per_reads ~init:[]
    ~f:(fun acc {Output.name; d} ->
          match d.most_likely with
          | None  -> invalid_argf "Odd %s has no most likely!" name
          | Some (l, allele) ->
              begin match d.aaps with
              | MPaired mpr_lst ->
                begin match List.Assoc.get l mpr_lst with
                | None  -> invalid_argf "What? %s is missing loci: %s"
                              name (Nomenclature.show_locus l)
                | Some r ->
                    begin match remove_and_assoc l acc with
                    | exception Not_found -> (l, [allele, name, Pr r]) :: acc
                    | (lacc, rest)        -> (l, ((allele, name, Pr r) :: lacc)) :: rest
                    end
                end
              | Single_or_incremental soi_lst ->
                begin match List.Assoc.get l soi_lst with
                | None   -> invalid_argf "What? %s is missing loci: %s"
                              name (Nomenclature.show_locus l)
                | Some r ->
                  begin match remove_and_assoc l acc with
                  | exception Not_found -> (l, [allele, name, Soi r]) :: acc
                  | (lacc, rest)        -> (l, ((allele, name, Soi r) :: lacc)) :: rest
                  end
                end
              end)

