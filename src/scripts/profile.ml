
(*
#require "csv";
*)

open Util

let string_to_locus s =
  match Nomenclature.parse_locus s with
  | Ok l    -> l
  | Error e -> invalid_arg e

let locus_to_string = Nomenclature.show_locus

let keys =
  [ "120021"; "120074"; "120984"; "122879"; "123768"; "125785"; "126119"
  ; "126966"; "127051"; "127079"; "127245"; "128964"; "130445"; "160426"
  ; "160452"; "160917"; "161549"; "161702"; "162569"; "162636"; "162781"
  ; "164054"; "164115"; "164543"; "164964"; "165266"; "200220"; "201858"
  ; "203288"; "203482"; "204378"; "204587"; "205302"; "205855"; "205856"
  ; "205944"; "206039"; "206130"; "206139"; "206704"; "206708"; "206709"
  ; "206961"; "207956"; "207978"; "120013"; "120989"; "125881"; "127064"
  ; "129670"; "160908"; "161712"; "163581"; "164779"; "201063"; "203778"
  ; "205673"; "206038"; "206240"; "206876"
  ]

let _hdr :: data = Csv.load res_fname

let to_correct_map data =
  List.fold_left data ~init:StringMap.empty ~f:(fun smap lst ->
    let k = List.nth_exn lst 4 in
    match Nomenclature.parse (List.nth_exn lst 1) with
    | Error e         -> invalid_arg e
    | Ok (locus, nc) ->
        let l = try StringMap.find k smap with Not_found -> [] in
        StringMap.add k ((locus, nc) :: l) smap)
  |> StringMap.map ~f:(fun l ->
      group_by_assoc l
      |> List.filter_map ~f:(fun (locus, allele_lst) ->
            match allele_lst with
            | a1 :: a2 :: [] -> Some (locus, Nomenclature.Diploid.init a1 a2)
            | _ -> eprintf "Odd! locus: %s has <2 or >2 alleles: %s\n"
                    (locus_to_string locus)
                    (string_of_list ~sep:" " ~f:Nomenclature.to_string allele_lst);
                   None))

let f_of_yojson =
  ParPHMM_drivers.(Output.of_yojson Multiple_loci.final_read_info_of_yojson)

let of_file f =
  Yojson.Safe.stream_from_file f
  |> Stream.next (* only one element per file. *)
  |> f_of_yojson
  |> unwrap_ok

type record =
  { ca : (Nomenclature.locus * Nomenclature.Diploid.tt) list
  ; po : ParPHMM_drivers.Multiple_loci.final_read_info ParPHMM_drivers.Output.t
  }

let verbose_ref = ref true

let load_data ~dir ~fname_to_key =
  let correct = to_correct_map data in
  Sys.readdir dir
  |> Array.to_list
  |> List.fold_left ~init:StringMap.empty ~f:(fun jmap file ->
      if !verbose_ref then printf "At %s\n%!" file;
      if Filename.check_suffix file ".json" then
        let k = fname_to_key file in
        let po = of_file (Filename.concat dir file) in
        try
          let ca = StringMap.find k correct in
          StringMap.add k {ca; po} jmap
        with Not_found ->
          eprintf "Didn't find key %s in correct results! Ignoring" k;
          jmap
      else
        jmap)

let diploid_of_zll { ParPHMM_drivers.Output.allele1; allele2; _} =
  match Nomenclature.parse allele1 with
  | Error e         -> invalid_arg e
  | Ok (_l, nt1)    ->
      match Nomenclature.parse allele2 with
      | Error e     -> invalid_arg e
      | Ok (_, nt2) -> Nomenclature.Diploid.init nt1 nt2

let most_likely_minimum_diploid_distance ref_diploid = function
  | []     -> invalid_arg "Empty prohlatype zygosity list!"
  | h :: t ->
    let dd = Nomenclature.Diploid.distance ~n:2 ref_diploid in
    let dh = dd (diploid_of_zll h) in
    List.fold_left t ~init:(h, dh) ~f:(fun (bz, bd) zll ->
        let nd = dd (diploid_of_zll zll) in
        let c = Nomenclature.Diploid.Distance.compare bd nd in
        if c < 0 then
          (bz, bd)
        else if c = 0 then begin
          (* Have the same distance -> choose the most likely zll! *)
          if bz.ParPHMM_drivers.Output.prob >= zll.ParPHMM_drivers.Output.prob then
            (bz, bd)
          else
            (zll, nd)
        end else
          (zll, nd))

(* How do we classify/categorize the resulting distances when our reference
   (or "validated") data has a resolution of two (ex. A*01:02) *)

type distance_category_vs_two_depth_spec =
  (* But prohlatypes output can have more resolution. *)
  | BothRight of
      { number_of_alternatives    : int
      (* For the designated validated 2 type, was is the sum of the 2 digit
       * resolution probability. *)
      ; matched_2_res_probability : float
      }
  (* One has (0,0,_,_) and the other has (a1, a2, _, _)
   * where not both a1 and a2 = 0 *)
  | OneRight of
      { number_of_alternatives  : int
      ; right_is_first          : bool
      ; right                   : string
      ; wrong                   : string
      ; versus                  : string
      ; wrong_homozygosity      : bool
      }
  (* (a1, a2, _, _), (b1, b2, _, _), (a1 <> 0 || a2 <> 0)
                         && (b1 <> 0 || b2 <> 0) *)
  | Different of
      { number_of_alternatives  : int
      ; most_likely_probability : float
      ; right1                  : string
      ; right2                  : string
      ; versus1                 : string
      ; versus2                 : string
      ; wrong_homozygosity      : bool
      }
  [@@deriving show]

let second_digit_mismatches d =
  let open Nomenclature.Diploid in
  let lf = Distance.to_lowest_form d in
  match lf with
  | (0, 0, _, _), (0, 0, _, _) -> false
  | _                          -> true

let categorize_distance dp mz zplst d =
  let open Nomenclature.Diploid in
  let open ParPHMM_drivers.Output in
  let number_of_alternatives_fst =
    List.map zplst ~f:(fun { allele1; prob; _ } -> allele1)
    |> List.dedup
    |> List.length
  in
  let number_of_alternatives_snd =
    List.map zplst ~f:(fun { allele2; prob; _ } -> allele2)
    |> List.dedup
    |> List.length
  in
  let lf = Distance.to_lowest_form d in
  match lf with
  | (0, 0, _, _), (0, 0, _, _)                              ->
      let number_of_alternatives = List.length zplst in
      let matched_2_res_probability =
        let dz2 = reduce_two (diploid_of_zll mz) in
        List.fold_left zplst ~init:0. ~f:(fun p zll ->
            if reduce_two (diploid_of_zll zll) = dz2 then
              p +. zll.prob
            else
              p)
      in
      BothRight
        { number_of_alternatives
        ; matched_2_res_probability
        }
  | (0, 0, _, _), (_c, _d, _, _) (* (_c <> 0 || _d <> 0) *) ->  (* Since in lowest form *)
      let wrong_homozygosity = dp.lower = dp.upper && mz.allele1 <> mz.allele2 in
      begin match d with
      | Distance.Straight (d1, d2) ->
          if (fst lf) = d1 then         (* Disagree on 2nd distance -> snd allele *)
            OneRight
              { number_of_alternatives = number_of_alternatives_snd
              ; right_is_first = true
              ; right = mz.allele1
              ; wrong = mz.allele2
              ; versus = Nomenclature.to_string dp.upper
              ; wrong_homozygosity
              }
          else (* (fst lf) = d2 *)      (* Disagree on 1st distance -> fst allele *)
            OneRight
              { number_of_alternatives = number_of_alternatives_fst
              ; right_is_first = false
              ; right = mz.allele2
              ; wrong = mz.allele1
              ; versus = Nomenclature.to_string dp.lower
              ; wrong_homozygosity
              }
      | Distance.Crossed (d1, d2) ->
          if (fst lf) = d1 then
            OneRight
              { number_of_alternatives = number_of_alternatives_fst
              ; right_is_first = false
              ; right = mz.allele2
              ; wrong = mz.allele1
              ; versus = Nomenclature.to_string dp.lower
              ; wrong_homozygosity
              }
          else  (* (fst lf) = d2 *)
            OneRight
              { number_of_alternatives = number_of_alternatives_snd
              ; right_is_first = true
              ; right = mz.allele1
              ; wrong = mz.allele2
              ; versus = Nomenclature.to_string dp.upper
              ; wrong_homozygosity
              }
      end
  | _                                                       ->
      let versus1, versus2 =
        match d with
        | Distance.Straight _ -> Nomenclature.to_string dp.lower, Nomenclature.to_string dp.upper
        | Distance.Crossed  _ -> Nomenclature.to_string dp.upper, Nomenclature.to_string dp.lower
      in
      let wrong_homozygosity = dp.lower = dp.upper && mz.allele1 <> mz.allele2 in
      Different
        { number_of_alternatives = List.length zplst
        ; most_likely_probability = mz.prob
        ; right1  = mz.allele1
        ; right2  = mz.allele2
        ; versus1
        ; versus2
        ; wrong_homozygosity
        }

type read_info =
  | P of ParPHMM_drivers.Alleles_and_positions.t ParPHMM_drivers.Multiple_loci.paired
  | S of ParPHMM_drivers.Alleles_and_positions.t ParPHMM_drivers.Multiple_loci.single_or_incremental
  [@@deriving show]

let read_position =
  let open ParPHMM in
  let open ParPHMM_drivers in
  let open ParPHMM_drivers in
  let of_aalp_pr = function
    | Filtered _    -> invalid_argf "read was filtered ?!?"
    | Completed alp -> (List.hd_exn alp).position
  in
  let take_regular r c = Alleles_and_positions.descending_cmp r c <= 0 in
  let mlo fp =
    Orientation.most_likely_between ~take_regular fp
    |> map_completed ~f:snd
  in
  function
  | S (Multiple_loci.SingleRead or_) ->
      true, of_aalp_pr (mlo or_)
  | S (Multiple_loci.PairedDependent pd) ->
      false, of_aalp_pr pd.Multiple_loci.second
  | P (Multiple_loci.FirstFiltered ff) ->
      true, of_aalp_pr (Multiple_loci.(mlo ff.ff_second))
  | P (Multiple_loci.FirstOrientedSecond fos) ->
      false, of_aalp_pr (fos.Multiple_loci.second)

type per_locus_per_sample =
  { claimed_homozygous  : bool
  ; actually_homozygous : bool
  ; category            : distance_category_vs_two_depth_spec
  ; reads1              : (string * string * read_info) list
  ; reads2              : (string * string * read_info) list
  ; reads1_pct          : float
  ; reads2_pct          : float
  }
  [@@deriving show]

let aggregate_read_positions ?(readsize=100) =
  List.fold_left ~init:[] ~f:(fun acc (_, read_info) ->
      let single, end_ = read_position read_info in
      if single then
        (end_ - readsize, end_) :: acc
      else
        (end_ - 2 * readsize, end_) :: acc)


(* Fail if they don't have the same keys *)
let merge_assoc l1 l2 ~f =
  let cmp_by_fst (f1, _) (f2, _) = compare f1 f2 in
  let l1 = List.sort ~cmp:cmp_by_fst l1 in
  let l2 = List.sort ~cmp:cmp_by_fst l2 in
  List.map2 l1 l2 ~f:(fun (l1, k1) (l2, k2) ->
      if l1 <> l2 then invalid_arg "mismatched keys!" else l1, (f k1 k2))

let zll_to_string {ParPHMM_drivers.Output.allele1; allele2; prob; _} =
  sprintf "%s\t%s\t%0.5f" allele1 allele2 prob

let zlls_to_string =
  string_of_list ~sep:"\n\t" ~f:zll_to_string

let by_loci_empty_assoc = Nomenclature.[ A, []; B, []; C, [] ]

let reads_by_loci po =
  let open ParPHMM_drivers.Multiple_loci in
  let init = by_loci_empty_assoc in
  List.fold_left po.ParPHMM_drivers.Output.per_reads ~init
    ~f:(fun acc {ParPHMM_drivers.Output.name; d} ->
          match d.most_likely with
          | None  -> invalid_argf "Odd %s has no most likely!" name
          | Some (l, allele) ->
              begin match d.aaps with
              | MPaired mpr_lst ->
                begin match List.Assoc.get l mpr_lst with
                | None  -> invalid_argf "What? %s is missing loci: %s" name (Nomenclature.show_locus l)
                | Some r ->
                  let lacc, rest = remove_and_assoc l acc in
                  (l, ((allele, name, P r) :: lacc)) :: rest
                end
              | Single_or_incremental soi_lst ->
                begin match List.Assoc.get l soi_lst with
                | None   -> invalid_argf "What? %s is missing loci: %s" name (Nomenclature.show_locus l)
                | Some r ->
                  let lacc, rest = remove_and_assoc l acc in
                  (l, ((allele, name, S r) :: lacc)) :: rest
                end
              end)

let assign_reads mdz reads =
  let open ParPHMM_drivers.Output in
  let z1 = Nomenclature.parse mdz.allele1 |> unwrap_ok |> snd in
  let z2 = Nomenclature.parse mdz.allele2 |> unwrap_ok |> snd in
  List.fold_left reads ~init:([],[]) ~f:(fun (r1, r2) r ->
      let allele, _name, _ri = r in
      let az = Nomenclature.parse allele |> unwrap_ok |> snd in
      let d1 = Nomenclature.distance az z1 in
      let d2 = Nomenclature.distance az z2 in
      let cd = Nomenclature.compare_distance d1 d2 in
      if cd < 0 then
        (r :: r1, r2)
      else if cd > 0 then
        (r1, r :: r2)
     else begin
        printf "allele: %s has same distance to %s and %s"
          allele mdz.allele1 mdz.allele2 ;
        (r1,r2)
      end)

let categorize_wrap ?(output=false) {ca; po} =
  let open ParPHMM_drivers.Output in
  let open Nomenclature.Diploid in
  let zps =
    List.map po.per_loci
      ~f:(fun { locus; zygosity; _} -> locus, zygosity)
  in
  let classical = List.filter ca ~f:(fun (l, _) -> l <> Nomenclature.DRB1) in
  let reads_by_l = reads_by_loci po in
  let total_number_of_reafs =
    List.fold_left reads_by_l ~init:0 ~f:(fun s (_l, lst) -> s + List.length lst)
    |> float_of_int
  in
  let mgd = merge_assoc classical zps ~f:(fun x y -> (x,y)) in
  let mgd2 = merge_assoc mgd reads_by_l ~f:(fun (x,y) z -> (x,y,z)) in
  List.map mgd2 ~f:(fun (l, (dp, zplst, reads)) ->
    let mdz, mdd = most_likely_minimum_diploid_distance dp zplst in
    let fulld = distance dp (diploid_of_zll mdz) in
    let c = categorize_distance dp mdz zplst fulld in
    let reads1, reads2 = assign_reads mdz reads in
    let actually_homozygous =
      if mdz.allele1 = mdz.allele2 then begin
        printf "found a homozygous pair with p %f" mdz.prob;
        mdz.prob > 0.95
      end else
        false
    in
    let plps =
      { claimed_homozygous  = dp.lower = dp.upper
      ; actually_homozygous
      ; category            = c
      ; reads1
      ; reads2
      ; reads1_pct = (float (List.length reads1)) /. total_number_of_reafs
      ; reads2_pct = (float (List.length reads2)) /. total_number_of_reafs
      }
    in
    if output then begin
    printf "%s\n\t%s\n\t%s\n\n\t\t%s\n\t\t%s\n\t\t%s\n\t\t%s\n"
      (Nomenclature.show_locus l)
      (to_string dp)
      (zlls_to_string zplst)
      (zll_to_string mdz)
      (Distance.to_string mdd)
      (Distance.to_string fulld)
      (show_distance_category_vs_two_depth_spec c)
    end;
    l, plps)

let to_categories ?(output=false) r =
  StringMap.mapi r ~f:(fun key data ->
    if output then printf "%s\n" key;
    categorize_wrap ~output data)

let aggregate_by_locus r =
  let init = by_loci_empty_assoc in
  StringMap.fold r ~init ~f:(fun ~key ~data acc ->
    merge_assoc data acc ~f:(fun a l -> a :: l))

let record_map_to_aggregate_lsts ?(output=true) r =
  to_categories ~output r
  |> aggregate_by_locus

type per_locus =
  { br                    : int
  ; or_                   : int
  ; di                    : int
  ; incorrect_homozygous  : int
  (* Sum of the matched 2-digit resolution P *)
  ; br_match_2d           : float
  ; br_r1                 : float
  ; br_r2                 : float
  ; br_ratio              : float
  ; or_right              : float
  ; or_wrong              : float
  ; or_ratio              : float
  ; dr_r1                 : float
  ; dr_r2                 : float
  }

let init_per_locus =
  { br                    = 0
  ; or_                   = 0
  ; di                    = 0
  ; incorrect_homozygous  = 0
  ; br_match_2d           = 0.
  ; br_r1                 = 0.
  ; br_r2                 = 0.
  ; br_ratio              = 0.
  ; or_right              = 0.
  ; or_wrong              = 0.
  ; or_ratio              = 0.
  ; dr_r1                 = 0.
  ; dr_r2                 = 0.
  }

let aggregate_plps l =
  let ifch b = if b then 1 else 0 in
  let fs =
    List.fold_left l ~init:init_per_locus ~f:(fun ipl plps ->
        let l1 = float (List.length plps.reads1) in
        let l2 = float (List.length plps.reads2) in
        match plps.category with
        | BothRight b ->
          { ipl with
              br = ipl.br + 1
              ; incorrect_homozygous =
                  ipl.incorrect_homozygous + ifch (plps.claimed_homozygous && (not plps.actually_homozygous))
              ; br_match_2d = ipl.br_match_2d +. b.matched_2_res_probability
              ; br_r1 = ipl.br_r1 +. l1
              ; br_r2 = ipl.br_r2 +. l2
              ; br_ratio = ipl.br_ratio +. (l1 /. l2)
          }
        | OneRight o  ->
          { ipl with
              or_ = ipl.or_ + 1
              ; incorrect_homozygous =
                  ipl.incorrect_homozygous + ifch (plps.claimed_homozygous && (not plps.actually_homozygous))
              ; or_right = ipl.or_right +. (if o.right_is_first then l1 else l2)
              ; or_wrong = ipl.or_wrong +. (if o.right_is_first then l2 else l1)
              ; or_ratio = ipl.or_ratio +. (if o.right_is_first then (l1/.l2) else (l2/.l1))
          }
        | Different _ ->
          { ipl with
              di = ipl.di + 1
              ; incorrect_homozygous =
                  ipl.incorrect_homozygous + ifch (plps.claimed_homozygous && (not plps.actually_homozygous))
              ; dr_r1 = ipl.dr_r1 +. l1
              ; dr_r2 = ipl.dr_r2 +. l2
          })
  in
  { fs with
      br_match_2d = fs.br_match_2d /. (float fs.br)
      ; br_r1 = fs.br_r1 /. (float fs.br)
      ; br_r2 = fs.br_r2 /. (float fs.br)
      ; br_ratio = fs.br_ratio /. (float fs.br)

      ; or_right = fs.or_right /. (float fs.or_)
      ; or_wrong = fs.or_wrong /. (float fs.or_)
      ; or_ratio = fs.or_ratio /. (float fs.or_)

      ; dr_r1 = fs.dr_r1 /. (float fs.di)
      ; dr_r2 = fs.dr_r2 /. (float fs.di)
  }

let compare_two_plps_maps pm1 pm2 =
  StringMap.merge pm1 pm2 ~f:(fun k p1 p2 ->
      match p1, p2 with
      | None,    None     -> assert false
      | Some _,  None     -> invalid_argf "Different key %s in first" k
      | None  ,  Some _   -> invalid_argf "Different key %s in second" k
      | Some p1, Some p2  ->
          merge_assoc p1 p2 ~f:(fun r1 r2 -> (r1, r2))
          |> List.filter ~f:(fun (_l, (r1, r2)) ->
              match r1.category, r2.category with
              | BothRight _, BothRight _ -> false
              | OneRight _ , OneRight _  -> false
              | Different _, Different _ -> false
              | _          , _           -> true)
          |> function
              | [] -> None
              | l  -> Some l)
