
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

let of_json_file f =
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
        let po = of_json_file (Filename.concat dir file) in
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
  | []     -> invalid_argf "Empty prohlatype zygosity list! ref: %s"
                              (Nomenclature.Diploid.to_string ref_diploid)
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

let most_likely = function
  | []     -> None
  | h :: t ->
      Some (List.fold_left t ~init:h ~f:(fun bz zll ->
        if bz.ParPHMM_drivers.Output.prob >= zll.ParPHMM_drivers.Output.prob then
          bz
        else
          zll))

(* How do we classify/categorize the resulting distances when our reference
   (or "validated") data has a resolution of two (ex. A*01:02) *)

type distance_category_vs_two_depth_spec =
  (* Nothing to compare against, to make it easier to tell what's going on I'm
   * going to add a None option to this type. *)
  | NoSpec
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


let merge_assoc l1 l2 ~f =
  let cmp_by_fst (f1, _) (f2, _) = compare f1 f2 in
  let l1 = List.sort ~cmp:cmp_by_fst l1 in
  let l2 = List.sort ~cmp:cmp_by_fst l2 in
  let cons acc k = function
    | None  -> acc
    | Some v -> (k, v) :: acc
  in
  let rec loop l1 l2 acc =
    match l1, l2 with
    | [],             []              -> List.rev acc
    | (k1, e1) :: t1, []              -> loop t1 [] (cons acc k1 (f k1 (`First e1)))
    | [],             (k2, e2) :: t2  -> loop [] t2 (cons acc k2 (f k2 (`Second e2)))
    | (k1, e1) :: t1, (k2, e2) :: t2  -> if k1 < k2 then
                                           loop t1 l2 (cons acc k1 (f k1 (`First e1)))
                                         else if k1 = k2 then
                                           loop t1 t2 (cons acc k1 (f k1 (`Both (e1, e2))))
                                         else (* k1 > k2 *)
                                           loop l1 t2 (cons acc k2 (f k2 (`Second e2)))
  in
  loop l1 l2 []

(*  List.map2 l1 l2 ~f:(fun (l1, k1) (l2, k2) ->
      if l1 <> l2 then invalid_arg "mismatched keys!" else l1, (f k1 k2)) *)

let zll_to_string {ParPHMM_drivers.Output.allele1; allele2; prob; _} =
  sprintf "%s\t%s\t%0.5f" allele1 allele2 prob

let zlls_to_string =
  string_of_list ~sep:"\n\t" ~f:zll_to_string

let reads_by_loci po =
  let open ParPHMM_drivers.Multiple_loci in
  List.fold_left po.ParPHMM_drivers.Output.per_reads ~init:[]
    ~f:(fun acc {ParPHMM_drivers.Output.name; d} ->
          match d.most_likely with
          | None  -> invalid_argf "Odd %s has no most likely!" name
          | Some (l, allele) ->
              begin match d.aaps with
              | MPaired mpr_lst ->
                begin match List.Assoc.get l mpr_lst with
                | None  -> invalid_argf "What? %s is missing loci: %s" name (locus_to_string l)
                | Some r ->
                    begin match remove_and_assoc l acc with
                    | exception Not_found -> (l, [allele, name, P r]) :: acc
                    | (lacc, rest)        -> (l, ((allele, name, P r) :: lacc)) :: rest
                    end
                end
              | Single_or_incremental soi_lst ->
                begin match List.Assoc.get l soi_lst with
                | None   -> invalid_argf "What? %s is missing loci: %s" name (locus_to_string l)
                | Some r ->
                  begin match remove_and_assoc l acc with
                  | exception Not_found -> (l, [allele, name, S r]) :: acc
                  | (lacc, rest)        -> (l, ((allele, name, S r) :: lacc)) :: rest
                  end
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
      else if mdz.allele1 = mdz.allele2 then
        (r :: r1, r2)
      else begin
        eprintf "allele: %s has same distance to %s and %s, NOT assigning!\n" allele mdz.allele1 mdz.allele2 ;
        (r1, r2)
      end)

let categorize_wrap ?(output=false) {ca; po} =
  let open ParPHMM_drivers.Output in
  let open Nomenclature.Diploid in
  let zps =
    List.map po.per_loci ~f:(fun { locus; zygosity; _} -> locus, zygosity)
  in
  (*let classical = List.filter ca ~f:(fun (l, _) -> l <> Nomenclature.DRB1) in *)
  let reads_by_l = reads_by_loci po in
  let total_number_of_reads =
    List.fold_left reads_by_l ~init:0 ~f:(fun s (_l, lst) -> s + List.length lst)
    |> float_of_int
  in
  let mgd =
    merge_assoc ca zps ~f:(fun locus m -> match m with
      | `First _      -> printf "Dropping typed results for %s locus\n"
                            (locus_to_string locus);
                         None
      | `Second z     -> Some (None, z)
      | `Both (r, z)  -> Some (Some r, z))
  in
  let mgd2 =
    merge_assoc mgd reads_by_l ~f:(fun locus m -> match m with
        | `First (ro, zo)     -> Some (ro, zo, [])      (* Empty list for no reads assigned to locus *)
        | `Second _           -> invalid_argf "We have reads assigned to %s locus but no diploid!\n"
                                    (locus_to_string locus)
        | `Both ((ro, zo), l) -> Some (ro, zo, l))
  in
  List.map mgd2 ~f:(fun (l, (dp_opt, zplst, reads)) ->
    let plps =
      match dp_opt with
      | None    ->  (* We don't have a validated type to compare against. *)
          begin match most_likely zplst with
          | None  ->
            { claimed_homozygous = false
            ; actually_homozygous = false
            ; category = NoSpec
            ; reads1 = reads
            ; reads2 = []
            ; reads1_pct = (float (List.length reads)) /. total_number_of_reads
            ; reads2_pct = 0.
            }
          | Some mdz ->
            let reads1, reads2 = assign_reads mdz reads in
            (* TODO, output? *)
            { claimed_homozygous = false
            ; actually_homozygous = mdz.allele1 = mdz.allele2
            ; category = NoSpec
            ; reads1
            ; reads2
            ; reads1_pct = (float (List.length reads1)) /. total_number_of_reads
            ; reads2_pct = (float (List.length reads2)) /. total_number_of_reads
            }
          end
      | Some dp ->
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
          if output then begin
            printf "%s\n\t%s\n\t%s\n\n\t\t%s\n\t\t%s\n\t\t%s\n\t\t%s\n"
              (Nomenclature.show_locus l)
              (Option.value_map dp_opt ~default:"" ~f:to_string)
              (zlls_to_string zplst)
              (zll_to_string mdz)
              (Distance.to_string mdd)
              (Distance.to_string fulld)
              (show_distance_category_vs_two_depth_spec c)
          end;
          { claimed_homozygous  = dp.lower = dp.upper
          ; actually_homozygous
          ; category            = c
          ; reads1
          ; reads2
          ; reads1_pct = (float (List.length reads1)) /. total_number_of_reads
          ; reads2_pct = (float (List.length reads2)) /. total_number_of_reads
          }
    in
    l, plps)

let to_categories ?(output=false) r =
  StringMap.mapi r ~f:(fun key data ->
    if output then printf "%s\n" key;
    categorize_wrap ~output data)

let aggregate_by_locus r =
  StringMap.fold r ~init:[] ~f:(fun ~key ~data acc ->
    merge_assoc acc data ~f:(fun locus m -> match m with
          | `First f      -> Some (f)
          | `Second s     -> Some ([s])
          | `Both (l, a)  -> Some (a :: l)))

let record_map_to_aggregate_lsts ?(output=true) r =
  to_categories ~output r
  |> aggregate_by_locus

type 'rs per_locus =
  { br                    : int
  ; or_                   : int
  ; di                    : int
  ; incorrect_homozygous  : int
  (* Sum of the matched 2-digit resolution P *)
  ; br_match_2d           : float
  ; ns                    : int
  ; read_stat             : 'rs
  }


let init_per_locus =
  { br                    = 0
  ; or_                   = 0
  ; di                    = 0
  ; incorrect_homozygous  = 0
  ; br_match_2d           = 0.
  ; ns                    = 0
  ; read_stat             = []
  }

let aggregate_plps l =
  let ifch b = if b then 1 else 0 in
  let fs =
    List.fold_left l ~init:init_per_locus ~f:(fun ipl plps ->
        match plps.category with
        | NoSpec      ->
            { ipl with ns = ipl.ns + 1
                     ; read_stat = (plps.reads1_pct +. plps.reads2_pct) :: ipl.read_stat
            }
        | BothRight b ->
          { ipl with
              br = ipl.br + 1
              ; incorrect_homozygous =
                  ipl.incorrect_homozygous + ifch (plps.claimed_homozygous && (not plps.actually_homozygous))
              ; br_match_2d = ipl.br_match_2d +. b.matched_2_res_probability
              ; read_stat = (plps.reads1_pct +. plps.reads2_pct) :: ipl.read_stat
          }
        | OneRight o  ->
          { ipl with
              or_ = ipl.or_ + 1
              ; incorrect_homozygous =
                  ipl.incorrect_homozygous + ifch (plps.claimed_homozygous && (not plps.actually_homozygous))
              ; read_stat = (plps.reads1_pct +. plps.reads2_pct) :: ipl.read_stat
          }
        | Different _ ->
          { ipl with
              di = ipl.di + 1
              ; incorrect_homozygous =
                  ipl.incorrect_homozygous + ifch (plps.claimed_homozygous && (not plps.actually_homozygous))
              ; read_stat = (plps.reads1_pct +. plps.reads2_pct) :: ipl.read_stat
          })
  in
  let t = float (List.length l) in
  { fs with br_match_2d = fs.br_match_2d /. (float fs.br)
          ; read_stat = let s = List.fold_left ~f:(+.) ~init:0.0 fs.read_stat in
                        s /. t
  }

let just_top_level_categories_match r1 r2 =
  match r1.category, r2.category with
 | BothRight _, BothRight _ -> true
 | OneRight _ , OneRight _  -> true
 | Different _, Different _ -> true
 | _          , _           -> false

let compare_two_plps_maps ?(strict=false) pm1 pm2 =
  StringMap.merge pm1 pm2 ~f:(fun k rec1 rec2 ->
      match rec1, rec2 with
      | None,    None     -> assert false
      | Some _,  None     -> invalid_argf "Different key %s in first" k
      | None  ,  Some _   -> invalid_argf "Different key %s in second" k
      | Some p1, Some p2  ->
          merge_assoc p1 p2 ~f:(fun locus m ->
            match m with
            | `First r1      -> Some (`First r1)
            | `Second r2     -> Some (`Second r2)
            | `Both (r1, r2) ->
                  if strict && r1.category <> r2.category then
                    Some (`Both (r1, r2))
                  else if not (just_top_level_categories_match r1 r2) then
                    Some (`Both (r1, r2))
                  else
                    None)
          |> function
              | [] -> None
              | l  -> Some l)

let zygosity_of_plps_lst msg locus lst =
  List.find_map lst ~f:(fun p ->
    if p.ParPHMM_drivers.Output.locus = locus then
      Some p.ParPHMM_drivers.Output.zygosity
    else
      None)
  |> Option.value_exn ~msg

let compare_two_r_maps ?(strict=false) r1 r2 =
  let pm1 = to_categories ~output:false r1 in
  let pm2 = to_categories ~output:false r2 in
  let cm = compare_two_plps_maps ~strict pm1 pm2 in
  StringMap.mapi cm ~f:(fun k lst ->
      let rec1 = StringMap.find k r1 in
      let rec2 = StringMap.find k r2 in
      List.map lst ~f:(fun (locus, _) ->
          let open ParPHMM_drivers.Output in
          let zyg1 = zygosity_of_plps_lst (sprintf "Missing zyg1 for %s" k) locus rec1.po.per_loci in
          let zyg2 = zygosity_of_plps_lst (sprintf "Missing zyg2 for %s" k) locus rec2.po.per_loci in
          locus, (zyg1, zyg2)))


let get_pl key l m =
  let {po; _} = StringMap.find key m in
  List.find_map po.ParPHMM_drivers.Output.per_loci ~f:(fun o ->
      if o.ParPHMM_drivers.Output.locus = l then
        Some o
      else
        None)

let get_z key l m =
  Option.map (get_pl key l m)
    ~f:(fun { ParPHMM_drivers.Output.zygosity; _} -> zygosity)

type times =
  { real  : float
  ; user  : float
  ; sys   : float
  }

let scan_line line =
 Scanf.sscanf line "%s@\t%fm%fs" (fun t m s -> t, m, s)

let load_times file =
  let ic = open_in file in
  try
    let l1 = input_line ic in
    assert (l1 = "");
    let r, rm, rs = scan_line (input_line ic) in
    let u, um, us = scan_line (input_line ic) in
    let s, sm, ss = scan_line (input_line ic) in
    assert (r = "real" && u = "user" && s = "sys");
    close_in ic;
    { real  = rm *. 60. +. rs
    ; user  = um *. 60. +. us
    ; sys   = sm *. 60. +. ss
    }
  with e ->
    close_in ic;
    raise e

let load_times dir =
  Sys.readdir dir
  |> Array.to_list
  |> List.filter ~f:(fun s ->
      match String.split ~on:(`Character '_') s with
      | [_; "1"] -> true
      | _        -> false)
  |> List.fold_left ~init:StringMap.empty ~f:(fun jmap file ->
      let times = load_times (Filename.concat dir file) in
      StringMap.add file times jmap)

let compare_times tm1 tm2 =
  StringMap.merge tm1 tm2 ~f:(fun k t1 t2 ->
      match t1, t2 with
      | None,    None     -> assert false
      | Some _,  None     -> invalid_argf "Different key %s in first" k
      | None  ,  Some _   -> invalid_argf "Different key %s in second" k
      | Some p1, Some p2  -> Some (p1, p2))
