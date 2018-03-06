(* Aggregation of the forward passes. *)

open Util

module Pm = Partition_map

let max_pm =
  Pm.fold_values ~init:ParPHMM.Lp.zero
    ~f:(fun m (l,_p) -> ParPHMM.Lp.max m l)

let higher_value_in_first e1 e2 =
  let r1 = max_pm e1 in
  let r2 = max_pm e2 in
  ParPHMM.Lp.(r2 <= r1)

let descending_lg_llhd_cmp f1 f2 =
  ParPHMM.Lp.compare f2 f1 (* Just reverse the order *)

(* What to do when we know the result of a previous run. *)
module Past_threshold = struct

  open ParPHMM

  type t =
    [ `Don't
    | `Start
    | `Set of Lp.t
    ]

  let to_string = function
    | `Don't  -> "Don't"
    | `Start  -> "Start"
    | `Set v  -> sprintf "Set %s" (Lp.to_string v)

  let init b =
    if b then `Start else `Don't

  let new_ prev_threshold proc =
    let mm = proc.maximum_positions_median_match () in
    Lp.max prev_threshold mm

  let new_value prev_threshold proc = function
    | Pass_result.Filtered _  -> prev_threshold
    | Pass_result.Completed _ -> new_ prev_threshold proc

  let update prev_threshold proc pr =
    match prev_threshold with
    | `Don't  -> `Don't
    | `Start  -> `Set (new_value Lp.zero proc pr)
    | `Set pv -> `Set (new_value pv proc pr)

end  (* Past_threshold *)

(* Consider a single orientation: either regular or reverse complement. *)
let specific_orientation ?prev_threshold proc pass_result_map reverse_complement read read_errors =
  let open ParPHMM in
  Pass_result.map (proc.single ?prev_threshold ~read ~read_errors reverse_complement)
    ~f:(fun () -> pass_result_map proc reverse_complement)

(* Sometimes (such as in the beginning) we do not know whether to prefer a
   regular or a reverse-complement orientation for read. *)
module Orientation = struct

  open ParPHMM

  type 'a t =
    { regular     : 'a Pass_result.t
    ; complement  : 'a Pass_result.t
    }
    [@@deriving show,yojson]

  let map ~f t =
    { regular    = f t.regular
    ; complement = f t.complement
    }

  (* The output is ordered so that the left most stat is the most likely.
    We add characters 'R', 'C', 'F' as little annotations. *)
  let to_string ?(sep='\t') ?take_regular completed_to_string t =
    match t.regular, t.complement with
    | Completed r, Completed c  ->
        begin match take_regular with
        | Some p when not (p r c) ->
            sprintf "C%c%s%cR%c%s"
              sep (completed_to_string c) sep sep (completed_to_string r)
        | Some _ | None           ->
            sprintf "R%c%s%cC%c%s"
              sep (completed_to_string r) sep sep (completed_to_string c)
        end
    | Completed r, Filtered m   ->
        sprintf "R%c%s%cF%c%s" sep (completed_to_string r) sep sep m
    | Filtered m, Completed c   ->
        sprintf "C%c%s%cF%c%s" sep (completed_to_string c) sep sep m
    | Filtered mr, Filtered mc  ->
        sprintf "F%c%s%cF%c%s" sep mr sep sep mc

  let most_likely_between ~take_regular t =
    let open ParPHMM in
    let open Pass_result in
    match t.regular, t.complement with
    | Filtered mr, Filtered mc -> Filtered (sprintf "Both! %s, %s" mr mc)
    | Filtered _r, Completed c -> Completed (true, c)
    | Completed r, Filtered _c -> Completed (false, r)
    | Completed r, Completed c ->
      if take_regular r c then
        Completed (false, r)
      else
        Completed (true, c)

  (* Compute a specific orientation and label the other as Filtered. *)
  let specific ?prev_threshold proc pass_result_map rc read read_errors =
    let open Pass_result in
    if rc then
      { regular     = Filtered "Selected reverse complement"
      ; complement  = specific_orientation ?prev_threshold proc pass_result_map
                        true read read_errors
      }
    else
      { regular     = specific_orientation ?prev_threshold proc pass_result_map
                        false read read_errors
      ; complement  = Filtered "Selected regular"
      }

  (* Consider both orientations. *)
  let check pt pass_result_map proc read read_errors =
    let do_work ?prev_threshold rc =
      specific_orientation ?prev_threshold proc pass_result_map rc read read_errors
    in
    match pt with
    | `Don't ->
        let regular = do_work false in
        let complement = do_work true in
        { regular; complement }, `Don't
    | `Start ->
        let regular = do_work false in
        let prev_threshold = Past_threshold.new_value Lp.zero proc regular in
        let complement = do_work ~prev_threshold true in
        let last_threshold = Past_threshold.new_value prev_threshold proc complement in
        { regular; complement }, `Set last_threshold
    | `Set v ->
        let regular = do_work ~prev_threshold:v false in
        let prev_threshold = Past_threshold.new_value v proc regular in
        let complement = do_work ~prev_threshold true in
        let last_threshold = Past_threshold.new_value prev_threshold proc complement in
        { regular; complement; }, `Set last_threshold

end (* Orientation *)

(* Represent diagnostic result of reads.  *)
module Alleles_and_positions = struct

  open ParPHMM

  type datum =
    { allele    : string
    ; llhd      : Lp.t
    ; position  : int
    }
  and t = datum list
  [@@deriving show,yojson]

  let to_string {allele; llhd; position} =
    sprintf "%s,%d,%s" allele position (Lp.to_string ~precision:5 llhd)

  let string_of_list =
    string_of_list ~sep:";" ~f:to_string

  (* Descending *)
  let descending_cmp l1 l2 =
    match l1, l2 with
    | [],  _            -> -1  (* Should error? *)
    | _ , []            -> 1
    | t1  :: _, t2 :: _ -> descending_lg_llhd_cmp t1.llhd t2.llhd

  let of_mt lookup_allele depth mt =
    let lg = ParPHMM.largest depth in
    Pm.fold_indices_and_values mt ~init:[]
        ~f:(fun acc i (lp, pos) -> lg lp (i, pos) acc)
    |> List.map ~f:(fun (llhd, (i, position)) ->
        { allele = lookup_allele i            (* 2nd usually has alterations. *)
        ; llhd
        ; position
        })

  let of_mt_pr la depth =
    ParPHMM.Pass_result.map ~f:(of_mt la depth)

  let of_mt_or la depth or_ =
    Orientation.map or_ ~f:(of_mt_pr la depth)

end (* Alleles_and_positions *)

module Per_read_result = struct

  type t = Alleles_and_positions.t Orientation.t Sp.t [@@deriving yojson]

  let to_string =
    let take_regular r c = Alleles_and_positions.descending_cmp r c <= 0 in
    function
    | Sp.Single aap ->
      sprintf "\n\t%s"
        (Orientation.to_string ~take_regular
            Alleles_and_positions.string_of_list aap)
    | Sp.Paired (ap1, ap2) ->
      sprintf "\n\t%s\n\t%s"
        (Orientation.to_string ~take_regular
            Alleles_and_positions.string_of_list ap1)
        (Orientation.to_string ~take_regular
            Alleles_and_positions.string_of_list ap2)

end (* Per_read_result *)

module Zygosity_best = struct

  type t =
    | NonZero of float
    | Spec of int
    | NoSpec

  let default_non_zero = 0.000001

end (* Zygosity_best *)

module Zygosity_mixed = struct

  open ParPHMM

  (* When storing Zygosity (per pair of allele) information it is important
     to encode the information in the most compressed form: there is quiet
     a lot of it and since most of the pairs have another that is almost
     identical, a lot of is redundant.

     For the purposes of storing coverage, a position of maximum PHMM emission,
     we'll store a list. But across most zygosity pairs, for a given read, that
     position is almost identical.  Therefore, we'll order these values by that
     position.

     Afterwards we'll track a set (in the Partition map sense) of the
     cross-paired index (see cpair) and the first and second "weight".
     A weight is 2x the number of reads that emit from that position.
     This way if the likelihood is the same for a given allele pair
     (ex. homozygous) we can assign a "half"-weight of 1.
  *)

  type emissions = (Pm.Set.t * int * int) list

  type t =
    (* I believe that it still makes sense to store this value expicitly
       as embedding it into a Pm will slow down the merging at the end.*)
    { ta : Lp.t Triangular.Array.t
    ; mutable et : (int * emissions) list
    }

  let init number_alleles =
    { ta = Triangular.Array.make true number_alleles Lp.one
    ; et = []
    }

  let insert_at pos ~new_ ~merge et =
    let rec loop = function
      | []  -> [pos, new_ ()]
      | (p, c) :: tl ->
          if pos = p then
            (p, merge c) :: tl
          else if pos < p then
            (pos, new_ ()) :: (p, c) :: tl
          else (* pos > p *)
            (p, c) :: loop tl
    in
    loop et

  let insert_into st f s st_lst =
    let rec loop acc remaining = function
      | []                 -> List.rev ((remaining, f, s) :: acc)
      | (ps, fn, sn) :: tl ->
          let inter, to_add, didnotm = Pm.Set.all_intersections remaining ps in
          (* Using length as a faster proxy for set... *)
          let li = Pm.Set.length inter in
          let ld = Pm.Set.length didnotm in
          if li = 0 then begin      (* No intersection *)
            if ld = 0 then          (* Nothing remaining in ps ? *)
              invalid_argf "Empty sets? %s with weights %d %d"
                (Pm.Set.to_string ps) fn sn
            else (* ld > 0 ... didnotm = ps & to_add = remaining *)
              loop ((ps, fn, sn) :: acc) to_add tl
          end else (* li > 0 *) begin
            let nacc =
              if ld = 0 then
                (inter, fn + f, sn + s) :: acc
              else if ld > li then
                (inter, fn + f, sn + s) :: (didnotm, fn, sn) :: acc
              else (* ld <= li *)
                (didnotm, fn, sn) :: (inter, fn + f, sn + s) :: acc
            in
            if Pm.Set.length to_add = 0 then                         (* done! *)
              List.rev_append nacc tl
            else
              loop nacc to_add tl
          end
    in
    loop [] st st_lst

  let sp_to_string = Sp.to_string (sprintf "%d")

  let fa (l1, p1) (l2, p2) =
    if Lp.(l1 < l2) then
      (p2, 0, 2, l2)
    else if Lp.(l2 < l1) then
      (p1, 2, 0, l1)
    else begin (* l2 = l1 *)
      (* It is possible to have l2 = l1 but p1 <> p2: I'm not eactly certain
         why this happens, ie. if one of the alleles traverses a gap? We'll
         return the first rather arbitrarily instead of breaking the
         cross-paired methodology.
      if p1 <> p2 then
        invalid_argf "Same likelihoods %s %s at different positions %s %s!"
          (Lp.to_string l1) (Lp.to_string l2)
          (sp_to_string p1) (sp_to_string p2); *)
      (p1, 1, 1, l1)
    end

  let equal_quadruple (a, b, c, l1) (d, e, f, l2) =
    a = d && b = e && c = f && Lp.close_enough l1 l2

  let update t llhd_and_pos =
    let cp = Pm.cpair ~f:fa equal_quadruple llhd_and_pos in
    let net =
      Pm.fold_set_and_values cp ~init:t.et
        ~f:(fun lst st (sp_pos, f, s, l) ->
              (* Update the Lp Ta, by side-effect. *)
              Pm.Set.iter st ~f:(fun k ->
                  Triangular.Array.update_k t.ta k
                    ~f:(fun lp -> Lp.(lp * l)));
              let i p l = insert_at p l ~new_:(fun () -> [st, f, s])
                ~merge:(insert_into st f s)
              in
              match sp_pos with
              | Sp.Single p1       -> i p1 lst
              | Sp.Paired (p1, p2) -> i p2 (i p1 lst))
    in
    t.et <- net

  (* Create a sorted, grouped, and constrained, list of the paired likelihoods.
   * The total size (4000^2 = 16_000_000) is a little bit unwieldly and at this
   * point we don't care about the tail, only about where there is actual
   * probability mass. But in-order to properly normalize we need the max
   * likelihood of all values (for the log-sum-exp trick) and a normalizing
   * constant.
   *
   * [sorted_insert] keeps an ascending list list of the [desired_size] highest
   * likelihoods (all indices of the inner list have the same likelihood:
   * grouped), such that only these values are afterwards exported.
   *)
  let sorted_insert desired_size l k (current_size, lst) =
    let insert lst =
      let rec loop acc lst =
        match lst with
        | []     -> List.rev acc @ [l, 1, [k]]
        | h :: t ->
            let hl, hn, hv = h in
            if Lp.(hl < l) then
              loop (h :: acc) t
            else if Lp.(close_enough l hl) then
              List.rev acc @ ((hl, hn + 1, k :: hv) :: t)
            else (* l < hl *)
              List.rev acc @ ((l, 1, [k]) :: lst)
      in
      loop [] lst
    in
    match lst with
    | []
    | _ :: [] -> current_size + 1, insert lst
    | h :: t  ->
      let hl, hn, hv = h in
      if Lp.(hl < l) then begin
          if current_size - hn < desired_size then
          current_size + 1, h :: insert t
        else
          current_size + 1 - hn, insert t
      end else if Lp.(close_enough l hl) then
        current_size + 1, ((hl, hn + 1, k :: hv) :: t)
      else (* l < hl *) begin
        if current_size >= desired_size then
          current_size, lst
        else
          current_size + 1, insert lst
      end

  let best size t =
    let desired_size, nz =
      let open Zygosity_best in
      match size with
      | NonZero v -> Triangular.Array.size t.ta, Some v
      | Spec n    -> max 1 n, None
      | NoSpec    -> Triangular.Array.size t.ta, None
    in
    let _cs, res =
      Triangular.Array.foldk_left t.ta ~init:(0, [])
        ~f:(fun cs_acc k l -> sorted_insert desired_size l k cs_acc)
    in
    let most_likely = List.rev_map ~f:(fun (l, _n, ijlist) -> (l, ijlist)) res in
    let maxl, _ = List.hd_exn most_likely in
    (* If we've never added any values then just return empty list. *)
    if maxl = Lp.one then
      []
    else
      let prob = Lp.probability ~maxl in
      (* Compute the normalizing constant. *)
      let ns =
        let f s l = s +. prob l in
        Triangular.Array.fold_left t.ta ~init:0. ~f
      in
      let kinv k =
        Triangular.(Indices.full_upper_inverse t.ta.Array.n k)
      in
      let lookup_emissions k =
        (* A leaky abstraction here. *)
        let i, j = kinv k in
        let e1, e2 =
          List.fold_left t.et ~init:([], [])
            ~f:(fun (a1, a2) (pos, lst) ->
                  let inside_opt =
                    List.find_map lst ~f:(fun (st, f, s) ->
                        if Pm.Set.inside k st then Some (f, s) else None)
                  in
                  match inside_opt with
                  | None        -> (a1, a2)
                  | Some (0, 0) -> invalid_argf "for %d -> %d, %d, double zero count at %d"
                                      k i j pos
                  | Some (n, 0) -> (pos, n) :: a1, a2
                  | Some (0, n) -> a1, (pos, n) :: a2
                  | Some (n, m) -> (pos, n) :: a1, (pos, m) :: a2)
        in
        i, j, List.rev e1, List.rev e2
      in
      match nz with
      | Some lb ->
        List.filter_map most_likely ~f:(fun (l, klst) ->
          let p = prob l /. ns in
          if p > lb then     (* This isn't necessarily the right notion of zero. *)
            let ijlst = List.map klst ~f:lookup_emissions in
            Some (l, p, ijlst)
          else
            None)
      | None ->
          List.map most_likely ~f:(fun (lp, klst) ->
            let ijlst = List.map klst ~f:lookup_emissions in
            (lp, prob lp /. ns, ijlst))

  let number_of_reads lst =
    let s = List.fold_left lst ~init:0 ~f:(fun a (_p, w) -> a + w) in
    (float s) /. 2.

  let single_emission_to_string lst =
    string_of_list lst ~show_empty:false ~sep:";"
      ~f:(fun (p,w) -> sprintf "%d,%f" p (float w /. 2.0))
    |> sprintf "[%s]"

end (* Zygosity_mixed *)

module Likelihoods_and_zygosity = struct

  (* TODO: We're opening Util which has likehood defined. *)
  open ParPHMM

  let add_ll state_llhd llhd =
    Pm.iter_set llhd ~f:(fun i (v, _p) ->
      state_llhd.(i) <- Lp.(state_llhd.(i) * v))

  let add_lz zygosity llhd =
    Zygosity_mixed.update zygosity llhd

  let add_ll_and_lz state_llhd zygosity llhd =
    add_ll state_llhd llhd;
    add_lz zygosity llhd

end (* Likelihoods_and_zygosity *)

(* Outputing results.

   TODO: Expose a separator argument besides tab. *)
module Output = struct

  open ParPHMM

  type header =
    { prohlatype_version  : string
    ; commandline         : string
    }
  and locus = Nomenclature.locus
  and per_allele_log_likelihood =
    { allele : string
    ; a_llhd : Lp.t
    ; alters : MSA.Alteration.t list
    }
  and zygosity_log_likelihood =
    { allele1           : string
    ; allele2           : string
    ; z_llhd            : Lp.t
    ; prob              : float
    ; read1_emissions   : (int * int) list
    ; read2_emissions   : (int * int) list
    }
  and per_locus =
    { imgt_release  : string
    ; locus         : locus
    ; zygosity      : zygosity_log_likelihood list
    ; per_allele    : per_allele_log_likelihood array
    }
  and 'p per_read =
    { name : string
    ; d    : 'p
    }
  and 'p t =
    { header    : header
    ; per_loci  : per_locus list
    ; per_reads : 'p per_read list
    }
  [@@deriving yojson]

  type format_spec =
    [ `TabSeparated
    | `Json
    ]

  type 'd format_ =
    | TabSeparated of ('d -> string)
    | Json of ('d -> Yojson.Safe.json)

  (* Control how much information to print out. If not specified (ie None) then
     print everything. This makes sense for num_likelihoods and num_per_read
     but not for num_zygosities (and maybe not so much for num_likelihoods).
     Set to Some 0 to not print anything. *)
  type depth =
    { num_likelihoods : int option
    ; num_zygosities  : Zygosity_best.t
    ; num_per_read    : int option
    }

  let default_depth =
    { num_likelihoods = None
    ; num_zygosities  = Zygosity_best.(NonZero default_non_zero)
    ; num_per_read    = None
    }

  let cmp_pall pl1 pl2 =
    let lc = descending_lg_llhd_cmp pl1.a_llhd pl2.a_llhd in
    if lc <>  0 then
      lc
    else
      (* To Schwartz transform or not ? *)
      Nomenclature.(compare_by_resolution
        (parse_to_resolution_exn pl1.allele)
        (parse_to_resolution_exn pl2.allele))

  let create_per_locus imgt_release locus allele_arr allele_llhd_arr zbest =
    let per_allele =
      Array.mapi allele_arr ~f:(fun i (allele, alters) ->
        { allele; a_llhd = allele_llhd_arr.(i); alters })
    in
    Array.sort ~cmp:cmp_pall per_allele;
    let zygosity =
      List.map zbest ~f:(fun (l, p, ijlist) ->
        List.map ijlist ~f:(fun (i, j, read1_emissions, read2_emissions) ->
            { allele1 = fst allele_arr.(i)
            ; allele2 = fst allele_arr.(j)
            ; z_llhd  = l
            ; prob    = p
            ; read1_emissions
            ; read2_emissions
            }))
        |> List.concat
    in
    { imgt_release
    ; locus
    ; per_allele
    ; zygosity
    }

  let create header per_loci per_reads =
    { header; per_loci; per_reads }

  let shrink d t =
    let shrink_list l =
      Option.value_map ~default:l ~f:(List.take l)
    in
    let shrink_array a =
      Option.value_map ~default:a ~f:(fun l ->
          Array.sub a ~pos:0 ~len:(min l (Array.length a)))
    in
    let shrink_pl { imgt_release; locus; per_allele; zygosity } =
      { imgt_release
      ; locus
      ; per_allele = shrink_array per_allele d.num_likelihoods
      ; zygosity = zygosity (* Shrunk on construction. *)
      }
    in
    { header    = t.header
    ; per_loci  = List.map t.per_loci ~f:shrink_pl
    ; per_reads = shrink_list t.per_reads d.num_per_read
    }

  let array_iter_on_non_empty f b = function
    | [||] -> ()
    | a    -> ignore (b ()); Array.iter ~f a

  let list_iter_on_non_empty f b = function
    | []   -> ()
    | lst  -> ignore (b ()); List.iter ~f lst

  let tabs oc d_to_string { header; per_loci; per_reads } =
    let open Nomenclature in
    let alters_to_string a =
      string_of_list ~show_empty:false ~sep:";"
        ~f:MSA.Alteration.to_string a
    in
    let fprint_per_allele_log_likelihood { allele; a_llhd; alters } =
      fprintf oc "%-16s\t%s\t%s\n"
        allele (Lp.to_string ~precision:20 a_llhd) (alters_to_string alters)
    in
    let fprint_zygosity_log_likelihood
      { allele1; allele2; z_llhd; prob; read1_emissions; read2_emissions } =
      fprintf oc "%-16s\t%-16s\t%s\t%0.6f\t%0.1f\t%0.1f\t%s\t%s\n"
        allele1 allele2 (Lp.to_string ~precision:10 z_llhd) prob
        (Zygosity_mixed.number_of_reads read1_emissions)
        (Zygosity_mixed.number_of_reads read2_emissions)
        (Zygosity_mixed.single_emission_to_string read1_emissions)
        (Zygosity_mixed.single_emission_to_string read2_emissions)
    in
    fprintf oc "Prohlatype version: %s\n" header.prohlatype_version;
    fprintf oc "Command Line: %s\n" header.commandline;
    List.iter per_loci ~f:(fun { locus; per_allele; zygosity} ->
      list_iter_on_non_empty
        fprint_zygosity_log_likelihood
        (fun () ->
          fprintf oc "Zygosity %s:\n%-16s\t%-16s\tlog likelihood\tProb    \t\
                      # of reads1\t# of reads2\t    E1\t    E2\n"
            (show_locus locus)
            "Allele 1" "Allele 2")
        zygosity;
      array_iter_on_non_empty
        fprint_per_allele_log_likelihood
        (fun () ->
           fprintf oc "Likelihood %s:\nAllele\tlog likelihood\tImputation Description\n"
             (show_locus locus))
        per_allele);
    list_iter_on_non_empty
      (fun { name; d} -> fprintf oc "%s%s\n" name (d_to_string d))
      (fun () -> fprintf oc "Per Read:\n")
      per_reads

  let json oc d_to_yojson t =
    let yj = to_yojson d_to_yojson t in
    Yojson.Safe.to_channel ~std:true oc yj

  let f oc t d = function
    | TabSeparated d_to_string -> tabs oc d_to_string (shrink d t)
    | Json d_to_yojson         -> json oc d_to_yojson (shrink d t)

end (* Output *)

(* Enclose the extraction of read errors from FASTQ reads. *)
module Fastq_items = struct

  exception Read_error_parsing of string

  let repf fmt =
    ksprintf (fun s -> raise (Read_error_parsing s)) fmt

  let single_utimed fqi ~k =
    let open Core_kernel.Std in
    let open Biocaml_unix.Fastq in
    match Fastq.phred_log_probs fqi.qualities with
    | Result.Error e        -> repf "%s" (Error.to_string_hum e)
    | Result.Ok read_errors -> k fqi.name fqi.sequence read_errors

  let single oc fqi ~k =
    let open Biocaml_unix.Fastq in
    time oc (sprintf "updating on single read %s" fqi.name)
      (fun () -> single_utimed fqi ~k)

  let paired_untimed fq1 fq2 ~k =
    let open Core_kernel.Std in
    let open Biocaml_unix.Fastq in
    match Fastq.phred_log_probs fq1.qualities with
    | Result.Error e    -> repf "%s" (Error.to_string_hum e)
    | Result.Ok re1     ->
      match Fastq.phred_log_probs fq2.qualities with
      | Result.Error e  -> repf "%s" (Error.to_string_hum e)
      | Result.Ok re2   -> k fq1.name fq1.sequence re1 fq2.sequence re2

  let paired oc fq1 fq2 ~k =
    let open Biocaml_unix.Fastq in
    assert (fq1.name = fq2.name);
    time oc (sprintf "updating on double read %s" fq1.name)
      (fun () -> paired_untimed fq1 fq2 ~k)

end (* Fastq_items *)

type output_opt =
  { allele_depth  : int
  ; output_format : Output.format_spec
  ; depth         : Output.depth
  }

type single_conf =
  { commandline           : string
  (** How the program was invoked.
   * We'll save it here explicitly because we want to preserve it in the final
   * output for reproducibility.*)

  ; prealigned_transition_model : bool
  (** Use a transition model (transition probability between hidden states of
      the Markov model that represent Start, End, Match, Insert and Delete),
      that assumes the entire read will fit completely inside the reference
      area. *)

  ; allele                : string option
  (** Perform the operation against only the specified allele instead *)

  ; insert_p              : float option
  (* Override the default insert emission probability. *)

  ; band                  : ParPHMM.band_config option
  (* Perform banded passes. *)

  ; max_number_mismatches : int option
  (* Max number of allowable mismatches to use a threshold filter for the
     forward pass. *)

  ; past_threshold_filter : bool
  (* Use the previous match likelihood, when available (ex. against reverse
     complement), as a threshold filter for the forward pass. *)

  ; check_rc              : bool
  (** Compare against the reverse complement of a read and take the best likelihood. *)

  ; split                 : int option
  (** Split the forward pass to operate over this many segments of the read.
      The value must divide the read_length evenly or an exception is thrown. *)

  ; output_opt            : output_opt
  }

let single_conf ?allele ?insert_p ?band ?max_number_mismatches ?split
  ~prealigned_transition_model
  ~past_threshold_filter ~check_rc
  ~output_opt commandline =
    { commandline
    ; prealigned_transition_model
    ; allele
    ; insert_p
    ; band
    ; max_number_mismatches
    ; past_threshold_filter
    ; check_rc
    ; split
    ; output_opt
    }

module type Worker = sig

  type conf
  type input
  type state
  type read_result
  type final_state

  val init      : out_channel
                  -> conf
                  -> read_length:int
                  -> input
                  -> state
  val single    : out_channel
                  -> state
                  -> Biocaml_unix.Fastq.item
                  -> read_result
  val paired    : out_channel
                  -> state
                  -> Biocaml_unix.Fastq.item
                  -> Biocaml_unix.Fastq.item
                  -> read_result
  val merge     : out_channel
                -> state
                -> read_result
                -> unit
  val finalize  : state -> final_state
  val output    : final_state -> out_channel -> unit

end (* Worker *)

(* Conforms to the Worker signature, but I'm not limiting it to it, because
   we reuse some of the logic in Multiple_loci. *)
module Forward (* : Worker *) = struct

  type conf = single_conf
  type input = ParPHMM.t

  open ParPHMM

  type stat = (Lp.t * int) mt

  let proc_to_stat proc _rc = proc.per_allele_llhd_and_pos ()

  let forward oc past_threshold proc read read_errors =
    time oc (sprintf "forward of %s, past_threshold: %s; read: %s"
      proc.ParPHMM.name (Past_threshold.to_string past_threshold)
        (String.sub_exn read ~index:0 ~length:10))
    (fun () ->
      Orientation.check past_threshold proc_to_stat
        proc read read_errors)


  type opt = output_opt

  type state =
    { commandline       : string
    ; imgt_release      : string
    ; locus             : Nomenclature.locus
    ; per_allele_lhood  : Lp.t array
    ; zygosity          : Zygosity_mixed.t
    ; mutable per_reads : Per_read_result.t Output.per_read list
    ; initial_pt        : Past_threshold.t
    ; proc              : ParPHMM.proc
    ; allele_arr        : (string * MSA.Alteration.t list) array
    ; opt               : output_opt
    }

  type final_state = state

  type read_result = stat Orientation.t Sp.t Output.per_read

  let to_proc_and_allele_arr
    ?allele ?insert_p ?max_number_mismatches ?band ?split
    ~prealigned_transition_model read_length parPHMM_t =
    match allele with
    | None ->
      let proc =
        match split with
        | None  ->
          setup_single_pass ?insert_p ?max_number_mismatches ?band
            ~prealigned_transition_model read_length parPHMM_t
        | Some n ->
          setup_splitting_pass ?insert_p ?max_number_mismatches ?band
            ~prealigned_transition_model read_length n parPHMM_t
      in
      proc, parPHMM_t.alleles
    | Some allele ->
      let proc =
        setup_single_allele_forward_pass ?insert_p ?max_number_mismatches
          ~prealigned_transition_model read_length allele parPHMM_t
      in
      let arr =
        let index_opt = Array.findi parPHMM_t.alleles ~f:(fun (a,_) -> a = allele) in
        match index_opt with
        | None  -> invalid_argf "Allele %s isn't part of gene!" allele
        | Some i -> [| parPHMM_t.alleles.(i) |]
      in
      proc, arr

  let init _oc conf ~read_length parPHMM_t =
    let { prealigned_transition_model; allele; insert_p; band
        ; max_number_mismatches; past_threshold_filter ; split; _} = conf
    in
    let open ParPHMM in
    let proc, allele_arr = to_proc_and_allele_arr
        ?allele ?insert_p ?max_number_mismatches ?band ?split
        ~prealigned_transition_model read_length parPHMM_t
    in
    let initial_pt = Past_threshold.init past_threshold_filter in
    (* Curry away arguments and discard final threshold; can't use it. *)
    { commandline      = conf.commandline
    ; imgt_release     = parPHMM_t.release
    ; locus            = parPHMM_t.locus
    ; per_allele_lhood = proc.init_global_state ()
    ; zygosity         = Zygosity_mixed.init (Array.length allele_arr)
    ; per_reads        = []
    ; initial_pt
    ; proc
    ; allele_arr
    ; opt              = conf.output_opt
    }

  let single oc { initial_pt; proc; opt; _} fqi =
    let forward pt r re = fst (forward oc pt proc r re) in
    Fastq_items.single oc fqi ~k:(fun name read read_errors ->
      { Output.name; d = Sp.Single (forward initial_pt read read_errors) })

  let take_regular = higher_value_in_first

  let paired oc { initial_pt; proc; opt } fq1 fq2 =
    let forward pt r re = fst (forward oc pt proc r re) in
    Fastq_items.paired oc fq1 fq2 ~k:(fun name rd1 re1 rd2 re2 ->
      let r1 = forward initial_pt rd1 re1 in
      match Orientation.most_likely_between r1 ~take_regular with
      | Filtered m  ->
          let r2 = forward initial_pt rd2 re2 in
          { Output.name; d = Sp.Paired (r1, r2)}
      | Completed (c, _)  ->
          let r2 = Orientation.specific proc proc_to_stat
                    (not c) rd2 re2 in
          { Output.name; d = Sp.Paired (r1, r2)})

  let cons_per_reads state name stat_or_sp =
    let la i = fst state.allele_arr.(i) in
    let f = Alleles_and_positions.of_mt_or la state.opt.allele_depth in
    let d = Sp.map f stat_or_sp in
    state.per_reads <- { Output.name; d } :: state.per_reads

  let merge_paired_pms s1 s2 =
    let m (l1,p1) (l2,p2) =
      let p =
        if p1 <= p2 then
          Sp.Paired (p1, p2)
        else
          Sp.Paired (p2, p1)
      in
      Lp.(l1 * l2), p
    in
    Pm.merge ~eq:Lpr.equal s1 s2 m

  let to_single_pms s =
    Pm.map s (fun _ _ -> false) ~f:(fun (l,p) -> l, Sp.Single p)

  let merge oc state rr =
    let pr s = Orientation.most_likely_between s ~take_regular in
    cons_per_reads state rr.Output.name rr.Output.d;
    match rr.Output.d with
    | Sp.Single stat_or  ->
        begin match pr stat_or with
        | Filtered m        ->
            fprintf oc "Both orientations filtered %s: %s\n" rr.Output.name m
        | Completed (_c, s) ->
            Likelihoods_and_zygosity.add_ll_and_lz
              state.per_allele_lhood state.zygosity (to_single_pms s)
        end
    | Sp.Paired (s1, s2) ->
        begin match pr s1 with
        | Filtered m        ->
            fprintf oc "Both orientations of first read filtered %s: %s\n" rr.Output.name m
        | Completed (_, sf1) ->
            begin match pr s2 with
            | Filtered m    ->
                fprintf oc "Both orientations of second read filtered %s: %s\n" rr.Output.name m
            | Completed (_, sf2) ->
                let cml = merge_paired_pms sf1 sf2 in
                Likelihoods_and_zygosity.add_ll_and_lz
                  state.per_allele_lhood state.zygosity cml
            end
        end

  let finalize s = s

  let output state oc =
    let ot =
      let zbest   = Zygosity_mixed.best state.opt.depth.Output.num_zygosities state.zygosity in
      let pl      = Output.create_per_locus state.imgt_release state.locus state.allele_arr
                      state.per_allele_lhood zbest in
      let header  =
        { Output.prohlatype_version = Version.version
        ; Output.commandline        = state.commandline } in
      Output.create header [pl] state.per_reads
    in
    let of_ =
      match state.opt.output_format with
      | `TabSeparated ->
          Output.TabSeparated Per_read_result.to_string
      | `Json ->
          Output.Json Per_read_result.to_yojson
    in
    Output.f oc ot state.opt.depth of_

end (* Forward *)

(* A map of reads to a Viterbi decoded path *)
module Viterbi :
    Worker with type input = ParPHMM.t
            and type conf = single_conf
  = struct

  type conf = single_conf
  type input = ParPHMM.t

  open ParPHMM
  open Path

  let center_justify str size =
    let n = String.length str in
    let rem = size - n in
    if rem <= 0 then
      str
    else
      let l = rem / 2 in
      let r = if (rem mod 2) = 1 then l + 1 else l in
      sprintf "%*s%s%*s" l "" str r ""

  let viterbi_result_to_string ?(width=250) ?(labels=("reference ", "read "))
    {reverse_complement; emission; path_list} =
    let separate ?(from="") rc { reference; read; start; end_ } =
      sprintf "Reverse complement: %b, final emission: %s, path length: %d from %d to %d %s\n%s"
        rc (Lp.to_string emission) (List.length path_list) start end_ from
        (manual_comp_display ~width ~labels reference read)
    in
    match to_strings path_list with
    | Sp.Single s        -> separate reverse_complement s
    | Sp.Paired (r1, r2) ->
        let s1, rc1, f1, s2, rc2, f2 =
          if r1.start < r2.start then
            r1, reverse_complement, "read 1", r2, not reverse_complement, "read 2"
          else
            r2, not reverse_complement, "read2", r1, reverse_complement, "read 1"
        in
        let inner_distance = s2.start - s1.end_ in
        if inner_distance < 0 then
          sprintf "Negative Inner distance:\n%s\n%s"
            (separate ~from:f1 rc1 s1) (separate ~from:f2 rc2 s2)
        else
          let idmsg = sprintf " inner: %d " inner_distance in
          let allele_i = center_justify idmsg inner_distance in
          let actual_size = String.length allele_i in
          let read_i = String.make actual_size '=' in
          let reference = s1.reference ^ allele_i ^ s2.reference in
          let read = s1.read ^ read_i ^ s2.read in
          let msm_offset = - actual_size in
          sprintf "Reverse complement: %b, final emission: %s, path length: %d from %d to %d %s->%s\n%s"
            rc1 (Lp.to_string emission) (List.length path_list) s1.start s2.end_ f1 f2
            (manual_comp_display ~msm_offset ~width ~labels reference read)

  type read_result =
    { name  : string
    ; res   : viterbi_result
    }

  type state =
    { mutable results : read_result list
    ; vfuncs          : viterbi_proc
    ; labels          : string * string
    }

  type final_state = state

  type opt = unit

  let init _oc conf ~read_length (parPHMM_t : ParPHMM.t) =
    let { prealigned_transition_model; allele; insert_p; _ } = conf in
    let open ParPHMM in
    (* Relying on reference being first in the specified ParPHMM.t allele list. *)
    let allele = Option.value ~default:(fst parPHMM_t.alleles.(0)) allele in
    { results = []
    ; vfuncs  = setup_single_allele_viterbi_pass ?insert_p ~allele
                  ~prealigned_transition_model read_length parPHMM_t
    ; labels = allele ^ " ", "read "
    }

  let single oc { vfuncs; _} =
    Fastq_items.single oc ~k:(fun name read read_errors ->
      { name; res = vfuncs.single ~read ~read_errors})

  let paired oc { vfuncs; _} fq1 fq2 =
    Fastq_items.paired oc fq1 fq2 ~k:(fun name read1 read_errors1 read2 read_errors2 ->
      { name; res = vfuncs.paired ~read1 ~read_errors1 ~read2 ~read_errors2})

  let merge _oc state pr =
    state.results <- pr :: state.results

  let finalize s = s

  let output { results; labels; _} oc =
      List.iter results ~f:(fun {name; res } ->
        fprintf oc "%s:\n%s\n" name (viterbi_result_to_string ~labels res))

end (* Viterbi *)


(* Combine results from multiple loci. *)

type multiple_conf =
  { commandline                 : string
  (* How the program was invoked. *)

  ; prealigned_transition_model : bool
  (** Use a transition model (transition probability between hidden states of
      the Markov model that represent Start, End, Match, Insert and Delete),
      that assumes the entire read will fit completely inside the reference
      area. *)

  ; insert_p                : float option
  (* Override the default insert emission probability. *)

  ; band                    : ParPHMM.band_config option
  (* Perform banded passes. *)

  ; max_number_mismatches   : int option
  (* Max number of allowable mismatches to use a threshold filter for the
     forward pass. *)

  ; past_threshold_filter   : bool
  (* Use the previous match likelihood, when available (ex. against reverse
     complement), as a threshold filter for the forward pass. *)

  ; split                   : int option
  (** Split the forward pass to operate over this many segments of the read.
      The value must divide the read_length evenly or an exception is thrown. *)

  ; output_opt              : output_opt

  (* Not robust see comment below
  Incremental logic has been disabled ... it still introduced errors over a
  small but discoverable set of reads. Maybe there is a way to revive it and
  make it robust, because it does make things faster but for now it is
  commented out because it violates the LogProbability.t abstraction

  ; incremental_pairs       : bool
  (* This option only applies to paired typing. Instead of naively modeling
     the two read pairs at the same time, use the first as guidance of which
     loci is the best, and then apply the second read to only the best loci.
     The default should is false as there seem to always be reads where the
     first is too ambiguous and the error profile shifts such that
     [first_read_best_log_gap] is not good enough. We need a better algorithm
     for this. *)

  One idea to make things better:
  ; first_read_best_log_gap : float option
  (** When performing paired read inference, for the sake of expediency we want
      to make a decision about the best (most likely aligned to read) gene
      based upon the first read; then we know the orientation of second read and
      where to measure. Unfortunately, sometimes, the first read aligns almost
      equally well to different alleles within 2 genes. This value controls a
      band within which we'll keep the best genes (based on the likelihood of the
      first read) and afterwards align the second. The likelihoods, for other
      genes, outside of this band are discarded. By default this value is set
      to |log_10 1/100|. *)
*)
  }

let multiple_conf ?insert_p ?band ?max_number_mismatches ?split
  ~prealigned_transition_model
  ~past_threshold_filter
  ~output_opt
  (*~incremental_pairs *)
  (*?first_read_best_log_gap *)
  commandline =
    { commandline
    ; prealigned_transition_model
    ; insert_p
    ; band
    ; max_number_mismatches
    ; past_threshold_filter
    ; split
    ; output_opt
    (*; incremental_pairs
      ; first_read_best_log_gap *)
    }

(* Conforms to this signature but I'm not imposing it to expose all of the
   yojson input/output logic. *)
module Multiple_loci (* :
  Worker with type input = ParPHMM.t list
          and type conf = multiple_conf
  *) = struct

  type conf = multiple_conf
  type input = ParPHMM.t list

  open ParPHMM

  type 'single_result paired_dependent =
      { first_o : bool
      ; first   : 'single_result
      ; second  : 'single_result Pass_result.t
      }
    (* Orientation and gene of 2nd read is determined by 1st. There is no
       {Pass_result.t} since we're not going to apply a filter to the 2nd read.
       Note that one could "naturally" combine the result of the first by
       passing the first's emission (via per_allele_hood) into the second
       read's forward pass (this is similar to thinking of the 2 reads as just
       one long read), but that will dramastically slow down the pass (it will
       start out already branched). It is better to keep the two likelihoods
       separate and then just multiply (log -> add) the results together. This
       is the desired outcome for paired reads as it is the most
       descriminative. *)
  and 'single_result single_or_incremental =
    | SingleRead of 'single_result Orientation.t
    (* We'll test both regular and complement options. *)
    | PairedDependent of 'single_result paired_dependent
    [@@deriving show,yojson]

  let orientation_char_of_bool reverse_complement =
    if reverse_complement then 'C' else 'R'

  let pass_result_to_short_string sep sr_to_string = function
    | Pass_result.Filtered m -> sprintf "F%c%s" sep m
    | Pass_result.Completed r -> sr_to_string r

  let soi_to_string ?take_regular ?(sep='\t') sr_to_string soi =
    let ots = Orientation.to_string ?take_regular ~sep sr_to_string in
    match soi with
    | SingleRead ors                             ->
        sprintf "SingleRead%c%s" sep (ots ors)
    | PairedDependent { first_o; first; second } ->
        sprintf "PairedDependent%c%c%c%s%c%s"
          sep (orientation_char_of_bool first_o)
          sep (sr_to_string first)
          sep (pass_result_to_short_string sep sr_to_string second)

  let take_regular = Forward.take_regular

  let most_likely_orientation or_ =
    Orientation.most_likely_between ~take_regular or_

  (* true if s1 is better than s2 *)
  let better_soi s1 s2 =
    let l s = max_pm s in
    match s1, s2 with
    | PairedDependent p1, PairedDependent p2  ->
        begin
          match p1.second, p2.second with
          | Completed ss1, Completed ss2 -> let s2 = Lp.(l p2.first * l ss2) in
                                            let s1 = Lp.(l p1.first * l ss1) in
                                            Lp.(s2 < s1)
          | Filtered _,    _             -> false
          | _,             Filtered _    -> true
        end
    | PairedDependent p1, SingleRead sr2      ->
        begin
          match most_likely_orientation sr2 with
          | Completed (_, s)  ->
              begin match p1.second with
              | Completed _               -> true
              | Filtered  _               -> Lp.(l s <= l p1.first)
              end
          | Filtered  _                   -> true
        end
    | SingleRead sr1,     PairedDependent p2  ->
        begin
          match most_likely_orientation sr1 with
          | Completed (_, s)  ->
              begin match p2.second with
              | Completed _               -> false
              | Filtered  _               -> Lp.(l p2.first <= l s)
              end
          | Filtered _                    -> false
        end
    | SingleRead sr1,     SingleRead sr2      ->
        begin
          match most_likely_orientation sr1 with
          | Completed (_, l1)      ->
            begin match most_likely_orientation sr2 with
            | Completed (_, l2)           -> Lp.(l l2 <= l l1)
            | Filtered  _                 -> true
            end
          | Filtered _                    -> false
        end

  type 'single_result first_filtered =
    { ff_first  : string
    ; ff_second : 'single_result Orientation.t
    }
    (* Both orientations of the first read were filtered, so the 2nd read is
       entirely independent of the first. This can occur if we have a filter
       (such as past threshold from another loci) that ignores that loci. *)
  and 'single_result first_oriented_second =
    { first_o : bool
    ; first   : 'single_result
    ; second  : 'single_result Pass_result.t
    }
    (* One of the orientations of the first wasn't filtered and we used
       it to orient the second read. There is a {Pass_result.t} for second
       because in the cases where we have a split, those may trigger a
       filter that stops computation. *)
  and 'single_result paired =
    | FirstFiltered of 'single_result first_filtered
    | FirstOrientedSecond of 'single_result first_oriented_second
    [@@deriving show,yojson]

  let pr_to_string ?(sep='\t') sr_to_string pr =
    let ots = Orientation.to_string ~sep sr_to_string in
    match pr with
    | FirstFiltered { ff_first; ff_second}                 ->
        sprintf "FirstFiltered%c%s%c%s" sep ff_first sep (ots ff_second)
    | FirstOrientedSecond { first_o; first; second } ->
        sprintf "FirstOrientedSecond%c%c%c%s%c%s"
          sep (orientation_char_of_bool first_o)
          sep (sr_to_string first)
          sep (pass_result_to_short_string sep sr_to_string second)

  (* Can be GADT, pp = potentially paired *)
  type 'sr pp =
    | Single_or_incremental of
        (Nomenclature.locus * 'sr single_or_incremental) list
    | MPaired of                  (* M for multi to not conflict with Paired *)
        ((Nomenclature.locus * 'sr paired) list)
    [@@deriving yojson]

  let aap_soi_to_string =
    let take_regular r c = Alleles_and_positions.descending_cmp r c <= 0 in
    soi_to_string ~take_regular

  let pp_to_string ?(sep='\t') sr_to_string = function
    | Single_or_incremental soi_lst ->
        string_of_list soi_lst ~sep:"\n" ~f:(fun (l, s) ->
          sprintf "%s%c%s" (Nomenclature.show_locus l) sep
            (aap_soi_to_string ~sep sr_to_string s))
    | MPaired pr_lst ->
        string_of_list pr_lst ~sep:"\n" ~f:(fun (l, p) ->
          sprintf "%s%c%s" (Nomenclature.show_locus l) sep
            (pr_to_string ~sep sr_to_string p))

  type read_result = Forward.stat pp Output.per_read

  (* The information about each read that we store during the processing of all
   * the read data.  *)
  type initial_read_info =
    { locus         : Nomenclature.locus
    ; llhd_and_pos  : (Lp.t * int Sp.t) mt
    ; aaps          : Alleles_and_positions.t pp
    }

  (* What we'll report. Once we know the final, most likely (P > ~0.0) set of
   * alleles. We'll use the orignal likelihoods to guess which chromosome,
   * specifically which of two alleles, this read comes from.
   *)
  type final_read_info =
    { most_likely   : (Nomenclature.locus * string) option
    (* most likely locus and allele. *)
    ; aaps          : Alleles_and_positions.t pp
    }
    [@@deriving yojson]

  let final_read_info_to_string { most_likely; aaps } =
    sprintf "%s\n%s"
      (Option.value_map most_likely
        ~default:"Unassigned! Is the non-zero threshold too high?"
        ~f:(fun (l, a) -> sprintf "%s %s" (Nomenclature.show_locus l) a))
      (pp_to_string Alleles_and_positions.string_of_list aaps)

  (* Easily convert to Output.per_locus *)
  type per_locus =
    { imgt_release  : string
    ; pl_locus      : Nomenclature.locus
    ; allele_arr    : (string * MSA.Alteration.t list) array
    ; likelihood    : Lp.t array
    ; zygosity      : Zygosity_mixed.t
    }

  type state =
    { commandline       : string
    ; per_loci          : per_locus list
    ; mutable per_reads : initial_read_info Output.per_read list
    ; single_read       : bytes -> float array -> Forward.stat pp
    ; ordinary_paired   : string -> bytes -> float array
                                 -> bytes -> float array
                                 -> Forward.stat pp
    ; opt               : Forward.opt
    }

  type final_state =
    { f_commandline : string
    ; f_opt         : Forward.opt
    ; f_pl_lst      : Output.per_locus list
    ; f_per_reads   : final_read_info Output.per_read list
    }

  (* A PairedDependent with completed second should be chosen over any SingleRead's.
     A PairedDependent with a Filtered second should be compared by first.
     Compare 2 PairedDependent by both.
   *)
  let reduce_soi soi_lst =
    let bsoi_opt =
      List.reduce soi_lst ~f:(fun (l1, s1) (l2, s2) ->
        if better_soi s1 s2 then
          l1, s1
        else
          l2, s2)
    in
    match bsoi_opt with
    | None  -> invalid_arg "Empty soi list!"
    | Some (loci, soi) ->
      let sp =
        match soi with
        | PairedDependent { first; second; _} ->
            begin match second with
            | Filtered _  -> Sp.Single first
            | Completed s -> Sp.Paired (first, s)
            end
        | SingleRead  or_ ->
            begin match most_likely_orientation or_ with
            | Filtered _      -> invalid_argf "Final single read filtered? %s"
                                  (Nomenclature.show_locus loci)
            | Completed (_,s) -> Sp.Single s
            end
      in
      loci, sp

  let allele_arr_from_pl per_loci loci =
    List.find_map per_loci (function
      | {pl_locus; allele_arr} when pl_locus = loci -> Some allele_arr
      | _                                           -> None)
    |> Option.value_exn
          ~msg:(sprintf "Missig %s loci!" (Nomenclature.show_locus loci))

  let map_soi_to_aap per_loci depth lst =
    List.map lst ~f:(fun (loci, soi) ->
      let alleles = allele_arr_from_pl per_loci loci in
      let la i = fst alleles.(i) in
      let of_mt = Alleles_and_positions.of_mt la depth in
      let of_mt_pr = Alleles_and_positions.of_mt_pr la depth in
      let of_mt_or = Alleles_and_positions.of_mt_or la depth in
      let n =
        match soi with
        | PairedDependent { first_o; first; second }  ->
            PairedDependent
              { first_o; first = of_mt first; second = of_mt_pr second }
        | SingleRead oi                               ->
           SingleRead (of_mt_or oi)
      in
      (loci, n))

  (* true if m1 is better than m2, ie higher likelihood. *)
  let better_mpr m1 m2 =
    let l s = max_pm s in
    let ml = most_likely_orientation in
    match m1, m2 with
    | FirstFiltered ff1, FirstFiltered ff2 ->
        begin match ml ff1.ff_second, ml ff2.ff_second with
        | Filtered _,        _                 -> false
        | _,                 Filtered _        -> true
        | Completed (_, s1), Completed (_, s2) -> Lp.(l s2 < l s1)
        end
    | FirstFiltered _,        FirstOrientedSecond _  -> false
    | FirstOrientedSecond _,  FirstFiltered _        -> true
    | FirstOrientedSecond f1, FirstOrientedSecond f2 ->
        begin
          match f1.second, f2.second with
          | Completed ss1, Completed ss2 -> let s1 = Lp.(l f1.first * l ss1) in
                                            let s2 = Lp.(l f2.first * l ss2) in
                                            Lp.(s2 < s1)
          | Filtered _,    _             -> false
          | _,             Filtered _    -> true
        end

  (* We can only compare based upon the result of the 2nd read, but if a
     FirstOrientedSecond scenario is best, we want to use the likelihood from
     BOTH reads in the inference. *)
  let reduce_mpr pr_lst =
    let bmpr_opt =
      List.reduce pr_lst ~f:(fun (l1, m1) (l2, m2) ->
        if better_mpr m1 m2 then
          l1, m1
        else
          l2, m2)
    in
    match bmpr_opt with
    | None            -> invalid_arg "Emmpty pr list!"
    | Some (loci, pr) ->
        let sp =
          match pr with
          | FirstFiltered { ff_second; _} ->
              begin match most_likely_orientation ff_second with
              | Filtered _       -> invalid_argf "Best paired filtered: %s"
                                      (Nomenclature.show_locus loci)
              | Completed (_, s) -> Sp.Single s
              end
          | FirstOrientedSecond {first; second; _} ->
              begin match second with
              | Filtered _  -> Sp.Single first
              | Completed s -> Sp.Paired (first, s)
              end
        in
        loci, sp

  let map_pr_to_aap per_loci depth lst =
    List.map lst ~f:(fun (loci, pr) ->
      let alleles = allele_arr_from_pl per_loci loci in
      let la i = fst alleles.(i) in
      let of_mt = Alleles_and_positions.of_mt la depth in
      let of_mt_pr = Alleles_and_positions.of_mt_pr la depth in
      let of_mt_or = Alleles_and_positions.of_mt_or la depth in
      let pr =
        match pr with
        | FirstFiltered { ff_first; ff_second } ->
            FirstFiltered { ff_first; ff_second = of_mt_or ff_second }
        | FirstOrientedSecond { first_o; first; second } ->
            FirstOrientedSecond
              { first_o; first = of_mt first; second = of_mt_pr second}
      in
      loci, pr)

  let init oc conf ~read_length parPHMM_t_lst =
    let { prealigned_transition_model; insert_p; band; max_number_mismatches
        ; split; output_opt; _ } = conf in
    let paa =
      List.map parPHMM_t_lst ~f:(fun parPHMM_t ->
        let proc, allele_arr = Forward.to_proc_and_allele_arr
            ?insert_p ?max_number_mismatches ?band ?split
            ~prealigned_transition_model read_length parPHMM_t
        in
        let n = Array.length allele_arr in
        let iml = parPHMM_t.ParPHMM.release,parPHMM_t.ParPHMM.locus in
        iml, proc, allele_arr, n)
    in
    let initial_pt = Past_threshold.init conf.past_threshold_filter in
    let forward = Forward.forward oc in
    let specific_orientation p rc r re =
      specific_orientation p Forward.proc_to_stat rc r re
    in
    let single_read read read_errors =
      List.fold_left paa ~init:(initial_pt, [])
            ~f:(fun (pt, acc) ((_, locus), p, _, _) ->
                  let r, npt = forward pt p read read_errors in
                  npt, (locus, SingleRead r) :: acc)
      |> snd                                 (* Discard final threshold *)
      |> List.rev                                (* Restore loci order. *)
      |> fun l -> Single_or_incremental l
    in
    let ordinary_paired name rd1 re1 rd2 re2 =
      let _fpt1, _fpt2, res =
        List.fold_left paa ~init:(initial_pt, initial_pt, [])
          ~f:(fun (pt1, pt2, acc) ((_, locus), p, _, _) ->
                let result1, npt1 = forward pt1 p rd1 re1 in
                match most_likely_orientation result1 with
                | Filtered  m ->
                    let ff_second, npt2 = forward pt2 p rd2 re2 in
                    npt1, npt2, (locus, FirstFiltered { ff_first = m; ff_second }) :: acc
                | Completed (first_o, first) ->
                    let second = specific_orientation p (not first_o) rd2 re2 in
                    let npt2 = Past_threshold.update pt2 p second in
                    npt1, npt2, (locus, FirstOrientedSecond { first_o; first; second }) :: acc)
      in
      MPaired (List.rev res)
    in
    (*let decreasing_likelihood = insert_sorted ( < ) in
      let incremental_paired name rd1 re1 rd2 re2 =
      let best, rest, _fpt1 =
        List.fold_left paa ~init:([], [], initial_pt)
          ~f:(fun (best, acc, pt) ((_, locus), p, _, _) ->
                let result1, npt = forward pt p rd1 re1 in
                match most_likely_orientation result1 with
                | Filtered m ->
                    let nacc = (locus, SingleRead result1) :: acc in
                    best, nacc, npt
                | Completed (c, first) ->
                    let l = max_pm first.Forward.likelihood in
                    let for_second = locus, p, result1, c, first in
                    match best with
                    | []           ->
                            [l, for_second], acc, npt
                    | (bl, (_,_,_,_,_)) :: _ ->
                        begin match classify_likelihood ~best:bl l with
                        | LessLikely ->
                            let nacc = (locus, SingleRead result1) :: acc in
                            best, nacc, npt
                        | MoreLikely ->
                            let nacc =
                              List.map best ~f:(fun (_, (locus, _, r1, _, _)) -> (locus, SingleRead r1))
                              @ acc
                            in
                            [l, for_second], nacc, npt
                        | Within     ->
                            let nbest = decreasing_likelihood l for_second best in
                            nbest, acc, npt
                        end)
      in
      match best with
      | []  ->
          eprintf "All loci were filtered for first read of %s, ignoring second.\n"
            name;
          Single_or_incremental (List.rev rest)
      | lst ->
          let seconds = List.map lst ~f:(fun (_l, (locus, p, _, first_o, first)) ->
            let second = specific_orientation p (not first_o) rd2 re2 in
            (locus, PairedDependent {first_o; first; second}))
          in
          Single_or_incremental (seconds @ (List.rev rest))
    in *)
    { commandline = conf.commandline
    ; per_loci =
      List.map paa ~f:(fun ((release, pl_locus), p, allele_arr, number_alleles) ->
        { imgt_release = release
        ; pl_locus
        ; allele_arr
        ; likelihood = p.init_global_state ()
        ; zygosity   = Zygosity_mixed.init number_alleles
        })
    ; per_reads = []
    ; single_read
    ; ordinary_paired
    ; opt       = output_opt
    }

  let single oc { single_read; _} fqi =
    Fastq_items.single oc fqi ~k:(fun name read read_errors ->
      { Output.name
      ; d = single_read read read_errors
      })

  let paired oc { ordinary_paired; _} fq1 fq2 =
    Fastq_items.paired oc fq1 fq2 ~k:(fun name rd1 re1 rd2 re2 ->
      { Output.name
      ; d = ordinary_paired name rd1 re1 rd2 re2
      })

  let cons_per_reads state pr =
    state.per_reads <- pr :: state.per_reads

  let assign_to_per_locus oc name { pl_locus; likelihood; zygosity; _ } llhd_and_pos =
    let msg =
      sprintf "For %s merging into %s zygosity array read result: %d"
        name
        (Nomenclature.show_locus pl_locus)
        (Pm.length llhd_and_pos)
    in
    time oc msg (fun () ->
      Likelihoods_and_zygosity.add_ll_and_lz likelihood zygosity llhd_and_pos)

  let assign_to_best oc name state best_locus llhd_and_pos =
    List.iter state.per_loci ~f:(fun pl ->
      if best_locus = pl.pl_locus then
        assign_to_per_locus oc name pl llhd_and_pos)

  let merge oc state pr =
    let locus, best_stat =
      match pr.Output.d with
      | Single_or_incremental l -> reduce_soi l
      | MPaired l               -> reduce_mpr l
    in
    let llhd_and_pos =
      match best_stat with
      | Sp.Single s        -> Forward.to_single_pms s
      (* Since a paired read is a physical thing there should only be 1
       * likelihood for each allele, but we'll separate the positions. *)
      | Sp.Paired (s1, s2) -> Forward.merge_paired_pms s1 s2
    in
    (* modify state by side effect *)
    assign_to_best oc pr.Output.name state locus llhd_and_pos;
    let aaps =
      match pr.Output.d with
      | Single_or_incremental l ->
          let nl = map_soi_to_aap state.per_loci state.opt.allele_depth l in
          Single_or_incremental nl
      | MPaired l               ->
          let nl = map_pr_to_aap state.per_loci state.opt.allele_depth l in
          MPaired nl
    in
    let iri = { locus; llhd_and_pos; aaps } in
    cons_per_reads state { pr with Output.d = iri }

  let per_locus_to_output_per_locus pl zb =
    Output.create_per_locus pl.imgt_release pl.pl_locus
      pl.allele_arr pl.likelihood zb

  let finalize { commandline; opt; per_reads; per_loci; _} =
    let f_pl_lst, zygosities_by_locus_assoc =
      List.map per_loci ~f:(fun pl ->
        let zb = Zygosity_mixed.best opt.depth.Output.num_zygosities pl.zygosity in
        let opl = per_locus_to_output_per_locus pl zb in
        opl, (pl.pl_locus, (zb, pl.allele_arr)))
      |> List.split
    in
    let la allele_arr i = fst (allele_arr.(i)) in
    let assign { locus; llhd_and_pos; _ } =
      let open Option in
      List.Assoc.get locus zygosities_by_locus_assoc >>= (fun (zlst, aa) ->
        let lopt =
          List.fold_left zlst ~init:None ~f:(fun s (_, _, ijlist) ->
            List.fold_left ijlist ~init:s ~f:(fun s (i, j, _nr1, _nr2) ->
              let li, _pi = Pm.get llhd_and_pos i in
              let lj, _pj = Pm.get llhd_and_pos j in
              let l, k = if Lp.(lj <= li) then li, i else lj, j in
              match s with
              | None         -> Some (l, (la aa k))
              | Some (bl, _) -> if Lp.(l <= bl) then s else Some (l, la aa k)))
              in
              Option.map lopt ~f:(fun (_likelihood, allele) ->
                  locus, allele))
    in
    let f_per_reads =
      List.map per_reads ~f:(fun pr ->
        let most_likely = assign pr.Output.d in
        { pr with Output.d = { most_likely; aaps = pr.Output.d.aaps }})
    in
    { f_commandline = commandline
    ; f_opt         = opt
    ; f_pl_lst
    ; f_per_reads
    }

  let output { f_commandline; f_opt; f_pl_lst; f_per_reads } oc =
    let header  =
      { Output.prohlatype_version = Version.version
      ; Output.commandline        = f_commandline }
    in
    let ot = Output.create header f_pl_lst f_per_reads in
    let of_ =
      match f_opt.output_format with
      | `TabSeparated ->
        let out_to_string rrp =
          sprintf "\n%s" (final_read_info_to_string rrp)
        in
        Output.TabSeparated out_to_string
      | `Json ->
        Output.Json final_read_info_to_yojson
    in
    Output.f oc ot f_opt.depth of_

end (* Multiple_loci *)

module type Execution = sig

  type conf
  type input
  type state

  val init : out_channel
           -> (int -> input)
           -> conf
           -> int
           -> state

  val across_fastq : log_oc:out_channel
                   -> data_oc:out_channel
                   -> conf
                   -> ?number_of_reads:int
                   -> specific_reads:string list
                   -> [ `Setup of int -> input | `Set of state ]
                   -> string
                   -> unit

  val across_paired : log_oc:out_channel
                    -> data_oc:out_channel
                    -> finish_singles:bool
                    -> conf
                    -> ?number_of_reads:int
                    -> specific_reads:string list
                    -> [ `Setup of int -> input | `Set of state ]
                    -> string
                    -> string
                    -> unit

end (* Execution *)

module Sequential (W : Worker) :
  Execution with type conf = W.conf
             and type input = W.input
             and type state = W.state
  = struct

  type conf = W.conf
  type input = W.input
  type state = W.state

  let init oc rp conf read_length =
    let pt =
      time oc (sprintf "Setting up ParPHMM transitions with %d read_length"
                read_length)
        (fun () -> rp read_length)
    in
    time oc "Allocating forward pass workspace"
      (fun () -> W.init oc conf ~read_length pt)

  let single oc conf =
    fun acc fqi ->
      match acc with
      | `Setup rp ->
          let read_length = String.length fqi.Biocaml_unix.Fastq.sequence in
          let s = init oc rp conf read_length in
          W.merge oc s (W.single oc s fqi);
          `Set s
      | `Set s ->
          W.merge oc s (W.single oc s fqi);
          `Set s

  let paired oc conf =
    fun acc fq1 fq2 ->
      match acc with
      | `Setup rp ->
          let read_length = String.length fq1.Biocaml_unix.Fastq.sequence in
          let s = init oc rp conf read_length in
          W.merge oc s (W.paired oc s fq1 fq2);
          `Set s
      | `Set s ->
          W.merge oc s (W.paired oc s fq1 fq2);
          `Set s

  let across_fastq ~log_oc ~data_oc conf ?number_of_reads ~specific_reads init file =
    let f = single log_oc conf in
    try
      Fastq.fold ?number_of_reads ~specific_reads file ~init ~f
      |> function
          | `Setup _ -> fprintf log_oc "Didn't find any reads."
          | `Set  s  -> W.output (W.finalize s) data_oc
    with Fastq_items.Read_error_parsing e ->
      fprintf log_oc "%s" e

  let across_paired ~log_oc ~data_oc ~finish_singles conf ?number_of_reads ~specific_reads
    init file1 file2 =
    let f = paired log_oc conf in
    try
      begin
        if finish_singles then
          let fs = single log_oc conf in
          let ff = fs in
          Fastq.fold_paired_both ?number_of_reads ~specific_reads file1 file2
            ~init ~f ~ff ~fs
        else
          Fastq.fold_paired ?number_of_reads ~specific_reads file1 file2 ~init ~f
      end
      |> function
          | `BothFinished o
          | `FinishedSingle o
          | `OneReadPairedFinished (_, o)
          | `StoppedByFilter o
          | `DesiredReads o ->
              match o with
              | `Setup _ -> fprintf log_oc "Didn't find any reads."
              | `Set s   -> W.output (W.finalize s) data_oc
    with Fastq_items.Read_error_parsing e ->
      fprintf log_oc "%s" e

end (* Sequential *)

module Parallel (W : Worker) = struct

  let init oc rp conf read_length =
    let pt =
      time oc (sprintf "Setting up ParPHMM transitions with %d read_length"
        read_length)
        (fun () -> rp read_length)
    in
    time oc "Allocating forward pass workspace"
      (fun () -> W.init oc conf ~read_length pt)

  let across_fastq ~log_oc ~data_oc conf ?number_of_reads ~specific_reads ~nprocs
      state file =
    try
      let map fqi = W.single log_oc state fqi in
      let mux rr = W.merge log_oc state rr in
      Fastq.fold_parany ?number_of_reads ~specific_reads ~nprocs ~map ~mux file;
      W.output (W.finalize state) data_oc
    with Fastq_items.Read_error_parsing e ->
      fprintf log_oc "%s" e

  let across_paired ~log_oc ~data_oc conf ?number_of_reads ~specific_reads ~nprocs
      state file1 file2 =
    let map = function
      | Sp.Single fqi      -> W.single log_oc state fqi
      | Sp.Paired (f1, f2) -> W.paired log_oc state f1 f2
    in
    let mux rr = W.merge log_oc state rr in
    try
      Fastq.fold_paired_parany ?number_of_reads ~specific_reads
        ~nprocs ~map ~mux file1 file2;
      W.output (W.finalize state) data_oc
    with Fastq_items.Read_error_parsing e ->
      fprintf log_oc "%s" e

end (* Parallel *)
