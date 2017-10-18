(* Aggregation of the forward passes. *)

open Util

module Pm = Partition_map

let max_pm = Pm.fold_values ~init:ParPHMM.Lp.zero ~f:ParPHMM.Lp.max

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
    Lp.max prev_threshold (proc.maximum_match ())

  let new_value prev_threshold proc = function
    | Filtered _  -> prev_threshold
    | Completed _ -> new_ prev_threshold proc

  let update prev_threshold proc pr =
    match prev_threshold with
    | `Don't  -> `Don't
    | `Start  -> `Set (new_value Lp.zero proc pr)
    | `Set pv -> `Set (new_value pv proc pr)

end  (* Past_threshold *)

(* Consider a single orientation: either regular or reverse complement. *)
let specific_orientation ?prev_threshold proc pass_result_map reverse_complement read read_errors =
  let open ParPHMM in
  match proc.single ?prev_threshold ~read ~read_errors reverse_complement with
  | Filtered m    -> Filtered m
  | Completed ()  -> Completed (pass_result_map proc reverse_complement)

(* Sometimes (such as in the beginning) we do not know whether to prefer a
   regular or a reverse-complement orientation for read. *)
module Orientation = struct

  open ParPHMM

  type 'a t =
    { regular     : 'a pass_result
    ; complement  : 'a pass_result
    }
    [@@deriving yojson]

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

(* Represent a diagnostc result of reads. *)
module Alleles_and_positions = struct

  open ParPHMM

  type datum = per_allele_datum
  and t = datum list
    [@@deriving yojson]

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

end (* Alleles_and_positions *)

(* Storing the log likelihoods of heterozygous pairs.
   Homozygous is equivalent to the individual pairs. *)
module Zygosity_array = struct

  open ParPHMM

  (* Not exactly the triangular number, T_(n-1);
     just have less subtractions this way. *)
  let triangular_number n =
    n * (n - 1) / 2

  let triangular_inverse t =
    let tf = float t in
    (truncate (sqrt (1.0 +. 8.0 *. tf)) - 1) / 2

  let k_n n i j =
    triangular_number n - triangular_number (n - i) + j - i - 1

  let kinv_n n k =
    let tn = triangular_number n in
    let i  = n - 1 - (triangular_inverse (tn - k - 1)) - 1 in
    let j  = k + i + 1 - triangular_number n + triangular_number (n - i) in
    i, j

  (* Test the k <-> i, j mapping logic. *)
  let test n =
    let kn = k_n n in
    let ki = kinv_n n in
    for i = 0 to n - 2 do
      for j = i + 1 to n - 1 do
        let k = kn i j in
        let ni, nj = ki k in
        printf "i: %d\tj: %d\t:k: %d\ ---- ni: %d\tnk: %d %b\n"
          i j k ni nj (i = ni && j = nj)
      done
    done

  (* Top triangle of pair wise comparison *)
  type t =
    { n : int
    ; l : Lp.t array
    }

  let storage_size number_alleles =
    triangular_number number_alleles                               (* T_(n-1) *)

  let init number_alleles =
    let s = storage_size number_alleles in
    { n = number_alleles
    ; l = Array.make s Lp.zero
    }

  let add t ~allele1 ~allele2 likelihood =
    let i = min allele1 allele2 in
    let j = max allele1 allele2 in
    let k = k_n t.n i j in
    t.l.(k) <- Lp.(t.l.(k) + likelihood)

  type best_size =
    | NonZero
    | Spec of int
    | NoSpec

  (* log likelihood, probability, index 1st, index 2nd *)
  let best size t likelihood_arr =
    let desired_size, nz =
      match size with
      | NonZero -> storage_size t.n, true
      | Spec n  -> max 1 n, false
      | NoSpec  -> storage_size t.n, false
    in
(* Create a sorted, grouped, and constrained, list of the paired likelihoods.
  The total size (4000^2 = 16_000_000) is a little bit unwieldly and at this
  point we don't care about the tail only about where there is actual
  probability mass. But in-order to properly normalize we need the max
  likelihood of all values (for the log-sum-exp trick) and a normalizing
  constant.

  [sorted_insert] keeps an ascending list list of the [desired_size] highest
  likelihoods (all indices of the inner list have the same likelihood: grouped),
  such that only these values are afterwards exported.
  *)
    let sorted_insert l i j current_size lst =
      let rec insert lst =
        match lst with
        | []     -> [l, 1, [i,j]]
        | h :: t ->
            let hl, hn, hv = h in
            if l > hl then
              h :: insert t
            else if l = hl then
              (hl, hn + 1, (i,j) :: hv) :: t
            else (* l < hl *)
              (l, 1, [i,j]) :: lst
      in
      match lst with
      | []
      | _ :: [] -> current_size + 1, insert lst
      | h :: t  ->
        let hl, hn, hv = h in
        if l > hl then begin
           if current_size - hn < desired_size then
            current_size + 1, h :: insert t
          else
            current_size + 1 - hn, insert t
        end else if l = hl then
          current_size + 1, ((hl, hn + 1, (i,j) :: hv) :: t)
        else (* l < hl *) begin
          if current_size >= desired_size then
            current_size, lst
          else
            current_size + 1, insert lst
        end
    in
    let _fk, fcs, fres =                                     (* Heterozygous. *)
      Array.fold_left t.l ~init:(0, 0, [])
        ~f:(fun (k, cs, acc) l ->
              let i, j = kinv_n t.n k in
              let ncs, nacc = sorted_insert l i j cs acc in
              (k + 1, ncs, nacc))
    in
    let _fk, _cs, res =                                        (* Homozygous. *)
      Array.fold_left likelihood_arr ~init:(0, fcs, fres)
        ~f:(fun (i, cs, acc) l ->
              let ncs, nacc = sorted_insert l i i cs acc in
              (i + 1, ncs, nacc))
    in
    let most_likely = List.rev_map ~f:(fun (l, _n, ijlst) -> (l, ijlst)) res in
    let maxl, _ = List.hd_exn most_likely in
    (*printf "maxl: %20.20f \n%!" maxl; *)
    let prob = Lp.probability ~maxl in
(* Compute the normalizing constant. *)
    let ns =
      let f s l = s +. prob l in
      let heterozygous_sum = Array.fold_left ~init:0. ~f t.l in
      Array.fold_left ~init:heterozygous_sum ~f likelihood_arr
    in
    if nz then begin
      (*printf "hoooray nonzero requested! %20.20f \n%!" ns; *)
      List.filter_map most_likely ~f:(fun (l, ij) ->
        (*printf "l: %0.20f, p:%0.20f, %s \n%!" l (prob l)
          (string_of_list ~sep:";" ~f:(fun (i,j) -> sprintf "%d,%d" i j) ij); *)
        let p = prob l /. ns in
        (* This isn't necessarily the right notion of zero. *)
        if p > !ParPHMM.dx then
          Some (l, prob l /. ns, ij)
        else
          None)
    end else
      List.map most_likely ~f:(fun (l, ij) -> (l, prob l /. ns, ij))

end (* Zygosity_array *)

module Likelihoods_and_zygosity = struct

  (* TODO: We're opening Util which has likehood defined. *)
  open ParPHMM

  let add_ll state_llhd llhd =
    Pm.iter_set llhd ~f:(fun i v ->
      state_llhd.(i) <- Lp.(state_llhd.(i) + v))

  let add_lz zygosity llhd =
    Pm.iter_set llhd ~f:(fun allele1 v1 ->
      Pm.iter_set llhd ~f:(fun allele2 v2 ->
        (* The Zygosity array only stores the upper triangle of the full
           pairwise likelihood matrix. Therefore, ignore the diagonal and
           lower triangular elements to avoid double counting (per read). *)
        if allele1 >= allele2 then
          ()
        else
          Zygosity_array.add zygosity ~allele1 ~allele2 (max v1 v2)))

  (* State is labled because it is modified. *)
  let add_ll_and_lz state_llhd zygosity llhd =
    add_ll state_llhd llhd;
    add_lz zygosity llhd

end (* Likelihoods_and_zygosity *)

(* Outputing results.

   TODO: Expose a separator argument. *)
module Output = struct

  open ParPHMM

  type locus = Nomenclature.locus
  and per_allele_log_likelihood =
    { allele : string
    ; a_llhd : Lp.t
    ; alters : MSA.Alteration.t list
    }
  and zygosity_log_likelihood =
    { allele1 : string
    ; allele2 : string
    ; z_llhd  : Lp.t
    ; prob    : float
    }
  and per_locus =
    { locus       : locus
    ; zygosity    : zygosity_log_likelihood list
    ; per_allele  : per_allele_log_likelihood array
    }
  and 'p per_read =
    { name : string
    ; d    : 'p
    }
  and 'p t =
    { per_loci  : per_locus list
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
    ; num_zygosities  : Zygosity_array.best_size
    ; num_per_read    : int option
    }

  let default_depth =
    { num_likelihoods = None
    ; num_zygosities  = Zygosity_array.NonZero
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

  let create_pl size locus allele_arr allele_llhd_arr zt =
    let per_allele =
      Array.mapi allele_arr ~f:(fun i (allele, alters) ->
        { allele; a_llhd = allele_llhd_arr.(i); alters })
    in
    Array.sort ~cmp:cmp_pall per_allele;
    let zygosity =
      let zb = Zygosity_array.best size zt allele_llhd_arr in
      List.map zb ~f:(fun (l, p, ijlist) ->
        List.map ijlist ~f:(fun (i, j) ->
          { allele1 = fst allele_arr.(i)
          ; allele2 = fst allele_arr.(j)
          ; z_llhd  = l
          ; prob    = p
          }))
      |> List.concat
    in
    { locus
    ; per_allele
    ; zygosity
    }

  let create per_loci per_reads =
    { per_loci; per_reads }

  let shrink d t =
    let shrink_list l =
      Option.value_map ~default:l ~f:(List.take l)
    in
    let shrink_array a =
      Option.value_map ~default:a ~f:(fun l ->
          Array.sub a ~pos:0 ~len:(min l (Array.length a)))
    in
    let shrink_pl { locus; per_allele; zygosity } =
      { locus
      ; per_allele = shrink_array per_allele d.num_likelihoods
      ; zygosity = zygosity (* Shrunk on construction. *)
      }
    in
    { per_loci = List.map t.per_loci ~f:shrink_pl
    ; per_reads = shrink_list t.per_reads d.num_per_read
    }

  let array_iter_on_non_empty f b = function
    | [||] -> ()
    | a    -> ignore (b ()); Array.iter ~f a

  let list_iter_on_non_empty f b = function
    | []   -> ()
    | lst  -> ignore (b ()); List.iter ~f lst

  let tabs oc d_to_string { per_loci; per_reads } =
    let open Nomenclature in
    let alters_to_string a =
      string_of_list ~show_empty:false ~sep:";"
        ~f:MSA.Alteration.to_string a
    in
    let fprint_per_allele_log_likelihood { allele; a_llhd; alters } =
      fprintf oc "%16s\t%s\t%s\n"
        allele (Lp.to_string ~precision:20 a_llhd) (alters_to_string alters)
    in
    let fprint_zygosity_log_likelihood { allele1; allele2; z_llhd; prob} =
      fprintf oc "%16s\t%16s\t%s\t%0.4f\n"
        allele1 allele2 (Lp.to_string ~precision:20 z_llhd) prob
    in
    List.iter per_loci ~f:(fun { locus; per_allele; zygosity} ->
      list_iter_on_non_empty
        fprint_zygosity_log_likelihood
        (fun () -> fprintf oc "Zygosity %s:\n" (show_locus locus))
        zygosity;
      array_iter_on_non_empty
        fprint_per_allele_log_likelihood
        (fun () -> fprintf oc "Likelihood %s:\n" (show_locus locus))
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

  let single fqi ~k =
    let open Core_kernel.Std in
    let open Biocaml_unix.Fastq in
    time (sprintf "updating on single read %s" fqi.name) (fun () ->
      match Fastq.phred_log_probs fqi.qualities with
      | Result.Error e        -> repf "%s" (Error.to_string_hum e)
      | Result.Ok read_errors -> k fqi.name fqi.sequence read_errors)

  let paired fq1 fq2 ~k =
    let open Biocaml_unix.Fastq in
    let open Core_kernel.Std in
    assert (fq1.name = fq2.name);
    time (sprintf "updating on double read %s" fq1.name) (fun () ->
      match Fastq.phred_log_probs fq1.Biocaml_unix.Fastq.qualities with
      | Result.Error e       -> repf "%s" (Error.to_string_hum e)
      | Result.Ok re1 ->
          match Fastq.phred_log_probs fq2.Biocaml_unix.Fastq.qualities with
          | Result.Error e       -> repf "%s" (Error.to_string_hum e)
          | Result.Ok re2 ->
              k fq1.name fq1.sequence re1 fq2.sequence re2)

end (* Fastq_items *)

type single_conf =
  { prealigned_transition_model : bool
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
  }

let single_conf ?allele ?insert_p ?band ?max_number_mismatches ?split
  ~prealigned_transition_model
  ~past_threshold_filter ~check_rc () =
    { prealigned_transition_model
    ; allele
    ; insert_p
    ; band
    ; max_number_mismatches
    ; past_threshold_filter
    ; check_rc
    ; split
    }

type ('state, 'read_result) worker =
  { single  : Biocaml_unix.Fastq.item -> 'read_result
  ; paired  : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> 'read_result
  ; merge   : 'state -> 'read_result -> unit
  ; output  : 'state -> out_channel -> unit
  }

type output_opt =
  { allele_depth  : int
  ; output_format : Output.format_spec
  ; depth         : Output.depth
  }

module Forward = struct

  open ParPHMM

  type stat =
    { aap        : Alleles_and_positions.t
    ; likelihood : Lp.t mt
    }

  let proc_to_stat allele_depth proc _comp =
    { aap        = proc.best_allele_pos allele_depth
    ; likelihood = proc.per_allele_llhd ()
    }

  let forward allele_depth past_threshold proc read read_errors =
    time (sprintf "forward of allele_depth: %d; past_threshold: %s; read: %s"
            allele_depth (Past_threshold.to_string past_threshold)
              (String.sub_exn read ~index:0 ~length:10))
    (fun () ->
      Orientation.check past_threshold (proc_to_stat allele_depth)
        proc read read_errors)

  type aapt = Alleles_and_positions.t Orientation.t [@@deriving yojson]

  type state =
    { locus             : Nomenclature.locus
    ; per_allele_lhood  : Lp.t array
    ; zygosity          : Zygosity_array.t
    ; mutable per_reads : aapt Output.per_read list
    }

  type read_result =
    { name  : string
    ; stat  : stat Orientation.t single_or_paired
    }

  let just_aap = function
    | Filtered m   -> Filtered m
    | Completed rm -> Completed rm.aap

  let just_aap_or = Orientation.map ~f:just_aap

  let take_regular_by_likelihood  r c =
    higher_value_in_first r.likelihood c.likelihood

  let cons_per_reads state name stat_or =
    let d = just_aap_or stat_or in
    state.per_reads <- { Output.name; d } :: state.per_reads

  let update_state_single state name stat_or =
    let take_regular = take_regular_by_likelihood in
    cons_per_reads state name stat_or;
    let pr = Orientation.most_likely_between stat_or ~take_regular in
    match pr with
    | Filtered m        ->
        printf "Both orientations filtered %s: %s\n" name m
    | Completed (_c, s) ->
        Likelihoods_and_zygosity.add_ll_and_lz
          state.per_allele_lhood state.zygosity s.likelihood

  let output_state oo allele_arr state oc =
    let take_regular r c = Alleles_and_positions.descending_cmp r c <= 0 in
    let size = oo.depth.Output.num_zygosities in
    let ot =
      let pl = Output.create_pl size state.locus allele_arr
                state.per_allele_lhood state.zygosity in
      Output.create [pl] state.per_reads
    in
    let of_ =
      match oo.output_format with
      | `TabSeparated ->
        let aapt_to_string aapt =
          sprintf "\n\t%s"
            (Orientation.to_string ~take_regular
                Alleles_and_positions.string_of_list aapt)
        in
        Output.TabSeparated aapt_to_string
      | `Json ->
        Output.Json aapt_to_yojson
    in
    Output.f oc ot oo.depth of_

  let to_proc_allele_arr ?allele ?insert_p ?max_number_mismatches ?band ?split
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

  type opt = output_opt

  let init conf ~read_length parPHMM_t oo =
    let { prealigned_transition_model; allele; insert_p; band
        ; max_number_mismatches; past_threshold_filter ; split; _} = conf
    in
    let open ParPHMM in
    let proc, allele_arr =
      to_proc_allele_arr ?allele ?insert_p ?max_number_mismatches ?band ?split
        ~prealigned_transition_model read_length parPHMM_t
    in
    let initial_pt = Past_threshold.init past_threshold_filter in
    (* Curry away arguments and discard final threshold; can't use it. *)
    let forward pt r re = fst (forward oo.allele_depth pt proc r re) in
    let rec state =
      { locus = parPHMM_t.locus; per_allele_lhood; zygosity; per_reads; }
    and per_allele_lhood = proc.init_global_state ()
    and zygosity = Zygosity_array.init (Array.length allele_arr)
    and per_reads = []
    in
    let rec single fqi =
      Fastq_items.single fqi ~k:(fun name read read_errors ->
        { name; stat = Single (forward initial_pt read read_errors) })
    and paired fq1 fq2 =
      let take_regular = take_regular_by_likelihood in
      Fastq_items.paired fq1 fq2 ~k:(fun name rd1 re1 rd2 re2 ->
        let r1 = forward initial_pt rd1 re1 in
        match Orientation.most_likely_between r1 ~take_regular with
        | Filtered m  ->
            let r2 = forward initial_pt rd2 re2 in
            { name; stat = Paired (r1, r2)}
        | Completed (c, _)  ->
            let r2 = Orientation.specific proc (proc_to_stat oo.allele_depth)
                      (not c) rd2 re2 in
            { name; stat = Paired (r1, r2)})
    and merge state rr =
      match rr.stat with
      | Single stat_or  -> update_state_single state rr.name  stat_or
      | Paired (s1, s2) -> update_state_single state (rr.name ^ " 1") s1;
                           update_state_single state (rr.name ^ " 2") s2
    and output t oc =
      output_state oo allele_arr t oc
    and t = { single; paired; merge; output }
    in
    t, state

end (* Forward *)

(* A map of reads to a Viterbi decoded path *)
module Viterbi = struct

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
    | Single s        -> separate reverse_complement s
    | Paired (r1, r2) ->
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
    }

  type opt = unit

  let init conf ~read_length parPHMM_t () =
    let { prealigned_transition_model; allele; insert_p; _ } = conf in
    let open ParPHMM in
    (* Relying on reference being first in the specified ParPHMM.t allele list. *)
    let allele = Option.value ~default:(fst parPHMM_t.alleles.(0)) allele in
    let labels = allele ^ " ", "read " in
    let s, p =
      setup_single_allele_viterbi_pass ?insert_p ~allele
        ~prealigned_transition_model read_length parPHMM_t
    in
    let state = { results = [] } in
    let rec t = { single; paired; merge; output }
    and single =
      Fastq_items.single ~k:(fun name read read_errors ->
        { name; res = s read read_errors})
    and paired fq1 fq2 =
      Fastq_items.paired fq1 fq2 ~k:(fun name rd1 re1 rd2 re2 ->
        { name; res =  p rd1 re1 rd2 re2})
    and merge state pr =
      state.results <- pr :: state.results
    and output state oc =
      List.iter state.results ~f:(fun {name; res } ->
        fprintf oc "%s:\n%s\n" name (viterbi_result_to_string ~labels res))
    in
    t, state

end (* Viterbi *)

(** Single Loci Drivers **)
module type S = sig

  type state
  type read_result

  type opt

  val init : single_conf
           -> read_length:int
           -> ParPHMM.t
           -> opt
           -> (state, read_result) worker * state

end (* S *)

module Make_single (S : S) = struct

  let init = S.init

  let single d fqi = d.single fqi
  let paired d fq1 fq2 = d.paired fq1 fq2
  let merge d s r = d.merge s r
  let output d s oc = d.output s oc

end (* Make_single *)

(* Combine results from multiple loci. *)

type multiple_conf =
  { prealigned_transition_model : bool
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
  (*~incremental_pairs *)
  (*?first_read_best_log_gap *)
  () =
    { prealigned_transition_model
    ; insert_p
    ; band
    ; max_number_mismatches
    ; past_threshold_filter
    ; split
    (*; incremental_pairs
      ; first_read_best_log_gap *)
    }

module Multiple_loci = struct

  open ParPHMM

  type 'single_result paired_dependent =
      { first_o : bool
      ; first   : 'single_result
      ; second  : 'single_result pass_result
      }
    (* Orientation and gene of 2nd read is determined by 1st. There is no
       {pass_result} since we're not going to apply a filter to the 2nd read.
       Note that one could "naturally" combine the result of the first by
       passing the first's emission (via per_allele_hood) into the second
       read's forward pass (this is similar to thinking of the 2 reads as just
       one long read), but that will dramatically slow down the pass (it will
       start out already branched). It is better to keep the two likelihoods
       separate and then just multiply (log -> add) the results together. This
       is the desired outcome for paired reads as it is the most
       descriminative. *)
  and 'single_result single_or_incremental =
    | SingleRead of 'single_result Orientation.t
    (* We'll test both regular and complement options. *)
    | PairedDependent of 'single_result paired_dependent
    [@@deriving yojson]

  let orientation_char_of_bool reverse_complement =
    if reverse_complement then 'C' else 'R'

  let pass_result_to_short_string sep sr_to_string = function
    | Filtered m -> sprintf "F%c%s" sep m
    | Completed r -> sr_to_string r

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

  let take_regular = Forward.take_regular_by_likelihood

  let most_likely_orientation or_ =
    Orientation.most_likely_between ~take_regular or_

  (* true if s1 is better than s2 *)
  let better_soi s1 s2 =
    let l s = max_pm s.Forward.likelihood in
    match s1, s2 with
    | PairedDependent p1, PairedDependent p2  ->
        begin
          match p1.second, p2.second with
          | Filtered _,    _             -> false
          | _,             Filtered _    -> true
          | Completed ss1, Completed ss2 -> Lp.(l p2.first + l ss2 < l p1.first + l ss1)
        end
    | PairedDependent p1, SingleRead sr2      ->
        begin
          match most_likely_orientation sr2 with
          | Filtered  _                   -> true
          | Completed (_, s)  ->
              begin match p1.second with
              | Filtered  _               -> Lp.(l s <= l p1.first)
              | Completed _               -> true
              end
        end
    | SingleRead sr1,     PairedDependent p2  ->
        begin
          match most_likely_orientation sr1 with
          | Filtered _                    -> false
          | Completed (_, s)  ->
              begin match p2.second with
              | Filtered  _               -> Lp.(l p2.first <= l s)
              | Completed _               -> false
              end
        end
    | SingleRead sr1,     SingleRead sr2      ->
        begin
          match most_likely_orientation sr1 with
          | Filtered _                    -> false
          | Completed (_, l1)      ->
            begin match most_likely_orientation sr2 with
            | Filtered  _                 -> true
            | Completed (_, l2)           -> Lp.(l l2 <= l l1)
            end
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
    ; second  : 'single_result pass_result
    }
    (* One of the orientations of the first wasn't filtered and we used
       it to orient the second read. There is a {pass_result} for second
       because in the cases where we have a split, those may trigger a
       filter that stops computation. *)
  and 'single_result paired =
    | FirstFiltered of 'single_result first_filtered
    | FirstOrientedSecond of 'single_result first_oriented_second
    [@@deriving yojson]

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

  type pd = Alleles_and_positions.t pp

  type 'sr per_read =
    { name  : string
    ; rrs   : 'sr pp
    }

  (* Easily convert to Output.per_locus *)
  type per_locus =
    { pl_locus    : Nomenclature.locus
    ; allele_arr  : (string * MSA.Alteration.t list) array
    ; likelihood  : Lp.t array
    ; zygosity    : Zygosity_array.t
    }
  type state =
    { per_loci          : per_locus list
    ; mutable per_reads : pd Output.per_read list
    }

  type read_result = Forward.stat per_read

  (* A PairedDependent with completed second should be chosen over any SingleRead's.
     A PairedDependent with a Filtered second should be compared by first.
     Compare 2 PairedDependent by both.
   *)
  let reduce_soi soi_lst =
    List.reduce soi_lst ~f:(fun (l1, s1) (l2, s2) ->
      if better_soi s1 s2 then
        l1, s1
      else
        l2, s2)
    |> Option.bind ~f:(fun (loci, soi) ->
        match soi with
        | PairedDependent { first; second; _} ->
            begin match second with
            | Filtered _  -> Some (loci, Single first)
            | Completed s -> Some (loci, Paired (first, s))
            end
        | SingleRead  or_ ->
            begin match most_likely_orientation or_ with
            | Filtered _      -> None             (* After all that! *)
            | Completed (_,s) -> Some (loci, Single s)
            end)

  let map_soi_to_aap lst =
    List.map lst ~f:(fun (loci, soi) ->
      let n =
        match soi with
        | PairedDependent { first_o; first; second } ->
            PairedDependent
              { first_o
              ; first = first.Forward.aap
              ; second = map_completed ~f:(fun s -> s.Forward.aap) second
              }
        | SingleRead oi                             ->
           SingleRead (Forward.just_aap_or oi)
      in
      (loci, n))

  (* We can only compare based upon the result of the 2nd read, but if a
     FirstOrientedSecond scenario is best, we want to use the likelihood from
     BOTH reads in the inference. *)
  let reduce_mpr pr_lst =
    let update loci stat pr best =
      let l = max_pm stat.Forward.likelihood in
      match best with
      | None                        -> Some (l, loci, pr)
      | Some (bl, _, _) when l > bl -> Some (l, loci, pr)
      | Some _                      -> best
    in
    let rec loop current_best = function
      | []               -> current_best
      | (loci, pr) :: tl ->
          begin match pr with
          | FirstFiltered { ff_second } ->
              begin match most_likely_orientation ff_second with
              | Filtered  _      -> loop current_best tl
              | Completed (_, s) -> loop (update loci s (Single s) current_best) tl
              end
          | FirstOrientedSecond { first; second } ->
              begin match second with
              | Filtered _  -> loop current_best tl
              | Completed s -> loop (update loci s (Paired (first, s)) current_best) tl
              end
          end
    in
    loop None pr_lst

  let map_pr_to_aap lst =
    List.map_snd lst ~f:(function
        | FirstFiltered { ff_first; ff_second } ->
            FirstFiltered { ff_first; ff_second = Forward.just_aap_or ff_second }
        | FirstOrientedSecond { first_o; first; second } ->
            FirstOrientedSecond
              { first_o
              ; first   = first.Forward.aap
              ; second  = map_completed ~f:(fun s -> s.Forward.aap) second
              })

  let cons_per_reads state name d =
    state.per_reads <- { Output.name; d } :: state.per_reads

  let merge state { name ; rrs } =
    let assign_to_per_locus { likelihood; zygosity; _ } = function
      | Single stat ->
          Likelihoods_and_zygosity.add_ll_and_lz
            likelihood zygosity stat.Forward.likelihood
      | Paired (s1, s2) ->
          Likelihoods_and_zygosity.add_ll_and_lz
            likelihood zygosity s1.Forward.likelihood;
          Likelihoods_and_zygosity.add_ll_and_lz
            likelihood zygosity s2.Forward.likelihood
    in
    let assign_to_best best_locus best_stat =
      List.iter state.per_loci ~f:(fun pl ->
        if best_locus = pl.pl_locus then
          assign_to_per_locus pl best_stat)
    in
    match rrs with
    | Single_or_incremental soi_lst ->
        begin match reduce_soi soi_lst with
        | None  -> eprintf "All read results filtered for %s\n%!" name
        | Some (best_locus, best_stat)  ->
            assign_to_best best_locus best_stat
        end;
        cons_per_reads state name (Single_or_incremental (map_soi_to_aap soi_lst))
    | MPaired pr_lst ->
        begin match reduce_mpr pr_lst with
        | None -> eprintf "All 2nd read results filtered for %s\n%!" name
        | Some (_bl, best_locus, best_stat)  ->
            assign_to_best best_locus best_stat
        end;
        cons_per_reads state name (MPaired (map_pr_to_aap pr_lst))

  let output_per_reads prlst oc =
    fprintf oc "Per Read:\n";
    List.iter prlst ~f:(fun {name; rrs} ->
      fprintf oc "%s\n%s\n" name
        (pp_to_string Alleles_and_positions.string_of_list rrs))

  (*   Easier to reason about this as a distance, hence >= 0
  let abs_logp_one_over_hundred = abs_float (log10 0.01)

  type incremental_classification =
    | LessLikely
    | Within
    | MoreLikely

  let classify_likelihood ?(log_gap=abs_logp_one_over_hundred) ~best newl =
    if newl <= best then begin
      if (best -. newl) < log_gap then
        Within
      else
        LessLikely
    end else begin (* newl > best *)
      if (newl -. best) < log_gap then
        Within
      else
        MoreLikely
    end *)

  let init conf opt read_length parPHMM_t_lst =
    let { prealigned_transition_model; insert_p; band; max_number_mismatches
        ; split; _ } = conf in
    let paa =
      List.map parPHMM_t_lst ~f:(fun parPHMM_t ->
        let proc, allele_arr =
          Forward.to_proc_allele_arr ?insert_p ?max_number_mismatches ?band
            ?split ~prealigned_transition_model read_length parPHMM_t
        in
        let n = Array.length allele_arr in
        parPHMM_t.locus, proc, allele_arr, n)
    in
    let initial_pt = Past_threshold.init conf.past_threshold_filter in
    let forward = Forward.forward opt.allele_depth in
    let specific_orientation p rc r re =
      specific_orientation p (Forward.proc_to_stat opt.allele_depth) rc r re
    in
    let ordinary_paired name rd1 re1 rd2 re2 =
      let _fpt1, _fpt2, res =
        List.fold_left paa ~init:(initial_pt, initial_pt, [])
          ~f:(fun (pt1, pt2, acc) (locus, p, _, _) ->
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
          ~f:(fun (best, acc, pt) (locus, p, _, _) ->
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
    let rec state = { per_loci; per_reads = []}
    and per_loci =
      List.map paa ~f:(fun (pl_locus, p, allele_arr, number_alleles) ->
        { pl_locus
        ; allele_arr
        ; likelihood = p.init_global_state ()
        ; zygosity   = Zygosity_array.init number_alleles
        })
    in
    let rec t = { single; paired; merge; output }
    and single fqi =
      Fastq_items.single fqi ~k:(fun name read read_errors ->
        { name
        ; rrs =
            List.fold_left paa ~init:(initial_pt, [])
                ~f:(fun (pt, acc) (locus, p, _, _) ->
                      let r, npt = forward pt p read read_errors in
                      npt, (locus, SingleRead r) :: acc)
            |> snd                                 (* Discard final threshold *)
            |> List.rev                                (* Restore loci order. *)
            |> fun l -> Single_or_incremental l
        })
    and paired fq1 fq2 =
     Fastq_items.paired fq1 fq2 ~k:(fun name rd1 re1 rd2 re2 ->
        { name
        ; rrs =
            (*if conf.incremental_pairs then
              incremental_paired name rd1 re1 rd2 re2
            else *)
              ordinary_paired name rd1 re1 rd2 re2
        })
    and output state oc =
      let pl_lst =
        List.map state.per_loci
          ~f:(fun { pl_locus; allele_arr; likelihood; zygosity } ->
                Output.create_pl opt.depth.Output.num_zygosities pl_locus
                  allele_arr likelihood zygosity)
      in
      let ot = Output.create pl_lst state.per_reads in
      let of_ =
        match opt.output_format with
        | `TabSeparated ->
          let out_to_string rrp =
            sprintf "\n%s" (pp_to_string Alleles_and_positions.string_of_list rrp)
          in
          Output.TabSeparated out_to_string
        | `Json ->
          Output.Json (pp_to_yojson Alleles_and_positions.to_yojson)
      in
      Output.f oc ot opt.depth of_
    in
    t, state

end (* Multiple_loci *)
