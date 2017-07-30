(* Aggregation of the forward passes. *)

open Util
open ParPHMM

type mapper_stats =
  { per_allele : (float * string) list
  ; positions  : (float * int) list
  }

let mapper_stats_to_string ?(sep='\t') s =
  let l_to_s fmt l =
    String.concat ~sep:";" (List.map l ~f:(fun (l,a) -> sprintf fmt  a l))
  in
  let al_to_s l = l_to_s "%s:%0.2f" l in
  let pl_to_s l = l_to_s "%d:%0.2f" l in
  sprintf "%s%c%s" (al_to_s s.per_allele) sep (pl_to_s s.positions)

type 'a rc_states =
  { regular         : 'a pass_result
  ; complement      : 'a pass_result
  }

(* The output is ordered so that the left most stat is the most likely.
   We add characters 'R', 'C', 'F' as little annotations. *)
let rc_stats_to_string ?(sep='\t') cmp comp_to_string ms =
  match ms.regular, ms.complement with
  | Completed r, Completed c ->
    if cmp r c then
      sprintf "R%c%s%cC%c%s"
        sep (comp_to_string r)
        sep sep (comp_to_string c)
    else
      sprintf "C%c%s%cR%c%s"
        sep (comp_to_string c)
        sep sep (comp_to_string r)
  | Completed r, Filtered m ->
      sprintf "R%c%s%cF%c%s" sep (comp_to_string r) sep sep m
  | Filtered m, Completed c ->
      sprintf "C%c%s%cF%c%s" sep (comp_to_string c) sep sep m
  | Filtered mr, Filtered mc ->
      sprintf "F%c%s%cF%c%s" sep mr sep sep mc

type mapped_complement_stats = mapper_stats rc_states

let mapped_complement_stats_to_string ?sep =
  rc_stats_to_string ?sep
    (fun r c -> fst (List.hd_exn r.positions) > fst (List.hd_exn c.positions))
    (fun r -> mapper_stats_to_string ?sep r )

type reducer_stats =
  { comp       : bool
  ; likelihood : float array
  }

let proc_to_reducer_stat proc comp =
  { comp ; likelihood = proc.per_allele_llhd ()}

(* Consider a single orientation: either regular or reverse complement. *)
let single ?prev_threshold proc pass_result_map reverse_complement read read_errors =
  match proc.single ?prev_threshold ~reverse_complement ~read ~read_errors with
  | Filtered m    -> Filtered m
  | Completed ()  -> Completed (pass_result_map proc reverse_complement)

let new_threshold prev_threshold proc =
  max prev_threshold (proc.maximum_match ())

let update_threshold prev_threshold proc = function
  | Filtered _  -> prev_threshold
  | Completed _ -> new_threshold prev_threshold proc

let init_past_threshold b =
  if b then `Start else `Don't

(* Consider both orientations. *)
let both pt pass_result_map proc read read_errors =
  let do_work ?prev_threshold rc =
    single ?prev_threshold proc pass_result_map rc read read_errors
  in
  match pt with
  | `Don't ->
      let regular = do_work false in
      let complement = do_work true in
      { regular; complement }, None
  | `Start ->
      let regular = do_work false in
      let prev_threshold = update_threshold neg_infinity proc regular in
      let complement = do_work ~prev_threshold true in
      let last_threshold = update_threshold prev_threshold proc complement in
      { regular; complement }, Some (last_threshold)
  | `Set v ->
      let regular = do_work ~prev_threshold:v false in
      let prev_threshold = update_threshold v proc regular in
      let complement = do_work ~prev_threshold true in
      let last_threshold = update_threshold prev_threshold proc complement in
      { regular; complement; }, Some (last_threshold)

let update_pt = function
  | None   -> `Start
  | Some v -> `Set v

let mapper ~n past_threshold =
  both past_threshold
    (fun proc _rc ->
      { per_allele = proc.best_alleles n
      ; positions  = proc.best_positions n
      })

let max_arr = Array.fold_left ~init:neg_infinity ~f:max

let compare_emissions e1 e2 =
  let r1 = max_arr e1 in
  let r2 = max_arr e2 in
  r1 >= r2

let most_likely_between_rc = function
  | Filtered mr, Filtered mc -> Filtered (sprintf "Both! %s, %s" mr mc)
  | Filtered _r, Completed c -> Completed c
  | Completed r, Filtered _c -> Completed r
  | Completed r, Completed c ->
    if compare_emissions r.likelihood c.likelihood then
      Completed r
    else
      Completed c

let reducer ?(check_rc=true) past_threshold proc read read_errors =
  if check_rc then begin
    let mp, pt = both past_threshold proc_to_reducer_stat proc read read_errors in
    most_likely_between_rc (mp.regular, mp.complement), pt
  end else
    (* Can't care about previous threshold! *)
    single proc proc_to_reducer_stat false read read_errors, None

type alleles_and_positions = (float * string * int) list

let aap_to_string (e,a,p) = sprintf "%s,%d,%0.5f" a p e

let aap_lst_to_string al =
  String.concat ~sep:";" (List.map al ~f:aap_to_string)

let fc_aap_to_string comp = function
  | Filtered m          ->
      sprintf "F %s" m
  | Completed aap ->
      sprintf "%c %s" (if comp then 'C' else 'R')
        (aap_lst_to_string aap)

type rammer_stats =
  { aap         : alleles_and_positions
  ; rlikelihood : float array
  }

let proc_to_rammer_stat report_size proc _comp =
  { aap         = proc.best_allele_pos report_size
  ; rlikelihood = proc.per_allele_llhd ()
  }

let rammer report_size past_threshold proc read read_errors =
  both past_threshold (proc_to_rammer_stat report_size) proc read read_errors

exception Read_error_parsing of string

let repf fmt =
  ksprintf (fun s -> raise (Read_error_parsing s)) fmt

let fqi_to_read_and_probs fqi ~k =
  let open Core_kernel.Std in
  let open Biocaml_unix.Fastq in
  time (sprintf "updating on single read %s" fqi.name) (fun () ->
    match Fastq.phred_log_probs fqi.qualities with
    | Result.Error e        -> repf "%s" (Error.to_string_hum e)
    | Result.Ok read_errors -> k fqi.name fqi.sequence read_errors)

let fqi2_to_read_and_probs fq1 fq2 ~k =
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

type single_conf =
  { allele                : string option
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
  }

module Zygosity_likelihood_array = struct

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
    ; l : float array
    }

  let storage_size number_alleles =
    triangular_number number_alleles                               (* T_(n-1) *)

  let init number_alleles =
    let s = storage_size number_alleles in
    { n = number_alleles
    ; l = Array.make s 0.
    }

  let add t ~allele1 ~allele2 likelihood =
    let i = min allele1 allele2 in
    let j = max allele1 allele2 in
    let k = k_n t.n i j in
    t.l.(k) <- t.l.(k) +. likelihood

  let best ?size t = (* TODO: Figure out the k -> i,j map *)
    let max_size = max 1 (Option.value size ~default:(storage_size t.n)) in
    let sorted_insert big_enough ((l, _, _) as p) lst =
      let rec insert lst = match lst with
        | []     -> [p]
        | h :: t ->
            let hl, _, _ = h in
            if hl < l then
              h :: insert t
            else
              p :: lst
      in
      if big_enough then begin
        match lst with
        | []            -> assert false        (* Be big enough and be empty? *)
        | (hl,_,_) :: t ->
          if l > hl then
            insert t
          else
            lst
      end else
        insert lst
    in
    let _fk, res =
      Array.fold_left t.l ~init:(0, [])
        ~f:(fun (k, acc) l ->
              let i, j = kinv_n t.n k in
              let p = l, i, j in
              let big_enough = k >= max_size in
              let nacc = sorted_insert big_enough p acc in
              (k + 1, nacc))
    in
    List.rev res

end

(* A reducer reduces the information from a set of reads into a per allele
   likelihood state. We can construct reducers that operate on just a single
   allele. *)
module Reducer = struct

  type t =
    { state     : float array  (* Need to track the per allele likelihoods *)
    ; zygosity  : Zygosity_likelihood_array.t
    ; apply     : Biocaml_unix.Fastq.item -> unit
    ; paired    : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> unit
    ; output    : out_channel -> unit
    }

  let higher_likelihood_compare (l1, r1, _i1, _a1) (l2, r2, _i2, _a2) =
    if l2 < l1 then
      -1
    else if l2 > l1 then
      1
    else
      Nomenclature.compare r1 r2

  let by_likelihood_arr allele_arr final_likelihoods =
    let o =
      Array.mapi allele_arr ~f:(fun i a ->
        final_likelihoods.(i)
        , Nomenclature.parse_to_resolution_exn a
        , i         (* Preserver the index for merging with zygosity results. *)
        , a)
    in
    Array.sort o ~cmp:higher_likelihood_compare;
    o

  let output_allele_first o oc =
    Array.iter o ~f:(fun (l,_,_,a) -> fprintf oc "%16s\t%0.20f\n" a l)

  let output_likelihood_first allele_arr final_likelihoods oc =
    let o =
      Array.mapi allele_arr ~f:(fun i a -> final_likelihoods.(i), a)
      |> Array.to_list
      |> group_by_assoc
      |> List.map ~f:(fun (l, alst) ->
            let clst = Alleles.CompressNames.f alst in
            l, String.concat ~sep:";" clst)
      |> List.sort ~cmp:(fun (l1, _) (l2, _) -> compare l2 l1) (* higher fst *)
    in
    List.iter o ~f:(fun (l, a) -> fprintf oc "%0.20f\t%16s\n" l a)

  let default_zygosity_report_size = 100

  let output_zygosity ?(size=default_zygosity_report_size) likelihood_first aa
    blarr zt oc =
    fprintf oc "Zygosity:\n";
    let zb = Zygosity_likelihood_array.best ~size zt in
    let nzb = (* Merge with blarr to weight homozygous pairs. *)
      match List.last zb with
      | None                -> eprintf "Didn't find any zygosity weights!\n"; zb
      | Some (lowest, _, _) ->
          match Array.findi blarr ~f:(fun (l, _, _, _) -> l < lowest) with
          | None            -> zb                     (* Nothing more likely. *)
          | Some len        ->
              let append_me =
                Array.sub blarr ~pos:0 ~len
                |> Array.map ~f:(fun (l, _, i, _) -> l, i, i)
                |> Array.to_list
              in
              zb @ append_me
              |> List.sort ~cmp:(fun (l1, _, _) (l2, _, _) -> compare l2 l1)
              |> fun l -> List.take l size
    in
    if likelihood_first then
      List.map nzb ~f:(fun (l, i, j) -> l, (i, j))
      |> group_by_assoc
      |> List.iter ~f:(fun (l, ijlist) ->
          let alleles_str =
            List.map ijlist ~f:(fun (i,j) -> sprintf "%s,%s" aa.(i) aa.(j))
            |> String.concat ~sep:"\t"
          in
          fprintf oc "%0.20f\t%16s\n" l alleles_str)
    else
      List.iter nzb ~f:(fun (l, i, j) -> fprintf oc "%16s\t%16s\t%0.20f\n" aa.(i) aa.(j) l)

  let output_i likelihood_first allele_arr final_likelihoods oc =
    let o = by_likelihood_arr allele_arr final_likelihoods in
    fprintf oc "Likelihood:\n";
    if likelihood_first then
      output_likelihood_first allele_arr final_likelihoods oc
    else
      output_allele_first o oc;
    o

  let output_zi likelihood_first zygosity_report_size allele_arr
    final_likelihoods zt oc =
    let blarr = output_i likelihood_first allele_arr final_likelihoods oc in
    output_zygosity ~size:zygosity_report_size likelihood_first allele_arr
      blarr zt oc

  let to_proc_allele_arr ?allele ?insert_p ?max_number_mismatches ?band
    read_length pt =
    match allele with
    | None ->
      let proc = setup_single_pass ?insert_p ?max_number_mismatches ?band
                  read_length pt in
      proc, (Array.map ~f:fst pt.alleles)  (* Drop merge info for now *)
    | Some allele ->
      let proc = setup_single_allele_forward_pass ?insert_p
                  ?max_number_mismatches read_length allele pt in
      proc, [| allele |]

  let add_log_likelihoods n state = function
    | Filtered _  -> ()              (* Should we warn about ignoring a read? *)
    | Completed r ->
        for i = 0 to n - 1 do
          state.(i) <- state.(i) +. r.likelihood.(i)
        done

  let add_log_likelihoods_from_llhd ~state zygosity likelihood =
    let n = Array.length state in
    for i = 0 to n - 1 do
      state.(i) <- state.(i) +. likelihood.(i)
    done;
    for allele1 = 0 to n - 2 do
      for allele2 = allele1 + 1 to n - 1 do
        Zygosity_likelihood_array.add zygosity ~allele1 ~allele2
          (max likelihood.(allele1) likelihood.(allele2)) ;
      done
    done

  let add_log_likelihoods_andz ~state zygosity = function
    | Filtered _  -> ()              (* Should we warn about ignoring a read? *)
    | Completed r -> add_log_likelihoods_from_llhd ~state zygosity r.likelihood

  let init conf likelihood_first zygosity_report_size read_length pt =
    let { allele; insert_p; band; max_number_mismatches; past_threshold_filter
        ; check_rc } = conf
    in
    let open ParPHMM in
    let proc, allele_arr =
      to_proc_allele_arr ?allele ?insert_p ?max_number_mismatches ?band
        read_length pt
    in
    let past_threshold = init_past_threshold past_threshold_filter in
    (* Discard resulting past threshold *)
    let r rd rp = fst (reducer ~check_rc past_threshold proc rd rp) in
    let rec t = { state; zygosity; apply; paired; output }
    and state = proc.init_global_state ()
    and zygosity = Zygosity_likelihood_array.init (Array.length allele_arr)
    and apply fqi =
      fqi_to_read_and_probs fqi ~k:(fun _name read read_errors ->
        r read read_errors
        |> add_log_likelihoods_andz ~state:t.state t.zygosity)
    and paired fq1 fq2 =
      fqi2_to_read_and_probs fq1 fq2 ~k:(fun _name rd1 re1 rd2 re2 ->
        let r1 = r rd1 re1 in
        let r2 =
          match r1 with
          | Filtered _  -> r rd2 re2 (* Filtered so check both conf *)
          | Completed c -> single proc proc_to_reducer_stat (not c.comp) rd2 re2
        in
        add_log_likelihoods_andz ~state:t.state t.zygosity r1;
        add_log_likelihoods_andz ~state:t.state t.zygosity r2)
    and output oc =
      output_zi likelihood_first zygosity_report_size allele_arr t.state
        t.zygosity oc
    in
    t

  let apply fqi t =
    t.apply fqi

  let paired fq1 fq2 t =
    t.paired fq1 fq2

  let output oc t =
    t.output oc

end (* Reducer *)

(* A mapper maps the information from each read into a summarization of how
   a specific HLA-loci aligns the read. *)
module Mapper = struct

  type mp = mapped_complement_stats single_or_paired

  let pass_result_to_comparable_likelihood = function
    | Completed s -> fst (List.hd_exn s.positions)
    | Filtered _  -> neg_infinity

  let best_stat_rc ms =
    max (pass_result_to_comparable_likelihood ms.regular)
        (pass_result_to_comparable_likelihood ms.complement)

  let best_stat = function
    | Single ms         -> best_stat_rc ms
    | Paired (ms1, ms2) -> max (best_stat_rc ms1) (best_stat_rc ms2)

  let mp_to_string = function
    | Single ms         -> mapped_complement_stats_to_string ~sep:'\t' ms
    | Paired (ms1, ms2) -> sprintf "%s\t%s"
                            (mapped_complement_stats_to_string ~sep:'\t' ms1)
                            (mapped_complement_stats_to_string ~sep:'\t' ms2)
  type t =
    { mutable state : (string * mp) list
    ; apply         : Biocaml_unix.Fastq.item -> unit
    ; paired        : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> unit
    ; output        : out_channel -> unit
    }

  let to_proc ?allele ?insert_p ?max_number_mismatches ?band read_length pt =
    match allele with
    | None ->
        setup_single_pass ?insert_p ?max_number_mismatches ?band read_length pt
    | Some allele ->
        setup_single_allele_forward_pass ?insert_p ?max_number_mismatches
          read_length allele pt

  let init conf read_length n pt =
    let { allele; insert_p; band; max_number_mismatches; past_threshold_filter } = conf in
    let open ParPHMM in
    let m =
      let proc =
        to_proc ?allele ?insert_p ?max_number_mismatches ?band read_length pt
      in
      match allele with
      | None ->
          let pt = init_past_threshold past_threshold_filter in
          fun rd rp -> fst (mapper ~n pt proc rd rp)
      | Some _ -> (* Ignore since we do want to map *)
          fun rd rp -> fst (mapper ~n `Don't proc rd rp)
    in
    let rec t = { state = []; apply; paired; output }
    and apply =
      fqi_to_read_and_probs ~k:(fun name read read_errors ->
        let rc = m read read_errors in
        t.state <- (name, Single rc) :: t.state)
    and paired fq1 fq2 =
      fqi2_to_read_and_probs fq1 fq2 ~k:(fun name rd1 re1 rd2 re2 ->
        let rc1 = m rd1 re1 in
        let rc2 = m rd2 re2 in
        t.state <- (name, (Paired (rc1, rc2))) :: t.state)
    and output oc =
      fprintf oc "Reads: %d\n" (List.length t.state);
      List.sort t.state ~cmp:(fun (_n1, ms1) (_n2, ms2) ->
        compare (best_stat ms1) (best_stat ms2))
      |> List.rev
      |> List.iter ~f:(fun (n, mp) ->
          fprintf oc "%s\t%s\n" n (mp_to_string mp))
    in
    t

  let apply fqi t =
    t.apply fqi

  let paired fq1 fq2 t =
    t.paired fq1 fq2

  let output oc t =
    t.output oc

end (* Mapper *)

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
      sprintf "Reverse complement: %b, final emission: %f, path length: %d from %d to %d %s\n%s"
        rc emission (List.length path_list) start end_ from
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
          sprintf "Reverse complement: %b, final emission: %f, path length: %d from %d to %d %s->%s\n%s"
            rc1 emission (List.length path_list) s1.start s2.end_ f1 f2
            (manual_comp_display ~msm_offset ~width ~labels reference read)

  type t =
    { mutable state : (string * viterbi_result) list
    ; apply         : Biocaml_unix.Fastq.item -> unit
    ; paired        : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> unit
    ; output        : out_channel -> unit
    }

  let init conf read_length pt =
    let { allele; insert_p; _ } = conf in
    let open ParPHMM in
    (* Relying on reference being first. *)
    let allele = Option.value ~default:(fst pt.alleles.(0)) allele in
    let labels = allele ^ " ", "read " in
    let s, p =
      setup_single_allele_viterbi_pass ?insert_p ~allele
        read_length pt
    in
    let rec t = { state = []; apply; paired; output }
    and apply =
      fqi_to_read_and_probs ~k:(fun read_name read read_errors ->
        t.state <- (read_name, s read read_errors) :: t.state)
    and paired fq1 fq2 =
      fqi2_to_read_and_probs fq1 fq2 ~k:(fun read_name rd1 re1 rd2 re2 ->
        t.state <- (read_name, p rd1 re1 rd2 re2) :: t.state)
    and output oc =
      List.iter t.state ~f:(fun (read_name, vr ) ->
        fprintf oc "%s:\n%s\n" read_name (viterbi_result_to_string ~labels vr))
    in
  t

  let apply fqi t =
    t.apply fqi

  let paired fq1 fq2 t =
    t.paired fq1 fq2

  let output oc t =
    t.output oc

end (* Viterbi *)

(* Reduce and map: Adding this to simplify matters. *)
module Rammer = struct

  open ParPHMM

  type per_read =
    { name  : string
    ; runr  : alleles_and_positions rc_states
    }

  type t =
    { per_allele_lhood  : float array
    ; zygosity          : Zygosity_likelihood_array.t
    ; mutable per_reads : per_read list
    ; apply             : Biocaml_unix.Fastq.item -> unit
    ; paired            : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> unit
    ; output            : out_channel -> unit
    }

  let just_aap = function
    | Filtered m   -> Filtered m
    | Completed rm -> Completed rm.aap

  let rc_just_aap rc =
    { regular    = just_aap rc.regular
    ; complement = just_aap rc.complement
    }

  let update_state t name rcs =
    let allars = Reducer.add_log_likelihoods_from_llhd in
    t.per_reads <- {name; runr = rc_just_aap rcs} :: t.per_reads;
    match rcs.regular, rcs.complement with
    | Filtered rm, Filtered cm  ->
        invalid_argf "Both were filtered!: %s %s" rm cm
    | Filtered _m, Completed c  ->
        allars ~state:t.per_allele_lhood t.zygosity c.rlikelihood
    | Completed r, Filtered _   ->
        allars ~state:t.per_allele_lhood t.zygosity r.rlikelihood
    | Completed r, Completed c  ->
        if compare_emissions r.rlikelihood c.rlikelihood then
          allars ~state:t.per_allele_lhood t.zygosity r.rlikelihood
        else
          allars ~state:t.per_allele_lhood t.zygosity c.rlikelihood

  let output_state likelihood_first zygosity_report_size allele_arr t oc =
    Reducer.output_zi likelihood_first zygosity_report_size allele_arr
      t.per_allele_lhood t.zygosity oc;
    fprintf oc "Per Read Positions:\n";
    List.iter t.per_reads ~f:(fun {name; runr} ->
      fprintf oc "%s\n\t%s\n\t%s\n"
        name
        (fc_aap_to_string false runr.regular)
        (fc_aap_to_string true runr.complement))

  let init conf likelihood_first zygosity_report_size ~report_size read_length pt =
    let { allele; insert_p; band; max_number_mismatches; past_threshold_filter
        ; _} = conf
    in
    let open ParPHMM in
    let proc, allele_arr =
      Reducer.to_proc_allele_arr ?allele ?insert_p ?max_number_mismatches ?band
        read_length pt
    in
    let past_threshold = init_past_threshold past_threshold_filter in
    (* Discard resulting past threshold *)
    let r rd rp = fst (rammer report_size past_threshold proc rd rp) in
    let rec t = { per_allele_lhood; zygosity; per_reads; apply; paired; output }
    and per_allele_lhood = proc.init_global_state ()
    and zygosity = Zygosity_likelihood_array.init (Array.length allele_arr)
    and per_reads = []
    and apply fqi =
      fqi_to_read_and_probs fqi ~k:(fun name read read_errors ->
        r read read_errors
        |> update_state t name)
    and paired fq1 fq2 =
      fqi2_to_read_and_probs fq1 fq2 ~k:(fun name rd1 re1 rd2 re2 ->
        let r1 = r rd1 re1 in
        let r2 = r rd2 re2 in               (* A bit wasteful. *)
        update_state t name r1;
        update_state t name r2)
    and output oc =
      output_state likelihood_first zygosity_report_size allele_arr t oc
    in
    t

  let apply fqi t =
    t.apply fqi

  let paired fq1 fq2 t =
    t.paired fq1 fq2

  let output oc t =
    t.output oc

end (* Rammer *)

(** Single Loci Drivers **)
module Single_loci = struct

  let conf ?allele ?insert_p ?band ?max_number_mismatches
    ~past_threshold_filter ~check_rc () =
    { allele; insert_p; band; max_number_mismatches
    ; past_threshold_filter ; check_rc
    }

  type t =
    | Reducer of Reducer.t
    | Mapper of Mapper.t
    | Viterbi of Viterbi.t
    | Rammer of Rammer.t

  let init conf read_length pt mode =
      match mode with
      | `Reducer (lfst, zrs)    ->
          Reducer (Reducer.init conf lfst zrs read_length pt)
      | `Mapper n               ->
          Mapper (Mapper.init conf read_length n pt)
      | `Viterbi                ->
          Viterbi (Viterbi.init conf read_length pt)
      | `Rammer (lfst, zrs, rs) ->
          Rammer (Rammer.init conf lfst zrs ~report_size:rs read_length pt)

  (* fqi = FastQ Item *)
  let apply t fqi = match t with
    | Reducer r -> Reducer.apply fqi r
    | Mapper m  -> Mapper.apply fqi m
    | Viterbi v -> Viterbi.apply fqi v
    | Rammer v  -> Rammer.apply fqi v

  let paired t fq1 fq2 = match t with
    | Reducer r -> Reducer.paired fq1 fq2 r
    | Mapper m  -> Mapper.paired fq1 fq2 m
    | Viterbi v -> Viterbi.paired fq1 fq2 v
    | Rammer v  -> Rammer.paired fq1 fq2 v

  let output t oc = match t with
    | Reducer r -> Reducer.output oc r
    | Mapper m  -> Mapper.output oc m
    | Viterbi v -> Viterbi.output oc v
    | Rammer v  -> Rammer.output oc v

end (* Single_loci *)

(* Combine results from multiple loci. *)

type multiple_conf =
  { insert_p              : float option
  (* Override the default insert emission probability. *)

  ; band                  : ParPHMM.band_config option
  (* Perform banded passes. *)

  ; max_number_mismatches : int option
  (* Max number of allowable mismatches to use a threshold filter for the
     forward pass. *)

  ; past_threshold_filter : bool
  (* Use the previous match likelihood, when available (ex. against reverse
     complement), as a threshold filter for the forward pass. *)

  ; incremental_pairs : bool
  (* This option only applies to paired typing. Instead of naively modeling
     the two read pairs at the same time, use the first as guidance of which
     loci is the best, and then apply the second read to only the best loci.
     The default should be true. *)
  }

(* Reduce to multiple loci. *)
module Multiple_reducer = struct

  type best_state =
    { name  : string
    ; rspr  : reducer_stats pass_result single_or_paired
    ; maxl  : float
    }

  let pr_to_llhd = function
    | Filtered _  -> neg_infinity
    | Completed c -> max_arr c.likelihood

  let sp_to_llhd = function
    | Single s        -> pr_to_llhd s
    | Paired (p1, p2) -> max (pr_to_llhd p1) (pr_to_llhd p2)

  let to_bs (name, rspr) =
    { name; rspr; maxl = sp_to_llhd rspr}

  let best_bs bs1 bs2 =
    (*printf "%s %f vs %s %f\n%!" bs1.name bs1.maxl bs2.name bs2.maxl; *)
    if bs1.maxl >= bs2.maxl then bs1 else bs2

  let best = function
    | []          -> invalid_argf "Can't select best from empty!"
    | h :: t ->
        let init = to_bs h in
        List.fold_left t ~init ~f:(fun bs p ->
          best_bs bs (to_bs p))

  type per_loci =
    { loci           : string
    ; likelihood     : float array
    ; zygosity       : Zygosity_likelihood_array.t
    }

  let add_log_likelihoods  { likelihood; zygosity; _} = function
    | Single r        ->
        Reducer.add_log_likelihoods_andz ~state:likelihood zygosity r
    | Paired (r1, r2) ->
        Reducer.add_log_likelihoods_andz ~state:likelihood zygosity r1;
        Reducer.add_log_likelihoods_andz ~state:likelihood zygosity r2

  let merge state current =
    let b = best current in
    List.iter state ~f:(fun pl ->
      if pl.loci = b.name then add_log_likelihoods pl b.rspr)

  let to_best_loci ~init ~f ~c = function
    | [] -> invalid_argf "Empty loci list"
    | (name, (proc, _, _)) :: tl ->
        let acc, m = f init proc in
        let rec loop acc bl best = function
          | []  -> best
          | (name, (proc, _, _)) :: t ->
              let nacc, m = f acc proc in
              let l = c m in
              if l > bl then
                loop nacc l (name, proc, m) t
              else
                loop nacc bl best t
        in
        loop acc (c m) (name, proc, m) tl

  type t =
    { state   : per_loci list
    ; apply   : Biocaml_unix.Fastq.item -> unit
    ; paired  : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> unit
    ; output  : out_channel -> unit
    }

  let init conf likelihood_first zygosity_report_size read_length pt_lst =
    let { insert_p; band; max_number_mismatches; _ } = conf in
    let open ParPHMM in
    let paa =
      List.map_snd pt_lst ~f:(fun pt ->
        let proc, allele_arr =
          Reducer.to_proc_allele_arr ?insert_p ?max_number_mismatches ?band
            read_length pt
        in
        let n = Array.length allele_arr in
        proc, allele_arr, n)
    in
    let pt = init_past_threshold conf.past_threshold_filter in
    let rec t = { state; apply; paired; output }
    and state =
      List.map paa ~f:(fun (loci, (p, _, number_alleles)) ->
        { loci
        ; likelihood = p.init_global_state ()
        ; zygosity = Zygosity_likelihood_array.init number_alleles
        })
    and output oc =
      List.iter2 paa t.state
        ~f:(fun (_, (_, allele_arr, _)) { loci; likelihood; zygosity; _} ->
              fprintf oc "%s\n" loci;
              Reducer.output_zi likelihood_first zygosity_report_size
                allele_arr likelihood zygosity oc)
    and apply fqi =
      match pt with
      | `Don't ->
          fqi_to_read_and_probs fqi ~k:(fun _name read read_errors ->
            List.map_snd paa ~f:(fun (p, _, _) -> Single (fst (reducer pt p read read_errors)))
            |> merge t.state)
      | `Start ->
          fqi_to_read_and_probs fqi ~k:(fun _name read read_errors ->
          List.fold_left paa ~init:(pt, [])
            ~f:(fun (pt, acc) (name, (p, _, _)) ->
                  let m, npto = reducer pt p read read_errors in
                  let npt = update_pt npto in
                  (npt, (name, Single m) :: acc))
          |> function _, acc ->
                merge t.state acc)
    and paired fq1 fq2 =
      let add_to_state best_loci m1 m2 =
        List.iter t.state ~f:(fun {loci; likelihood; zygosity; _ } ->
          if loci = best_loci then begin
            Reducer.add_log_likelihoods_andz likelihood zygosity m1;
            Reducer.add_log_likelihoods_andz likelihood zygosity m2
         end)
      in
      match pt with
      | `Don't ->
          fqi2_to_read_and_probs fq1 fq2 ~k:(fun _name rd1 re1 rd2 re2 ->
            if conf.incremental_pairs then
              let best_loci, best_proc, m1 =
                to_best_loci paa ~init:()
                  ~f:(fun () p -> (), fst (reducer `Don't p rd1 re1))
                  ~c:pr_to_llhd
              in
              let m2 =
                match m1 with
                | Filtered  _ -> fst (reducer `Don't best_proc rd2 re2)
                | Completed c -> single best_proc proc_to_reducer_stat (not c.comp) rd2 re2
              in
              add_to_state best_loci m1 m2
            else
              List.map_snd paa ~f:(fun (p, _, _) ->
                let r1 = fst (reducer `Don't p rd1 re1) in
                let r2 =
                  match r1 with
                  | Filtered _  -> fst (reducer `Don't p rd2 re2) (* Filtered so check both conf *)
                  | Completed c -> single p proc_to_reducer_stat (not c.comp) rd2 re2
                in
                Paired (r1, r2))
              |> merge t.state)
      | `Start ->
          fqi2_to_read_and_probs fq1 fq2 ~k:(fun _name rd1 re1 rd2 re2 ->
            if conf.incremental_pairs then
              let best_loci, best_proc, m1 =
                to_best_loci paa ~init:`Start
                  ~f:(fun pt proc ->
                        let m1, npto = reducer pt proc rd1 re1 in
                        update_pt npto, m1)
                  ~c:pr_to_llhd
              in
              let m2 =
                match m1 with
                | Filtered  _ -> fst (reducer `Don't best_proc rd2 re2)
                | Completed c -> single best_proc proc_to_reducer_stat (not c.comp) rd2 re2
              in
              (*let m2 = fst (reducer `Don't best_proc rd2 re2) in *)
              add_to_state best_loci m1 m2
            else
              List.fold_left paa ~init:(`Start, `Start, [])
                ~f:(fun (pt1, pt2, acc) (name, (p, _, _)) ->
                      let m1, npto1 = reducer pt1 p rd1 re1 in
                      let npt1 = update_pt npto1 in
                      let m2, npto2 = reducer pt2 p rd2 re2 in
                      let npt2 = update_pt npto2 in
                      let nacc = (name, Paired (m1, m2)) :: acc in
                      (npt1, npt2, nacc))
              |> function _, _, nacc -> merge t.state nacc)
    in
    t

  let apply fqi t =
    t.apply fqi

  let paired fq1 fq2 t =
    t.paired fq1 fq2

  let output oc t =
    t.output oc

end (* Multiple_reducer *)

module Multiple_mapper = struct

  type t =
    { mutable state : (string * (string * Mapper.mp) list) list
    ; apply         : Biocaml_unix.Fastq.item -> unit
    ; paired        : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> unit
    ; output        : out_channel -> unit
    }

  let sort_mapped_output lst =
    let highest_llhd l =
      List.map l ~f:(fun (_name, mp) -> Mapper.best_stat mp)
      |> List.reduce ~f:max
    in
    List.sort lst ~cmp:(fun (_rn1, nlst1) (_rn2, nlst2) ->
      let sl1 = highest_llhd nlst1 in
      let sl2 = highest_llhd nlst2 in
      compare sl2 sl1 (* higher is better *))

  let init conf read_length n pt_lst =
    let { insert_p; band; max_number_mismatches; past_threshold_filter } = conf in
    let open ParPHMM in
    let procs =
      List.map_snd pt_lst ~f:(fun pt ->
        Mapper.to_proc ?insert_p ?max_number_mismatches ?band read_length pt)
    in
    let pt = init_past_threshold past_threshold_filter in
    let rec t = { state; apply; paired; output }
    and state = []
    and output oc =
      fprintf oc "Reads: %d\n" (List.length t.state);
      sort_mapped_output t.state
      |> List.iter ~f:(fun (read, lst) ->
        fprintf oc "%s\n\t%s\n"
          read (String.concat ~sep:"\n\t"
            (List.map (List.sort ~cmp:compare lst) (* Sort by name too *)
              ~f:(fun (name, mp) -> sprintf "%s\t%s"
                      name (Mapper.mp_to_string mp)))))
    and apply fqi =
      match pt with
      | `Don't  ->
          fqi_to_read_and_probs fqi ~k:(fun read_name read read_errors ->
            let r = List.map_snd procs ~f:(fun p ->
                      Single (fst (mapper ~n `Don't p read read_errors)))
            in
            t.state <- (read_name, r) :: t.state)
      | `Start  ->
          (* Now we have to take into account previous filter results. *)
          fqi_to_read_and_probs fqi ~k:(fun read_name read read_errors ->
            List.fold_left procs ~init:(`Start, [])
              ~f:(fun (pt, acc) (name, p) ->
                    let m, npto = mapper ~n pt p read read_errors in
                    let npt = update_pt npto in
                    (npt, (name, Single m) :: acc))
            |> function _, acc ->
                  t.state <- (read_name, acc) :: t.state)
    and paired fq1 fq2 =
      match pt with
      | `Don't ->
          fqi2_to_read_and_probs fq1 fq2 ~k:(fun read_name rd1 re1 rd2 re2 ->
            let r = List.map_snd procs ~f:(fun p ->
              let m1, _ = mapper `Don't ~n p rd1 re1 in
              let m2, _ = mapper `Don't ~n p rd2 re2 in
              Paired (m1, m2))
            in
            t.state <- (read_name, r) :: t.state)
      | `Start ->
          (* Now we have to take into account previous filter results. *)
          fqi2_to_read_and_probs fq1 fq2 ~k:(fun read_name rd1 re1 rd2 re2 ->
            List.fold_left procs ~init:(`Start, `Start, [])
              ~f:(fun (pt1, pt2, acc) (name, p) ->
                    let m1, npto1 = mapper ~n pt1 p rd1 re1 in
                    let npt1 = update_pt npto1 in
                    let m2, npto2 = mapper ~n pt2 p rd2 re2 in
                    let npt2 = update_pt npto2 in
                    let nacc = (name, Paired (m1, m2)) :: acc in
                    (npt1, npt2, nacc))
            |> function _, _, acc ->
                  t.state <- (read_name, acc) :: t.state)
    in
    t

  let apply fqi t =
    t.apply fqi

  let paired fq1 fq2 t =
    t.paired fq1 fq2

  let output oc t =
    t.output oc

end (* Multiple_mapper *)

module Multiple_rammer = struct

  type read_result =
    | SingleRead of alleles_and_positions rc_states
    (* We'll test both regular and complement options. *)

    | FirstFiltered of
      { first   : string
      ; second  : alleles_and_positions rc_states
      }
    (* Both orientations filered *)

    | FirstOrientedSecond of
      { first   : alleles_and_positions rc_states
      ; second  : alleles_and_positions pass_result
      }
    (* At least one orientation of the first read isn't Filtered. *)

    | PairedDependent of alleles_and_positions rc_states
    (* Orientation and loci of 2nd is determined by 1st. *)

  let rrpr_to_likelihood_opt = function
    | Filtered  _ -> None
    | Completed r -> Some r.rlikelihood

  let rrcs_to_likelihood_opt rrcs =
    match rrcs.regular, rrcs.complement with
    | Filtered mr, Filtered mc -> None
    | Filtered _r, Completed c -> Some c.rlikelihood
    | Completed r, Filtered _c -> Some r.rlikelihood
    | Completed r, Completed c ->
      if compare_emissions r.rlikelihood c.rlikelihood then
        Some r.rlikelihood
      else
        Some c.rlikelihood

  let aap_rc_to_string ?sep aaprc =
    rc_stats_to_string ?sep
      (fun al1 al2 ->
        let e1, _, _ = List.hd_exn al1 in
        let e2, _, _ = List.hd_exn al2 in
        e1 > e2)
      aap_lst_to_string aaprc

  let read_result_to_string = function
    | SingleRead aaprc                      ->
        sprintf "SR\t%s" (aap_rc_to_string aaprc)
    | FirstFiltered { first; second}        ->
        sprintf "FF\t%s\t%s" first (aap_rc_to_string second)
    | FirstOrientedSecond { first; second } ->
        sprintf "FO\t%s\t%s" (aap_rc_to_string first)
          (pass_result_to_string aap_lst_to_string second)
    | PairedDependent aaprc ->
        sprintf "PD\t%s" (aap_rc_to_string aaprc)

  let single_filtered_later comp s =
    if comp then
      { regular = Filtered "Not best loci"
      ; complement = Completed s
      }
    else
      { regular = Completed s
      ; complement = Filtered "Not best loci"
      }

  type per_read =
    { name  : string
    ; rrs   : (string * read_result) list
    }

  let output_per_reads prlst oc =
    List.iter prlst ~f:(fun {name ; rrs} ->
      fprintf oc "%s\n" name;
      List.iter rrs ~f:(fun (loci, rr) ->
        fprintf oc "\t%s\t%s\n" loci (read_result_to_string rr)))

  type per_loci =
    { loci        : string
    ; likelihood  : float array
    ; zygosity    : Zygosity_likelihood_array.t
    }

  type t =
    { per_loci          : per_loci list
    ; mutable per_read  : per_read list
    ; apply             : Biocaml_unix.Fastq.item -> unit
    ; paired            : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> unit
    ; output            : out_channel -> unit
    }

  (* I'm doing an aweful lot of tagging to delay adding likelihoods. *)
  let max_opt_sp = function
    | Single None               -> None
    | Single (Some la)          -> Some (max_arr la, Single la)
    | Paired (None,    None)    -> None
    | Paired (None,    Some lb) -> Some (max_arr lb, Single lb)
    | Paired (Some la, None)    -> Some (max_arr la, Single la)
    | Paired (Some la, Some lb) -> Some (max (max_arr la) (max_arr lb)
                                        , Paired (la, lb))

  let merge t = function
    | []                 -> invalid_argf "Empty loci list!"
    | likelihood_opt_lst ->
        List.fold_left likelihood_opt_lst ~init:None
          ~f:(fun best_opt (loci, lldho) ->
                match max_opt_sp lldho with
                | None            -> best_opt
                | Some (nla, nsp) ->
                    match best_opt with
                    | None -> Some (loci, nsp, nla)
                    | Some (_current_best_loci, _current_best_sp, best_l) ->
                            if nla > best_l then
                              Some (loci, nsp, nla)
                            else
                              best_opt)
        |> function
            | None  -> eprintf "All loci filtered out\n"
            | Some (best_loci, best_sp_llhd, _bl) ->
                List.iter t.per_loci ~f:(fun {loci; likelihood; zygosity} ->
                  if loci = best_loci then begin
                    match best_sp_llhd with
                    | Single best_llhd ->
                        Reducer.add_log_likelihoods_from_llhd
                            ~state:likelihood zygosity best_llhd
                    | Paired (b1, b2)  ->
                        Reducer.add_log_likelihoods_from_llhd
                          ~state:likelihood zygosity b1;
                        Reducer.add_log_likelihoods_from_llhd
                          ~state:likelihood zygosity b2
                  end)

  let informative comp rc =
    match rc.regular, rc.complement with
    | Filtered mr, Filtered mc -> `FilteredBoth (sprintf "Both! %s, %s" mr mc)
    | Filtered _r, Completed c -> `OneGood (true, c)
    | Completed r, Filtered _c -> `OneGood (false, r)
    | Completed r, Completed c ->
        if comp r c then
          `OneGood (false, r)
        else
          `OneGood (true, c)

  let informative_rammer_res =
    informative (fun r c -> compare_emissions r.rlikelihood c.rlikelihood)

  let init conf likelihood_first zygosity_report_size ~report_size read_length
    pt_lst =
    let { insert_p; band; max_number_mismatches; _ } = conf in
    let open ParPHMM in
    let module Mr = Multiple_reducer in
    let paa =
      List.map_snd pt_lst ~f:(fun pt ->
        let proc, allele_arr =
          Reducer.to_proc_allele_arr ?insert_p ?max_number_mismatches ?band
            read_length pt
        in
        let n = Array.length allele_arr in
        proc, allele_arr, n)
    in
    let pt = init_past_threshold conf.past_threshold_filter in
    let rammer = rammer report_size in
    let single p = single p (proc_to_rammer_stat report_size) in
    let split_rammer pt p read read_errors =
      let rammer_rc_stats, npt = rammer pt p read read_errors in
      let rlikelihood_opt = rrcs_to_likelihood_opt rammer_rc_stats in
      let read_result = Rammer.rc_just_aap rammer_rc_stats in
      npt, rlikelihood_opt, read_result
    in
    let rec t = { per_loci; per_read; apply; paired; output }
    and per_loci =
      List.map paa ~f:(fun (loci, (p, _, number_alleles)) ->
        { loci
        ; likelihood = p.init_global_state ()
        ; zygosity = Zygosity_likelihood_array.init number_alleles
        })
    and per_read = []
    and output oc =
      List.iter2 paa t.per_loci
        ~f:(fun (_, (_, allele_arr, _)) { loci; likelihood; zygosity; _} ->
              fprintf oc "%s\n" loci;
              Reducer.output_zi likelihood_first zygosity_report_size
                allele_arr likelihood zygosity oc);
      output_per_reads t.per_read oc
    and apply fqi =
      match pt with
      | `Don't ->
          fqi_to_read_and_probs fqi ~k:(fun name read read_errors ->
            let rlikelihoods, rrs =
              List.map paa ~f:(fun (loci, (p, _, _)) ->
                let _npt, rlikelihood_opt, read_result =
                  split_rammer pt p read read_errors
                in
                (loci, Single rlikelihood_opt)
                , (loci, SingleRead read_result))
              |> List.split
            in
            t.per_read <- { name; rrs } :: t.per_read;
            merge t rlikelihoods)
      | `Start ->
          fqi_to_read_and_probs fqi ~k:(fun name read read_errors ->
            let _final_pt, rlikelihoods, rrslst =
              List.fold_left paa ~init:(pt, [], [])
                ~f:(fun (pt, acc1, acc2) (loci, (p, _, _)) ->
                    let npto, rlikelihood_opt, read_result =
                      split_rammer pt p read read_errors
                    in
                    let npt = update_pt npto in
                    (npt
                    , (loci, Single rlikelihood_opt) :: acc1
                    , (loci, SingleRead read_result) :: acc2))
            in
            t.per_read <- {name; rrs = List.rev rrslst} :: t.per_read;
            merge t rlikelihoods)
    and paired fq1 fq2 =
      fqi2_to_read_and_probs fq1 fq2 ~k:(fun name rd1 re1 rd2 re2 ->
        match pt with
        | `Don't ->                                           (* No previous threshold. *)
            if conf.incremental_pairs then begin
              let rlikelihoods, rrs, best_opt =
                List.fold_left paa ~init:([], [], None)
                  ~f:(fun (acc1, acc2, best) (loci, (p, _, _)) ->
                        let rammer_rc_stats1, _npt = rammer `Don't p rd1 re1 in
                        match informative_rammer_res rammer_rc_stats1 with
                        | `FilteredBoth _  ->
                            (loci, Single None) :: acc1
                            , (loci, SingleRead (Rammer.rc_just_aap rammer_rc_stats1)) :: acc2
                            , best
                        | `OneGood (comp, rammer_rc) ->
                            let l = max_arr rammer_rc.rlikelihood in
                            begin match best with
                            | None  ->
                                let nbest = Some (loci, p, comp, rammer_rc_stats1, l) in
                                acc1, acc2, nbest
                            | Some (bloci, bp, bcomp, brammer_rc_stat, bl) ->
                                if l > bl then    (* replace *)
                                  let nbest = Some (loci, p, comp, rammer_rc_stats1, l) in
                                  (loci, Single (rrcs_to_likelihood_opt brammer_rc_stat)) :: acc1
                                  , (loci, SingleRead (Rammer.rc_just_aap brammer_rc_stat)) :: acc2
                                  , nbest
                                else
                                  (loci, Single (Some rammer_rc.rlikelihood)) :: acc1
                                  , (loci, SingleRead (Rammer.rc_just_aap rammer_rc_stats1)) :: acc2
                                  , best
                            end)
              in
              let rlikelihoods, rrs =
                match best_opt with
                | None  ->
                    rlikelihoods, rrs (* Possible if all filtered. *)
                | Some (bloci, bp, bcomp, brammer_rc_stat, bl) ->
                    let second = single bp (not bcomp) rd2 re2 in
                    let combined_likelihood =
                      Paired (rrcs_to_likelihood_opt brammer_rc_stat
                             , rrpr_to_likelihood_opt second)
                    in
                    (bloci, combined_likelihood) :: rlikelihoods
                   , (bloci, FirstOrientedSecond
                              { first = Rammer.rc_just_aap brammer_rc_stat
                              ; second = Rammer.just_aap second }) :: rrs
              in
              t.per_read <- { name; rrs } :: t.per_read;
              merge t rlikelihoods
            end else begin
              let rlikelihoods, rrs =
                List.map paa ~f:(fun (loci, (p, _, _)) ->
                  let rammer_rc_stats1, _npt = rammer `Don't p rd1 re1 in
                  match informative_rammer_res rammer_rc_stats1 with
                  | `FilteredBoth first  ->
                      let _npt, rlikelihood_opt2, second = split_rammer `Don't p rd2 re2 in
                      (loci, Single rlikelihood_opt2)
                      , (loci, FirstFiltered {first; second})
                  | `OneGood (comp, rammer_rc) ->
                      let second = single p (not comp) rd2 re2 in
                      let combined_likelihood =
                        Paired (Some rammer_rc.rlikelihood
                               , rrpr_to_likelihood_opt second)
                      in
                      (loci, combined_likelihood)
                      , (loci, FirstOrientedSecond
                                { first = Rammer.rc_just_aap rammer_rc_stats1
                                ; second = Rammer.just_aap second }))
                |> List.split
              in
              t.per_read <- { name; rrs } :: t.per_read;
              merge t rlikelihoods
            end
      | `Start ->
          if conf.incremental_pairs then begin
            let npt, rlikelihoods, rrs, best_opt =
              List.fold_left paa ~init:(`Start, [], [], None)
                ~f:(fun (pt, acc1, acc2, best) (loci, (p, _, _)) ->
                      let rammer_rc_stats1, npt = rammer pt p rd1 re1 in
                      let npto = update_pt npt in
                      match informative_rammer_res rammer_rc_stats1 with
                      | `FilteredBoth _  ->
                          npto
                          , (loci, Single None) :: acc1
                          , (loci, SingleRead (Rammer.rc_just_aap rammer_rc_stats1)) :: acc2
                          , best
                      | `OneGood (comp, rammer_rc) ->
                          let l = max_arr rammer_rc.rlikelihood in
                          begin match best with
                          | None  ->
                              let nbest = Some (loci, p, comp, rammer_rc_stats1, l) in
                              npto, acc1, acc2, nbest
                          | Some (bloci, bp, bcomp, brammer_rc_stat, bl) ->
                              if l > bl then    (* replace *)
                                let nbest = Some (loci, p, comp, rammer_rc_stats1, l) in
                                npto
                                , (loci, Single (rrcs_to_likelihood_opt brammer_rc_stat)) :: acc1
                                , (loci, SingleRead (Rammer.rc_just_aap brammer_rc_stat)) :: acc2
                                , nbest
                              else
                                npto
                                , (loci, Single (Some rammer_rc.rlikelihood)) :: acc1
                                , (loci, SingleRead (Rammer.rc_just_aap rammer_rc_stats1)) :: acc2
                                , best
                          end)
            in
            let rlikelihoods, rrs =
              match best_opt with
              | None  ->
                  rlikelihoods, rrs (* Possible if all filtered. *)
              | Some (bloci, bp, bcomp, brammer_rc_stat, bl) ->
                  let second = single bp (not bcomp) rd2 re2 in
                  let combined_likelihood =
                    Paired (rrcs_to_likelihood_opt brammer_rc_stat
                            , rrpr_to_likelihood_opt second)
                  in
                  (bloci, combined_likelihood) :: rlikelihoods
                  , (bloci, FirstOrientedSecond
                            { first = Rammer.rc_just_aap brammer_rc_stat
                            ; second = Rammer.just_aap second }) :: rrs
            in
            t.per_read <- { name; rrs } :: t.per_read;
            merge t rlikelihoods
          end else begin
            failwith "Not Implemented: Do you REALLY want filtered non-increment pairs?"
          end)
    in
    t

  let apply fqi t =
    t.apply fqi

  let paired fq1 fq2 t =
    t.paired fq1 fq2

  let output oc t =
    t.output oc


end (* Multiple_rammer *)

module Mulitple_loci = struct

  let conf ?insert_p ?band ?max_number_mismatches ~past_threshold_filter
    ~incremental_pairs () =
    { insert_p; band; max_number_mismatches; past_threshold_filter
    ; incremental_pairs }

  type t =
    | Reducer of Multiple_reducer.t
    | Mapper of Multiple_mapper.t
    | Rammer of Multiple_rammer.t

  let init conf read_length pt mode =
      match mode with
      | `Reducer (lfst, zrs) ->
          Reducer (Multiple_reducer.init conf lfst zrs read_length pt)
      | `Mapper n ->
          Mapper (Multiple_mapper.init conf read_length n pt)
      | `Rammer (lfst, zrs, rs) ->
          Rammer (Multiple_rammer.init conf lfst zrs ~report_size:rs
                    read_length pt)

  (* fqi = FastQ Item *)
  let apply t fqi = match t with
    | Reducer r -> Multiple_reducer.apply fqi r
    | Mapper m  -> Multiple_mapper.apply fqi m
    | Rammer m  -> Multiple_rammer.apply fqi m

  let paired t fq1 fq2 = match t with
    | Reducer r -> Multiple_reducer.paired fq1 fq2 r
    | Mapper m  -> Multiple_mapper.paired fq1 fq2 m
    | Rammer m  -> Multiple_rammer.paired fq1 fq2 m

  let output t oc = match t with
    | Reducer r -> Multiple_reducer.output oc r
    | Mapper m  -> Multiple_mapper.output oc m
    | Rammer m  -> Multiple_rammer.output oc m

end (* Mulitple_loci *)
