(* Aggregation of the forward passes. *)

open Util

module Pm = Partition_map

let max_arr = Array.fold_left ~init:neg_infinity ~f:max
let max_pm = Pm.fold_values ~init:neg_infinity ~f:max

let higher_value_in_first e1 e2 =
  let r1 = max_pm e1 in
  let r2 = max_pm e2 in
  r1 >= r2

let descending_float_cmp f1 f2 =
  if f1 > f2 then -1
  else if f1 = f2 then 0
  else 1

(* What to do when we know the result of a previous run. *)
module Past_threshold = struct

  open ParPHMM

  type t =
    [ `Don't
    | `Start
    | `Set of float
    ]

  let to_string = function
    | `Don't  -> "Don't"
    | `Start  -> "Start"
    | `Set v  -> sprintf "Set %f" v

  let init b =
    if b then `Start else `Don't

  let new_ prev_threshold proc =
    max prev_threshold (proc.maximum_match ())

  let new_value prev_threshold proc = function
    | Filtered _  -> prev_threshold
    | Completed _ -> new_ prev_threshold proc

  let update prev_threshold proc pr =
    match prev_threshold with
    | `Don't  -> `Don't
    | `Start  -> `Set (new_value neg_infinity proc pr)
    | `Set pv -> `Set (new_value pv proc pr)

end  (* Past_threshold *)

(* Consider a single orientation: either regular or reverse complement. *)
let single ?base_p ?prev_threshold proc pass_result_map reverse_complement read read_errors =
  let open ParPHMM in
  match proc.single ?base_p ?prev_threshold ~reverse_complement ~read ~read_errors with
  | Filtered m    -> Filtered m
  | Completed ()  -> Completed (pass_result_map proc reverse_complement)

let most_likely_between_two_pr ~take_regular pr1 pr2 =
  let open ParPHMM in
  match pr1, pr2 with
  | Filtered mr, Filtered mc -> Filtered (sprintf "Both! %s, %s" mr mc)
  | Filtered _r, Completed c -> Completed (`Snd, c)
  | Completed r, Filtered _c -> Completed (`Fst, r)
  | Completed r, Completed c ->
    if take_regular r c then
      Completed (`Fst, r)
    else
      Completed (`Snd, c)

(* Sometimes (such as in the beginning) we do not know whether to prefer a
   regular or a reverse-complement orientation for read. *)
module Orientation = struct

  open ParPHMM

  type 'a t =
    { regular     : 'a pass_result
    ; complement  : 'a pass_result
    }

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
    most_likely_between_two_pr ~take_regular t.regular t.complement
    |> function
        | Filtered m          -> Filtered m
        | Completed (`Snd, c) -> Completed (true, c)
        | Completed (`Fst, r) -> Completed (false, r)

  (* Consider both orientations. *)
  let paired pt pass_result_map proc read read_errors =
    let do_work ?prev_threshold rc =
      single ?prev_threshold proc pass_result_map rc read read_errors
    in
    match pt with
    | `Don't ->
        let regular = do_work false in
        let complement = do_work true in
        { regular; complement }, `Don't
    | `Start ->
        let regular = do_work false in
        let prev_threshold = Past_threshold.new_value neg_infinity proc regular in
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

  (* Likelyhood emission, allele, position of highest emission. *)
  type t = (float * string * int) list

  let to_string (e,a,p) = sprintf "%s,%d,%0.5f" a p e

  let lst_to_string al =
    String.concat ~sep:";" (List.map al ~f:to_string)

  (* Descending *)
  let descending_cmp l1 l2 =
    match l1, l2 with
    | [],  _                           -> -1  (* Should error? *)
    | _ , []                           -> 1
    | (e1, _, _) :: _, (e2, _, _) :: _ -> descending_float_cmp e1 e2

end (* Alleles_and_positions *)

(* Storing the log likelihoods of heterozygous pairs.
   Homozygous is equivalent to the individual pairs. *)
module Zygosity_array = struct

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

  let best ?size t =
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

end (* Zygosity_array *)

module Likelihoods_and_zygosity = struct

  (* TODO: We're opening Util which has likehood defined. *)

  let add_ll state_llhd llhd =
    Pm.iter_set llhd ~f:(fun i v -> state_llhd.(i) <- state_llhd.(i) +. v)

  let add_lz zygosity llhd =
    Pm.iter_set llhd ~f:(fun allele1 v1 ->
      Pm.iter_set llhd ~f:(fun allele2 v2 ->
        if allele1 = allele2 then
          ()
        else
          Zygosity_array.add zygosity ~allele1 ~allele2 (max v1 v2)))

  (* State is labled because it is modified. *)
  let add_ll_and_lz state_llhd zygosity llhd =
    add_ll state_llhd llhd;
    add_lz zygosity llhd

end (* Likelihoods_and_zygosity *)

(* Outputing results *)
module Output = struct

  let compare2 (l1, _) (l2, _) =
    descending_float_cmp l1 l2

  let compare3 (l1, _, _) (l2, _, _) =
    descending_float_cmp l1 l2

  let compare4 (l1, r1, _i1, _a1) (l2, r2, _i2, _a2) =
    let lc = descending_float_cmp l1 l2 in
    if lc <> 0 then
      lc
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
    Array.sort o ~cmp:compare4;
    o

  let allele_first o oc =
    Array.iter o ~f:(fun (l,_,_,a) -> fprintf oc "%16s\t%0.20f\n" a l)

  let likelihood_first allele_arr final_likelihoods oc =
    let o =
      Array.mapi allele_arr ~f:(fun i a -> final_likelihoods.(i), a)
      |> Array.to_list
      |> group_by_assoc
      |> List.map ~f:(fun (l, alst) ->
            let clst = Alleles.CompressNames.f alst in
            l, String.concat ~sep:";" clst)
      |> List.sort ~cmp:compare2
    in
    List.iter o ~f:(fun (l, a) -> fprintf oc "%0.20f\t%16s\n" l a)

  let default_zygosity_report_size = 100

  let zygosity ?(size=default_zygosity_report_size) likelihood_first aa
    blarr zt oc =
    fprintf oc "Zygosity:\n";
    let zb = Zygosity_array.best ~size zt in
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
              |> List.sort ~cmp:compare3
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

  let output_likelihood lfirst allele_arr final_likelihoods oc =
    let o = by_likelihood_arr allele_arr final_likelihoods in
    fprintf oc "Likelihood:\n";
    if lfirst then
      likelihood_first allele_arr final_likelihoods oc
    else
      allele_first o oc;
    o

  let f lfirst zygosity_report_size allele_arr
    final_likelihoods zt oc =
    let blarr = output_likelihood lfirst allele_arr final_likelihoods oc in
    zygosity ~size:zygosity_report_size lfirst allele_arr
      blarr zt oc

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

let single_conf ?allele ?insert_p ?band ?max_number_mismatches
  ~past_threshold_filter ~check_rc () =
    { allele
    ; insert_p
    ; band
    ; max_number_mismatches
    ; past_threshold_filter
    ; check_rc
    }

type ('state, 'read_result) t =
  { single  : Biocaml_unix.Fastq.item -> 'read_result
  ; paired  : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> 'read_result
  ; merge   : 'state -> 'read_result -> unit
  ; output  : 'state -> out_channel -> unit
  }

type forward_opt =
  { likelihood_first     : bool
  ; zygosity_report_size : int
  ; report_size          : int
  }

module Forward = struct

  open ParPHMM

  type stat =
    { aap        : Alleles_and_positions.t
    ; likelihood : float mt
    }

  let proc_to_stat report_size proc _comp =
    { aap        = proc.best_allele_pos report_size
    ; likelihood = proc.per_allele_llhd ()
    }

  let forward report_size past_threshold proc read read_errors =
    Orientation.paired past_threshold (proc_to_stat report_size) proc read read_errors

  type per_read =
    { name  : string
    ; aapt  : Alleles_and_positions.t Orientation.t
    }

  type state =
    { per_allele_lhood  : float array
    ; zygosity          : Zygosity_array.t
    ; mutable per_reads : per_read list
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

  let update_state_single t name stat_or =
    let take_regular = take_regular_by_likelihood in
    t.per_reads <- { name; aapt = just_aap_or stat_or } :: t.per_reads;
    let pr = Orientation.most_likely_between stat_or ~take_regular in
    match pr with
    | Filtered m        -> printf "For %s: %s\n" name m
    | Completed (_c, s) ->
        Likelihoods_and_zygosity.add_ll_and_lz
              t.per_allele_lhood t.zygosity s.likelihood

  let output_state likelihood_first zygosity_report_size allele_arr t oc =
    let take_regular r c = Alleles_and_positions.descending_cmp r c <= 0 in
    Output.f likelihood_first zygosity_report_size allele_arr
      t.per_allele_lhood t.zygosity oc;
    fprintf oc "Per Read:\n";
    List.iter t.per_reads ~f:(fun {name; aapt} ->
      fprintf oc "%s\n\t%s\n" name
        (Orientation.to_string ~take_regular
            Alleles_and_positions.lst_to_string aapt))

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

  type opt = forward_opt

  let init conf ~read_length pt opt =
    let { likelihood_first; zygosity_report_size; report_size } = opt in
    let { allele; insert_p; band; max_number_mismatches; past_threshold_filter
        ; _} = conf
    in
    let open ParPHMM in
    let proc, allele_arr =
      to_proc_allele_arr ?allele ?insert_p ?max_number_mismatches ?band
        read_length pt
    in
    let initial_pt = Past_threshold.init past_threshold_filter in
    (* Curry away arguments and discard final threshold; can't use it. *)
    let forward pt r re = fst (forward report_size pt proc r re) in
    let rec state = { per_allele_lhood; zygosity; per_reads; }
    and per_allele_lhood = proc.init_global_state ()
    and zygosity = Zygosity_array.init (Array.length allele_arr)
    and per_reads = []
    in
    let rec single fqi =
      Fastq_items.single fqi ~k:(fun name read read_errors ->
        { name; stat = Single (forward initial_pt read read_errors) })
    and paired fq1 fq2 =
      Fastq_items.paired fq1 fq2 ~k:(fun name rd1 re1 rd2 re2 ->
        let r1 = forward initial_pt rd1 re1 in
        (* A little bit wasteful as we can check the best orientation for the
           2nd read; not implementing this here because it is much better, more
           corner cases are considered, in the Multiple_loci implementation. *)
        let r2 = forward initial_pt rd2 re2 in
        { name; stat = Paired (r1, r2)})
    and merge state rr =
      match rr.stat with
      | Single stat_or  -> update_state_single state rr.name  stat_or
      | Paired (s1, s2) -> update_state_single state (rr.name ^ " 1") s1;
                           update_state_single state (rr.name ^ " 2") s2
    and output t oc =
      output_state likelihood_first zygosity_report_size allele_arr t oc
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

  type read_result =
    { name  : string
    ; res   : viterbi_result
    }

  type state =
    { mutable results : read_result list
    }

  type opt = unit

  let init conf ~read_length pt () =
    let { allele; insert_p; _ } = conf in
    let open ParPHMM in
    (* Relying on reference being first. *)
    let allele = Option.value ~default:(fst pt.alleles.(0)) allele in
    let labels = allele ^ " ", "read " in
    let s, p =
      setup_single_allele_viterbi_pass ?insert_p ~allele
        read_length pt
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
           -> (state, read_result) t * state

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

let multiple_conf ?insert_p ?band ?max_number_mismatches
  ~past_threshold_filter ~incremental_pairs () =
    { insert_p
    ; band
    ; max_number_mismatches
    ; past_threshold_filter
    ; incremental_pairs
    }

module Multiple_loci = struct

  type 'single_result single_or_incremental =
    | SingleRead of 'single_result Orientation.t
    (* We'll test both regular and complement options. *)

    | PairedDependent of
      { first_o : bool
      ; first   : 'single_result
      ; second  : 'single_result
      }
    (* Orientation and loci of 2nd is determined by 1st There is no
       {pass_result} since we're not going to apply a filter to the 2nd read.
       Note that one could "naturally" combine the result of the first by
       passing the first's emission (via per_allele_hood) into the second
       read's forward pass via [base_p] (this is similar to thinking of the 2
       reads as just one long read), but that will dramatically slow down the
       pass (it will start out already branched). It is better to keep the two
       likelihoods separate and then just multiply (log -> add) the results
       together. This is the desired outcome for paired reads as it is the most
       descriminative. *)

  let orientation_char_of_bool reverse_complement =
    if reverse_complement then 'C' else 'R'

  let soi_to_string ?take_regular ?(sep='\t') sr_to_string soi =
    let ots = Orientation.to_string ?take_regular ~sep sr_to_string in
    match soi with
    | SingleRead ors                             ->
        sprintf "SingleRead%c%s" sep (ots ors)
    | PairedDependent { first_o; first; second } ->
        sprintf "PairedDependent%c%c%c%s%c%s"
          sep (orientation_char_of_bool first_o)
          sep (sr_to_string first)
          sep (sr_to_string second)

  type 'single_result paired =
    | FirstFiltered of
      { first   : string
      ; second  : 'single_result Orientation.t
      }
    (* Both orientations of the first read were filtered, so the 2nd read is
       entirely independent of the first. This can occur if we have a filter
       (such as past threshold from another loci) that ignores that loci. *)
    | FirstOrientedSecond of
      { first_o : bool
      ; first   : 'single_result
      ; second  : 'single_result
      }
    (* One of the orientations of the first wasn't filtered and we used
       it to orient the second read. Similar to PairedDependent there is no
       {pass_result} for second since we don't have a comparable filter. *)

  let pr_to_string ?(sep='\t') sr_to_string pr =
    let ots = Orientation.to_string ~sep sr_to_string in
    match pr with
    | FirstFiltered { first; second}                 ->
        sprintf "FirstFiltered%c%s%c%s" sep first sep (ots second)
    | FirstOrientedSecond { first_o; first; second } ->
        sprintf "FirstOrientedSecond%c%c%c%s%c%s"
          sep (orientation_char_of_bool first_o) sep
          (sr_to_string first) sep (sr_to_string second)

  type 'sr potentially_paired =
    | Single_or_incremental of (string * 'sr single_or_incremental) list
    | MPaired of ((string * 'sr paired) list)   (* M for multi to not conflict with Paired *)

  let aap_soi_to_string =
    let take_regular r c = Alleles_and_positions.descending_cmp r c <= 0 in
    soi_to_string ~take_regular

  let pp_to_string ?(sep='\t') sr_to_string = function
    | Single_or_incremental soi_lst ->
        List.map soi_lst ~f:(fun (l, s) ->
          sprintf "%s%c%s" l sep (aap_soi_to_string ~sep sr_to_string s))
        |> String.concat ~sep:"\n"
    | MPaired pr_lst ->
        List.map pr_lst ~f:(fun (l, p) ->
          sprintf "%s%c%s" l sep (pr_to_string ~sep sr_to_string p))
        |> String.concat ~sep:"\n"

  type 'sr per_read =
    { name  : string
    ; rrs   : 'sr potentially_paired
    }

  type per_loci =
    { loci        : string
    ; likelihood  : float array
    ; zygosity    : Zygosity_array.t
    }

  type state =
    { per_loci          : per_loci list

    ; mutable per_read  : Alleles_and_positions.t per_read list
    (* At the final state we'll discard the individual likelihood float array's
       and just store the alleles and positions. *)
    }

  type read_result = Forward.stat per_read

  let take_regular = Forward.take_regular_by_likelihood

  let most_likely_orientation or_ =
    Orientation.most_likely_between ~take_regular or_

  (* A PairedDependent should be chosen over any SingleRead's. *)
  let reduce_soi soi_lst =
    let open ParPHMM in
    let rec update loci stat best tl =
      let l = max_pm stat.Forward.likelihood in
      match best with
      | None                        -> loop (Some (l, loci, Single stat)) tl
      | Some (bl, _, _) when l > bl -> loop (Some (l, loci, Single stat)) tl
      | Some _                      -> loop best tl
    and loop best = function
      | []                -> best
      | (loci, soi) :: tl ->
          begin match soi with
              (* ignore previous bests and remaining loci. *)
          | PairedDependent { first; second; _ } ->
              Some (0.0, loci, Paired (first, second))
          | SingleRead or_                       ->
              begin match most_likely_orientation or_ with
              | Filtered _      -> loop best tl
              | Completed (_,s) -> update loci s best tl
              end
          end
      in
      loop None soi_lst

  let map_soi_to_aap lst =
    List.map lst ~f:(fun (loci, soi) ->
      let n =
        match soi with
        | PairedDependent { first_o; first; second } ->
            PairedDependent { first_o
                            ; first = first.Forward.aap
                            ; second = second.Forward.aap
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
    let open ParPHMM in
    let rec loop best = function
      | []               -> best
      | (loci, pr) :: tl ->
          begin match pr with
          | FirstFiltered { second } ->
              begin match most_likely_orientation second with
              | Filtered  _         -> loop best tl
              | Completed (_, stat) -> loop (update loci stat (Single stat) best) tl
              end
          | FirstOrientedSecond { first; second } ->
              loop (update loci second (Paired (first, second)) best) tl
          end
    in
    loop None pr_lst

  let map_pr_to_aap lst =
    List.map_snd lst ~f:(function
        | FirstFiltered { first; second } ->
            FirstFiltered { first; second = Forward.just_aap_or second }
        | FirstOrientedSecond { first_o; first; second } ->
            FirstOrientedSecond { first_o
                                ; first = first.Forward.aap
                                ; second = second.Forward.aap })

  let merge state { name ; rrs } =
    let open ParPHMM in
    let assign_to_per_loci { loci; likelihood; zygosity } = function
      | Single stat ->
          Likelihoods_and_zygosity.add_ll_and_lz
            likelihood zygosity stat.Forward.likelihood
      | Paired (s1, s2) ->
          Likelihoods_and_zygosity.add_ll_and_lz
            likelihood zygosity s1.Forward.likelihood;
          Likelihoods_and_zygosity.add_ll_and_lz
            likelihood zygosity s2.Forward.likelihood
    in
    let assign_to_best best_loci best_stat =
      List.iter state.per_loci ~f:(fun pl ->
        if best_loci = pl.loci then
          assign_to_per_loci pl best_stat)
    in
    match rrs with
    | Single_or_incremental soi_lst ->
        begin match reduce_soi soi_lst with
        | None  -> eprintf "All read results filtered for %s\n%!" name
        | Some (_bl, best_loci, best_stat)  ->
            assign_to_best best_loci best_stat
        end;
        let mrs = Single_or_incremental (map_soi_to_aap soi_lst) in
        state.per_read <- { name; rrs = mrs } :: state.per_read
    | MPaired pr_lst ->
        begin match reduce_mpr pr_lst with
        | None -> eprintf "All 2nd read results filtered for %s\n%!" name
        | Some (_bl, best_loci, best_stat)  ->
            assign_to_best best_loci best_stat
        end;
        let mrs = MPaired (map_pr_to_aap pr_lst) in
        state.per_read <- { name; rrs = mrs } :: state.per_read

  let output_per_reads prlst oc =
    fprintf oc "Per Read:\n";
    List.iter prlst ~f:(fun {name; rrs} ->
      fprintf oc "%s\n%s\n" name
        (pp_to_string Alleles_and_positions.lst_to_string rrs))

  let init conf opt read_length parPHMM_t_lst =
    let { likelihood_first; zygosity_report_size; report_size } = opt in
    let { insert_p; band; max_number_mismatches; _ } = conf in
    let paa =
      List.map_snd parPHMM_t_lst ~f:(fun parPHMM_t ->
        let proc, allele_arr =
          Forward.to_proc_allele_arr ?insert_p ?max_number_mismatches ?band
            read_length parPHMM_t
        in
        let n = Array.length allele_arr in
        proc, allele_arr, n)
    in
    let initial_pt = Past_threshold.init conf.past_threshold_filter in
    let forward = Forward.forward report_size in
    let filterless_single p rc r re =
      match single p (Forward.proc_to_stat report_size) rc r re with
      | ParPHMM.Filtered m  -> invalid_argf "Filtered result %s without filter!" m
      | ParPHMM.Completed s -> s
    in
    let ordinary_paired name rd1 re1 rd2 re2 =
      let _fpt1, _fpt2, res =
        List.fold_left paa ~init:(initial_pt, initial_pt, [])
          ~f:(fun (pt1, pt2, acc) (loci, (p, _, _)) ->
                let result1, npt1 = forward pt1 p rd1 re1 in
                match most_likely_orientation result1 with
                | ParPHMM.Filtered  m ->
                    let second, npt2 = forward pt2 p rd2 re2 in
                    npt1, npt2, (loci, FirstFiltered { first = m; second }) :: acc
                | ParPHMM.Completed (first_o, first) ->
                    let second = filterless_single p (not first_o) rd2 re2 in
                    let npt2 = Past_threshold.update pt2 p (ParPHMM.Completed second) in
                    npt1, npt2, (loci, FirstOrientedSecond { first_o; first; second }) :: acc)
      in
      MPaired (List.rev res)
    in
    let incremental_paired name rd1 re1 rd2 re2 =
      let best, rest, _fpt1 =
        List.fold_left paa ~init:(None, [], initial_pt)
          ~f:(fun (best, acc, pt) (loci, (p, _, _)) ->
                let result1, npt = forward pt p rd1 re1 in
                match most_likely_orientation result1 with
                | ParPHMM.Filtered m ->
                    let nacc = (loci, SingleRead result1) :: acc in
                    best, nacc, npt
                | ParPHMM.Completed (c, first) ->
                    let l = max_pm first.Forward.likelihood in
                    begin match best with
                    | None ->
                        Some (loci, p, result1, c, first, l), acc, npt
                    | Some (bloci, _bp, bresult1, _bcomp, _bfirst, bl) ->
                        if l > bl then  (* new best *)
                          let nacc = (bloci, SingleRead bresult1) :: acc in
                          Some (loci, p, result1, c, first, l), nacc, npt
                        else
                          let nacc = (loci, SingleRead result1) :: acc in
                          best, nacc, npt
                    end)
      in
      match best with
      | None  ->
          eprintf "All loci were filtered for first read of %s, ignoring second.\n"
            name;
          Single_or_incremental (List.rev rest)
      | Some (bloci, bp, _bresult1, first_o, first, _bl) ->
          let second = filterless_single bp (not first_o) rd2 re2 in
          let nrest = (bloci, PairedDependent {first_o; first; second}) :: (List.rev rest) in
          Single_or_incremental nrest
    in
    let rec state = { per_loci; per_read = []}
    and per_loci =
      List.map paa ~f:(fun (loci, (p, _, number_alleles)) ->
        { loci
        ; likelihood = p.ParPHMM.init_global_state ()
        ; zygosity = Zygosity_array.init number_alleles
        })
    in
    let rec t = { single; paired; merge; output }
    and single fqi =
      Fastq_items.single fqi ~k:(fun name read read_errors ->
        { name
        ; rrs =
            List.fold_left paa ~init:(initial_pt, [])
                ~f:(fun (pt, acc) (loci, (p, _, _)) ->
                      let r, npt = forward pt p read read_errors in
                      npt, (loci, SingleRead r) :: acc)
            |> snd                                 (* Discard final threshold *)
            |> List.rev                                (* Restore loci order. *)
            |> fun l -> Single_or_incremental l
        })
    and paired fq1 fq2 =
     Fastq_items.paired fq1 fq2 ~k:(fun name rd1 re1 rd2 re2 ->
        { name
        ; rrs =
            if conf.incremental_pairs then
              incremental_paired name rd1 re1 rd2 re2
            else
              ordinary_paired name rd1 re1 rd2 re2
        })
    and output state oc =
      List.iter2 paa state.per_loci
        ~f:(fun (_, (_, allele_arr, _)) { loci; likelihood; zygosity; _} ->
              fprintf oc "%s\n" loci;
              Output.f likelihood_first zygosity_report_size
                allele_arr likelihood zygosity oc);
      output_per_reads state.per_read oc
    in
    t, state

end (* Multiple_loci *)
