(* Aggregation of the forward passes. *)

open Util
open ParPHMM

type stats =
  { per_allele : (float * string) list
  ; positions  : (float * int) list
  }

let stats_to_strings ?(sep='\t') s =
  let l_to_s fmt l =
    String.concat ~sep:";" (List.map l ~f:(fun (l,a) -> sprintf fmt  a l))
  in
  let al_to_s l = l_to_s "%s:%0.2f" l in
  let pl_to_s l = l_to_s "%d:%0.2f" l in
  sprintf "%s%c%s" (al_to_s s.per_allele) sep (pl_to_s s.positions)

type 'a rc_states =
  { regular         : 'a pass_result
  ; complement      : 'a pass_result
  ; last_threshold  : float
  }

type mapped_stats = stats rc_states

(* The output is ordered so that the left most stat is the most likely.
   We add characters 'R', 'C', 'F' as little annotations. *)
let mapped_stats_to_string ?(sep='\t') ms =
  match ms.regular, ms.complement with
  | Completed r, Completed c ->
    if fst (List.hd_exn r.positions) > fst (List.hd_exn c.positions) then
      sprintf "R%c%s%cC%c%s"
        sep (stats_to_strings ~sep r)
        sep sep (stats_to_strings ~sep c)
    else
      sprintf "C%c%s%cR%c%s"
        sep (stats_to_strings ~sep c)
        sep sep (stats_to_strings ~sep r)
  | Completed r, Filtered m ->
      sprintf "R%c%s%cF%c%s" sep (stats_to_strings ~sep r) sep sep m
  | Filtered m, Completed c ->
      sprintf "C%c%s%cF%c%s" sep (stats_to_strings ~sep c) sep sep m
  | Filtered mr, Filtered mc ->
      sprintf "F%c%s%cF%c%s" sep mr sep sep mc

let pass_result_to_stat = function
  | Completed s -> fst (List.hd_exn s.positions)
  | Filtered _  -> neg_infinity

let best_stat ms =
  max (pass_result_to_stat ms.regular)
      (pass_result_to_stat ms.complement)

let or_neg_inf =
  Option.value ~default:neg_infinity

let new_threshold proc prev_threshold =
  let r = proc.maximum_match () in
  let p = or_neg_inf prev_threshold in
  max r p

let latest_threshold_of_pass_result prev_threshold proc = function
  | Filtered _  -> or_neg_inf prev_threshold
  | Completed _ -> new_threshold proc prev_threshold

let both ?past_threshold_filter pass_res_to_stat proc read read_prob =
  let do_work prev_threshold rc =
    match proc.doit ?prev_threshold rc read read_prob with
    | Filtered m   -> Filtered m
    | Completed () -> Completed (pass_res_to_stat rc proc)
  in
  let reg_prev_threshold, comp_prev_threshold =
    match past_threshold_filter with
    | None          -> None,    false
    | Some `Start   -> None,    true
    | Some (`Set v) -> Some v,  true
  in
  let regular = do_work reg_prev_threshold false in
  let comp_prev_threshold =
    if comp_prev_threshold then
      Some (latest_threshold_of_pass_result reg_prev_threshold
              proc regular)
    else
      None
  in
  let complement = do_work comp_prev_threshold true in
  let last_threshold =
    latest_threshold_of_pass_result comp_prev_threshold
      proc complement
  in
  { regular; complement; last_threshold }

let mapper ~n ?past_threshold_filter =
  both ?past_threshold_filter
    (fun _rc proc ->
      { per_allele = proc.best_alleles n
      ; positions  = proc.best_positions n
      })

let max_arr = Array.fold_left ~init:neg_infinity ~f:max

let compare_emissions e1 e2 =
  let r1 = max_arr e1 in
  let r2 = max_arr e2 in
  r1 >= r2

type reducer_stat =
  { comp       : bool
  ; likelihood : float array
  }

let most_likely_between_rc = function
  | Filtered mr, Filtered mc -> Filtered (sprintf "Both! %s, %s" mr mc)
  | Filtered _r, Completed c -> Completed c
  | Completed r, Filtered _c -> Completed r
  | Completed r, Completed c ->
    if compare_emissions r.likelihood c.likelihood then
      Completed r
    else
      Completed c

let reducer proc ?(check_rc=true) ?past_threshold_filter read read_prob =
  if check_rc then begin
    let mp =
      both ?past_threshold_filter
        (fun comp proc ->
            let likelihood = proc.per_allele_llhd () in
            { comp ; likelihood })
        proc read read_prob
    in
    most_likely_between_rc (mp.regular, mp.complement)
  end else begin
    match proc.doit false read read_prob with
    | Filtered m    ->
        Filtered m
    | Completed ()  ->
        let likelihood = proc.per_allele_llhd () in
        Completed { comp = false; likelihood }
  end

exception Read_error_parsing of string

let repf fmt =
  ksprintf (fun s -> raise (Read_error_parsing s)) fmt

let fqi_to_read_and_probs fqi ~k =
  let name = fqi.Biocaml_unix.Fastq.name in
  time (sprintf "updating on %s" name) (fun () ->
    let open Core_kernel.Std in
    match Fastq.phred_log_probs fqi.Biocaml_unix.Fastq.qualities with
    | Result.Error e       -> repf "%s" (Error.to_string_hum e)
    | Result.Ok read_probs -> k name fqi.Biocaml_unix.Fastq.sequence read_probs)

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

  ; past_threshold_filter : bool option
  (* Use the previous match likelihood, when available (ex. against reverse
     complement), as a threshold filter for the forward pass. *)

  ; check_rc              : bool option
  (** Compare against the reverse complement of a read and take the best likelihood. *)
  }


(* A reducer reduces the information from a set of reads into a per allele
   likelihood state. We can construct reducers that operate on just a single
   allele. *)
module Reducer = struct

  type t =
    { state   : float array  (* Need to track the per allele likelihoods *)
    ; apply   : Biocaml_unix.Fastq.item -> unit
    ; output  : out_channel -> unit
    }

  let output_i allele_arr final_likelihoods =
    let o = Array.mapi allele_arr ~f:(fun i a -> final_likelihoods.(i), a) in
    let cmp (l1, a1) (l2, a2) =
      if l2 < l1 then
        -1
      else if l2 > l1 then
        1
      else
        compare a1 a2 (* Sort by the nomenclature order, not strings. *)
    in
    Array.sort o ~cmp;
    fun oc ->
      Array.iter o ~f:(fun (l,a) -> fprintf oc "%10s\t%0.20f\n" a l)

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
    | Filtered _  -> ()                     (* Should we warn about ignoring a read? *)
    | Completed r ->
        for i = 0 to n - 1 do
          state.(i) <- state.(i) +. r.likelihood.(i)
        done

  let init conf read_length pt =
    let { allele; insert_p; band; max_number_mismatches; past_threshold_filter
        ; check_rc } = conf
    in
    let open ParPHMM in
    let proc, allele_arr =
      to_proc_allele_arr ?allele ?insert_p ?max_number_mismatches ?band
        read_length pt
    in
    let past_threshold_filter =
      Option.bind past_threshold_filter
        ~f:(function | true -> Some `Start | false -> None)
    in
    let r = reducer ?check_rc ?past_threshold_filter proc in
    let rec t = { state ; apply ; output }
    and state = proc.init_global_state ()
    and apply fqi =
      let n = Array.length allele_arr in
      fqi_to_read_and_probs fqi ~k:(fun _name read read_probs ->
        r read read_probs
        |> add_log_likelihoods n t.state)
    and output oc = output_i allele_arr t.state oc in
    t

  let apply fqi t =
    t.apply fqi

  let output oc t =
    t.output oc

end (* Reducer *)

(* A mapper maps the information from each read into a summarization of how
   a specific HLA-loci aligns the read. *)
module Mapper = struct

  type t =
    { mutable state : (string * mapped_stats) list
    ; apply         : Biocaml_unix.Fastq.item -> unit
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
          let past_threshold_filter =
            Option.bind past_threshold_filter
              ~f:(function | true -> Some `Start | false -> None)
          in
          mapper ~n ?past_threshold_filter proc
      | Some _ -> (* Ignore since we do want to map *)
          mapper ~n proc
    in
    let rec t = { state = [] ; apply ; output }
    and apply =
      fqi_to_read_and_probs ~k:(fun name read read_probs ->
        let rc = m read read_probs in
        t.state <- (name, rc) :: t.state)
    and output oc =
      fprintf oc "Reads: %d\n" (List.length t.state);
      List.sort t.state ~cmp:(fun (_n1, ms1) (_n2, ms2) ->
        compare (best_stat ms1) (best_stat ms2))
      |> List.rev
      |> List.iter ~f:(fun (n, s) ->
          fprintf oc "%s\t%s\n" n (mapped_stats_to_string ~sep:'\t' s))
    in
    t

  let apply fqi t =
    t.apply fqi

  let output oc t =
    t.output oc

end (* Mapper *)

(** Single Loci Drivers **)
module Single_loci = struct

  let conf ?allele ?insert_p ?band ?max_number_mismatches ?past_threshold_filter
    ?check_rc () =
    { allele; insert_p; band; max_number_mismatches; past_threshold_filter
    ; check_rc }

  type t =
    | Reducer of Reducer.t
    | Mapper of Mapper.t

  let init conf read_length pt mode =
      match mode with
      | `Reducer ->
          Reducer (Reducer.init conf read_length pt)
      | `Mapper n ->
          Mapper (Mapper.init conf read_length n pt)

  (* fqi = FastQ Item *)
  let apply t fqi = match t with
    | Reducer r -> Reducer.apply fqi r
    | Mapper m  -> Mapper.apply fqi m

  let output t oc = match t with
    | Reducer r -> Reducer.output oc r
    | Mapper m  -> Mapper.output oc m

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

  ; past_threshold_filter : bool option
  (* Use the previous match likelihood, when available (ex. against reverse
     complement), as a threshold filter for the forward pass. *)

  }

(* Reduce to multiple loci. *)
module Multiple_reducer = struct

  type best_state =
    { name  : string
    ; rspr  : reducer_stat pass_result
    ; maxl  : float
    }

  let pr_to_llhd = function
    | Filtered _  -> neg_infinity
    | Completed c -> max_arr c.likelihood

  let to_bs (name, rspr) =
    { name; rspr; maxl = pr_to_llhd rspr}

  let best_bs bs1 bs2 =
    (*printf "%s %f vs %s %f\n%!" bs1.name bs1.maxl bs2.name bs2.maxl; *)
    if bs1.maxl >= bs2.maxl then bs1 else bs2

  let best = function
    | []          -> invalid_argf "Can't select best from empty!"
    | h :: t ->
        let init = to_bs h in
        List.fold_left t ~init ~f:(fun bs p ->
          best_bs bs (to_bs p))

  let merge state current =
    let b = best current in
    List.iter state ~f:(fun (name, n, lklhd_arr) ->
      if name = b.name then Reducer.add_log_likelihoods n lklhd_arr b.rspr)

  type t =
    { state   : (string * int * float array) list
    ; apply   : Biocaml_unix.Fastq.item -> unit
    ; output  : out_channel -> unit
    }

  let init conf read_length pt_lst =
    let { insert_p; band; max_number_mismatches; past_threshold_filter } = conf in
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
    let rec t = { state; apply; output }
    and state =
      List.map paa ~f:(fun (name, (p, _, n)) -> name, n, p.init_global_state ())
    and output oc =
      List.iter2 paa t.state
        ~f:(fun (_, (_, allele_arr, _)) (name, _, lkld_arr) ->
              fprintf oc "%s\n" name;
              Reducer.output_i allele_arr lkld_arr oc)
    and apply fqi =
      match past_threshold_filter with
      | None
      | Some false  ->
          fqi_to_read_and_probs fqi ~k:(fun _name read read_probs ->
            List.map_snd paa ~f:(fun (p, _, n) -> reducer p read read_probs)
            |> merge t.state)
      | Some true   ->
          (* Now we have to take into account previous filter results. *)
          fqi_to_read_and_probs fqi ~k:(fun _name read read_probs ->
            List.fold_left paa ~init:(`Start, neg_infinity, [])
              ~f:(fun (past_threshold_filter, pt, acc) (name, (p, _, _)) ->
                    let m = reducer ~past_threshold_filter p read read_probs in
                    let nt = max pt (p.maximum_match ()) in
                    (`Set nt, nt, (name, m) :: acc))
            |> function _, _, acc ->
                  merge t.state acc)
    in
    t

  let apply fqi t =
    t.apply fqi

  let output oc t =
    t.output oc

end (* Multiple_reducer *)

module Multiple_mapper = struct

  type t =
    { mutable state : (string * (string * mapped_stats) list) list
    ; apply         : Biocaml_unix.Fastq.item -> unit
    ; output        : out_channel -> unit
    }

  let sort_mapped_output lst =
    let highest_llhd l =
      List.map l ~f:(fun (_name, s) -> best_stat s)
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
    let rec t = { state; apply; output }
    and state = []
    and output oc =
      fprintf oc "Reads: %d\n" (List.length t.state);
      sort_mapped_output t.state
      |> List.iter ~f:(fun (read, lst) ->
        fprintf oc "%s\n\t%s\n"
          read (String.concat ~sep:"\n\t"
                  (List.map lst ~f:(fun (name, ms) ->
                    sprintf "%s\t%s"
                      name (mapped_stats_to_string ~sep:'\t' ms)))))
    and apply fqi =
      match past_threshold_filter with
      | None
      | Some false  ->
          fqi_to_read_and_probs fqi ~k:(fun read_name read read_probs ->
            let r = List.map_snd procs ~f:(fun p -> mapper ~n p read read_probs) in
            t.state <- (read_name, r) :: t.state)
      | Some true   ->
          (* Now we have to take into account previous filter results. *)
          fqi_to_read_and_probs fqi ~k:(fun read_name read read_probs ->
            List.fold_left procs ~init:(`Start, neg_infinity, [])
              ~f:(fun (past_threshold_filter, pt, acc) (name, p) ->
                    let m = mapper ~n ~past_threshold_filter p read read_probs in
                    let nt = max pt (p.maximum_match ()) in
                    (`Set nt, nt, (name, m) :: acc))
            |> function _, _, acc ->
                  t.state <- (read_name, acc) :: t.state)
    in
    t

  let apply fqi t =
    t.apply fqi

  let output oc t =
    t.output oc

end (* Multiple_mapper *)

module Mulitple_loci = struct

  let conf ?insert_p ?band ?max_number_mismatches ?past_threshold_filter () =
    { insert_p; band; max_number_mismatches; past_threshold_filter }

  type t =
    | Reducer of Multiple_reducer.t
    | Mapper of Multiple_mapper.t

  let init conf read_length pt mode =
      match mode with
      | `Reducer ->
          Reducer (Multiple_reducer.init conf read_length pt)
      | `Mapper n ->
          Mapper (Multiple_mapper.init conf read_length n pt)

  (* fqi = FastQ Item *)
  let apply t fqi = match t with
    | Reducer r -> Multiple_reducer.apply fqi r
    | Mapper m  -> Multiple_mapper.apply fqi m

  let output t oc = match t with
    | Reducer r -> Multiple_reducer.output oc r
    | Mapper m  -> Multiple_mapper.output oc m

end (* Mulitple_loci *)

