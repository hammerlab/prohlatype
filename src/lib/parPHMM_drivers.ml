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

type mapped_complement_stats = mapper_stats rc_states

(* The output is ordered so that the left most stat is the most likely.
   We add characters 'R', 'C', 'F' as little annotations. *)
let mapped_complement_stats_to_string ?(sep='\t') ms =
  match ms.regular, ms.complement with
  | Completed r, Completed c ->
    if fst (List.hd_exn r.positions) > fst (List.hd_exn c.positions) then
      sprintf "R%c%s%cC%c%s"
        sep (mapper_stats_to_string ~sep r)
        sep sep (mapper_stats_to_string ~sep c)
    else
      sprintf "C%c%s%cR%c%s"
        sep (mapper_stats_to_string ~sep c)
        sep sep (mapper_stats_to_string ~sep r)
  | Completed r, Filtered m ->
      sprintf "R%c%s%cF%c%s" sep (mapper_stats_to_string ~sep r) sep sep m
  | Filtered m, Completed c ->
      sprintf "C%c%s%cF%c%s" sep (mapper_stats_to_string ~sep c) sep sep m
  | Filtered mr, Filtered mc ->
      sprintf "F%c%s%cF%c%s" sep mr sep sep mc

type reducer_stats =
  { comp       : bool
  ; likelihood : float array
  }

let proc_to_reducer_stat proc comp =
  { comp ; likelihood = proc.per_allele_llhd ()}

(* Consider a single orientation: either regular or reverse complement. *)
let single ?prev_threshold proc pass_result_map reverse_complement read read_prob =
  match proc.doit ?prev_threshold reverse_complement read read_prob with
  | Filtered m    -> Filtered m
  | Completed ()  -> Completed (pass_result_map proc reverse_complement)

let new_threshold prev_threshold proc =
  max prev_threshold (proc.maximum_match ())

let update_threshold prev_threshold proc = function
  | Filtered _  -> prev_threshold
  | Completed _ -> new_threshold prev_threshold proc

let init_past_threshold = function
  | None
  | Some false  -> `Don't
  | Some true   -> `Start

(* Consider both orientations. *)
let both pt pass_result_map proc read read_prob =
  let do_work ?prev_threshold rc =
    single ?prev_threshold proc pass_result_map rc read read_prob
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

let reducer ?(check_rc=true) past_threshold proc read read_prob =
  if check_rc then begin
    let mp, pt =
      both past_threshold
        (fun proc comp ->
            let likelihood = proc.per_allele_llhd () in
            { comp ; likelihood })
        proc read read_prob
    in
    most_likely_between_rc (mp.regular, mp.complement), pt
  end else
    (* Can't care about previous threshold! *)
    single proc proc_to_reducer_stat false read read_prob, None

exception Read_error_parsing of string

let repf fmt =
  ksprintf (fun s -> raise (Read_error_parsing s)) fmt

let fqi_to_read_and_probs fqi ~k =
  let open Core_kernel.Std in
  let open Biocaml_unix.Fastq in
  time (sprintf "updating on %s" fqi.name) (fun () ->
    match Fastq.phred_log_probs fqi.qualities with
    | Result.Error e       -> repf "%s" (Error.to_string_hum e)
    | Result.Ok read_probs -> k fqi.name fqi.sequence read_probs)

let fqi2_to_read_and_probs fq1 fq2 ~k =
  let open Biocaml_unix.Fastq in
  let open Core_kernel.Std in
  assert (fq1.name = fq2.name);
  time (sprintf "updating on %s" fq1.name) (fun () ->
    match Fastq.phred_log_probs fq1.Biocaml_unix.Fastq.qualities with
    | Result.Error e       -> repf "%s" (Error.to_string_hum e)
    | Result.Ok rp1 ->
        match Fastq.phred_log_probs fq2.Biocaml_unix.Fastq.qualities with
        | Result.Error e       -> repf "%s" (Error.to_string_hum e)
        | Result.Ok rp2 ->
            k fq1.name fq1.sequence rp1 fq2.sequence rp2)

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

  ; joined_pairs          : bool option

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
    ; paired  : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> unit
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
    let past_threshold = init_past_threshold past_threshold_filter in
    let n = Array.length allele_arr in
    (* Discard resulting past threshold *)
    let r rd rp = fst (reducer ?check_rc past_threshold proc rd rp) in
    let rec t = { state; apply; paired; output }
    and state = proc.init_global_state ()
    and apply fqi =
      fqi_to_read_and_probs fqi ~k:(fun _name read read_probs ->
        r read read_probs
        |> add_log_likelihoods n t.state)
    and paired fq1 fq2 =
      fqi2_to_read_and_probs fq1 fq2 ~k:(fun _name rd1 rp1 rd2 rp2 ->
        let r1 = r rd1 rp1 in
        let r2 =
          match r1 with
          | Filtered _  -> r rd2 rp2 (* Filtered so check both conf *)
          | Completed c -> single proc proc_to_reducer_stat (not c.comp) rd2 rp2
        in
        add_log_likelihoods n t.state r1;
        add_log_likelihoods n t.state r2)
    and output oc = output_i allele_arr t.state oc in
    t

  let apply fqi t =
    t.apply fqi

  let paired fq1 fq2 t =
    t.paired fq1 fq2

  let output oc t =
    t.output oc

end (* Reducer *)

type 'a single_or_paired =
  | Single of 'a
  | Paired of ('a * 'a)

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
      fqi_to_read_and_probs ~k:(fun name read read_probs ->
        let rc = m read read_probs in
        t.state <- (name, Single rc) :: t.state)
    and paired fq1 fq2 =
      fqi2_to_read_and_probs fq1 fq2 ~k:(fun name rd1 rp1 rd2 rp2 ->
        let rc1 = m rd1 rp1 in
        let rc2 = m rd2 rp2 in
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

  type viterbi_result =
    { reverse_complement : bool       (* Was the reverse complement best? *)
    ; emission           : float      (* Final emission likelihood *)
    ; path_list          : path list
    }

  let viterbi_result_to_string ?(width=120) ?(labels=("reference ", "read "))
    {reverse_complement; emission; path_list} =
    let { reference; read } = path_list_to_strings path_list in
    sprintf "Reverse complement: %b, final emission: %f, path length: %d\n%s"
      reverse_complement emission (List.length path_list)
        (manual_comp_display ~width ~labels reference read)

  type vr = viterbi_result single_or_paired

  let vr_to_string ~labels = function
    | Single v -> viterbi_result_to_string ~labels v
    | Paired (v1, v2) -> sprintf "%s\n%s" (viterbi_result_to_string ~labels v1)
                                          (viterbi_result_to_string ~labels v2)

  type t =
    { mutable state : (string * vr) list
    ; apply         : Biocaml_unix.Fastq.item -> unit
    ; paired        : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> unit
    ; output        : out_channel -> unit
    }

  let init conf read_length pt =
    let { allele; insert_p; joined_pairs; _ } = conf in
    let open ParPHMM in
    (* Relying on reference being first. *)
    let allele = Option.value ~default:(fst pt.alleles.(0)) allele in
    let labels = allele ^ " ", "read " in
    let of_triple v read read_probs =
      let (reverse_complement, emission, path_list) = v read read_probs in
      { reverse_complement; emission; path_list}
    in
    match joined_pairs with
    | Some true ->
      begin
        let length = read_length in
        let read_length = read_length * 2 in
        let v = setup_single_allele_viterbi_pass ?insert_p ~allele read_length pt in
        let vr = of_triple v in
        let rec t = { state = []; apply; paired; output }
        and apply =
          fqi_to_read_and_probs ~k:(fun read_name read read_probs ->
            t.state <- (read_name, Single (vr read read_probs)) :: t.state)
        and paired fq1 fq2 =
          fqi2_to_read_and_probs fq1 fq2 ~k:(fun read_name rd1 rp1 rd2 rp2 ->
            let rd = String.make read_length 'A' in
            String.blit ~src:rd1 ~src_index:0 ~dst:rd ~dst_index:0 ~length;
            let rc = reverse_complement rd2 in
            String.blit ~src:rc ~src_index:0 ~dst:rd ~dst_index:length ~length;
            let rp = Array.create_float read_length in
            Array.blit ~src:rp1 ~src_pos:0 ~dst:rp ~dst_pos:0 ~len:length;
            for i = 0 to length - 1 do
              rp.(length + i) <- rp2.(length - 1 - i)
            done;
            let v = vr rd rp in
            t.state <- (read_name, Single v) :: t.state)
        and output oc =
          List.iter t.state ~f:(fun (read_name, vr ) ->
            fprintf oc "%s:\n%s\n" read_name (vr_to_string ~labels vr))
        in
      t
      end
    | _ ->
      let v = setup_single_allele_viterbi_pass ?insert_p ~allele
                read_length pt in
      let vr = of_triple v in
      let rec t = { state = []; apply; paired; output }
      and apply =
        fqi_to_read_and_probs ~k:(fun read_name read read_probs ->
          t.state <- (read_name, Single (vr read read_probs)) :: t.state)
      and paired fq1 fq2 =
        fqi2_to_read_and_probs fq1 fq2 ~k:(fun read_name rd1 rp1 rd2 rp2 ->
          let v1 = vr rd1 rp1 in
          let v2 = vr rd2 rp2 in
          t.state <- (read_name, Paired (v1, v2)) :: t.state)
      and output oc =
        List.iter t.state ~f:(fun (read_name, vr ) ->
          fprintf oc "%s:\n%s\n" read_name (vr_to_string ~labels vr))
      in
      t

  let apply fqi t =
    t.apply fqi

  let paired fq1 fq2 t =
    t.paired fq1 fq2

  let output oc t =
    t.output oc

end (* Viterbi *)

(** Single Loci Drivers **)
module Single_loci = struct

  let conf ?allele ?insert_p ?band ?max_number_mismatches
    ?joined_pairs ?past_threshold_filter ?check_rc () =
    { allele; insert_p; band; max_number_mismatches
    ; joined_pairs ; past_threshold_filter ; check_rc
    }

  type t =
    | Reducer of Reducer.t
    | Mapper of Mapper.t
    | Viterbi of Viterbi.t

  let init conf read_length pt mode =
      match mode with
      | `Reducer  -> Reducer (Reducer.init conf read_length pt)
      | `Mapper n -> Mapper (Mapper.init conf read_length n pt)
      | `Viterbi  -> Viterbi (Viterbi.init conf read_length pt)

  (* fqi = FastQ Item *)
  let apply t fqi = match t with
    | Reducer r -> Reducer.apply fqi r
    | Mapper m  -> Mapper.apply fqi m
    | Viterbi v -> Viterbi.apply fqi v

  let paired t fq1 fq2 = match t with
    | Reducer r -> Reducer.paired fq1 fq2 r
    | Mapper m  -> Mapper.paired fq1 fq2 m
    | Viterbi v -> Viterbi.paired fq1 fq2 v

  let output t oc = match t with
    | Reducer r -> Reducer.output oc r
    | Mapper m  -> Mapper.output oc m
    | Viterbi v -> Viterbi.output oc v

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

  let add_log_likelihoods n lklhd_arr = function
    | Single r        -> Reducer.add_log_likelihoods n lklhd_arr r
    | Paired (r1, r2) -> Reducer.add_log_likelihoods n lklhd_arr r1;
                         Reducer.add_log_likelihoods n lklhd_arr r2

  let merge state current =
    let b = best current in
    List.iter state ~f:(fun (name, n, lklhd_arr) ->
      if name = b.name then add_log_likelihoods n lklhd_arr b.rspr)

  type t =      (* loci, # alleles, likelihoods *)
    { state   : (string * int * float array) list
    ; apply   : Biocaml_unix.Fastq.item -> unit
    ; paired  : Biocaml_unix.Fastq.item -> Biocaml_unix.Fastq.item -> unit
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
    let pt = init_past_threshold past_threshold_filter in
    let rec t = { state; apply; paired; output }
    and state =
      List.map paa ~f:(fun (name, (p, _, n)) -> name, n, p.init_global_state ())
    and output oc =
      List.iter2 paa t.state
        ~f:(fun (_, (_, allele_arr, _)) (name, _, lkld_arr) ->
              fprintf oc "%s\n" name;
              Reducer.output_i allele_arr lkld_arr oc)
    and apply fqi =
      match pt with
      | `Don't ->
          fqi_to_read_and_probs fqi ~k:(fun _name read read_probs ->
            List.map_snd paa ~f:(fun (p, _, _) -> Single (fst (reducer pt p read read_probs)))
            |> merge t.state)
      | `Start ->
          fqi_to_read_and_probs fqi ~k:(fun _name read read_probs ->
          List.fold_left paa ~init:(pt, [])
            ~f:(fun (pt, acc) (name, (p, _, _)) ->
                  let m, npto = reducer pt p read read_probs in
                  let npt = update_pt npto in
                  (npt, (name, Single m) :: acc))
          |> function _, acc ->
                merge t.state acc)
    and paired fq1 fq2 =
      match pt with
      | `Don't ->
          fqi2_to_read_and_probs fq1 fq2 ~k:(fun _name rd1 rp1 rd2 rp2 ->
            List.map_snd paa ~f:(fun (p, _, _) ->
              let r1 = fst (reducer `Don't p rd1 rp1) in
              let r2 =
                match r1 with
                | Filtered _  -> fst (reducer `Don't p rd2 rp2) (* Filtered so check both conf *)
                | Completed c -> single p proc_to_reducer_stat (not c.comp) rd2 rp2
              in
              Paired (r1, r2))
            |> merge t.state)
      | `Start ->
          fqi2_to_read_and_probs fq1 fq2 ~k:(fun _name rd1 rp1 rd2 rp2 ->
            List.fold_left paa ~init:(pt, pt, [])
              ~f:(fun (pt1, pt2, acc) (name, (p, _, _)) ->
                    let m1, npto1 = reducer pt1 p rd1 rp1 in
                    let npt1 = update_pt npto1 in
                    let m2, npto2 = reducer pt2 p rd2 rp2 in
                    let npt2 = update_pt npto2 in
                    let nacc = (name, Paired (m1, m2)) :: acc in
                    (npt1, npt2, nacc)))
            |> function _, _, nacc -> merge t.state nacc
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
          fqi_to_read_and_probs fqi ~k:(fun read_name read read_probs ->
            let r = List.map_snd procs ~f:(fun p ->
                      Single (fst (mapper ~n `Don't p read read_probs)))
            in
            t.state <- (read_name, r) :: t.state)
      | `Start  ->
          (* Now we have to take into account previous filter results. *)
          fqi_to_read_and_probs fqi ~k:(fun read_name read read_probs ->
            List.fold_left procs ~init:(`Start, [])
              ~f:(fun (pt, acc) (name, p) ->
                    let m, npto = mapper ~n pt p read read_probs in
                    let npt = update_pt npto in
                    (npt, (name, Single m) :: acc))
            |> function _, acc ->
                  t.state <- (read_name, acc) :: t.state)
    and paired fq1 fq2 =
      match pt with
      | `Don't ->
          fqi2_to_read_and_probs fq1 fq2 ~k:(fun read_name rd1 rp1 rd2 rp2 ->
            let r = List.map_snd procs ~f:(fun p ->
              let m1, _ = mapper `Don't ~n p rd1 rp1 in
              let m2, _ = mapper `Don't ~n p rd2 rp2 in
              Paired (m1, m2))
            in
            t.state <- (read_name, r) :: t.state)
      | `Start ->
          (* Now we have to take into account previous filter results. *)
          fqi2_to_read_and_probs fq1 fq2 ~k:(fun read_name rd1 rp1 rd2 rp2 ->
            List.fold_left procs ~init:(`Start, `Start, [])
              ~f:(fun (pt1, pt2, acc) (name, p) ->
                    let m1, npto1 = mapper ~n pt1 p rd1 rp1 in
                    let npt1 = update_pt npto1 in
                    let m2, npto2 = mapper ~n pt2 p rd2 rp2 in
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

  let paired t fq1 fq2 = match t with
    | Reducer r -> Multiple_reducer.paired fq1 fq2 r
    | Mapper m  -> Multiple_mapper.paired fq1 fq2 m

  let output t oc = match t with
    | Reducer r -> Multiple_reducer.output oc r
    | Mapper m  -> Multiple_mapper.output oc m

end (* Mulitple_loci *)

