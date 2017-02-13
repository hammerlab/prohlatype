(** Path Inference : Given a read (a sequence of characters return a PDF over
    the edges (alleles) in the graph. *)

open Util

type 'a product = 'a Alignment.product =
  | Finished of 'a Alleles.Map.t
  | Stopped of 'a Alleles.Map.t

(* The configuration of how we compute a metric of a single nucleotide sequence
    (abstracted by the [thread] type) versus the alleles in our variant graphs.
    This interface acts as a configuration for the [AgainstSequence] functor.

    In practice, this interface is not used, but acts as a stepping stone to
    understanding the elements necessary for [Multiple_config]. *)
module type Single_config = sig

  (* The metric of a single sequence against a single allele.
     For example:
        - an integer (probably non-negative) of the number of mismatches
        - a float for the product log-likelihood's across all the positions
        - a list of all the mismatch locations! *)
  type m

  (* For display purposes. *)
  val metric_to_string : m -> string

  (* An abstract nucleotide sequence. This might be:
      - just a sequence of letters.
      - the sequence paired with the Phred quality scores of that sequence. *)
  type thread

  (* To perform [Index.lookup]'s we need to expose the nucleotide sequence. *)
  val thread_to_seq : thread -> string

  (* To consider the reverse complement. *)
  val reverse_complement : thread -> thread

  (* Because [compute]ing the metric might be costly we expose a type to
     support an early stopping criterion. This split is modeled is also modeled
     via the [product] type. *)
  type stop_parameter

  (* [compute ~early_stop g thread position] Measure the metric for [thread]
     against graph [g] at [position] and optionally stopping because of
     [early_stop]. *)
  val compute : ?early_stop:stop_parameter -> Ref_graph.t -> thread ->
    Index.position -> m product

  (* Sometimes a read may have a valid lookup to multiple positions in the
     graph, in order to choose the best position we need to order based
     upon a metric. We need a similar way to choose when comparing whether a
     regular or the reverse complement orientation are better for a read.

     To this end we use [compare] to find the lowest (best) value of the
     metric across all alleles and then reuse that value across different
     positions or orientations. Note that this double use might not be ideal,
     but only metrics that have been computed and NOT stopped (via
     [stop_parameter]) are compared. *)
  val compare : m -> m -> int

  (* For paired end sequencing we want to compute a single metric by combining
     the result for both pairs. We force this reduction because if the user
     doesn't want to combine the results, they can always treat each read as
     a single (soloed) read. *)
  val combine_pairs : m Alleles.Map.t -> m Alleles.Map.t -> m Alleles.Map.t

end (* Single_config *)

type sequence_alignment_error =
  | NoIndexedPositions
  | AllComputesStopped of int             (* number of positions. *)
  | SequenceTooShort of too_short
  | InconsistentPairOrientation of bool   (* true = both regular
                                          , false = both complement *)
  | ToThread of string
  [@@deriving show]

type kmer_neighborhood =
  | UpTo of int             (* Scan from 0 up to, not including, boundary,
                               stop if we have good results. *)
  | AllTo of int            (* Add all neighbors, not including, boundary into
                               one scan. *)
  [@@deriving show]

(* Given a Single_config return functions that process a single (refered to as
   soloed from now on) or paired read. *)
module AgainstSequence (C : Single_config) = struct

  let metric_min m1 m2 =
    let r = C.compare m1 m2 in
    if r <= 0 then m1 else m2

  let min_el (_pos, ma) =
    Alleles.Map.(fold_wa ~init:(choose ma) ~f:metric_min ma)

  let reduce_across_positions s = function
    | [] -> s
    | ts -> List.fold_left ts ~init:(min_el s, s) ~f:(fun (b, s) m ->
              let bm = min_el m in
              if bm < b then (bm, m) else (b, s))
            |> snd

  let lookup ~distance idx thread =
    Index.lookup ~distance idx (C.thread_to_seq thread)

  let compute_and_keep_position =
    (* TODO: 1. Figure out better initialization number.
             2. Move to functoral interface.
             3. Currently ignores early_stop and g as hash parameters, but
                those shouldn't change between individual runs of a program
                like type, but caveat emptor.  *)
    let ht = Hashtbl.create 500 in
    fun ?early_stop g seq p ->
      let key = (p, seq) in
      if Hashtbl.mem ht key then
        p, Hashtbl.find ht key
      else
        let data = C.compute ?early_stop g seq p in
        Hashtbl.add ht ~key ~data;
        p, data

  let against_positions ?early_stop g seq ps =
    let res = List.map ps ~f:(compute_and_keep_position ?early_stop g seq) in
    let f = function
            | (p, Finished f) -> `Fst (p, f)
            | (p, Stopped s)  -> `Snd (p, s)
    in
    match List.partition_map res ~f with
    | [], []      -> Error NoIndexedPositions     (* assert false, call only on non-empty list! *)
    | [], als     -> Error (AllComputesStopped (List.length als))
    | r :: rs, _  -> Ok (reduce_across_positions r rs)

  (* Separate the lookup, aggregating positions, logic from the evaluating
     position logic. *)
  let lookup_positions ~distance idx seq k =
    match lookup ~distance idx seq with
    | Error e -> Error (SequenceTooShort e)
    | Ok ps   -> k ps

  let soloed_no_rc ?early_stop g seq = function
    | []  -> Error NoIndexedPositions
    | ps  -> against_positions ?early_stop g seq ps

  let regular_or_complement ~regular ~complement =
    let rm = min_el regular in
    let cm = min_el complement in
    (C.compare rm cm) <= 0

  let lookup_with_complement_positions ~distance idx seq rseq k =
    lookup_positions ~distance idx seq begin fun ps ->
      lookup_positions ~distance idx rseq begin fun rps ->
        k (ps, rps)
      end
    end

  let soloed_with_rc_choice ?early_stop g seq rseq = function
    | [], []  -> Error NoIndexedPositions
    | ps, []  -> against_positions ?early_stop g seq ps >>| fun e -> (e, `Reg)
    | [], rps -> against_positions ?early_stop g rseq rps >>| fun e -> (e, `Comp)
    | ps, rps ->
        let reg_e = against_positions ?early_stop g seq ps >>| fun e -> (e, `Reg) in
        let rcp_e = against_positions ?early_stop g rseq rps >>| fun e -> (e, `Comp) in
        begin
          match reg_e, rcp_e with
          | Error e,          Error _             -> reg_e
          | Ok _,             Error _             -> reg_e
          | Error _,          Ok _                -> rcp_e
          | Ok (regular, _),  Ok (complement, _)  ->
            if regular_or_complement ~regular ~complement then
              reg_e
            else
              rcp_e
        end

  let soloed_with_rc ?early_stop g seq rseq ps_pair =
    soloed_with_rc_choice ?early_stop g seq rseq ps_pair >>| fun (r, _) -> r

  let up_to max_distance ~f =
    let rec loop distance = function
        | Error NoIndexedPositions
        | Error (AllComputesStopped _) when distance < max_distance ->
            loop (distance + 1) (f ~distance)
        | res ->
            res
    in
    loop 1 (f ~distance:0)

  let all_to max_distance ~f ~k ~comb ~init =
    let rec loop distance acc =
      if distance >= max_distance then
        k acc
      else
        f ~distance (fun lst -> loop (distance + 1) (comb lst acc))
    in
    loop 0 init

  let unique_positions lst =
    List.dedup ~compare:Index.compare_position (List.concat lst)

  let soloed ?(check_rc=true) ?kmer_hood ?early_stop g idx seq =
    let upto_f distance =
      if check_rc then
        let rseq = C.reverse_complement seq in
        up_to distance
          ~f:(lookup_with_complement_positions idx seq rseq
                (soloed_with_rc ?early_stop g seq rseq))
      else
        up_to distance
          ~f:(lookup_positions idx seq (soloed_no_rc ?early_stop g seq))
    in
    match kmer_hood with
    | None                  -> upto_f 0
    | Some (UpTo distance)  -> upto_f distance
    | Some (AllTo distance) ->
        if check_rc then
          let rseq = C.reverse_complement seq in
          all_to distance ~f:(lookup_with_complement_positions idx seq rseq)
              ~k:(fun (acc, racc) ->
                    let uacc = unique_positions acc in
                    let uracc = unique_positions racc in
                    soloed_with_rc ?early_stop g seq rseq (uacc, uracc))
              ~comb:(fun (ps, rps) (psa, rpsa) -> (ps :: psa, rps :: rpsa))
              ~init:([], [])
        else
          all_to distance ~f:(lookup_positions idx seq)
            ~k:(fun acc ->
                  unique_positions acc
                  |> soloed_no_rc ?early_stop g seq)
            ~comb:(fun c a -> c :: a) ~init:[]

  let combine_pairs first_rc (p1, m1) (p2, m2) =
    (first_rc, p1, p2, C.combine_pairs m1 m2)

  let paired ?kmer_hood ?early_stop g idx seq1 seq2 =
    let f ~distance =
      let rseq1 = C.reverse_complement seq1 in
      lookup_with_complement_positions ~distance idx seq1 rseq1
        (soloed_with_rc_choice ?early_stop g seq1 rseq1) >>= begin fun (m1, c1) ->
        let rseq2 = C.reverse_complement seq2 in
        lookup_with_complement_positions ~distance idx seq2 rseq2
          (soloed_with_rc_choice ?early_stop g seq2 rseq2) >>= begin fun (m2, c2) ->
            match c1, c2 with
            | `Reg,  `Reg   -> Error (InconsistentPairOrientation true)
            | `Comp, `Comp  -> Error (InconsistentPairOrientation false)
            | `Reg,  `Comp  -> Ok (combine_pairs false m1 m2)
            | `Comp, `Reg   -> Ok (combine_pairs true m1 m2)
          end
        end
    in
    match kmer_hood with
    | None                  -> f ~distance:0
    | Some (UpTo distance)  -> up_to distance ~f
    | Some (AllTo distance) ->
      let rseq1 = C.reverse_complement seq1 in
      let rseq2 = C.reverse_complement seq2 in
      all_to distance
        ~f:(fun ~distance k ->
              lookup_with_complement_positions ~distance idx seq1 rseq1 (fun (ps1, rps1) ->
                lookup_with_complement_positions ~distance idx seq2 rseq2 (fun (ps2, rps2) ->
                  k (ps1, rps1, ps2, rps2))))
        ~k:(fun (acc1, racc1, acc2, racc2) ->
              let uacc1 = unique_positions acc1 in
              let uracc1 = unique_positions racc1 in
              let uacc2 = unique_positions acc2 in
              let uracc2 = unique_positions racc2 in
              soloed_with_rc_choice ?early_stop g seq1 rseq1 (uacc1, uracc1) >>=
                begin fun (m1, c1) ->
                  soloed_with_rc_choice ?early_stop g seq2 rseq2 (uacc2, uracc2) >>=
                    begin fun (m2, c2) ->
                      match c1, c2 with
                      | `Reg,  `Reg   -> Error (InconsistentPairOrientation true)
                      | `Comp, `Comp  -> Error (InconsistentPairOrientation false)
                      | `Reg,  `Comp  -> Ok (combine_pairs false m1 m2)
                      | `Comp, `Reg   -> Ok (combine_pairs true m1 m2)
                    end
                end)
        ~comb:(fun (a,b,c,d) (ac,bc,cc,dc) -> (a::ac, b::bc, c::cc, d::dc))
        ~init:([],[],[],[])

  type stop_parameter = C.stop_parameter
  type thread = C.thread

end (* AgainstSequence *)

(* Configurations of different metrics to measure against a soloed or
   paired read. *)
module ThreadIsJustSequence = struct
  let thread_to_seq s = s
  let reverse_complement = Util.reverse_complement
end

(* Keep track of all the positions where we have mismatches. *)
module List_mismatches_config = struct

  include Alignment.Position_mismatches
  include ThreadIsJustSequence

  let sum_mismatches =
    List.fold_left ~init:0 ~f:(fun s (_, v) -> s + v)

  (* Lowest mismatch across all alleles. *)
  let compare l1 l2 =
    (sum_mismatches l1) - (sum_mismatches l2)

  let combine_pairs =
    Alleles.Map.map2_wa ~f:(@)

end (* List_mismatches_config *)

(* Count the number of mismatches when the read is aligned against the graph. *)
module Count_mismatches_config = struct

  include Alignment.Mismatches
  include ThreadIsJustSequence

  (* Lowest mismatch across all alleles. *)
  let compare m1 m2 =
    m1 - m2

  let combine_pairs =
    Alleles.Map.map2_wa ~f:(+)

end (* Count_mismatches_config *)

(* Compute the 'Phred' likelihood of the read aligning aginst the graph. *)
module Phred_likelihood_config = struct

  include Alignment.Phred_likelihood

  let reverse_complement (s, a) =
    let n = Array.length a in
    Util.reverse_complement s
    , Array.init n ~f:(fun i -> a.(n - i - 1))

  let thread_to_seq (s, _) = s

  (* There are 2 ways to compute Best in this case:
     1. lowest number of mismatches as in the other mismatch counting algorithms
     2. highest sum of log likelihoods -> highest probability
      We use #2 and reverse the notion of comparison so that highest likelihoods
      are folded to the bottom. *)
  let compare p1 p2 =
    let open Alignment in (* TODO: Can this be simplified? *)
    if p1.sum_llhd < p2.sum_llhd then
      1
    else if p1.sum_llhd > p2.sum_llhd then
      -1
    else
      0

  let combine_pairs =
    let open Alignment in
    Alleles.Map.map2_wa ~f:(fun m1 m2 ->
      { mismatches = m1.mismatches +. m2.mismatches
      ; sum_llhd   = m1.sum_llhd +. m2.sum_llhd
      })

end (* Phred_likelihood_config *)

module List_mismatches_of_sequence = AgainstSequence (List_mismatches_config)
module Count_mismatches_of_sequence = AgainstSequence (Count_mismatches_config)
module Phred_likelihood_of_sequence = AgainstSequence (Phred_likelihood_config)

(* How we compute a metric across a set of reads .*)
module type Multiple_config = sig

  (* The same types from the Single_config, how we measure for one read.*)
  type stop_parameter
  type thread

  val to_thread : Biocaml_unix.Fastq.item -> (thread, string) result

  (* The metric for a read (or read pair) across a set of alleles. This type
     can frequently be the same as Single_config.m, but doesn't need to be
     there can be instances where we want to mutate between the two; this is
     accomplished be how [map], [map_paired] work with [AgainstSequence.soloed]
     and [AgainstSequence.paired]. *)
  type me

  (* Finally, a metric across all (soloed or paired) reads *)
  type re
  val empty : re

  val map :
    ?check_rc:bool ->
    ?kmer_hood:kmer_neighborhood ->
    ?early_stop:stop_parameter ->
    Ref_graph.t -> Index.t ->
    thread ->
    (Index.position * me Alleles.Map.t, sequence_alignment_error) result

  val map_paired :
    ?kmer_hood:kmer_neighborhood ->
    ?early_stop:stop_parameter ->
    Ref_graph.t -> Index.t ->
    thread ->
    thread ->
    (bool * Index.position * Index.position * me Alleles.Map.t, sequence_alignment_error) result

  (* Aggregate the mapped [me]'s into a single [re]. *)
  val reduce : me -> re -> re

end (* Multiple_config *)

(* TODO: Consider adding early_stop to this data structure if it proves that
  we end up never varying this parameter. We haven't developed a different
  (better?) stopping logic then when the min mismatches goes above some
  boundary. *)

type fastq_fold_args =
  { number_of_reads : int option
  ; specific_reads  : string list
  ; check_rc        : bool option
  ; kmer_hood       : kmer_neighborhood option
  } [@@deriving show]

type 'a paired_datum =
  { name      : string
  ; read1     : string
  ; first_rc  : bool              (* First or second is reverse complement *)
  ; read2     : string
  ; position1 : Index.position
  ; position2 : Index.position
  ; amap      : 'a Alleles.Map.t
  }

module Multiple (C : Multiple_config) = struct

  let fold_over_fastq fastq_file { number_of_reads; specific_reads; check_rc; kmer_hood }
      ?early_stop g idx =
    let amap = Alleles.Map.make g.Ref_graph.aindex C.empty in
    let f (errors, amap) fqi =
      match C.to_thread fqi with
      | Error e     -> ((ToThread e, fqi) :: errors), amap
      | Ok seq  ->
          match C.map ?kmer_hood ?check_rc ?early_stop g idx seq with
          | Error e     -> ((e, fqi) :: errors), amap
          | Ok (_p, a)  -> Alleles.Map.update2 ~source:a ~dest:amap C.reduce;
                           errors, amap
    in
    Fastq.fold ?number_of_reads ~specific_reads ~init:([], amap) ~f fastq_file

  let fold_over_paired file1 file2 { number_of_reads; specific_reads; kmer_hood; _}
      ?early_stop g idx =
    let amap = Alleles.Map.make g.Ref_graph.aindex C.empty in
    let f (errors, amap) fqi1 fqi2 =
      match C.to_thread fqi1 with
      | Error e     -> ((ToThread e, fqi1) :: errors), amap
      | Ok seq1     ->
          match C.to_thread fqi2 with
          | Error e -> ((ToThread e, fqi2) :: errors), amap
          | Ok seq2 ->
            match C.map_paired ?kmer_hood ?early_stop g idx seq1 seq2 with
            | Error e               -> ((e, fqi1) :: errors), amap      (* TODO: return both fqi1 and fqi2 ? *)
            | Ok (_rc, _p1, _p2, a) -> Alleles.Map.update2 ~source:a ~dest:amap C.reduce;
                                       errors, amap
    in
    let res =
      Fastq.fold_paired ?number_of_reads ~specific_reads ~init:([], amap) ~f file1 file2
        (fun fqi -> fqi.Biocaml_unix.Fastq.name)
    in
    match res with
    | `StoppedByFilter (errors, amap)            -> errors, amap
    | `DesiredReads (errors, amap)               -> errors, amap
    | `BothFinished (errors, amap)               -> errors, amap
    (* Continue typing on the soloed? *)
    | `OneReadPairedFinished (_, (errors, amap)) -> errors, amap

  let map_over_paired file1 file2 { number_of_reads; specific_reads; kmer_hood; _}
    ?early_stop g idx =
    let f (errors, amap_accum) fqi1 fqi2 =
      match C.to_thread fqi1 with
      | Error e     -> ((ToThread e, fqi1) :: errors), amap_accum
      | Ok seq1     ->
          match C.to_thread fqi2 with
          | Error e -> ((ToThread e, fqi2) :: errors), amap_accum
          | Ok seq2 ->
            match C.map_paired ?kmer_hood ?early_stop g idx seq1 seq2 with
            | Error e             -> ((e, fqi1) :: errors), amap_accum
            | Ok (rc, p1, p2, a)  -> errors,
                                    ( { name      = fqi1.Biocaml_unix.Fastq.name
                                      ; read1     = fqi1.Biocaml_unix.Fastq.sequence
                                      ; first_rc  = rc
                                      ; read2     = fqi2.Biocaml_unix.Fastq.sequence
                                      ; position1 = p1
                                      ; position2 = p2
                                      ; amap      = a
                                      } :: amap_accum)
    in
    let res =
      Fastq.fold_paired ?number_of_reads ~specific_reads ~init:([], []) ~f file1 file2
        (fun fqi -> fqi.Biocaml_unix.Fastq.name)
    in
    match res with
    | `StoppedByFilter (errors, acc)            -> errors, acc
    | `DesiredReads (errors, acc)               -> errors, acc
    | `BothFinished (errors, acc)               -> errors, acc
    (* Continue typing on the soloed? *)
    | `OneReadPairedFinished (_, (errors, acc)) -> errors, acc

end (* Multiple *)

module List_mismatches_of_reads = Multiple (struct
  include List_mismatches_of_sequence

  type me = (int * int) list    (* Same as m *)
  type re = me list
  let empty = []

  let to_thread fqi = Ok fqi.Biocaml_unix.Fastq.sequence
  let map = soloed
  let map_paired = paired
  let reduce v l = v :: l
end)  (* List_mismatches_of_reads *)

module Number_of_mismatches_of_reads = Multiple (struct

  include Count_mismatches_of_sequence
  type me = int                 (* Same as m *)
  type re = int
  let empty = 0

  let to_thread fqi = Ok fqi.Biocaml_unix.Fastq.sequence
  let map = soloed
  let map_paired = paired
  let reduce = (+)
end)  (* Number_of_mismatches_of_reads *)

module Likelihood_of_multiple_config = struct

  include Count_mismatches_of_sequence

  type me = float   (* Not the same as m, # mismatches -> likelihood *)
  type re = float

  let to_thread fqi = Ok fqi.Biocaml_unix.Fastq.sequence
  let map likelihood ?check_rc ?kmer_hood ?early_stop g idx th =
    let len = String.length th in
    soloed ?check_rc ?kmer_hood ?early_stop g idx th >>= fun (p, m) ->
      Ok (p, Alleles.Map.map_wa ~f:(likelihood ~len) m)

  let map_paired likelihood ?kmer_hood ?early_stop g idx th1 th2 =
    let len = String.length th1 in
    paired ?kmer_hood ?early_stop g idx th1 th2 >>= fun (rc, p1, p2, m) ->
      Ok (rc, p1, p2, Alleles.Map.map_wa ~f:(likelihood ~len) m)

  let reduce = (+.)

end (* Likelihood_of_multiple_config *)

module Phred_likelihood_of_reads = Multiple (struct

  include Phred_likelihood_of_sequence

  type me = float
  type re = float

  let empty = 0.

  let to_thread fqi =
    let module CE = Core_kernel.Error in
    let module CR = Core_kernel.Std.Result in
    match Fastq.phred_probabilities fqi.Biocaml_unix.Fastq.qualities with
    | CR.Error e -> Error (CE.to_string_hum e)
    | CR.Ok qarr -> Ok (fqi.Biocaml_unix.Fastq.sequence, qarr)
    (*match Fastq.phred_log_probs fqi.Biocaml_unix.Fastq.qualities with
    | CR.Error e -> Error (CE.to_string_hum e)
    | CR.Ok qarr -> Ok (fqi.Biocaml_unix.Fastq.sequence, qarr) *)

  let map ?check_rc ?kmer_hood ?early_stop g idx th =
    soloed ?kmer_hood ?check_rc ?early_stop g idx th >>= fun (p1, amap) ->
      Ok (p1, Alleles.Map.map_wa amap ~f:(fun pt -> pt.Alignment.sum_llhd))

  let map_paired ?kmer_hood ?early_stop g idx th1 th2 =
    paired ?kmer_hood ?early_stop g idx th1 th2 >>= fun (rc, p1, p2, amap) ->
      Ok (rc, p1, p2, Alleles.Map.map_wa amap ~f:(fun pt -> pt.Alignment.sum_llhd))

  (* log likelihoods -> + *)
  let reduce = (+.)

end) (* Phred_likelihood_of_reads *)


let type_soloed ?filter ?(as_=`Phred_likelihood_of_reads) fold_args g idx ~fastq_file =
  let early_stop =
    Option.map filter ~f:(fun n -> Ref_graph.number_of_alleles g, float n)
  in
  match as_ with
  | `List_mismatches_of_reads       ->
      `List_mismatches_of_reads (List_mismatches_of_reads.fold_over_fastq
          fastq_file fold_args ?early_stop g idx)
  | `Number_of_mismatches_of_reads  ->
      `Number_of_mismatches_of_reads (Number_of_mismatches_of_reads.fold_over_fastq
          fastq_file fold_args ?early_stop g idx)
  | `Likelihood_of_reads error      ->
      (* These modules are created dynamically because they are parameterized
         the error that we assign to matching a base incorrectly. *)
      let module Likelihood_of_reads = Multiple (struct
          include Likelihood_of_multiple_config
          let empty = 1.
          let map = map (fun ~len m -> likelihood ~er:error ~len (float m))
          let map_paired = map_paired (fun ~len m -> likelihood ~er:error ~len (float m))
          let reduce l a = a *. l
        end) (* Likelihood_of_reads *)
      in
      `Likelihood_of_reads (Likelihood_of_reads.fold_over_fastq
          fastq_file fold_args ?early_stop g idx)
  | `LogLikelihood_of_reads error   ->
      let module LogLikelihood_of_reads = Multiple (struct
          include Likelihood_of_multiple_config
          let empty = 0.
          let map = map (fun ~len m -> log_likelihood ~er:error ~len (float m))
          let map_paired = map_paired (fun ~len m -> log_likelihood ~er:error ~len (float m))
          let reduce l a = a +. l
        end)  (* LogLikelihood_of_reads *)
      in
      `LogLikelihood_of_reads (LogLikelihood_of_reads.fold_over_fastq
          fastq_file fold_args ?early_stop g idx)
  | `Phred_likelihood_of_reads      ->
      `Phred_likelihood_of_reads (Phred_likelihood_of_reads.fold_over_fastq
          fastq_file fold_args ?early_stop g idx)

let type_paired ?filter ?(as_=`Phred_likelihood_of_reads) fold_args g idx f1 f2 =
  let early_stop =
    Option.map filter ~f:(fun n -> Ref_graph.number_of_alleles g, float n)
  in
  match as_ with
  | `List_mismatches_of_reads       ->
      `List_mismatches_of_reads (List_mismatches_of_reads.fold_over_paired
          f1 f2 fold_args ?early_stop g idx)
  | `Number_of_mismatches_of_reads  ->
      `Number_of_mismatches_of_reads (Number_of_mismatches_of_reads.fold_over_paired
          f1 f2 fold_args ?early_stop g idx)
  | `Likelihood_of_reads error      ->
      let module Likelihood_of_reads = Multiple (struct
          include Likelihood_of_multiple_config
          let empty = 1.
          let map = map (fun ~len m -> likelihood ~er:error ~len (float m))
          let map_paired = map_paired (fun ~len m -> likelihood ~er:error ~len (float m))
          let reduce l a = a *. l
        end) (* Likelihood_of_reads *)
      in
      `Likelihood_of_reads (Likelihood_of_reads.fold_over_paired
          f1 f2 fold_args ?early_stop g idx)
  | `LogLikelihood_of_reads error   ->
      let module LogLikelihood_of_reads = Multiple (struct
          include Likelihood_of_multiple_config
          let empty = 0.
          let map = map (fun ~len m -> log_likelihood ~er:error ~len (float m))
          let map_paired = map_paired (fun ~len m -> log_likelihood ~er:error ~len (float m))
          let reduce l a = a +. l
        end)  (* LogLikelihood_of_reads *)
      in
      `LogLikelihood_of_reads (LogLikelihood_of_reads.fold_over_paired
          f1 f2 fold_args ?early_stop g idx)
  | `Phred_likelihood_of_reads      ->
      `Phred_likelihood_of_reads (Phred_likelihood_of_reads.fold_over_paired
          f1 f2 fold_args ?early_stop g idx)

let type_map_paired ?filter ?(as_=`Phred_likelihood_of_reads) fold_args g idx f1 f2 =
  let early_stop =
    Option.map filter ~f:(fun n -> Ref_graph.number_of_alleles g, float n)
  in
  match as_ with
  | `List_mismatches_of_reads       ->
      `List_mismatches_of_reads (List_mismatches_of_reads.map_over_paired
          f1 f2 fold_args ?early_stop g idx)
  | `Number_of_mismatches_of_reads  ->
      `Number_of_mismatches_of_reads (Number_of_mismatches_of_reads.map_over_paired
          f1 f2 fold_args ?early_stop g idx)
  | `Likelihood_of_reads error      ->
      let module Likelihood_of_reads = Multiple (struct
          include Likelihood_of_multiple_config
          let empty = 1.
          let map = map (fun ~len m -> likelihood ~er:error ~len (float m))
          let map_paired = map_paired (fun ~len m -> likelihood ~er:error ~len (float m))
          let reduce l a = a *. l
        end) (* Likelihood_of_reads *)
      in
      `Likelihood_of_reads (Likelihood_of_reads.map_over_paired
          f1 f2 fold_args ?early_stop g idx)
  | `LogLikelihood_of_reads error   ->
      let module LogLikelihood_of_reads = Multiple (struct
          include Likelihood_of_multiple_config
          let empty = 0.
          let map = map (fun ~len m -> log_likelihood ~er:error ~len (float m))
          let map_paired = map_paired (fun ~len m -> log_likelihood ~er:error ~len (float m))
          let reduce l a = a +. l
        end)  (* LogLikelihood_of_reads *)
      in
      `LogLikelihood_of_reads (LogLikelihood_of_reads.map_over_paired
          f1 f2 fold_args ?early_stop g idx)
  | `Phred_likelihood_of_reads      ->
      `Phred_likelihood_of_reads (Phred_likelihood_of_reads.map_over_paired
          f1 f2 fold_args ?early_stop g idx)
