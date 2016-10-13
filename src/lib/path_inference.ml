(** Path Inference : Given a read (a sequence of characters return a PDF over
    the edges (alleles) in the graph. *)

open Util

let log_likelihood ?(alph_size=4) ?(er=0.01) ~len mismatches =
  let lmp = log (er /. (float (alph_size - 1))) in
  let lcp = log (1. -. er) in
  let c = (float len) -. mismatches in
  c *. lcp +. mismatches *. lmp

let likelihood ?alph_size ?er ~len m =
  exp (log_likelihood ?alph_size ?er ~len m)

let report_mismatches el_to_string state aindex mm =
  printf "reporting %s mismatches\n" state;
  let () =
    Alleles.Map.values_assoc aindex mm
    |> List.sort ~cmp:(fun (i1,_) (i2,_) -> compare i1 i2)
    |> List.iter ~f:(fun (w, a) ->
        printf "%s\t%s\n"
          (el_to_string w)
          (insert_chars ['\t'; '\t'; '\n']
            (Alleles.Set.to_human_readable aindex ~max_length:10000 ~complement:`No a)))
  in
  printf "finished reporting mismatches\n%!"

module type Single_config = sig

  type t  (* How we measure the alignment of a sequence against an allele. *)

  val to_string : t -> string           (* Display *)

  type stop_parameter
  type thread

  val thread_to_seq : thread -> string
  val reverse_complement : thread -> thread

  val compute : ?early_stop:stop_parameter -> Ref_graph.t -> thread -> Index.position ->
    ([ `Finished of t Alleles.Map.t | `Stopped of t Alleles.Map.t ], string) result

  val combine_pairs : t Alleles.Map.t -> t Alleles.Map.t -> t Alleles.Map.t

  val reduce_across_positions : t Alleles.Map.t -> t Alleles.Map.t list -> t Alleles.Map.t

end (* Single_config *)

type sequence_alignment_error =
  | NoPositions
  | AllStopped of int
  | Other of string
  | ToThread of string

module AgainstSequence (C : Single_config) = struct

  let against_positions ?early_stop g ps seq =
    let rec loop acc = function
      | []     -> `Fin acc
      | h :: t ->
          match C.compute ?early_stop g seq h with
          | Error e -> `Errored e
          | Ok r    -> loop (r :: acc) t
    in
    match loop [] ps with
    | `Errored e -> Error (Other e)
    | `Fin res   ->
        let f = function | `Finished f -> `Fst f | `Stopped s -> `Snd s in
        match List.partition_map res ~f with
        | [], []      -> Error NoPositions
        | [], als     -> Error (AllStopped (List.length als))
        | r :: rs, _  -> Ok (C.reduce_across_positions r rs)

  let single ?(check_rc=true) ?early_stop g idx seq =
    match Index.lookup idx (C.thread_to_seq seq) with
    | Error e -> Error (Other e)
    | Ok []   ->
        if check_rc then
          let revt = C.reverse_complement seq in
            match Index.lookup idx (C.thread_to_seq revt) with
            | Error e -> Error (Other e)
            | Ok ps   -> against_positions ?early_stop g ps revt
        else
          Error NoPositions
    | Ok ps   -> against_positions ?early_stop g ps seq

  let paired ?early_stop g idx seq1 seq2 =
    let where t = Index.lookup idx (C.thread_to_seq t) in
    let ap = against_positions ?early_stop g in
    match Index.lookup idx (C.thread_to_seq seq1) with
    | Error e -> Error (Other e)
    | Ok []   ->
        let rev1 = C.reverse_complement seq1 in
        begin match where rev1, where seq2 with
        | Error e, _
        | _,      Error e -> Error (Other e)
        | Ok ps1, Ok ps2  ->
            match ap ps1 rev1, ap ps2 seq2 with
            | Ok r, Ok l -> Ok (C.combine_pairs r l)
            | Ok o, _
            | _,    Ok o -> Ok o
            | Error e, _ -> Error e
        end
    | Ok ps ->
        let rev2 = C.reverse_complement seq2 in
        begin match where seq1, where rev2 with
        | Error e, _
        | _,      Error e -> Error (Other e)
        | Ok ps1, Ok ps2  ->
            match ap ps1 seq1, ap ps2 rev2 with
            | Ok r, Ok l -> Ok (C.combine_pairs r l)
            | Ok o, _
            | _,    Ok o -> Ok o
            | Error e, _ -> Error e
        end

  type stop_parameter = C.stop_parameter
  type thread = C.thread

end (* AgainstSequence *)

module ThreadIsJustSequence = struct
  type thread = string
  let thread_to_seq s = s
  let reverse_complement = Util.reverse_complement
end

module ListMismatches_config = struct

  type t = (int * int) list
  let to_string l =
    String.concat ~sep:"; "
      (List.map l ~f:(fun (p,v) -> sprintf "(%d,%d)" p v))

  type stop_parameter = int * float

  include ThreadIsJustSequence

  let compute = Alignment.compute_mismatches_lst

  (* Lowest mismatch across all alleles. *)
  let to_min =
    Alleles.Map.fold_wa ~init:max_int ~f:(fun mx lst ->
      let s = List.fold_left lst ~init:0 ~f:(fun s (_, v) -> s + v) in
      min mx s)

  let reduce_across_positions s = function
    | [] -> s
    | ts ->
        let b = to_min s in
        List.fold_left ts ~init:(b, s) ~f:(fun (b, s) m ->
          let bm = to_min m in
          if bm < b then
            (bm, m)
          else
            (b, s))
        |> snd

  let combine_pairs = Alleles.Map.map2_wa ~f:(@)

end (* ListMismatches_config *)

module ListMismatches = AgainstSequence (ListMismatches_config)

module SequenceMismatches = AgainstSequence (struct

  type t = int
  let to_string = sprintf "%d"

  type stop_parameter = int * float

  include ThreadIsJustSequence

  let compute = Alignment.compute_mismatches

  (* Lowest mismatch across all alleles. *)
  let to_min = Alleles.Map.fold_wa ~init:max_int ~f:min
  let reduce_across_positions s = function
    | [] -> s
    | ts ->
        let b = to_min s in
        List.fold_left ts ~init:(b, s) ~f:(fun (b, s) m ->
          let bm = to_min m in
          if bm < b then
            (bm, m)
          else
            (b, s))
        |> snd


  let combine_pairs = Alleles.Map.map2_wa ~f:(+)

end)

module PhredLlhdMismatches = AgainstSequence ( struct

  open Alignment

  type t = PhredLikelihood_config.t
  let to_string = PhredLikelihood_config.to_string

  type stop_parameter = int * float
  type thread = string * float array

  let reverse_complement (s, a) =
    let n = Array.length a in
    Util.reverse_complement s
    , Array.init n ~f:(fun i -> a.(n - i - 1))

  let thread_to_seq (s, _) = s

  let compute = Alignment.compute_plhd

  (* There are 2 ways to compute Best in this case:
     1. lowest number of mismatches as in the other mismatch counting algorithms
     2. highest sum of log likelihoods -> highest probability
     we use 2. *)
  let to_max =
    Alleles.Map.fold_wa ~init:neg_infinity
      ~f:(fun a t -> max a t.PhredLikelihood_config.sum_llhd)

  let reduce_across_positions s = function
    | [] -> s
    | ts ->
        let b = to_max s in
        List.fold_left ts ~init:(b, s) ~f:(fun (b, s) m ->
          let bm = to_max m in
          if bm > b then
            (bm, m)
          else
            (b, s))
        |> snd

  let combine_pairs =
    let open PhredLikelihood_config in
    Alleles.Map.map2_wa ~f:(fun m1 m2 ->
      { mismatches = m1.mismatches +. m2.mismatches
      ; sum_llhd   = m1.sum_llhd +. m2.sum_llhd
      })

end) (* PhredLlhdMismatches *)

module type Multiple_config = sig

  type mp     (* map step *)
  type re     (* reduce across alleles *)

  val empty : re

  type stop_parameter
  type thread

  val to_thread : Biocaml_unix.Fastq.item -> (thread, string) result

  val map :
    ?check_rc:bool ->
    ?early_stop:stop_parameter ->
    Ref_graph.t -> Index.t ->
    thread ->
    (mp Alleles.Map.t, sequence_alignment_error) result

  val map_paired :
    ?early_stop:stop_parameter ->
    Ref_graph.t -> Index.t ->
    thread ->
    thread ->
    (mp Alleles.Map.t, sequence_alignment_error) result

  val reduce : mp -> re -> re

end (* Multiple_config *)

module Multiple (C : Multiple_config) = struct

  let fold_over_fastq ?number_of_reads fastq_file ?check_rc ?early_stop g idx =
    let amap = Alleles.Map.make g.Ref_graph.aindex C.empty in
    let f (errors, amap) fqi =
      match C.to_thread fqi with
      | Error e     -> ((ToThread e, fqi) :: errors), amap
      | Ok seq  ->
          match C.map ?check_rc ?early_stop g idx seq with
          | Error e -> ((e, fqi) :: errors), amap
          | Ok a    -> Alleles.Map.update2 ~source:a ~dest:amap C.reduce;
                       errors, amap
    in
    Fastq.fold ?number_of_reads ~init:([], amap) ~f fastq_file

  let map_over_fastq ?number_of_reads fastq_file ?check_rc ?early_stop g idx =
    let f fqi =
      match C.to_thread fqi with
      | Error e     -> Error (ToThread e, fqi)
      | Ok seq  ->
          match C.map ?check_rc ?early_stop g idx seq with
          | Error e -> Error (e, fqi)
          | Ok a    -> Ok a
    in
    Fastq.fold ?number_of_reads ~init:[] fastq_file
      ~f:(fun acc fqi -> (f fqi) :: acc)
    |> List.rev

  let fold_over_paired ?number_of_reads file1 file2 ?early_stop g idx =
    let amap = Alleles.Map.make g.Ref_graph.aindex C.empty in
    let f (errors, amap) fqi1 fqi2 =
      match C.to_thread fqi1 with
      | Error e     -> ((ToThread e, fqi1) :: errors), amap
      | Ok seq1     ->
          match C.to_thread fqi2 with
          | Error e -> ((ToThread e, fqi2) :: errors), amap
          | Ok seq2 ->
            match C.map_paired ?early_stop g idx seq1 seq2 with
            | Error e -> ((e, fqi1) :: errors), amap      (* TODO: return both fqi1 and fqi2 ? *)
            | Ok a    -> Alleles.Map.update2 ~source:a ~dest:amap C.reduce;
                         errors, amap
    in
    let res =
      Fastq.fold_paired ?number_of_reads ~init:([], amap) ~f file1 file2
        (fun fqi -> fqi.Biocaml_unix.Fastq.name)
    in
    match res with
    | `DesiredReads (errors, amap)               -> errors, amap
    | `BothFinished (errors, amap)               -> errors, amap
    (* TODO, go to single. *)
    | `OneReadPairedFinished (_, (errors, amap)) -> errors, amap

end (* Multiple *)

(** Typing. *)
let sequence_alignment_error_to_string = function
  | NoPositions  -> "No positions found"
  | AllStopped n -> sprintf "Stopped %d positions" n
  | Other m      -> sprintf "Error: %s" m
  | ToThread m   -> sprintf "ToThread: %s" m

module MismatchesList = Multiple (struct
  type mp = (int * int) list
  type re = (int * int) list list
  let empty = []

  type stop_parameter = ListMismatches.stop_parameter
  type thread = ListMismatches.thread

  let to_thread fqi = Ok fqi.Biocaml_unix.Fastq.sequence
  let map = ListMismatches.single
  let map_paired = ListMismatches.paired
  let reduce v l = v :: l
end)

module Mismatches = Multiple (struct
  type mp = int
  type re = int
  let empty = 0

  type stop_parameter = SequenceMismatches.stop_parameter
  type thread = SequenceMismatches.thread

  let to_thread fqi = Ok fqi.Biocaml_unix.Fastq.sequence
  let map = SequenceMismatches.single
  let map_paired = SequenceMismatches.paired
  let reduce = (+)
end)

module Llhd_config = struct
  type mp = float
  type re = float

  type stop_parameter = SequenceMismatches.stop_parameter
  type thread = SequenceMismatches.thread

  let to_thread fqi = Ok fqi.Biocaml_unix.Fastq.sequence
  let map l ?check_rc ?early_stop g idx th =
    let len = String.length th in
    SequenceMismatches.single ?check_rc ?early_stop g idx th >>= fun m ->
      Ok (Alleles.Map.map_wa ~f:(l ~len) m)

  let map_paired l ?early_stop g idx th1 th2 =
    let len = String.length th1 in
    SequenceMismatches.paired ?early_stop g idx th1 th2 >>= fun m ->
      Ok (Alleles.Map.map_wa ~f:(l ~len) m)

end (* Llhd_config *)

module Phred_lhd = Multiple (struct
  type mp = float
  type re = float

  let empty = 0.

  type stop_parameter = PhredLlhdMismatches.stop_parameter
  type thread = PhredLlhdMismatches.thread

  let to_thread fqi =
    let module CE = Core_kernel.Error in
    let module CR = Core_kernel.Std.Result in
    match Fastq.phred_probabilities fqi.Biocaml_unix.Fastq.qualities with
    | CR.Error e -> Error (CE.to_string_hum e)
    | CR.Ok qarr -> Ok (fqi.Biocaml_unix.Fastq.sequence, qarr)

  let map ?check_rc ?early_stop g idx th =
    let open Alignment in
    PhredLlhdMismatches.single ?check_rc ?early_stop g idx th >>= fun amap ->
      Ok (Alleles.Map.map_wa amap ~f:(fun pt -> pt.PhredLikelihood_config.sum_llhd))

  let map_paired ?early_stop g idx th1 th2 =
    let open Alignment in
    PhredLlhdMismatches.paired ?early_stop g idx th1 th2 >>= fun amap ->
      Ok (Alleles.Map.map_wa amap ~f:(fun pt -> pt.PhredLikelihood_config.sum_llhd))

  let reduce = (+.)

end) (* Phred_lhd *)


let map ?filter ?(as_=`PhredLikelihood) g idx ?number_of_reads ~fastq_file =
  let early_stop =
    Option.map filter ~f:(fun n -> Ref_graph.number_of_alleles g, float n)
  in
  match as_ with
  | `MismatchesList       ->
      `MismatchesList (MismatchesList.map_over_fastq
          ?number_of_reads fastq_file ?early_stop g idx)
  | `Mismatches           ->
      `Mismatches (Mismatches.map_over_fastq
          ?number_of_reads fastq_file ?early_stop g idx)
  | `Likelihood error     ->
      let module Likelihood = Multiple (struct
        include Llhd_config
        let empty = 1.
        let map = map (fun ~len m -> likelihood ~er:error ~len (float m))
        let map_paired = map_paired (fun ~len m -> likelihood ~er:error ~len (float m))
        let reduce l a = a *. l
      end) in
      `Likelihood (Likelihood.map_over_fastq
          ?number_of_reads fastq_file ?early_stop g idx)
  | `LogLikelihood error  ->
      let module LogLikelihood = Multiple (struct
        include Llhd_config
        let empty = 0.
        let map = map (fun ~len m -> log_likelihood ~er:error ~len (float m))
        let map_paired = map_paired (fun ~len m -> log_likelihood ~er:error ~len (float m))
        let reduce l a = a +. l
      end) in
     `LogLikelihood (LogLikelihood.map_over_fastq
          ?number_of_reads fastq_file ?early_stop g idx)
  | `PhredLikelihood ->
      `PhredLikelihood (Phred_lhd.map_over_fastq
          ?number_of_reads fastq_file ?early_stop g idx)

let type_ ?filter ?check_rc ?(as_=`PhredLikelihood) g idx ?number_of_reads ~fastq_file =
  let early_stop =
    Option.map filter ~f:(fun n -> Ref_graph.number_of_alleles g, float n)
  in
  match as_ with
  | `MismatchesList       ->
      `MismatchesList (MismatchesList.fold_over_fastq
          ?number_of_reads fastq_file ?check_rc ?early_stop g idx)
  | `Mismatches           ->
      `Mismatches (Mismatches.fold_over_fastq
          ?number_of_reads fastq_file ?check_rc ?early_stop g idx)
  | `Likelihood error     ->
      let module Likelihood = Multiple (struct
        include Llhd_config
        let empty = 1.
        let map = map (fun ~len m -> likelihood ~er:error ~len (float m))
        let map_paired = map_paired (fun ~len m -> likelihood ~er:error ~len (float m))
        let reduce l a = a *. l
      end) in
      `Likelihood (Likelihood.fold_over_fastq
          ?number_of_reads fastq_file ?check_rc ?early_stop g idx)
  | `LogLikelihood error  ->
      let module LogLikelihood = Multiple (struct
        include Llhd_config
        let empty = 0.
        let map = map (fun ~len m -> log_likelihood ~er:error ~len (float m))
        let map_paired = map_paired (fun ~len m -> log_likelihood ~er:error ~len (float m))
        let reduce l a = a +. l
      end) in
      `LogLikelihood (LogLikelihood.fold_over_fastq
          ?number_of_reads fastq_file ?early_stop g idx)
  | `PhredLikelihood ->
      `PhredLikelihood (Phred_lhd.fold_over_fastq
          ?number_of_reads fastq_file ?check_rc ?early_stop g idx)

let type_paired ?filter ?(as_=`PhredLikelihood) g idx ?number_of_reads f1 f2 =
  let early_stop =
    Option.map filter ~f:(fun n -> Ref_graph.number_of_alleles g, float n)
  in
  match as_ with
  | `MismatchesList       ->
      `MismatchesList (MismatchesList.fold_over_paired
          ?number_of_reads f1 f2 ?early_stop g idx)
  | `Mismatches           ->
      `Mismatches (Mismatches.fold_over_paired
          ?number_of_reads f1 f2 ?early_stop g idx)
  | `Likelihood error     ->
      let module Likelihood = Multiple (struct
        include Llhd_config
        let empty = 1.
        let map = map (fun ~len m -> likelihood ~er:error ~len (float m))
        let map_paired = map_paired (fun ~len m -> likelihood ~er:error ~len (float m))
        let reduce l a = a *. l
      end) in
      `Likelihood (Likelihood.fold_over_paired
          ?number_of_reads f1 f2 ?early_stop g idx)
  | `LogLikelihood error  ->
      let module LogLikelihood = Multiple (struct
        include Llhd_config
        let empty = 0.
        let map = map (fun ~len m -> log_likelihood ~er:error ~len (float m))
        let map_paired = map_paired (fun ~len m -> log_likelihood ~er:error ~len (float m))
        let reduce l a = a +. l
      end) in
      `LogLikelihood (LogLikelihood.fold_over_paired
          ?number_of_reads f1 f2 ?early_stop g idx)
  | `PhredLikelihood ->
      `PhredLikelihood (Phred_lhd.fold_over_paired
          ?number_of_reads f1 f2 ?early_stop g idx)
