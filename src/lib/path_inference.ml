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

module type Config = sig

  type t  (* How we measure the alignment of a sequence against an allele. *)

  val to_string : t -> string           (* Display *)

  val to_final_float : t -> float       (* For projecting to floats *)

  type stop_parameter
  type thread

  val thread_to_seq : thread -> string

  val compute : ?early_stop:stop_parameter -> Ref_graph.t -> thread -> Index.position ->
    ([ `Finished of t Alleles.Map.t | `Stopped of t Alleles.Map.t ], string) result

  type alignment_measure

  val to_am : t Alleles.Map.t -> alignment_measure

  val compare_am : alignment_measure -> alignment_measure -> int

end

module AgainstOneSeq (C : Config) = struct

  let one ?(verbose=false) ?(multi_pos=`Best) ?early_stop g idx seq =
    let report state mm =
      if verbose then
        report_mismatches C.to_string state g.Ref_graph.aindex mm
    in
    let to_final_map = Alleles.Map.map_wa ~f:C.to_final_float in
    let rec across_positions error_msg = function
      | []     -> error "%s" error_msg
      | h :: t ->
        C.compute ?early_stop g seq h >>= function
          | `Stopped mm  ->
              report "stopped" mm;
              across_positions "stopped on previous indexed positions" t
          | `Finished mm ->
              report "finished" mm;
              match t with
              | [] -> Ok (`Finished (to_final_map mm))
              | tl ->
                  match multi_pos with
                  | `TakeFirst    -> Ok (`Finished (to_final_map mm))
                  | `Average      ->
                      let lm = to_final_map mm in
                      let rec loop n = function
                        | []     -> Ok (`Finished lm)
                        | p :: t ->
                            C.compute ?early_stop g seq p >>= function
                              | `Stopped mm  ->
                                  report "stopped during averaging" mm;
                                  loop n t
                              | `Finished mm ->
                                  report "finished during averaging" mm;
                                  (* Online mean update formula *)
                                  Alleles.Map.update2 ~source:mm ~dest:lm
                                    (fun v m -> m +. n *. ((C.to_final_float v) -. m) /. (n +. 1.));
                                  loop (n +. 1.) t
                      in
                      loop 1.0 tl
                  | `Best ->
                      let current_best = C.to_am mm in
                      (*let current_best = Alleles.Map.fold_wa mm ~init:max_int ~f:min in *)
                      let rec loop b bm = function
                        | []     -> Ok (`Finished (to_final_map bm))
                        | p :: t ->
                            C.compute ?early_stop g seq p >>= function
                              | `Stopped mm ->
                                  report "stopped during best search" mm;
                                  loop b bm t
                              | `Finished mm ->
                                  report "finished during best search" mm;
                                  let nb = C.to_am mm in
                                  (*let nb = Alleles.Map.fold_wa mm ~init:max_int ~f:min in *)
                                  if C.compare_am nb b < 0 then
                                    loop nb mm t
                                  else
                                    loop b bm t
                      in
                      loop current_best mm tl
    in
    Index.lookup idx (C.thread_to_seq seq) >>= (across_positions "no index found!")

end (* AgainstOneSeq *)

module SequenceMismatches = AgainstOneSeq ( struct

  type t = int
  let to_string = sprintf "%d"
  let to_final_float m = float m

  type stop_parameter = int * float
  type thread = string

  let thread_to_seq t = t

  let compute = Alignment.compute_mismatches

  (* Lowest mismatch across all alleles. *)
  type alignment_measure = int

  (* Keep the map with the lowest mismatch across all alleles *)
  let to_am = Alleles.Map.fold_wa ~init:max_int ~f:min
  let compare_am am1 am2 = am1 - am2

end)

let multiple_fold_lst ?verbose ?multi_pos ?early_stop g idx =
  let list_map = Alleles.Map.make g.Ref_graph.aindex [] in
  let v = Option.value verbose ~default:false in
  let update v l = v :: l in
  let fold amap seq =
    SequenceMismatches.one ?verbose ?multi_pos ?early_stop g idx seq >>= begin function
      | `Stopped _m ->
          if v then printf "stopped result\n";
          Ok amap
      | `Finished m ->
          if v then printf "finished result\n";
          Alleles.Map.update2 ~source:m ~dest:amap update;
          Ok amap
      end
  in
  list_map, fold

let multiple_fold ?verbose ?multi_pos ?as_ ?early_stop ?er g idx =
  let zero_map = Alleles.Map.make g.Ref_graph.aindex 0. in
  let one_map = Alleles.Map.make g.Ref_graph.aindex 1. in
  let v = Option.value verbose ~default:false in
  let init, update =
    match as_ with
    | None             (* The order of arguments is not a "fold_left", it is 'a -> 'b -> 'b *)
    | Some `Mismatches    -> zero_map, (fun ?er ~len m s -> m +. s)
    | Some `LogLikelihood -> zero_map, (fun ?er ~len m s -> s +. log_likelihood ~len ?er m)
    | Some `Likelihood    -> one_map,  (fun ?er ~len m s -> s *. likelihood ~len ?er m)
  in
  let fold amap seq =
    let len = String.length seq in
    SequenceMismatches.one ?verbose ?multi_pos ?early_stop g idx seq >>= begin function
      | `Stopped _m ->
          if v then printf "stopped result\n";
          Ok amap
      | `Finished m ->
          if v then printf "finished result\n";
          Alleles.Map.update2 ~source:m ~dest:amap (update ?er ~len);
          Ok amap
      end
  in
  init, fold

module PhredLikelihood = AgainstOneSeq ( struct

  open Alignment

  type t = PhredLikelihood_config.t
  let to_string = PhredLikelihood_config.to_string
  let to_final_float t = t.PhredLikelihood_config.sum_llhd

  type stop_parameter = int * float
  type thread = string * float array

  let thread_to_seq (s, _) = s

  let compute = Alignment.compute_plhd

  (* There are 2 ways to compute Best in this case:
     1. lowest number of mismatches as in the other mismatch counting algorithms
     2. highest sum of log likelihoods -> highest probability
     we use 2. *)
  type alignment_measure = float

  let to_am =
    Alleles.Map.fold_wa ~init:neg_infinity
      ~f:(fun a t -> max a t.PhredLikelihood_config.sum_llhd)

  let compare_am am1 am2 =
    (* Higher is better. *)
    let d = am2 -. am1 in
    if d = 0. then 0 else if d < 0. then -1 else 1

end)

let multiple_phred ?verbose ?multi_pos as_ ?early_stop g idx =
  let init = Alleles.Map.make g.Ref_graph.aindex 0. in
  let v = Option.value verbose ~default:false in
  let fold amap seq =
    PhredLikelihood.one ?verbose ?multi_pos ?early_stop g idx seq >>= begin function
      | `Stopped _m ->
          if v then printf "stopped result\n";
          Ok amap
      | `Finished m ->
          if v then printf "finished result\n";
          Alleles.Map.update2 ~source:m ~dest:amap (+.);
          Ok amap
      end
  in
  init, fold
