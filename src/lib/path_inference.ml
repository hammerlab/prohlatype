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

let report_mismatches state aindex mm =
  printf "reporting %s mismatches\n" state;
  let () =
    Alleles.Map.values_assoc aindex mm
    |> List.sort ~cmp:(fun (i1,_) (i2,_) -> compare i1 i2)
    |> List.iter ~f:(fun (w, a) ->
        printf "%02d\t%s\n" w
          (insert_chars ['\t'; '\t'; '\n']
            (Alleles.Set.to_human_readable aindex ~max_length:10000 ~complement:`No a)))
  in
  printf "finished reporting mismatches\n%!"

let one ?(verbose=false) ?(multi_pos=`Best) ?early_stop g idx seq =
  let report state mm = if verbose then report_mismatches state g.Ref_graph.aindex mm in
  (*let len = String.length seq in *)
  let to_float_map = Alleles.Map.map_wa ~f:float in
  let early_stop = Option.value early_stop ~default:(String.length seq) in
  let rec across_positions error_msg = function
    | []     -> error "%s" error_msg
    | h :: t ->
      Alignment.compute_mismatches g early_stop seq h >>= function
        | `Stopped mm  ->
            report "stopped" mm;
            across_positions "stopped on previous indexed positions" t
        | `Finished mm ->
            report "finished" mm;
            match t with
            | [] -> Ok (to_float_map mm)
            | tl ->
                match multi_pos with
                | `TakeFirst    -> Ok (to_float_map mm)
                | `Average      ->
                    let lm = to_float_map mm in
                    let rec loop n = function
                      | []     -> Ok lm
                      | p :: t ->
                          Alignment.compute_mismatches g early_stop seq p >>= function
                            | `Stopped mm  ->
                                report "stopped during averaging" mm;
                                loop n t
                            | `Finished mm ->
                                report "finished during averaging" mm;
                                (* Online mean update formula *)
                                Alleles.Map.update2 ~source:mm ~dest:lm
                                  (fun v m -> m +. n *. ((float v) -. m) /. (n +. 1.));
                                loop (n +. 1.) t
                    in
                    loop 1.0 tl
                | `Best ->
                    let current_best = Alleles.Map.fold_wa mm ~init:max_int ~f:min in
                    let rec loop b bm = function
                      | []     -> Ok (to_float_map bm)
                      | p :: t ->
                          Alignment.compute_mismatches g early_stop seq p >>= function
                            | `Stopped mm ->
                                report "stopped during best search" mm;
                                loop b bm t
                            | `Finished mm ->
                                report "finished during best search" mm;
                                let nb = Alleles.Map.fold_wa mm ~init:max_int ~f:min in
                                if nb < b then
                                  loop nb mm t
                                else
                                  loop b bm t
                    in
                    loop current_best mm tl
  in
  Index.lookup idx seq >>= (across_positions "no index found!")

let mx b amap =
  let n = Alleles.Map.cardinal amap in
  let s = Alleles.Map.fold_wa ~init:0. ~f:(+.) amap in
  s /. (float n) < b

let yes _ = true

let multiple_fold_lst ?verbose ?multi_pos ?filter ?early_stop g idx =
  let list_map = Alleles.Map.make g.Ref_graph.aindex [] in
  let filter = match filter with | None -> yes | Some n -> mx n in
  let update v l = v :: l in
  let fold amap seq =
    one ?verbose ?multi_pos ?early_stop g idx seq >>= fun m ->
      if filter m then begin
        Alleles.Map.update2 ~source:m ~dest:amap update;
        match verbose with | Some true -> printf "passed filter\n" | _ -> ()
      end else begin
        match verbose with | Some true -> printf "did not pass filter\n" | _ -> ()
      end;
      Ok amap
  in
  list_map, fold

let multiple_fold ?verbose ?multi_pos ?as_ ?filter ?early_stop ?er g idx =
  let zero_map = Alleles.Map.make g.Ref_graph.aindex 0. in
  let one_map = Alleles.Map.make g.Ref_graph.aindex 1. in
  let filter = match filter with | None -> yes | Some n -> mx n in
  let init, update =
    match as_ with
    | None             (* The order of arguments is not a "fold_left", it is 'a -> 'b -> 'b *)
    | Some `Mismatches    -> zero_map, (fun ?er ~len m s -> m +. s)
    | Some `LogLikelihood -> zero_map, (fun ?er ~len m s -> s +. log_likelihood ~len ?er m)
    | Some `Likelihood    -> one_map,  (fun ?er ~len m s -> s *. likelihood ~len ?er m)
  in
  let fold amap seq =
    let len = String.length seq in
    one ?verbose ?multi_pos ?early_stop g idx seq >>= fun m ->
      if filter m then begin
        Alleles.Map.update2 ~source:m ~dest:amap (update ?er ~len);
        match verbose with | Some true -> printf "passed filter\n" | _ -> ()
      end else begin
        match verbose with | Some true -> printf "did not pass filter\n" | _ -> ()
      end;
      Ok amap
  in
  init, fold

let multiple ?multi_pos ?as_ ?er g idx =
  let init, f = multiple_fold ?er g idx in
  let rec loop a = function
    | []     -> Ok a
    | h :: t -> f a h >>= fun m -> loop m t
  in
  loop init
