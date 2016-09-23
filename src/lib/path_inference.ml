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

let report_mismatches aindex mm =
  printf "reporting mismatches\n";
  let () =
    Alleles.Map.values_assoc aindex mm
    |> List.sort ~cmp:(fun (i1,_) (i2,_) -> compare i1 i2)
    |> List.iter ~f:(fun (w, a) ->
        printf "%02d\t%s\n" w
          (insert_chars ['\t'; '\t'; '\n']
            (Alleles.Set.to_human_readable aindex ~max_length:10000 ~complement:`No a)))
  in
  printf "finished reporting mismatches\n%!"

let one ?(verbose=false) ?(multi_pos=`Best) g idx seq =
  let report mm =
    if verbose then report_mismatches g.Ref_graph.aindex mm else ()
  in
  (*let len = String.length seq in *)
  let to_float_map = Alleles.Map.map_wa ~f:float in
  Index.lookup idx seq >>= function
    | []      -> error "no index found!"
    | h :: t  ->
      Alignment.compute_mismatches g seq h >>= fun mm ->
        report mm;
        match t with
        | [] -> Ok (to_float_map mm)
        | tl ->
            match multi_pos with
            | `TakeFirst        -> Ok (to_float_map mm)
            | `Average          ->
                let n = 1 + List.length t in
                let weight = 1. /. (float n) in
                let lm = Alleles.Map.map_wa (to_float_map mm) ~f:(fun m -> m /. weight) in
                let rec loop = function
                  | []     -> Ok lm
                  | p :: t ->
                      Alignment.compute_mismatches g seq p >>= fun mm ->
                          report mm;
                          Alleles.Map.update2 mm lm (fun m c -> c +. weight *. (float m));
                          loop t
                in
                loop tl
            (* Find the best alignment, the one with the lowest number of
               mismatches, over all indexed positions! *)
            | `Best             ->
                let current_best = Alleles.Map.fold_wa mm ~init:max_int ~f:min in
                let rec loop b bm = function
                  | []     -> Ok (to_float_map bm)
                  | p :: t ->
                      (* TODO: If this becomes a bottleneck, we can stop [compute_mismatches]
                         early if the min mismatch > b *)
                      Alignment.compute_mismatches g seq p >>= fun nbm ->
                          report nbm;
                          let nb = Alleles.Map.fold_wa nbm ~init:max_int ~f:min in
                          if nb < b then
                            loop nb nbm t
                          else
                            loop b bm t
                in
                loop current_best mm tl

let mx b amap =
  let n, s = Alleles.Map.fold_wa ~f:(fun (n,s) v -> (n+1, s+.v)) ~init:(0,0.) amap in
  (float n) /. s < b

let yes _ = true

let multiple_fold_lst ?verbose ?multi_pos ?filter g idx =
  let list_map = Alleles.Map.make g.Ref_graph.aindex [] in
  let filter = match filter with | None -> yes | Some n -> mx n in
  let update v l = v :: l in
  let fold amap seq =
    one ?verbose ?multi_pos g idx seq >>= fun m ->
      if filter m then begin
        Alleles.Map.update2 ~source:m ~dest:amap update;
        match verbose with | Some true -> printf "passed filter\n" | _ -> ()
      end else begin
        match verbose with | Some true -> printf "did not pass filter\n" | _ -> ()
      end;
      Ok amap
  in
  list_map, fold

let multiple_fold ?verbose ?multi_pos ?as_ ?filter g idx =
  let zero_map = Alleles.Map.make g.Ref_graph.aindex 0. in
  let one_map = Alleles.Map.make g.Ref_graph.aindex 1. in
  let filter = match filter with | None -> yes | Some n -> mx n in
  let init, update =
    match as_ with
    | None
    | Some `Mismatches    -> zero_map, (+.)
    | Some `LogLikelihood -> zero_map, (fun s m -> s +. log_likelihood ~len:100 m) (* TODO: dep on seq *)
    | Some `Likelihood    -> one_map,  (fun s m -> s *. likelihood ~len:100 m)     (* TODO: dep on seq *)
  in
  let fold amap seq =
    one ?verbose ?multi_pos g idx seq >>= fun m ->
      if filter m then begin
        Alleles.Map.update2 ~source:m ~dest:amap update;
        match verbose with | Some true -> printf "passed filter\n" | _ -> ()
      end else begin
        match verbose with | Some true -> printf "did not pass filter\n" | _ -> ()
      end;
      Ok amap
  in
  init, fold

let multiple ?multi_pos ?as_ g idx =
  let init, f = multiple_fold g idx in
  let rec loop a = function
    | []     -> Ok a
    | h :: t -> f a h >>= fun m -> loop m t
  in
  loop init
