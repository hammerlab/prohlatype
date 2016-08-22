(** Path Inference : Given a read (a sequence of characters return a PDF over
    the edges (alleles) in the graph. *)

open Util

let likelihood ?(alph_size=4) ?(er=0.01) ~len mismatches =
  let lmp = log (er /. (float (alph_size - 1))) in
  let lcp = log (1. -. er) in
  let c = len - mismatches in
  exp ((float c) *. lcp +. (float mismatches) *. lmp)

let report_mismatches aindex mm =
  Alleles.Map.values_assoc aindex mm
  |> List.sort ~cmp:(fun (i1,_) (i2,_) -> compare i1 i2)
  |> List.iter ~f:(fun (w, a) ->
      printf "%02d\t%s\n" w
        (insert_chars ['\t'; '\t'; '\n']
          (Alleles.Set.to_human_readable aindex ~max_length:1000 ~complement:`No a)))

let one ?(verbose=false) ?(multi_pos=`TakeFirst) g idx seq =
  let report mm =
    if verbose then report_mismatches g.Ref_graph.aindex mm else ()
  in
  let len = String.length seq in
  Index.lookup idx seq >>= function
    | []      -> error "no index found!"
    | h :: t  ->
      Alignment.compute_mismatches g seq h >>= fun mm ->
        report mm;
        let lm = Alleles.Map.map_wa mm ~f:(likelihood ~len) in
        match t with
        | [] -> Ok lm
        | tl ->
            match multi_pos with
            | `TakeFirst        -> Ok lm
            | `Average          ->
                let n = 1 + List.length t in
                let weight = 1. /. (float n) in
                let lm = Alleles.Map.map_wa lm ~f:(fun m -> m /. weight) in
                let rec loop = function
                  | []     -> Ok lm
                  | p :: t ->
                      Alignment.compute_mismatches g seq p >>= fun mm ->
                          report mm;
                          Alleles.Map.update2 mm lm (fun m c -> c +. weight *. (likelihood ~len m));
                          loop t
                in
                loop tl
            (* Find the best alignment over all positions! *)
            | `Best             ->
                let current_best = Alleles.Map.fold_wa mm ~init:max_int ~f:min in
                let rec loop b bm = function
                  | []     -> Ok (Alleles.Map.map_wa mm ~f:(likelihood ~len))
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

let multiple_fold ?verbose ?multi_pos g idx =
  let init = Alleles.Map.make g.Ref_graph.aindex 0. in
  let fold amap seq =
    one ?verbose ?multi_pos g idx seq >>= fun m -> Alleles.Map.update2 m amap (+.); Ok amap
  in
  init, fold

let multiple ?multi_pos g idx =
  let init, f = multiple_fold g idx in
  let rec loop a = function
    | []     -> Ok a
    | h :: t -> f a h >>= fun m -> loop m t
  in
  loop init

