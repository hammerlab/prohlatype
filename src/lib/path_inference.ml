(** Path Inference : Given a read (a sequence of characters return a PDF over
    the edges (alleles) in the graph. *)

open Util

let likelihood ?(alph_size=4) ?(er=0.01) ~len mismatches =
  let lmp = log (er /. (float (alph_size - 1))) in
  let lcp = log (1. -. er) in
  let c = len - mismatches in
  exp ((float c) *. lcp +. (float mismatches) *. lmp)

let one ?(multi_pos=`TakeFirst) g idx seq =
  let len = String.length seq in
  Index.lookup idx seq >>= function
    | []      -> error "no index found!"
    | h :: t  ->
      Alignment.compute_mismatches g seq h >>= fun mm ->
        let lm = Alleles.Map.map_wa mm ~f:(likelihood ~len) in
        match t with
        | [] -> Ok lm
        | tl ->
            match multi_pos with
            | `TakeFirst -> Ok lm
            | `Average ->
                let n = 1 + List.length t in
                let weight = 1. /. (float n) in
                let lm = Alleles.Map.map_wa lm ~f:(fun m -> m /. weight) in
                let rec loop = function
                  | []     -> Ok lm
                  | p :: t ->
                      Alignment.compute_mismatches g seq p >>= fun mm ->
                          Alleles.Map.update2 mm lm (fun m c -> c +. weight *. (likelihood ~len m));
                          loop t
                in
                loop tl
            (* Find the best alignment over all positions! *)
            | `Best ->
                let current_best = Alleles.Map.fold_wa mm ~init:max_int ~f:min in
                let rec loop b bm = function
                  | []     -> Ok (Alleles.Map.map_wa mm ~f:(likelihood ~len))
                  | p :: t ->
                      (* TODO: If this becomes a bottleneck, we can stop [compute_mismatches]
                         early if the min mismatch > b *)
                      Alignment.compute_mismatches g seq p >>= fun nbm ->
                          let nb = Alleles.Map.fold_wa nbm ~init:max_int ~f:min in
                          if nb < b then
                            loop nb nbm t
                          else
                            loop b bm t
                in
                loop current_best mm tl

