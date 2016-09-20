(* How good is the alignment? *)

open Util
open Common

let k = 10
let n = None
let file = Some (to_alignment_file "A_nuc")

let rs = Fastq_reader.reads_from_fastq (List.hd_exn Reads.reads1)

let gall, idx = Cache.graph_and_two_index { Cache.k = k; Cache.g = cache_arg ?file ?n () }

let mismatches_from seq =
  Index.lookup idx seq >>= fun lst ->
    List.fold_left lst ~init:(Ok [])
      ~f:(fun ac_err p -> ac_err >>= fun acc ->
            Alignment.compute_mismatches gall seq p >>= fun al ->
              Ok ((p, al) :: acc))
    >>= fun lst -> Ok (seq, lst)

let mismatches_from_lst seq =
  Index.lookup idx seq >>= fun lst ->
      List.fold_left lst ~init:(Ok [])
        ~f:(fun ac_err p -> ac_err >>= fun acc ->
              Alignment.compute_mismatches_lst gall seq p >>= fun al ->
                Ok ((p, al) :: acc))
      >>= fun lst -> Ok (seq, lst)

let mismatches n =
  match List.nth rs n with
  | None     -> error "n %d out of bounds" n
  | Some seq -> mismatches_from seq

let mismatches_lst n =
  match List.nth rs n with
  | None     -> error "n %d out of bounds" n
  | Some seq -> mismatches_from_lst seq

let report_msm (seq, reslst) =
  List.iter reslst ~f:(fun (p, res) ->
    printf "At position: %s\n" (Index.show_position p);
    Alleles.Map.values_assoc gall.Ref_graph.aindex res
    |> List.sort ~cmp:compare
    |> List.iter ~f:(fun (n, als) ->
        printf "\t%d\t%s\n" n (Alleles.Set.to_human_readable gall.Ref_graph.aindex als)))

let sum_mismatches = List.fold_left ~init:0 ~f:(fun s (_, m) -> s + m)

let report_msm_lst (seq, reslst) =
  List.iter reslst ~f:(fun (p, res) ->
    printf "At position: %s\n" (Index.show_position p);
    Alleles.Map.values_assoc gall.Ref_graph.aindex res
    |> List.sort ~cmp:(fun (m1,_) (m2,_) -> compare (sum_mismatches m1) (sum_mismatches m2))
    |> List.iter ~f:(fun (nlst, als) ->
        let ma = Alleles.Set.min_elt gall.Ref_graph.aindex als in
        let as_seq = Ref_graph.sequence ~start:(`AtPos (p.Index.alignment + p.Index.offset)) ~stop:(`Length 100) gall ma |> unwrap_ok in
        printf "\t%d %s:\n%s\n\t%s\n\t\t%s\n"
          (sum_mismatches nlst)
          ma
          (manual_comp_display seq as_seq)
          (String.concat ~sep:";" (List.map nlst ~f:(fun (p,d) -> sprintf "(%d,%d)" p d)))
          (Alleles.Set.to_human_readable gall.Ref_graph.aindex als)))

let test n =
  match mismatches_lst n with
  | Error msg -> printf "error: %s\n" msg
  | Ok o -> report_msm_lst o

(* Takes a couple of minutes in toploop *)
let mismatches_for_reads ?n () =
  let mrs =
    match n with
    | None  ->
        printf "computing mismatches for every read!\n";
        rs
    | Some nn -> List.take rs nn
  in
  List.mapi mrs ~f:(fun i seq -> printf "%d, %!" i; (mismatches_from seq))

let best_match_allele_map = Alleles.Map.fold_wa ~init:max_int ~f:min

let best_matches ars =
  List.filter_map ars ~f:(function
    | Error e   -> printf "failed to match %s\n" e; None
    | Ok (s,pl) ->
        Some (List.map pl ~f:(fun (p,a) ->
          let bm = best_match_allele_map a in
          (bm, s, p, a))))
  |> List.concat
  |> List.sort ~cmp:(fun (bm1, _, _, _) (bm2, _, _, _) ->
      compare bm1 bm2)

let () =
  if !Sys.interactive then () else
    let ars =
      if Array.length Sys.argv >= 2 then
        mismatches_for_reads ~n:(int_of_string Sys.argv.(1)) ()
      else
        mismatches_for_reads ()
    in
    let foo = best_matches ars in
    let histo_me = List.map foo ~f:(fun (i, _, _, _) -> float i) |> Array.of_list in
    let hist = Oml.Statistics.Descriptive.histogram (`Width 10.) histo_me in
    print_newline ();
    Array.iter hist ~f:(fun (b,v) -> printf "%f \t %d\n" b v)
