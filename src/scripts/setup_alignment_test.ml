
#use "src/scripts/create_index.ml";;
#use "src/scripts/fastq_reader.ml";;
#use "src/scripts/load_test_case.ml";;

let pos = Index.lookup idxall greads.(0) |> unwrap_ok |>fun l -> List.nth_exn l 1 ;;

let r = Graph_alignment.compute_mismatches gall (String.sub_exn greads.(0) ~index:0 ~length:20) (*~search_pos_start:9*) pos  ;;
let rn length = Graph_alignment.compute_mismatches gall (String.sub_exn greads.(0) ~index:0 ~length) (*~search_pos_start:9*) pos  ;;

let al_to_list idx r = Alleles.Map.fold idx ~f:(fun acc c s -> (s, c) :: acc ) ~init:[]  r |> List.sort ~cmp:compare  ;;

let foo ?k ~gi ~length ~read =
  let _, (g, idx) = g_and_idx ?k gi in
  let sub_read = String.sub_exn ~index:0 ~length read in
  let pos = Index.lookup idx sub_read |> unwrap_ok |> List.hd_exn in
  (*let search_pos_start = (Option.value ~default:10 k) - 1 in*)
  let al = Graph_alignment.compute_mismatches g sub_read pos in
  let lal = al_to_list g.Ref_graph.aindex al in
  g, idx, pos, sub_read, (List.rev lal)

let foo_u ?k ~gi ~length ~read =
  let _, (g, idx) = g_and_idx ?k gi in
  let sub_read = String.sub_exn ~index:0 ~length read in
  Index.lookup idx sub_read >>= function
    | [] -> Error "missing position"
    | pos :: _ -> 
    (*let search_pos_start = (Option.value ~default:10 k) - 1 in*)
    let al = Graph_alignment.compute_mismatches g sub_read pos in
    let lal = al_to_list g.Ref_graph.aindex al in
    Ok (g, idx, pos, sub_read, (List.rev lal))



let reads_with_kmers ?k ~gi =
  let _, (g, idx) = g_and_idx ?k gi in
  Array.to_list greads
  |> List.filter_map ~f:(fun s -> Index.lookup idx s |> unwrap_ok |> function | [] -> None | _ -> Some s)

let just_lal ?k ~gi ~length read =
  let _g, _idx, _pos, _sub_read, lal = foo ?k ~gi ~length ~read in
  lal
 
let find_bad ?(length=100) ?(k=10) start_size =
  let start_reads = reads_with_kmers ~k ~gi:start_size in
  let start_lals =
    List.map ~f:(fun read -> read, just_lal ~k ~length ~gi:start_size read)
      start_reads
  in
  let diff_lals new_size prev_lals = 
    List.fold_left prev_lals ~init:([], []) ~f:(fun (nacc, wacc) (read, prev_lal) ->
      let lal_new = just_lal ~k ~length ~gi:new_size read in
      let diff_opt =
        List.fold_left prev_lal ~init:[] ~f:(fun acc (a, c) ->
          let new_c = List.assoc a lal_new in
          if new_c <> c then (a, c) :: acc else acc)
      in
      match diff_opt with
      | [] -> ((read, lal_new) :: nacc), wacc
      | dl -> ((read, lal_new) :: nacc), (read, dl) :: wacc)
  in
  let rec loop prev_size old_lals =
    let gi = prev_size + 1 in
    printf "new gi: %d\n%!" gi;
    match diff_lals gi old_lals with
    | lst, [] -> loop gi lst
    | lst, bd -> bd
  in
  loop start_size start_lals

let bad_read = "CGAGGACCTGCACTCCTGGACCGCCGCGAACACAGCGGCTCAGATCTCCCAGCACAAGTGGGAAGCGGACAAATACTCAGAGAAGGTCAGGGCCTACCTG" ;;
let g27, i27, p27, s27, al27 = foo ~k:10 ~gi:27 ~length:20 ~read:bad_read ;;
let g28, i28, p28, s28, al28 = foo ~k:10 ~gi:28 ~length:20 ~read:bad_read ;;

let bad_read = "CTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTTGGGCGGAGCAGTTGAGAGCCTACCTGGAGGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGA" ;;
let g2, i2, p2, s2, al2 = foo ~k:10 ~gi:2 ~length:20 ~read:bad_read ;;
let g3, i3, p3, s3, al3 = foo ~k:10 ~gi:3 ~length:20 ~read:bad_read ;;


