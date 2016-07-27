
open Util

let cargs n =
  { Cache.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
  ; Cache.which = Some (Ref_graph.NumberOfAlts n)
  ; Cache.join_same_sequence = true
  }

let g_and_idx ?(k=10) n = (n, Cache.graph_and_two_index { Cache.k = k; Cache.g = cargs n })

let reads_from_fastq file =
  let ic = open_in file in
  let li = ref [] in
  try
    let rec loop i =
      let line = input_line ic in
      if i mod 4 = 1 then li := line :: !li;
      loop (i + 1)
    in
    loop 0
  with End_of_file ->
    close_in ic;
    !li

let reads =
  reads_from_fastq
    "/Users/leonidrozenberg/Documents/projects/hlatyping/upenn/opti/merged/120013_TGACCA/120013_TGACCA_2fin.fastq"

let greads =
  List.filter reads ~f:(fun r -> match String.index_of_character r 'N' with | Some _ -> false | _ -> true)
  |> Array.of_list

let al_to_list idx r = Alleles.Map.fold idx ~f:(fun acc c s -> (s, c) :: acc ) ~init:[]  r |> List.sort ~cmp:compare

let test_case ?compare_pos ?k ~gi ~length read =
  let _, (g, idx) = g_and_idx ?k gi in
  let sub_read = String.sub_exn ~index:0 ~length read in
  let pos =
    let open Index in
    match lookup idx sub_read with
    | Error m     -> invalid_argf "error looking up %s in index: %s" sub_read m
    | Ok []       -> invalid_argf "empty position returned looking up %s in index" sub_read
    | Ok (h :: t) ->
          match compare_pos with
          | None   -> h
          | Some p ->
              let dp = p.alignment + p.offset in
              match List.find (h :: t) ~f:(fun pn -> pn.alignment + pn.offset = dp) with
              | None -> invalid_argf "Couldn't find desired pos %d in second graph index! for sub_read: %s" dp sub_read
              | Some p -> p
  in
  let al = Alignment.compute_mismatches g sub_read pos |> unwrap_ok in
  let lal = al_to_list g.Ref_graph.aindex al in
  g, idx, pos, sub_read, (List.rev lal)

let reads_with_kmers ?k ~gi =
  let _, (g, idx) = g_and_idx ?k gi in
  Array.to_list greads
  |> List.filter_map ~f:(fun s -> Index.lookup idx s |> unwrap_ok |> function | [] -> None | _ -> Some s)

let just_lal ?k ?compare_pos ~gi ~length read =
  let _g, _idx, pos, _sub_read, lal = test_case ?compare_pos ?k ~gi ~length read in
  pos, lal

let find_bad ?(length=100) ?(k=10) ?(stop=3000) start_size =
  let start_reads = reads_with_kmers ~k ~gi:start_size in
  let start_lals =
    List.map ~f:(fun read -> read, just_lal ~k ~length ~gi:start_size read)
      start_reads
  in
  let diff_lals new_size prev_lals =
    List.fold_left prev_lals ~init:([], [])
      ~f:(fun (nacc, wacc) (read, (compare_pos, prev_lal)) ->
            let (pos_new, lal_new) = just_lal ~compare_pos ~k ~length ~gi:new_size read in
            let diff_opt =
              List.fold_left prev_lal ~init:[] ~f:(fun acc (a, c) ->
                let new_c = List.assoc a lal_new in
                if new_c <> c then (a, c) :: acc else acc)
            in
            match diff_opt with
            | [] -> ((read, (pos_new, lal_new)) :: nacc), wacc
            | dl -> ((read, (pos_new, lal_new)) :: nacc), (read, dl) :: wacc)
  in
  let rec loop prev_size old_lals =
    let gi = prev_size + 1 in
    if gi > stop then (gi, []) else begin
      printf "new gi: %d\n%!" gi;
      match diff_lals gi old_lals with
      | lst, [] -> loop gi lst
      | lst, bd -> (gi, bd)
    end
  in
  loop start_size start_lals

let describe_error ?(length=20) read gi =
  let gim1 = gi - 1 in
  let cur = !Alignment.debug_ref in
  Alignment.debug_ref := true;
  let gnm1, i_nm1, pnm1, snm1, alnm1 = test_case ~k:10 ~gi:gim1 ~length read in
  let gn, i_n, pn, sn, aln = test_case ~compare_pos:pnm1 ~k:10 ~gi ~length read in
  Alignment.debug_ref := cur;
  (gnm1, i_nm1, pnm1, snm1, alnm1), (gn, i_n, pn, sn, aln)


let () =
  if !Sys.interactive then () else
    let n = Array.length Sys.argv in
    let length = if n <= 1 then 20 else int_of_string Sys.argv.(1) in
    let start = if n <= 2 then 2 else int_of_string Sys.argv.(2) in
    let bad_size, bad_elems = find_bad ~length start in
    printf "found bad alignments %d with graph size: %d and read length: %d\n"
      (List.length bad_elems) bad_size length

