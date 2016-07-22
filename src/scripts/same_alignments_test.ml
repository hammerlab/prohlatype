
open Util

let cargs n =
  { Cache.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
  ; Cache.which = Some (Ref_graph.NumberOfAlts n)
  ; Cache.normalize = true
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
  reads_from_fastq "/Users/leonidrozenberg/Documents/projects/hlatyping/upenn/opti/merged/120013_TGACCA/120013_TGACCA_2fin.fastq"

let greads =
  List.filter reads ~f:(fun r -> match String.index_of_character r 'N' with | Some _ -> false | _ -> true)
  |> Array.of_list

let al_to_list idx r = Alleles.Map.fold idx ~f:(fun acc c s -> (s, c) :: acc ) ~init:[]  r |> List.sort ~cmp:compare

let test_case ?k ~gi ~length ~read =
  let _, (g, idx) = g_and_idx ?k gi in
  let sub_read = String.sub_exn ~index:0 ~length read in
  let pos = Index.lookup idx sub_read |> unwrap_ok |> List.hd_exn in
  (*let search_pos_start = (Option.value ~default:10 k) - 1 in *)
  let al = Alignment.compute_mismatches g sub_read (*~search_pos_start*) pos in
  let lal = al_to_list g.Ref_graph.aindex al in
  g, idx, pos, sub_read, (List.rev lal)

let reads_with_kmers ?k ~gi =
  let _, (g, idx) = g_and_idx ?k gi in
  Array.to_list greads
  |> List.filter_map ~f:(fun s -> Index.lookup idx s |> unwrap_ok |> function | [] -> None | _ -> Some s)

let just_lal ?k ~gi ~length read =
  let _g, _idx, _pos, _sub_read, lal = test_case ?k ~gi ~length ~read in
  lal

let find_bad ?(length=100) ?(k=10) ?(stop=3000) start_size =
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
    if gi > stop then (gi, []) else begin
      printf "new gi: %d\n%!" gi;
      match diff_lals gi old_lals with
      | lst, [] -> loop gi lst
      | lst, bd -> (gi, bd)
    end
  in
  loop start_size start_lals


let () =
  if !Sys.interactive then () else
    let n = Array.length Sys.argv in
    let length = if n <= 1 then 20 else int_of_string Sys.argv.(1) in
    let bad_size, bad_elems = find_bad ~length 2 in
    printf "found bad alignments %d with graph size: %d\n"
      (List.length bad_elems) bad_size

