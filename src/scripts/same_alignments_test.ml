(* Tests of the alignment algorithm.
  There are two tests, Comparison (C below), default without dropping any
  reads from the front:
    Given a file with reads (fastq) and an alignment file (ex A_nuc), test that
    we get the same alignment of the read against the graph (using the
    multi-allele algorithm) versus the manual algorithm that first just reads
    each alleles sequence from the graph and then does a pairwise string
    comparison.

    ex:
    $ ./same_alignments_test.native /path/to/foo.fastq /path/to/alignment/file.txt
      [length of reads; default: 100] [C] [number of reads to drop; default:0]

    alignment file ex: IMGTHLA/alignments/A_nuc.txt
    If everything works, all of the output will look like:
    comparing alignments for 1569 TT...........GGAAG everything matched!

  and Stability (G below):
    Given a file with reads (fastq) and an alignment file (ex A_nuc), test that
    we get the same alignment of the read irrespective of the number of alleles
    we use (we should find all alignments of an n-1 graph in the n-graph, n=#of
    alleles to use in constructing graph). Defaults to starting with a graph with
    just the reference allele, 0 alternates.

    ex:
    $ ./same_alignments_test.native /path/to/foo.fastq /path/to//alignment/file.txt \
      [length of reads; default: 100] G [graph size start; default:0]

    alignment file ex: IMGTHLA/alignments/A_nuc.txt

    If everything works, the output will look like:
    Testing on 478 reads
    new gi: 1
    new gi: 2
    ..
    new gi: [# of alleles in alignment file]

ex:
  /same_alignments_test.native /path/to/foo.fastq ../foreign/IMGTHLA/alignments/A_nuc.txt \
    [length of reads; default: 100] [G|C] [int]

  C n -> Compare dropping (ignoring) the first n reads from fastq.
  G n -> Stability test, starting with graphs with n alleles and going up.
*)

open Util
open Common
open Alignment

let g_and_idx ?(cache=true) ?(k=10) ~file ?gi () =
  let input = Ref_graph.AlignmentFile file in
  let n = Option.map gi (fun n -> Ref_graph.NumberOfAlts n) in
  if cache then
    Cache.(graph_and_two_index { k = k; graph_args = graph_args ~input ?n () })
  else
    Cache.((invalid_arg_on_error "construct graph and index"
            graph_and_two_index_no_cache) { k = k; graph_args = graph_args ~input ?n () })

let unwrap_sf = function | Stopped r | Finished r -> r

let al_to_list idx r =
  (* Don't worry about stopping in these tests. *)
  let r = unwrap_sf r in
  Alleles.Map.fold idx ~f:(fun acc c s -> (s, c) :: acc ) ~init:[] r
  |> List.sort ~cmp:compare

let test_case ?compare_pos ~length (g, idx) read =
  let sub_read = String.sub_exn ~index:0 ~length read in
  let pos =
    let open Index in
    match lookup idx sub_read with
    | Error m     -> invalid_argf "error looking up %s in index: %s" sub_read (Kmer_table.too_short_to_string m)
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
  let stop = float (String.length read + 100) in (* Don't stop ever *)
  let size = Ref_graph.number_of_alleles g in
  let al = Alignment.compute_mismatches ~early_stop:(size, stop) g sub_read pos |> unwrap_ok in
  let lal = al_to_list g.Ref_graph.aindex al in
  pos, sub_read, (List.rev lal)

let reads_with_kmers reads_file (g, idx) =
  let reads = Fastq.all reads_file in
  let greads =
    List.filter_map reads ~f:(fun r ->
      let s = r.Biocaml_unix.Fastq.sequence in
      match String.index_of_character s 'N' with
      | Some _ -> None | _ -> Some s)
    |> Array.of_list
  in
  Array.to_list greads
  |> List.filter_map ~f:(fun s ->
      match Index.lookup idx s with
      | Ok [] | Error _ -> None
      | Ok ls           -> Some s)

let just_lal ?compare_pos ~length gidxp read =
  let pos, _sub_read, lal = test_case ?compare_pos ~length gidxp read in
  pos, lal

let find_bad ?cache ?(length=100) ?(k=10) ?stop reads_file ~file start_size =
  let stop =
    match stop with
    | Some s -> s
    | None -> let gall, _ = g_and_idx ?cache ~file () in
              Alleles.Map.cardinal gall.Ref_graph.bounds
  in
  let gsidx = g_and_idx ?cache ~k ~file ~gi:start_size () in
  let start_reads = reads_with_kmers reads_file gsidx in
  printf "Testing on %d reads\n" (List.length start_reads);
  let start_lals =
    List.map start_reads ~f:(fun read -> read, just_lal ~length gsidx read)
  in
  let diff_lals new_size prev_lals =
    List.fold_left prev_lals ~init:([], [])
      ~f:(fun (nacc, wacc) (read, (compare_pos, prev_lal)) ->
            let gsidx = g_and_idx ?cache ~k ~file ~gi:new_size () in
            let (pos_new, lal_new) = just_lal ~compare_pos ~length gsidx read in
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
    if gi > stop then
      Ok ("reached stop!")
    else begin
      printf "new gi: %d\n%!" gi;
      match diff_lals gi old_lals with
      | lst, [] -> loop gi lst
      | lst, bd -> Error (gi, bd)
    end
  in
  loop start_size start_lals

let describe_error ?(length=100) ?(k=10) file read gi =
  let gim1 = gi - 1 in
  let cur = !Alignment.debug_ref in
  Alignment.debug_ref := true;
  let (gnm1, i_nm1) as gsidx = g_and_idx ~k ~file ~gi:gim1 () in
  let pnm1, snm1, alnm1 = test_case ~length gsidx read in
  let (gn, i_n) as gsidx2 = g_and_idx ~k ~file ~gi () in
  let pn, sn, aln = test_case ~compare_pos:pnm1 ~length gsidx2 read in
  Alignment.debug_ref := cur;
  (gnm1, i_nm1, pnm1, snm1, alnm1), (gn, i_n, pn, sn, aln)

let manual g idx read =
  Index.lookup idx read
  |> error_map ~f:Kmer_table.too_short_to_string >>=
      function
      | h :: _ -> Ok (h, Alignment.manual_mismatches g read h)
      | []     -> error "read not in index"

let compare_manual g m fm read_len =
  let open Alignment in
  let aindx = g.Ref_graph.aindex in
  Alleles.Map.fold aindx m ~init:[]
    ~f:(fun acc cm allele ->
          let with_all = Alleles.Map.get aindx fm allele in
          match cm, with_all with
          | Ok (Mismatches.Fin sa)   , wa when wa = sa.MismatchesCounts.mismatches ->
            acc
          | Ok (Mismatches.GoOn (sa, sp)) , wa when wa = sa.MismatchesCounts.mismatches + (read_len - sp)  ->
            acc
          | Ok cmo             , wa                                  ->
            (allele, Ok (cmo, wa)) :: acc
          | Error ec           , wa                                  ->
            (allele, Error (ec, wa)) :: acc)


let compare_reads ?length ?(k=10) ?(drop=0) ?num_comp reads_file ~file =
  let g, idx = g_and_idx ~k ~file () in
  let reads = reads_with_kmers reads_file (g, idx) in
  let reads =
    match num_comp with
    | None -> List.drop reads drop
    | Some n -> List.take (List.drop reads drop) n
  in
  let n = List.length reads in
  let rec over_reads i = function
    | []          -> None
    | read :: tl  ->
      let sub_read, sub_read_len =
        match length with
        | None       -> read, String.length read
        | Some index -> String.take read index, index
      in
      printf "comparing alignments for %d %s" i read;
      match manual g idx sub_read with
      | Error em        ->
          eprintf "Skipping %s because wasn't able to map because of %s\n"
            sub_read em;
          over_reads (i + 1) tl
      | Ok (pos, manm)  ->
          let stop = float (String.length sub_read + 100) in
          let size = Ref_graph.number_of_alleles g in
          match Alignment.compute_mismatches ~early_stop:(size, stop) g sub_read pos with
          | Error mes ->
              eprintf "Wasn't able to compute mismatches for %s at %s because of %s"
                sub_read (Index.show_position pos) mes;
              over_reads (i + 1) tl
          | Ok m2     ->
              let m2 = unwrap_sf m2 in
              match compare_manual g manm m2 (String.length sub_read) with
              | [] -> printf " everything matched!\n%!"; over_reads (i + 1) tl
              | ba -> printf " see differences.\n%!"; Some (read, ba)
  in
  (n, over_reads 0 reads)

let () =
  if !Sys.interactive then () else
    let n = Array.length Sys.argv in
    let reads_file =
      if n <= 1 then
        invalid_argf
          "%s [reads_file] [alignment_file] [length] \
            'G' [graph_size for stability test] or \
            'C' [drop number (optional) for comparison test]"
          Sys.argv.(0)
      else
        Sys.argv.(1)
    in
    let file = if n <= 2 then "A_nuc" else Sys.argv.(2) in
    let length = if n <= 3 then 100 else int_of_string Sys.argv.(3) in
    let test =
      if n <= 4 then `Comparison None else
        begin match Sys.argv.(4) with
        | "G" -> `Stability (if n <= 5 then 0 else (int_of_string Sys.argv.(5)))
        | "C" -> `Comparison (Some (int_of_string Sys.argv.(5)))
        | x   -> invalid_argf "Unrecognized arg: %s" x
        end
    in
    match test with
    | `Stability start ->
        begin match find_bad reads_file ~length ~file start with
        | Ok s  -> print_endline s
        | Error (bad_size, bad_elems) ->
            printf "found bad alignments %d with graph size: %d and read length: %d\n"
              (List.length bad_elems) bad_size length;
            exit 1
        end
    | `Comparison drop ->
        begin match compare_reads ?drop reads_file ~length ~file with
        | n, None              ->
            printf "all %d reads match!\n" n
        | n, Some (read, blst) ->
            let open Alignment in
            printf "out of %d reads encountered the following errors:\n" n;
            printf "read: %s\n" read;
            List.iter blst ~f:(fun (allele, oe) ->
              printf "\t%s: %s\n" allele
                (match oe with
                 | Ok ((Mismatches.Fin m), m2)        -> sprintf "Fin %d vs %d" m.Mismatches_config.mismatches m2
                 | Ok ((Mismatches.GoOn (m, p)), m2)  -> sprintf "GoOn %d vs %d, sp: %d" m.Mismatches_config.mismatches m2 p
                 | Error (msg, d)          -> sprintf "Error %s %d" msg d));
            exit 1
        end

