(* Compare graph construction, from the alignment files, versus the sequence
   data represented in the IMGTHLA fasta files for the alleles. *)
open Util
open Common

let (//) = Filename.concat
let is_char = function 'A' | 'C' | 'G' | 'T' -> true | _ -> false

let apply_changes ~r d =
  let lr = String.to_character_list r in
  let ld = String.to_character_list d in
  let rec loop acc lr ld =
    match lr, ld with
    | [], []                              -> acc
    | ' ' :: tr, ' ' :: td                -> loop acc tr td
    | '|' :: tr, '|' :: td                -> loop acc tr td
    |  _  :: tr, '*' :: td                -> loop acc tr td
    |  c  :: tr, '-' :: td when is_char c -> loop (c :: acc) tr td
    |  c  :: tr, '.' :: td when is_char c -> loop acc tr td
    | '.' :: tr, '.' :: td                -> loop acc tr td
    |  c  :: tr,  d  :: td when is_char c && is_char d
                                          -> loop (d :: acc) tr td
    | '.' :: tr,  c  :: td when is_char c -> loop (c :: acc) tr td
    |  x  :: _,  []      -> invalid_argf "delta too short! last ref %c" x
    |  [],       y  :: _ -> invalid_argf "ref too short! last d %c" y
    |  x  :: _,  y  :: _ -> invalid_argf "not implemented characters: %c vs %c"
                                x y
  in
  loop [] lr ld
  |> List.rev
  |> String.of_character_list

let manual_diff ~reference ~allele ~file () =
  let cmd =
    sprintf "grep -E \"%s|%s\" %s | tr -s ' ' | cut -d ' ' -f 3-"
      reference allele (to_alignment_file file)
  in
  let ic = Unix.open_process_in cmd in
  let rec loop acc =
    try loop (input_line ic :: acc)
    with End_of_file ->
      close_in ic;
      acc
  in
  let lst = List.rev (loop []) in
  let parsed =
    let rec loop acc = function
      | []           -> List.rev acc
      | r :: d :: tl -> loop ((apply_changes ~r d) :: acc)  tl
      | x :: []      -> eprintf "Unbalanced grep for %s %s, adding last sequence.\n"
                          allele file;
                        let without_spaces = String.filter_map ~f:(function ' ' -> None | x -> Some x) x in
                        List.rev (without_spaces :: acc)
    in
    try
      loop [] lst
    with Invalid_argument m ->
      eprintf "%s couldn't manually apply changes for %s vs %s\n" m allele reference;
      []
  in
  String.concat parsed

let test_sequences file =
  let all_args = Cache.graph_args ~file:(to_alignment_file file) () in
  let gall = Cache.graph all_args in
  let a_fasta = Fasta.all (to_fasta_file file) in
  List.fold_left a_fasta ~init:[] ~f:(fun acc (header, seq) ->
    let allele = List.nth_exn (String.split ~on:(`Character ' ') header) 1 in
    match Ref_graph.sequence gall allele with
    (* TODO: This should be an Error not an exception! *)
    | exception Not_found -> Printf.printf "Couldn't find sequence for %s\n" allele;
                             (allele, (seq, "")) :: acc
    | Error msg           -> Printf.printf "missing sequence for %s because %s!\n" allele msg;
                             (allele, (seq, "")) :: acc
    | Ok graph_seq        ->
        if seq <> graph_seq then
          (allele, (seq, graph_seq)) :: acc
        else
          acc)
  |> fun lst -> List.length a_fasta, lst

let test known_differences ~reference file =
  let number_sequences, different = test_sequences file in
  if different = [] then begin
    Printf.printf "Nothing different in %d sequences!\n" number_sequences;
    0
  end else
    List.fold_left different ~init:0 ~f:(fun error_code (allele, (fasta, graph_seq)) ->
      let md = manual_diff ~reference ~allele ~file () in
      let known = Array.mem allele known_differences in
      let () =
        if md = fasta then
          Printf.printf "Known: %b\t%s different and manual matches fasta:\ngraph:\n%s\nmanual:\n%s\n"
            known allele graph_seq md
        else if md = graph_seq then
          Printf.printf "Known: %b\t%s different but manual matches graph seq.\nfasta:\n%s\ngraph:\n%s\n"
            known allele fasta graph_seq
        else
          Printf.printf
            "Known: %b\t%s different but manual doesn't match any sequence.\nfasta:\n%s\ngraph:\n%s\nmanual:\n%s\n"
              known allele fasta graph_seq md
      in
      if known then error_code else error_code + 1)

let () =
  let n = Array.length Sys.argv in
  if !Sys.interactive then
    ()
  else
    let file, reference, known_differences =
      if n <= 1 then "A_nuc", "A\\*01:01:01:01", [||]
      else if n = 2 then
        failwith "Specify file (ex 'A_nuc'),reference (A*01:01:01:01) and any know differences (A*31:02:24)"
      else
        Sys.argv.(1),  Sys.argv.(2), Array.sub Sys.argv ~pos:3 ~len:(n - 3)
    in
    test known_differences ~reference file |> exit
