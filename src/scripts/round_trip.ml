
open Util

let fasta_reader file =
  let ic = open_in file in
  let allele_from_line line =
    List.nth_exn (String.split ~on:(`Character ' ') line) 1
  in
  let rec read_read allele sacc acc =
    try
      let line = input_line ic in
      match String.get line ~index:0 with
      | None      -> raise (Invalid_argument "empty line!")
      | Some '>'  ->
          let nallele = allele_from_line line in
          let nacc    = (allele, (String.concat (List.rev sacc))) :: acc in
          read_read nallele [] nacc
      | Some _    ->
          read_read allele (line :: sacc) acc 
    with End_of_file ->
      close_in ic;
      (allele, (String.concat (List.rev sacc))) :: acc
  in
  let first_line = input_line ic in
  match String.get first_line ~index:0 with
  | Some '>'  -> read_read (allele_from_line first_line) [] [] 
  | None      -> eprintf "empty first line!"; []
  | Some c    -> eprintf "first line doesn't start with '>' %c" c ; []

let (//) = Filename.concat
let root_dir = "../foreign/IMGTHLA/"
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

let manual_diff ?(reference="A\\*01:01:01:01") ~allele ~file () =
  let cmd =
    sprintf "grep -E \"%s|%s\" %s | tr -s ' ' | cut -d ' ' -f 3-"
      reference allele (root_dir // "alignments" // (file ^ ".txt"))
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
      | x :: []      -> invalid_argf "Unbalanced grep for %s %s, last: %s"
                          allele file x
    in
    loop [] lst
  in
  String.concat parsed

let test_sequences file =
  let all_args =
    { Cache.alignment_file = root_dir // "alignments" // (file ^ ".txt")
    ; Cache.which = None
    } 
  in
  let gall = Cache.graph all_args in
  let a_fasta = fasta_reader (root_dir // "fasta" // (file ^ ".fasta")) in 
  List.fold_left a_fasta ~init:[] ~f:(fun acc (allele, seq) ->
    let gseq = 
      try Ref_graph.sequence gall allele 
      with Not_found ->
        Printf.eprintf "missing %s!\n" allele;
        ""
    in
    if seq <> gseq then
      (allele, (seq, gseq)) :: acc
    else
      acc)
  |> fun lst -> List.length a_fasta, lst

let test ~reference file =
  let number_sequences, different = test_sequences file in
  if different = [] then
    Printf.printf "Nothing different in %d sequences!" number_sequences
  else
    List.iter different ~f:(fun (allele, (fasta, gseq)) ->
      let md = manual_diff ~reference ~allele ~file () in
      if md = fasta then
        Printf.printf "%s different and manual matches fasta.\n" allele
      else if md = gseq then
        Printf.printf "%s different but manual matches graph seq.\n" allele
      else
        Printf.printf
          "%s different but manual doesn't match any sequence.\nfasta:\n%s\ngraph:\n%s\nmanual:\n%s\n"
            allele fasta gseq md)

let () =
  let n = Array.length Sys.argv in
  if !Sys.interactive then
    ()
  else
    let file, reference =
      if n <= 1 then "A_nuc", "A\\*01:01:01:01"
      else if n = 2 then
        failwith "Specify file (ex 'A_nuc') and reference (A*01:01:01:01)"
      else
      Sys.argv.(1),  Sys.argv.(2)
    in
    test ~reference file
