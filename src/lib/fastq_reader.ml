
open Util

let to_stop = function
  | None   -> fun _ -> false
  | Some n -> fun r -> r >= n

let fold ?number_of_reads file ~f ~init =
  let ic = open_in file in
  let li = ref init in
  let stop = to_stop number_of_reads in
  try
    let rec loop i =
      let line = input_line ic in
      if i mod 4 = 1 then li := f !li line;
      if stop ((i + 1) / 4) then begin
        close_in ic;
        !li
      end else
        loop (i + 1)
    in
    loop 0
  with End_of_file ->
    close_in ic;
    !li

let all ?number_of_reads =
  fold ?number_of_reads
    ~f:(fun acc el -> el :: acc)  ~init:[]

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
