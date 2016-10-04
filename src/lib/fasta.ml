
(*
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
  *)

open Core_kernel.Std
open CFStream

let to_stop = function
  | None   -> fun _ -> false
  | Some n -> fun r -> r >= n

let fold (type a) ?number_of_reads ~f ~(init:a) file =
  let stop = to_stop number_of_reads in
  let module M = struct exception Fin of a end in
  try
    Biocaml_unix.Fasta.with_file file ~f:(fun _hdr str ->
        Ok (Stream.foldi str ~init ~f:(fun i acc v ->
          if stop i then raise (M.Fin acc) else
            match v with
            | Error e -> eprintf "%s\n" (Error.to_string_hum e); acc
            | Ok fai  -> f acc fai)))
    |> Or_error.ok_exn
  with M.Fin a -> a

let all ?number_of_reads file =
  fold ?number_of_reads file ~f:(fun l fi ->
    (fi.Biocaml_unix.Fasta.description, fi.Biocaml_unix.Fasta.sequence) :: l)
    ~init:[]

