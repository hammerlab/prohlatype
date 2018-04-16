
open Core_kernel
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
