
open Core_kernel.Std

let to_stop = function
  | None   -> fun _ -> false
  | Some n -> fun r -> r >= n

let fold (type a) ?number_of_reads file ~f ~(init : a) =
  let open Biocaml_unix in
  let stop = to_stop number_of_reads in
  let module M = struct exception Fin of a end in
  try
    Future_unix.Reader.with_file file ~f:(fun rdr ->
      (* Avoid circular dep *)
      let fr = Biocaml_unix.Fastq.read rdr in
      Future_unix.Pipe.fold fr ~init:(0,init) ~f:(fun (c, acc) oe ->
        if stop c then raise (M.Fin acc) else
        match oe with
        | Error e -> eprintf "%s\n" (Error.to_string_hum e); (c, acc)
        | Ok fqi  -> (c + 1, f acc fqi)))
    |> snd
  with M.Fin m -> m

let all ?number_of_reads file =
  fold ?number_of_reads
    ~f:(fun acc el -> el :: acc)  ~init:[] file
  |> List.rev
