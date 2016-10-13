
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

let phred_probabilities s =
  let open Biocaml_unix in
  let n = String.length s in
  let a = Array.create ~len:n 0. in
  let rec loop i =
    if i = n then Ok a else
      Result.bind (Phred_score.of_char s.[i])
        (fun c -> a.(i) <- Phred_score.to_probability c; loop (i + 1))
  in
  loop 0

let rec same cmp l1 l2 =
  let rec single_loop search acc = function
    | []                   -> None, acc
    | h :: t when search h -> (Some h), (List.rev acc) @ t
    | h :: t               -> single_loop search (h :: acc) t
  in
  let rec double_loop l2 acc = function
    | []     -> None
    | h :: t ->
        match single_loop (cmp h) [] l2 with
        | None, _     -> double_loop l2 (h :: acc) t
        | Some f, lst -> Some (h, f, (List.rev acc) @ t, lst)
  in
  double_loop l2 [] l1

let fold_paired ?number_of_reads file1 file2 ~f ~init to_key =
  let open Biocaml_unix in
  let stop = to_stop number_of_reads in
  let cmp r1 r2 = (to_key r1) = (to_key r2) in
  let same = same cmp in
  let open Future_unix in
  let open Deferred in
  Reader.with_file file1 ~f:(fun rdr1 ->
    Reader.with_file file2 ~f:(fun rdr2 ->
      (* Avoid circular dep *)
      let fr1 = Biocaml_unix.Fastq.read rdr1 in
      let fr2 = Biocaml_unix.Fastq.read rdr2 in
      let rec two_pipe_fold c s1 s2 acc =
        if stop c then
          return (`DesiredReads acc)
        else
          Pipe.read fr1 >>= fun r1 ->
            Pipe.read fr2 >>= fun r2 ->
              match r1, r2 with
              | `Eof, `Eof                      -> return (`BothFinished acc)
              (* TODO: Add continuations to finish off the rest of the sequence as
                single read. *)
              | `Eof, _                         -> return (`OneReadPairedFinished (1, acc))
              | _   , `Eof                      -> return (`OneReadPairedFinished (2, acc))
              | `Ok (Error ea), `Ok (Error eb)  ->
                  eprintf "%s\n" (Error.to_string_hum ea);
                  eprintf "%s\n" (Error.to_string_hum eb);
                  two_pipe_fold c s1 s2 acc
              | `Ok (Error ea), `Ok (Ok b)      ->
                  eprintf "%s\n" (Error.to_string_hum ea);
                  two_pipe_fold c s1 (b::s2) acc
              | `Ok (Ok a), `Ok (Error eb)      ->
                  eprintf "%s\n" (Error.to_string_hum eb);
                  two_pipe_fold c (a::s1) s2 acc
              | `Ok (Ok a), `Ok (Ok b)          ->
                  if cmp a b then
                    two_pipe_fold (c + 1) s1 s2 (f acc a b)
                  else
                    let ns1 = a :: s1 in
                    let ns2 = b :: s2 in
                    match same ns1 ns2 with
                    | Some (r1, r2, ns1, ns2) ->
                        two_pipe_fold (c + 1) ns1 ns2 (f acc r1 r2)
                    | None ->
                        two_pipe_fold c ns1 ns2 acc
      in
      two_pipe_fold 0 [] [] init))
