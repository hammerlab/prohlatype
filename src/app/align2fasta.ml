
open Util
open Common_options

let app_name = "align2fasta"

let gap_char = '.'

let no_gaps s = String.filter ~f:(fun c -> c <> gap_char) s

let mutate (t1, s1) (t2, s2) =
  let off = t2 - t1 in
  for i = 0 to String.length s2 - 1 do
    s1.(i + off) <- s2.(i)
  done

let add_gaps (t1, s1) (t2, n) =
  let off = t2 - t1 in
  for i = 0 to n - 1 do
    s1.(i + off) <- '.'
  done

let zip ~ref ~allele =
  let open Mas_parser in
  let append (_, s) acc = (no_gaps s) :: acc
  and new_cur n s = function
    | Start _
    | End _
    | Boundary _  -> n ()
    | Sequence s  -> s (s.start, s.s)
    | Gap g       -> s (g.start, String.make g.length gap_char)
  and start acc r a =
    match r with
    | []     -> invalid_argf "Couldn't start!"
    | h :: t -> new_cur (fun () -> start [] t) (fun cp -> loop cp acc t a) h
  and end_ cur acc = function
    | []     -> List.rev (append cur acc)
    | h :: t -> new_cur (fun () -> end_ cur acc t)
                  (fun cp -> end_ cp (append cur acc)) h
  and mutate_current cp =
    | Start _
    | End _
    | Boundary _ as e -> invalid_argf "Found %s in sequence at %d"
                           (al_el_to_string e) (fst cp)
    | Sequence s      -> mutate cp (s.start, s.s)
    | Gap g           -> add_gaps cp (g.start, g.length)
  and loop cp acc r a =
    match a with
    | []        -> end_ current acc r
    | ah :: at  ->
        let sp, ss = cp in
        let sa = start_position ah in
        if sa < sp then
          loop current acc r at
        else if sa >= sp + String.length ss then
          start (append cp acc) r a
        else begin
          mutate_current cp ah;
          loop cp acc r at
        end
  in
  start [] ref allele

let convert file =
  let mp = Mas_parser.from_file file in


let () =
  let open Cmdliner in
  let output_fname_arg =
    let docv = "FILE" in
    let doc  = "Output file name, defaults to a (input file)_.[fasta]."
    in
    Arg.(value & opt (some string) None & info ~doc ~docv ["o"; "output"])
  in
  let construct =
    let version = "0.0.0" in
    let doc = "Transform MHC IMGT alignments to pdf graphs." in
    let bug =
      sprintf "Browse and report new issues at <https://github.com/hammerlab/%s"
        repo
    in
    let man =
      [ `S "AUTHORS"
      ; `P "Leonid Rozenberg <leonidr@gmail.com>"
      ; `Noblank
      ; `S "BUGS"
      ; `P bug
      ]
    in
    Term.(const construct
            (* input files *)
            $ file_arg $ merge_arg
            (* output file *)
            $ output_fname_arg
            $ num_alt_arg
            $ allele_arg
            $ allele_regex_arg
            $ remove_reference_flag
            $ not_short_flag $ no_pdf_flag $ no_open_flag
            $ no_cache_flag
            $ max_edge_char_length_flag
            $ not_human_edges_flag
            $ do_not_join_same_sequence_paths_flag
            $ do_not_compress_edges_flag
            $ do_not_compress_start_flag
            $ do_not_insert_newlines_flag
        , info app_name ~version ~doc ~man)
  in
  match Term.eval construct with
  | `Ok n -> exit n
  | `Error _ -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
