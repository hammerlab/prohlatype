
open Util
open Common_options

let app_name = "align2fasta"

let against_mp ?merge_assoc mp out =
  let open Mas_parser in
  let r = reference_sequence mp in
  let reference = mp.ref_elems in
  let a =
    List.map mp.alt_elems ~f:(fun (a, allele) ->
      a, apply ~reference ~allele)
  in
  let all = (mp.reference, r) :: a in
  let oc = open_out out in
  begin match merge_assoc with
  | None ->
    List.iter all ~f:(fun (r,s) ->
      fprintf oc ">%s %d bp\n%s\n" r (String.length s) s)
  | Some ma ->
    List.iter all ~f:(fun (r,s) ->
      let mg =
        try
          let i = List.assoc r ma in
          if i = r then "" else sprintf " %s introns" i
        with Not_found -> ""
      in
      fprintf oc ">%s%s %d bp\n%s\n" r mg (String.length s) s);
  end;
  close_out oc

let convert ofile ifile merge_file =
  let open Mas_parser in
  begin match merge_file, ifile with
    | None, None   -> None
    | None, Some f ->
        let mp = from_file f in
        let ofiledefault = Filename.(chop_extension (basename f)) in
        let out = sprintf "%s.fasta" (Option.value ofile ~default:ofiledefault) in
        Some (against_mp mp out)
    | Some p, _    ->
        match Merge_mas.and_check p with
        | Error e             -> failwith (sprintf "%s\n" e)
        | Ok (mp, merge_assoc) ->
          let ofiledefault = Filename.(basename p) in
          let out = sprintf "%s.fasta" (Option.value ofile ~default:ofiledefault) in
          Some (against_mp ~merge_assoc mp out)
  end

let () =
  let open Cmdliner in
  let output_fname_arg =
    let docv = "FILE" in
    let doc  = "Output file name, defaults to a (input file)_.[fasta]."
    in
    Arg.(value & opt (some string) None & info ~doc ~docv ["o"; "output"])
  in
  let convert =
    let version = "0.0.0" in
    let doc = "Transform MHC IMGT alignments to fasta." in
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
    Term.(const convert $ output_fname_arg $ file_arg $ merge_arg
        , info app_name ~version ~doc ~man)
  in
  match Term.eval convert with
  | `Ok (Some ())    -> exit 0
  | `Ok (None)       -> eprintf "Neither input alignment nor merge prefix specified!";
                        exit 1
  | `Error _         -> failwith "error"
  | `Version | `Help -> exit 0