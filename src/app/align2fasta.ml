
open Util
open Common_options

let app_name = "align2fasta"

let convert ofile ifile =
  let open Mas_parser in
  let file =
    match ifile with
    | Some i -> i
    | None -> sprintf "A_gen.txt"
  in
  let ofile =
    match ofile with
    | Some o -> o
    | None -> sprintf "%s.fasta" Filename.(chop_extension (basename file))
  in
  let mp = from_file file in
  let r = reference_sequence mp in
  let reference = mp.ref_elems in
  let a =
    List.map mp.alt_elems ~f:(fun (a, allele) ->
      a, apply ~reference ~allele)
  in
  let oc = open_out ofile in
  List.iter ((mp.reference, r) :: a)
    ~f:(fun (r,s) -> fprintf oc "%s\n%s\n" r s);
  close_out oc

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
    Term.(const convert $ output_fname_arg $ file_arg
        , info app_name ~version ~doc ~man)
  in
  match Term.eval convert with
  | `Ok ()   -> exit 0
  | `Error _ -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
