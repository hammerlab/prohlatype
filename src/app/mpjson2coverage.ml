
open Prohlatype
open Cmdline_options

let app_name = "mpjson2coverage"

let f_of_yojson =
  ParPHMM_drivers.(Output.of_yojson Multiple_loci.final_read_info_of_yojson)

let of_json_file f =
  Yojson.Safe.stream_from_file f
  |> Stream.next (* only one element per file. *)
  |> f_of_yojson
  |> unwrap_ok

let convert output_opt json_file =
  let open ParPHMM_drivers in
  let f = of_json_file json_file in
  let output_filename =
    Option.value output_opt
      ~default:((Filename.remove_extension json_file) ^ ".tsv")
  in
  let oc = open_out output_filename in
  let out_to_string rrp =
    sprintf "\n%s" (Multiple_loci.final_read_info_to_string rrp)
  in
  Output.f oc f Output.default_depth (Output.TabSeparated out_to_string);
  close_out oc

let () =
  let open Cmdliner in
  let output_fname_arg =
    let docv = "FILE" in
    let doc  = "Output file name, defaults to \"(input file).tsv\"." in
    Arg.(value
        & opt (some string) None
        & info ~doc ~docv ["o"; "output"]
        )
  in
  let input_fname_arg =
    let docv = "FILE" in
    let doc  = "Input json filename." in
    Arg.(required
        & pos 0 (some file) None
        & info ~doc ~docv []
        )
  in
  let convert =
    let doc =
      "Transform multi_par's or par_type's json output to a Tab Separated File.\
       \
       JSON provides a compact and machine useful format for communicating HLA \
       typing results. Unfortunately, it is not as useful for human \
       interpretation. This utility converts output as generated with the \
       \"-json\" flag into Tab Separated format."
    in
    let bugs =
      sprintf "Browse and report new issues at <https://github.com/hammerlab/%s>"
        repo
    in
    let man =
      [ `S "AUTHORS"
      ; `P "Leonid Rozenberg <leonidr@gmail.com>"
      ; `Noblank
      ; `S "BUGS"
      ; `P bugs
      ]
    in
    Term.(const convert $ output_fname_arg $ input_fname_arg
        , info app_name ~version ~doc ~man)
  in
  match Term.eval convert with
  | `Ok ()           -> exit 0
  | `Error _         -> failwith "Is this JSON input?"
  | `Version | `Help -> exit 0
