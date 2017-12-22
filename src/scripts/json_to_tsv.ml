
open Util

let f_of_yojson =
  ParPHMM_drivers.(Output.of_yojson Multiple_loci.final_read_info_of_yojson)


let of_json_file f =
  Yojson.Safe.stream_from_file f
  |> Stream.next (* only one element per file. *)
  |> f_of_yojson
  |> unwrap_ok

let () =
  let open ParPHMM_drivers in
  let json = Sys.argv.(1) in
  if Filename.check_suffix json "json" then
    let f = of_json_file json in
    let t = (Filename.chop_extension json) ^ ".tsv" in
    let oc = open_out t in
    let out_to_string rrp =
      sprintf "\n%s" (Multiple_loci.final_read_info_to_string rrp)
    in
    Output.f oc f Output.default_depth (Output.TabSeparated out_to_string);
    close_out oc
  else
    printf "not json: %s" json
