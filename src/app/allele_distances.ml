(** Compute distances between MHC alleles. *)
open Util
open MoreLabels

let app_name = "allele_distances"

let dt mp n =
  let open Mas_parser in
  let al, allele = List.nth_exn mp.alt_elems n in
  allele_distances ~reference:mp.ref_elems ~allele, al

let bt mp n m =
  let open Mas_parser in
  let al1, allele1 = List.nth_exn mp.alt_elems n in
  let al2, allele2 = List.nth_exn mp.alt_elems m in
  al1, al2, allele_distances_between ~reference:mp.ref_elems ~allele1 ~allele2

module Aset = Set.Make (struct
  type t = string
  let compare = compare
end)

let options_to_targets_and_candidates alignment_file_opt merge_opt =
  let open Mas_parser in
  let open Distances in
  match alignment_file_opt, merge_opt with
  | _, (Some prefix) ->
      let gen = from_file (prefix ^ "_gen.txt") in
      let nuc = from_file (prefix ^ "_nuc.txt") in
      let candidate_s =
        List.fold_left gen.alt_elems ~init:Aset.empty
          ~f:(fun s (allele, alst) -> Aset.add allele s)
      in
      (* we want the sequences from _nuc *)
      let targets, candidates =
        List.fold_left nuc.alt_elems ~init:(Amap.empty, Amap.empty)
          ~f:(fun (tm, cm) (allele, alst) ->
                Amap.add ~key:allele ~data:alst tm,
                if Aset.mem allele candidate_s then
                  Amap.add ~key:allele ~data:alst cm
                else
                  cm)
      in
      Ok (targets, candidates)
  | Some af, None ->
      let mp = from_file af in
      let targets =
        List.fold_left mp.alt_elems ~init:Amap.empty
          ~f:(fun m (allele, alst) -> Amap.add ~key:allele ~data:alst m)
      in
      Ok (targets, targets)
  | None, None  ->
      Error "Either a file or merge argument must be specified"

let allele_distances alignment_file_opt merge_opt mapping =
  options_to_targets_and_candidates alignment_file_opt merge_opt >>=
    begin fun (t,c) ->
      Distances.compute t c mapping >>= fun dmap ->
        printf "allele, closests alleles \n";
        Ok (Distances.Amap.iter dmap ~f:(fun ~key ~data ->
              printf "%s,%s\n" key (String.concat ~sep:"," data)))
  end
  |> function
      | Error e -> eprintf "%s" e; 1
      | Ok ()   -> 0

let () =
  let open Cmdliner in
  let open Common_options in
  let mapping_flag =
    let d = "How to compute the distance between alleles: " in
    Arg.(value & vflag `Trie
      [ `Trie,   info ~doc:(d ^ "trie based off of allele names.") ["trie"]
      ])
  in
  let allele_distances =
    let version = "0.0.0" in
    let doc = "Compute distances between HLA alleles." in
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
    Term.(const allele_distances
            (* Allele set construction args *)
            $ file_arg $ merge_arg
            (* How to measure distances. *)
            $ mapping_flag
        , info app_name ~version ~doc ~man)
  in
  match Term.eval allele_distances with
  | `Ok n            -> exit n
  | `Error _         -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
