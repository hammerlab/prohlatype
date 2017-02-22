

let tr =
  [ "120021"; "120074"; "120984"; "122879"; "123768"; "125785"; "126119"
  ; "126966"; "127051"; "127079"; "127245"; "128964"; "130445"; "160426"
  ; "160452"; "160917"; "161549"; "161702"; "162569"; "162636"; "162781"
  ; "164054"; "164115"; "164543"; "164964"; "165266"; "200220"; "201858"
  ; "203288"; "203482"; "204378"; "204587"; "205302"; "205855"; "205856"
  ; "205944"; "206039"; "206130"; "206139"; "206704"; "206708"; "206709"
  ; "206961"; "207956"; "207978"
  ]

let vd =
  [ "120013"; "120989"; "125881"; "127064"; "129670"; "160908"; "161712"
  ; "163581"; "164779"; "201063"; "203778"; "205673"; "206038"; "206240"
  ; "206876"
  ]

let _hdr :: data = Csv.load res_fname

let to_res g t =
  let els = match t with | `Train -> tr | `Valid -> vd in
  List.filter_map data ~f:(fun lst ->
    let k = List.nth_exn lst 4 in
    if List.mem k ~set:els then
      match Nomenclature.parse (List.nth_exn lst 1) with
      | Ok (s, (nc, _)) when s = g -> Some (k, nc)
      | Error e    -> invalid_arg e
      | _          -> None
    else
        None)

let load_prohlatype fname =
  Csv.load ~separator:'\t' fname
  |> List.map ~f:(function
    | [_; a] -> Nomenclature.parse a
    | _ -> error "not 2")

let lookup_position tsv = function
  | Nomenclature.Two (_, _) as nmt ->
      List.findi tsv ~f:(fun _i l ->
        match l with
        | Ok (_, (nmt2, _)) when nmt2 = nmt -> true
        | _ -> false)
      |> Option.map ~f:(fun (i, _) -> i)
      |> Option.value ~default:max_int
  | _ -> invalid_arg "Not 2 in mptt"

(*
let position to_prohlatype_fname (key, nmt) =
  match nmt with
  | Nomenclature.Two (_, _) ->
      let fname = to_prohlatype_fname key in
      let as_tv = load_prohlatype fname in
      List.findi as_tv ~f:(fun _i l ->
        match l with
        | Ok (_, (nmt2, _)) when nmt2 = nmt -> true
        | _ -> false)
      |> Option.map ~f:(fun (i, _) -> (i, key, nmt))
      |> Option.value ~default:(3000, key, nmt)
  | _ -> invalid_arg "Not 2 in mptt"
  *)

let prof_dir_example = "profc/type_out_%"

(*let open Oml.Statistics in *)
let analyze ?(genes=["A"; "B"; "C"]) ?(reads=[1;2]) ?(filters=[2;4;6;8;10;12]) ~prof_dir
  ?(suffix="_nuc.txt") () =
  List.fold_left genes ~init:[] ~f:(fun acc gene ->
    let train_samples = to_res gene `Train in
    List.fold_left reads ~init:acc ~f:(fun acc read ->
      List.fold_left filters ~init:acc ~f:(fun acc filter ->
        let to_filename sample = sprintf "%s_%d/%s_%d_%s%s" prof_dir filter sample read gene suffix in
        let positions =
          List.map train_samples ~f:(fun (sample, res) ->
            let filename = to_filename sample in
            sample, res, lookup_position (load_prohlatype filename) res)
        in
        (gene, read, filter, positions) :: acc)))

let display_analysis mp lst =
  printf "gene\tfilter\n";
  List.iter lst ~f:(fun (gene, _read, filter, m) ->
        printf "%s\t%d\t%s\n%!" gene filter (mp m))

let in_top ?(n=10.0) arr =
  Array.fold_left ~init:0.0 ~f:(fun s v -> if v < n then s +. 1. else s) arr
