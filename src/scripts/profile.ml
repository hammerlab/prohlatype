

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
      | Ok (s, nc) when s = g -> Some (k, nc)
      | Error e    -> invalid_arg e
      | _          -> None
    else
        None)

let position to_prohlatype_fname (key, nmt) =
  match nmt with
  | Nomenclature.Two (_, _) ->
      let fname = to_prohlatype_fname key in
      let as_tv =
        Csv.load ~separator:'\t' fname
        |> List.map ~f:(function
              | [_; a] -> Nomenclature.parse a
              | _ -> error "not 2")
      in
      List.findi as_tv ~f:(fun _i l ->
        match l with
        | Ok (_, nmt2) when nmt2 = nmt -> true
        | _ -> false)
      |> Option.map ~f:(fun (i, _) -> (i, key, nmt))
      |> Option.value ~default:(3000, key, nmt)
  | _ -> invalid_arg "Not 2 in mptt"

let prof_dir_example = "profc/type_out_%"

(*let open Oml.Statistics in *) *) *) *)
let analyze ?(genes=["A"; "B"; "C"]) ?(reads=[1;2]) ?(filters=[2;4;6;8;10;12]) ~prof_dir ~map () =
  List.fold_left genes ~init:[] ~f:(fun acc gene ->
    let train_samples = to_res gene `Train in
    List.fold_left reads ~init:acc ~f:(fun acc read ->
      List.fold_left filters ~init:acc ~f:(fun acc filter ->
        let to_filename k = sprintf "%s_%d/%s_%d_%s_nuc.txt" prof_dir filter k read gene in
        let rr = List.map train_samples ~f:(position to_filename) in
        let dt = List.map rr ~f:(fun (p,_,_) -> float p) |> Array.of_list in
        (gene, read, filter, (map dt)) :: acc)))

let display_analysis mp lst =
  printf "gene\tread\tfilter\n";
  List.iter lst ~f:(fun (gene, read, filter, m) ->
        printf "%s\t%d\t%d\t%s\n%!" gene read filter (mp m))

