

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

let both g = to_res g `Train @ to_res g `Valid

let load_prohlatype fname =
  Csv.load ~separator:'\t' fname
  |> List.map ~f:(function
    | [a; _] -> Nomenclature.parse a
    | _      -> error "not 2")
  |> fun l ->
      match List.find_map l ~f:(function | Error e -> Some e | _ -> None) with
      | Some e -> invalid_argf "Error parsing: %s" e
      | None -> List.map l ~f:unwrap_ok

let lookup_position tsv = function
  | Nomenclature.Two (_, _) as nmt ->
      List.findi tsv ~f:(fun _i l ->
        match l with
        | Ok (_, (nmt2, _)) when nmt2 = nmt -> true
        | _ -> false)
      |> Option.map ~f:(fun (i, _) -> i)
      |> Option.value ~default:max_int
  | _ -> invalid_arg "Not 2 in mptt"

let to_sample_and_locus s =
  match String.split s ~on:(`Character '_') with
  | sa :: _ :: l :: [] -> sa, l
  | _   -> invalid_arg ("to_sample_and_locus: " ^ s)

let dir = "performance/2017_04_22/filtered/"

let as_good_as_first f s =
  let open Nomenclature in
  match f, s with
  | One o          , One oo
  | One o          , Two (oo,_)
  | One o          , Three (oo,_,_)
  | One o          , Four (oo,_,_,_)    -> o = oo
  | Two (o,t)      , Two (oo,tt)
  | Two (o,t)      , Three (oo,tt,_)
  | Two (o,t)      , Four (oo,tt,_,_)   -> o = oo && t = tt
  | Three (o,t,r)  , Three (oo,tt,rr)
  | Three (o,t,r)  , Four (oo,tt,rr,_)  -> o = oo && t = tt && r = rr
  | Four (o,t,r,f) , Four (oo,tt,rr,ff) -> o = oo && t = tt && r = rr && f = ff
  | _              , _                  -> false

let do_it dir gene =
  let correct = both gene in
  Sys.readdir dir
  |> Array.to_list
  |> List.filter_map ~f:(fun f ->
      let sample, sample_gene = to_sample_and_locus f in
      if gene <> sample_gene then
        None
      else
        let res =
          load_prohlatype (Filename.concat dir f) 
          |> List.map ~f:(fun (_a, (v, _)) -> v)
        in
        List.filter correct ~f:(fun (s,_) -> s = sample)
        |> List.map ~f:(fun (_, nm) ->
              nm, List.findi res ~f:(fun _ v -> as_good_as_first nm v))
        |> fun l -> Some (sample, l))

let disp ?gene res =
  List.map res ~f:(fun (a,l) ->
    (a, List.map l 
      ~f:(fun (nm1, o) ->
          Nomenclature.resolution_to_string ?gene nm1
          , Option.map o ~f:(fun (i, nm2) ->
              i, Nomenclature.resolution_to_string ?gene nm2))))

let count_below k = 
  List.fold_left ~init:0 ~f:(fun c (_, l) ->
    List.fold_left l ~init:c ~f:(fun c (_, o) ->
      match o with | Some (v, _) when v < k -> c + 1  | _ -> c)) ;; 

  (*
(*let open Oml.Statistics in *)
let analyze ?(genes=["A"; "B"; "C"]) ?(reads=[1;2]) ~prof_dir ?(suffix="_nuc.txt") () =
  List.fold_left genes ~init:[] ~f:(fun acc gene ->
    let train_samples = to_res gene `Train in
    List.fold_left reads ~init:acc ~f:(fun acc read ->
        let to_filename sample = sprintf "%s_%d/%s_%d_%s%s" prof_dir sample read gene suffix in
        let positions =
          List.map train_samples ~f:(fun (sample, res) ->
            let filename = to_filename sample in
            sample, res, lookup_position (load_prohlatype filename) res)
        in
        (gene, read, filter, positions) :: acc))

let display_analysis mp lst =
  printf "gene\tfilter\n";
  List.iter lst ~f:(fun (gene, _read, filter, m) ->
        printf "%s\t%d\t%s\n%!" gene filter (mp m))

let in_top ?(n=10.0) arr =
  Array.fold_left ~init:0.0 ~f:(fun s v -> if v < n then s +. 1. else s) arr

  *)
