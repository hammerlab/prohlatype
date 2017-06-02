
(*
#require "csv";
let res_fname = ...;;
*)

let keys =
  [ "120021"; "120074"; "120984"; "122879"; "123768"; "125785"; "126119"
  ; "126966"; "127051"; "127079"; "127245"; "128964"; "130445"; "160426"
  ; "160452"; "160917"; "161549"; "161702"; "162569"; "162636"; "162781"
  ; "164054"; "164115"; "164543"; "164964"; "165266"; "200220"; "201858"
  ; "203288"; "203482"; "204378"; "204587"; "205302"; "205855"; "205856"
  ; "205944"; "206039"; "206130"; "206139"; "206704"; "206708"; "206709"
  ; "206961"; "207956"; "207978"; "120013"; "120989"; "125881"; "127064"
  ; "129670"; "160908"; "161712"; "163581"; "164779"; "201063"; "203778"
  ; "205673"; "206038"; "206240"; "206876"
  ]

let _hdr :: data = Csv.load res_fname

module Smap = Map.Make(struct type t = string let compare = compare end)

let to_correct_map data =
  List.fold_left data ~init:Smap.empty ~f:(fun smap lst ->
    let k = List.nth_exn lst 4 in
    match Nomenclature.parse (List.nth_exn lst 1) with
    | Error e       -> invalid_arg e
    | Ok (gene, nc) ->
        let l = try Smap.find k smap with Not_found -> [] in
       Smap.add k ((gene, nc) :: l) smap)

let load_prohlatype fname =
  let separator = '\t' in
  let sep = sprintf "%c" separator in
  Csv.load ~separator:'\t' fname
  |> List.filter_map ~f:(function
    | a :: r -> begin match Nomenclature.parse a with
                | Ok (gene, res) -> Some (gene, (res, r))
                | Error e        -> eprintf "%s: %s, ignoring\n" fname e;
                                    None
                end
    | s      -> eprintf "in %s ignoring: %s\n" fname (String.concat ~sep s);
                None)
  |> List.fold_left ~init:[] ~f:(fun acc (gene, r) ->
        match List.Assoc.remove_and_get gene acc with
        | None -> (gene, [r]) :: acc
        | Some (nr, nacc) -> (gene, r :: nr) :: nacc)

let hla_likelihoods lst =
  List.map lst ~f:(fun (gene, llhd) ->
    gene
    , (List.map llhd ~f:(fun (allele, l) -> (allele, List.hd_exn l |> float_of_string))
      |> List.sort ~cmp:(fun (_a1, l1) (_a2, l2) -> compare l2 l1) (* descending *)))

type nr = Nomenclature.resolution * Nomenclature.suffix option
type gene = string

type record =
  { correct_alleles     : (gene * nr) list
  ; allele_likelihoods  : (gene * (nr * float) list) list
  }

let load_data ~dir ~fname_to_key =
  let correct = to_correct_map data in
  Sys.readdir dir
  |> Array.to_list
  |> List.fold_left ~init:Smap.empty ~f:(fun jmap file ->
      let k = fname_to_key file in
      let allele_likelihoods =
        load_prohlatype (Filename.concat dir file)
        |> hla_likelihoods
      in
      try
        let correct_alleles = Smap.find k correct in
        Smap.add k {correct_alleles; allele_likelihoods} jmap
      with Not_found ->
        eprintf "Didn't find key %s in correct results! Ignoring" k;
        jmap)

let summarize r ~gene ~metric ~f ~init =
  Smap.fold r ~init ~f:(fun ~key ~data acc ->
    match List.filter data.correct_alleles ~f:(fun (g,_) -> g = gene) with
    | []  -> eprintf "Didn't find gene %s in %s" gene key;
             acc
    | gls ->
        begin match List.Assoc.get gene data.allele_likelihoods with
        | None      -> eprintf "Didn't find likelihoods for gene %s in %s\n" gene key;
                       acc
        | Some llhd -> f acc key (List.map gls ~f:(fun (_gene, nr) -> metric nr llhd))
        end)

let lookup_position ?wrong nr nr_assoc =
  let wrong = Option.value wrong ~default:(List.length nr_assoc) in
  match (fst nr) with                     (* Match on just the name, ignore suffix *)
  | Nomenclature.Two (r1, r2) ->
      List.findi nr_assoc ~f:(fun _i (nr2, _v) ->
        match (fst nr2) with                               (* Also, ignore suffix. *)
        | Nomenclature.One _ -> invalid_argf "What! %s only 1 in Prohlatype output!"
                                      (Nomenclature.resolution_and_suffix_opt_to_string nr2)
        | Nomenclature.Two (s1, s2)
        | Nomenclature.Three (s1, s2, _)
        | Nomenclature.Four (s1, s2, _, _) -> r1 = s1 && r2 = s2)
      |> Option.value_map ~default:wrong ~f:fst
  | fnr ->
    invalid_argf "Unsupported resolution: %s "
      (Nomenclature.resolution_to_string fnr)

let best_positions r ~gene =
  summarize r ~gene ~metric:(lookup_position ~wrong:10000)
    ~init:[]
    ~f:(fun l k il -> (k, (List.sort ~cmp:compare il)) :: l)

(*
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
*)

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
