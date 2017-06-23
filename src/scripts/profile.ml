
(*
#require "csv";
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

let group_by_assoc l =
  let insert assoc (k, v) =
    match List.Assoc.remove_and_get k assoc with
    | None              -> (k,[v]) :: assoc
    | Some (cv, rassoc) -> (k, v ::cv) :: rassoc
  in
  List.fold_left ~init:[] ~f:insert l

type nr = Nomenclature.resolution * Nomenclature.suffix option
type gene = string
(* Since prohlatype can report the same likelihood for multiple alleles,
   we need to group them. We'll assign all of the values in a group the
   best, ie the index of the lowest allele in the same group.  *)
type result_group =
  { alleles     : nr list
  ; pos         : int
  ; likelihood  : float
  }


let hla_likelihoods lst =
  List.map lst ~f:(fun (gene, llhd) ->
    let likelihoods =
      List.map llhd ~f:(fun (allele, l) -> (List.hd_exn l |> float_of_string, allele))  (* parse likelihood value. *)
      |> group_by_assoc                                                           (* group by the same likelihoods *)
      |> List.sort ~cmp:(fun (l1,_a1) (l2, _a2) -> compare l2 l1)                               (* sort descending *)
      |> List.fold_left ~init:(0, []) ~f:(fun (pos, acc) (likelihood, alleles) ->
          (pos + List.length alleles, { alleles; pos; likelihood } :: acc))
      |> snd                                                                               (* discard index count. *)
      |> List.rev                                                                (* reverse to make lookup easier. *)
    in
    gene, likelihoods)

type record =
  { correct_alleles     : (gene * nr) list
  ; allele_likelihoods  : (gene * result_group list) list
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

let alleles_match nr1 nr2 =
  match (fst nr1) with                     (* Match on just the name, ignore suffix *)
  | Nomenclature.Two (r1, r2) ->
      begin match (fst nr2) with                               (* Also, ignore suffix. *)
        | Nomenclature.One _ -> invalid_argf "What! %s only 1 in Prohlatype output!"
                                        (Nomenclature.resolution_and_suffix_opt_to_string nr2)
        | Nomenclature.Two (s1, s2)
        | Nomenclature.Three (s1, s2, _)
        | Nomenclature.Four (s1, s2, _, _) -> r1 = s1 && r2 = s2
      end
  | fnr ->
      invalid_argf "Unsupported resolution: %s "
        (Nomenclature.resolution_to_string fnr)

let lookup_position ~wrong nr group_list =
  List.find group_list ~f:(fun rg ->
    List.exists rg.alleles ~f:(alleles_match nr))
  |> Option.value_map ~default:wrong ~f:(fun rg -> rg.pos)

let best_positions r ~gene =
  summarize r ~gene ~metric:(lookup_position ~wrong:10000)
    ~init:[]
    ~f:(fun l k il -> (k, (List.sort ~cmp:compare il)) :: l)

let best ~gene r k =
  let cr = best_positions r ~gene in
  let d = List.length cr in
  let n = List.length (List.filter cr ~f:(fun (_, l) -> List.hd_exn l < k)) in
  (n, d, (float n) /. (float d))
