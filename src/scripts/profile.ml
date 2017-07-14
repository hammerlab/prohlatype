
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

let tab_character = '\t'

let load_prohlatype fname =
  let separator = tab_character in
  let sep = sprintf "%c" separator in
  Csv.load ~separator fname
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

let get_times fname =
  let ic = open_in fname in
  let rec loop acc =
    try
      loop (input_line ic :: acc)
    with End_of_file ->
      List.take acc 3
  in
  let last_three = loop [] in
  close_in ic;
  last_three

let load_times ~dir ~fname_to_key =
  Sys.readdir dir
  |> Array.fold_left ~init:[] ~f:(fun acc file ->
      let k = fname_to_key file in
      let t = get_times (Filename.concat dir file) in
      (k, t) :: acc)

let user_times =
  List.map ~f:(fun (sample_name, lst) ->
    sample_name,  Scanf.sscanf (List.nth_exn lst 1) "user\t%dm%fs" (fun m s -> (float m *. 60. +. s))) ;;

(*** Parsing mapped files *)
let key_to_gene = function
  | "AF_A_gen_true" -> "A"
  | "AF_B_gen_true" -> "B"
  | "AF_C_gen_true" -> "C"
  | s               -> invalid_argf "unrecognized gene: %s" s

type map_res =
  { alleles   : (string * float) list
  ; positions : (int * float) list
  }

let split_on_semicolon s =
  String.split ~on:(`Character ';') s

let split_on_colon s =
  String.split ~on:(`Character ':') s

let parse_allele_and_likelihood s =
  split_on_colon s
  |> List.rev
  |> function
      | lklhd :: reverse_me_for_allele ->
          (List.rev reverse_me_for_allele |> String.concat ~sep:":"),
          (float_of_string lklhd)
      | _   ->
          invalid_argf "cant parse alleles: %s" s

let parse_positions s =
  Scanf.sscanf s "%d:%f" (fun p l -> p, l)

let parse_res lst =
  match lst with
  | al :: pl :: tl ->
      { alleles = split_on_semicolon al |> List.map ~f:parse_allele_and_likelihood
      ; positions = split_on_semicolon pl |> List.map ~f:parse_positions
      } , tl
  | _ -> invalid_arg "Not enough elements"

let rec parse_gene lst =
  match lst with
  | "C" :: tl ->
      let map, tl = parse_res tl in
      (`C map), tl
  | "R" :: tl ->
      let map, tl = parse_res tl in
      (`R map), tl
  | "F" :: msg :: tl ->
      (`F msg), tl
  | s :: _ ->
      invalid_argf "unrecongized output: %s" s
  | []      ->
      invalid_arg "empty"


let parse_map_line_result l =
  try
    let split = String.split ~on:(`Character tab_character) l in
    match split with
    | "" :: gene_s :: rest ->
        let gene = key_to_gene gene_s in
        let r1, rest = parse_gene rest in
        let r2, rest = parse_gene rest in
        if rest = [] then
          Ok (gene, (r1, r2))
        else
          error "extra stuff"
    | _ -> error "unrecognized map line: %s" l
  with Invalid_argument e -> error "%s" e

let load_maps fname =
  let ic = open_in fname in
  let rec until_reads () =
    let s = input_line ic in
    if String.compare_substring ("Reads", 0, 5) (s, 0, 5) = 0 then
      grab_reads []
    else
      until_reads ()
  and grab_reads acc =
    try
      let read_line = input_line ic in
      let line1 = input_line ic in
      let line2 = input_line ic in
      let line3 = input_line ic in
      match
        parse_map_line_result line1 >>= fun l1 ->
          parse_map_line_result line2 >>= fun l2 ->
            parse_map_line_result line3 >>= fun l3 -> Ok [l1; l2; l3]
      with
      | Ok l -> grab_reads ((read_line, l) :: acc)
      | Error e ->
          eprintf "for read %s: %s" read_line e;
          acc
    with End_of_file ->
      List.rev acc
  in
  until_reads ()

let tcl (o, _) =
  match o with
  | `C mp -> List.hd_exn mp.alleles |> snd
  | `R mp -> List.hd_exn mp.alleles |> snd
  | `F _  -> neg_infinity

let get_gene gene glst =
  List.Assoc.get gene glst
  |> Option.value_exn ~msg:(sprintf "missing: %s" gene)

let mbest (_, l) =
  let la = get_gene "A" l in
  let lb = get_gene "B" l in
  let lc = get_gene "C" l in
  let ta = tcl la in
  let tb = tcl lb in
  let tc = tcl lc in
  if ta > tb then begin
    if ta > tc then
      `A la
    else
      `C lc
  end else begin
    if tb > tc then
      `B lb
    else
      `C lc
  end

let mjust a m =
  List.filter_map m ~f:(fun q ->
    match a, mbest q with
    | `A, `A (`R mp, _) -> Some (fst q, mp, true)
    | `A, `A (`C mp, _) -> Some (fst q, mp, false)
    | `A, `A (`F _, _)  -> invalid_arg "mapped to filter!"
    | `A,  _            -> None

    | `B, `B (`R mp, _) -> Some (fst q, mp, true)
    | `B, `B (`C mp, _) -> Some (fst q, mp, false)
    | `B, `B (`F _, _)  -> invalid_arg "mapped to filter!"
    | `B,  _            -> None

    | `C, `C (`R mp, _) -> Some (fst q, mp, true)
    | `C, `C (`C mp, _) -> Some (fst q, mp, false)
    | `C, `C (`F _, _)  -> invalid_arg "mapped to filter!"
    | `C,  _            -> None)

let msorted l =
  List.sort ~cmp:(fun (_, mp1, _) (_, mp2, _) ->
      compare (List.hd_exn mp2.alleles |> snd) (List.hd_exn mp1.alleles |> snd)) l
