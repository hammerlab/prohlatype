
(*
#require "csv";
*)

open Util

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

module Diploid = struct

  type t = Nomenclature.t * Nomenclature.t

  let init a1 a2 =
    if Nomenclature.compare a1 a2 <= 0 then
      (a1, a2)
    else
      (a2, a1)

  let lower_res_matches (l1, l2) (h1, h2) =
    Nomenclature.two_matches_full l1 h1 &&
    Nomenclature.two_matches_full l2 h2

  let one_lower_res_matches l (h1, h2) =
    Nomenclature.two_matches_full l h1 ||
    Nomenclature.two_matches_full l h2

end (* Diploid *)


let to_correct_map data =
  List.fold_left data ~init:StringMap.empty ~f:(fun smap lst ->
    let k = List.nth_exn lst 4 in
    match Nomenclature.parse (List.nth_exn lst 1) with
    | Error e       -> invalid_arg e
    | Ok (gene, nc) ->
        let l = try StringMap.find k smap with Not_found -> [] in
       StringMap.add k ((gene, nc) :: l) smap)
  |> StringMap.map ~f:(fun l ->
      group_by_assoc l
      |> List.filter_map ~f:(fun (gene, allele_lst) ->
            match allele_lst with
            | a1 :: a2 :: [] -> Some (gene, Diploid.init a1 a2)
            | _ -> eprintf "Odd! gene: %s has <2 or >2 alleles: %s\n"
                    gene
                    (string_of_list ~sep:" " ~f:Nomenclature.resolution_and_suffix_opt_to_string allele_lst);
                   None))

let tab_character = '\t'

let headers = ["Likelihood:"; "Zygosity:"; "Per Read:"]

let to_parse_state = function
  | "Likelihood:" :: [] -> `Likelihood
  | "Zygosity:" :: []   -> `Zygosity
  | "Per Read:" :: []   -> `PerRead
  | _                   -> `Unknown

type gene = string

type z =
  { d : Diploid.t
  ; l : float
  ; p : float
  }

(* Grouped zygosities. *)
type gz =
  { dl  : Diploid.t list          (* All pairs have the same likelihood/probability *)
  ; idx : int                     (* Index of the group. *)
  ; li  : float                   (* Likelihood *)
  ; pr  : float                   (* Probability *)
  }

let summarize_gz gz =
  let open Nomenclature in
  match gz.dl with
  | (a1, a2) :: [] -> sprintf "1 %s %s" (resolution_and_suffix_opt_to_string a1)
                          (resolution_and_suffix_opt_to_string a2)
  | (a1, a2) :: _  -> sprintf "%f (1/%d) %s\t%s" gz.pr (List.length gz.dl)
                          (resolution_and_suffix_opt_to_string a1)
                          (resolution_and_suffix_opt_to_string a2)
  | []             -> invalid_arg "can't summarize empty gz!"

let group_zygosities zlst =
  List.map zlst ~f:(fun z -> (z.l, z))
  |> group_by_assoc
  |> List.sort ~cmp:(fun (l1, _) (l2, _) -> compare l2 l1)  (* Make sure descending *)
  |> List.mapi ~f:(fun i (l, zlst) ->
      { dl  = List.map ~f:(fun z -> z.d) zlst
      ; idx = i
      ; li  = l
      ; pr  = (List.hd_exn zlst).p
      })

type prohlatype_output =
  { likelihoods : (gene * (Nomenclature.t * float) list) list
  ; zygosities  : (gene * gz list) list
  ; per_reads   : string list
  }

let load_prohlatype fname =
  let separator = tab_character in
  let sep = sprintf "%c" tab_character in
  let line l = String.concat ~sep l in
  let rec parse state lst =
    let searching_for, acc = state in
    match searching_for with
    | `Unknown  ->
        begin match to_parse_state lst with
        | `Unknown    -> eprintf "ignoring: %s\n" (line lst);
                         state
        | `Likelihood -> (`Likelihood ("", []), acc)
        | `Zygosity   -> (`Zygosity ("", []), acc)
        | `PerRead    -> (`PerRead [], acc)
        end
    | `Likelihood (gene, lacc) ->
        begin
          match lst with
          | a :: l :: [] ->
              begin
                match Nomenclature.parse a with
                | Ok (g, res) ->
                    if gene = "" || g = gene then
                      `Likelihood (g, (res, float_of_string l) :: lacc), acc
                    else (* g <> gene *)
                      invalid_argf "new gene %s in %s likelihood section!" g gene
                | Error e ->
                    eprintf "Parse error in likelihood: %s" e;
                    let nstate = `Unknown, (`Likelihood (gene, (List.rev lacc)) :: acc) in
                    parse nstate lst
              end
          | []  -> eprintf "Odd empty line\n";
                  `Unknown, (`Likelihood (gene, (List.rev lacc)) :: acc)
          | _   -> let nstate = `Unknown, (`Likelihood (gene, List.rev lacc)) :: acc in
                   parse nstate lst
        end
    | `Zygosity (gene, zacc) ->
        begin match lst with
        | a :: b :: l :: p :: [] ->
            begin match Nomenclature.parse a, Nomenclature.parse b with
            | Ok (ga, a1), Ok (gb, a2) ->
                if gene = "" || (ga = gene && gb = gene) then
                  let z =
                    { d = Diploid.init a1 a2
                    ; l = float_of_string l
                    ; p = float_of_string p
                    }
                  in
                  `Zygosity (ga, z :: zacc), acc
                else
                  invalid_argf "new gene %s, %s in %s zygosity section!" ga gb gene
            | Error e, _
            | _ , Error e ->
                eprintf "Parse error in zygosity: %s" e;
                let nstate = `Unknown, (`Zygosity (gene, List.rev zacc)) :: acc in
                parse nstate lst
            end
        | [] -> eprintf "Odd empty line\n";
                `Unknown, (`Zygosity (gene, List.rev zacc) :: acc)
        | _  -> let nstate = `Unknown, (`Zygosity (gene, List.rev zacc) :: acc) in
                parse nstate lst
        end
    | `PerRead pacc ->
        `PerRead ((line lst) :: pacc), acc
  in
  let final_state, acc =
    Csv.load ~separator fname
    |> List.fold_left ~init:(`Unknown, []) ~f:parse
  in
  let fs =
    match final_state with
    | `Unknown      -> eprintf "Finished as Unknown ?"; acc
    | `Likelihood _ -> eprintf "Finished as Likelihood ?"; (final_state :: acc)
    | `Zygosity _   -> eprintf "Finished as Zygosity ?"; (final_state :: acc)
    | `PerRead _    -> final_state :: acc
  in
  { likelihoods =
      List.filter_map fs ~f:(function
        | `Likelihood (gene, lst) -> Some (gene, lst)
        | _ -> None)
  ; zygosities  =
      List.filter_map fs ~f:(function
        | `Zygosity (gene, lst) -> Some (gene, (group_zygosities lst))
        | _ -> None)
  ; per_reads   =
    List.concat (List.filter_map fs ~f:(function
        | `PerRead lst -> Some lst
        | _ -> None))
  }

type record =
  { ca : (gene * Diploid.t) list
  ; po : prohlatype_output
  }

let load_data ~dir ~fname_to_key =
  let correct = to_correct_map data in
  Sys.readdir dir
  |> Array.to_list
  |> List.fold_left ~init:StringMap.empty ~f:(fun jmap file ->
      let k = fname_to_key file in
      let po = load_prohlatype (Filename.concat dir file) in
      try
        let ca = StringMap.find k correct in
        StringMap.add k {ca; po} jmap
      with Not_found ->
        eprintf "Didn't find key %s in correct results! Ignoring" k;
        jmap)


let classify correct_d glz =
  let first = summarize_gz (List.hd_exn glz) in
  let rec loop = function
    | []       -> `Didn'tFindIt
    | gz :: tl ->
        if List.exists gz.dl ~f:(Diploid.lower_res_matches correct_d) then begin
          if gz.idx = 0 then begin
            let n = List.length gz.dl in
            if n = 1 then
              `PerfectFirst
            else
              `OneAmong (n, gz.pr)
          end else
            `Found gz
        end else
          loop tl
  in
  loop glz

let classify1 a1 glz =
  let rec loop = function
    | []       -> `Didn'tFindIt
    | gz :: tl ->
        if List.exists gz.dl ~f:(Diploid.one_lower_res_matches a1) then
          `Single gz
        else
          loop tl
  in
  loop glz

let lookup_result ~gene r =
  match List.Assoc.get gene r.ca with
  | None  -> invalid_argf "Missing correct gene: %s type" gene
  | Some d  ->
      begin match List.Assoc.get gene r.po.zygosities with
      | None      -> invalid_argf "Missing zygosities output for gene %s" gene
      | Some zlst ->
          begin match classify d zlst with
          | `Didn'tFindIt  ->
              let d1, d2 = d in
              if d1 = d2 then
                match classify1 d1 zlst with
                | `Single gz1   -> d, `NotHomozygous gz1
                | `Didn'tFindIt -> d, `Didn'tFindIt
              else
                let m1 = classify1 d1 zlst in
                let m2 = classify1 d2 zlst in
                begin match m1, m2 with
                | `Single gz1,   `Single gz2    -> d, `Mixed (gz1, gz2)
                | `Single _,     `Didn'tFindIt  -> d, m1
                | `Didn'tFindIt, `Single _      -> d, m2
                | `Didn'tFindIt, _              -> d, m1
              end
          | s -> d, s
          end
      end

type sample_summary =
  { perfect   : int * (string list)
  ; one_among : int * (string list)
  ; found     : int * (string list)
  ; not_hmz   : int * (string list)
  ; single    : int * (string list)
  ; mixed     : int * (string list)
  ; didn't    : int * (string list)
  }

let to_ss l =
  let p = List.filter_map l ~f:(function | (s, (_, `Perfect)) -> Some s | _ -> None) in
  let o = List.filter_map l ~f:(function | (s, (_, `OneAmong _)) -> Some s | _ -> None) in
  let f = List.filter_map l ~f:(function | (s, (_, `Found _)) -> Some s | _ -> None) in
  let n = List.filter_map l ~f:(function | (s, (_, `NotHomozygous _)) -> Some s | _ -> None) in
  let s = List.filter_map l ~f:(function | (s, (_, `Single _)) -> Some s | _ -> None) in
  let m = List.filter_map l ~f:(function | (s, (_, `Mixed _)) -> Some s | _ -> None) in
  let d = List.filter_map l ~f:(function | (s, (_, `Didn'tFindIt)) -> Some s | _ -> None) in
  { perfect   = List.length p, p
  ; one_among = List.length o, o
  ; found     = List.length f, f
  ; not_hmz   = List.length n, n
  ; single    = List.length s, s
  ; mixed     = List.length m, m
  ; didn't    = List.length d, d
  }


(*
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
      |> Util.group_by_assoc                                                      (* group by the same likelihoods *)
      |> List.sort ~cmp:(fun (l1,_a1) (l2, _a2) -> compare l2 l1)                               (* sort descending *)
      |> List.fold_left ~init:(0, []) ~f:(fun (pos, acc) (likelihood, alleles) ->
          (pos + List.length alleles, { alleles; pos; likelihood } :: acc))
      |> snd                                                                               (* discard index count. *)
      |> List.rev                                                                (* reverse to make lookup easier. *)
    in
    gene, likelihoods)

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

let join_runs r1 r2 =
  let flatten_to_alleles lst =
    List.fold_left lst ~init:[] ~f:(fun acc (rg : result_group) ->
      acc @ List.map rg.alleles ~f:(fun nr -> nr, rg.likelihood))
      |> List.sort ~cmp:compare   (* sort by allele *)
  in
  let regroup acc =
    Util.group_by_assoc acc                                                     (* group by the same likelihoods *)
    |> List.sort ~cmp:(fun (l1,_a1) (l2, _a2) -> compare l2 l1)                               (* sort descending *)
    |> List.fold_left ~init:(0, []) ~f:(fun (pos, acc) (likelihood, alleles) ->
      (pos + List.length alleles, { alleles; pos; likelihood } :: acc))
    |> snd                                                                               (* discard index count. *)
    |> List.rev                                                                (* reverse to make lookup easier. *)
  in
  Smap.merge r1 r2 ~f:(fun k v1 v2 ->
    match v1, v2 with
    | None, None  -> assert false
    | None, _     -> printf "%s missing first record\n" k; v2
    | _,    None  -> printf "%s missing second record\n" k; v1
    | Some a1, Some a2 ->
        assert (a1.correct_alleles = a2.correct_alleles);
        let l =
          List.map2 (List.sort ~cmp:compare a1.allele_likelihoods)
                    (List.sort ~cmp:compare a2.allele_likelihoods) (* align genes *)
            ~f:(fun (g1, g1list) (g2, g2list) ->
                  assert (g1 = g2);
                  let al1 = flatten_to_alleles g1list in
                  let al2 = flatten_to_alleles g2list in
                  let al =
                    List.map2 al1 al2 ~f:(fun (a1, l1) (a2, l2) ->
                      assert (a1 = a2);
                      (l1 +. l2), a1)
                    |> regroup
                  in
                  g1, al)
        in
        Some {a1 with allele_likelihoods = l})

(* Running time analysis *)
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
      *)
