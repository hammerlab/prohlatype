(* Search for similarity (as measured by Levenshtein distance) across HLA
 * alleles. For the purpose of this work, we'll only search through Exons,
 * ( we have more alleles vs constraining this search somehow).
 *)

open Prohlatype
open Biology

let loci_to_check =
  let open Nomenclature in
  [ A ; B ; C ; DMA ; DMB ; DOA ; DOB ; DPA1 ; DPA2 ; DPB1 ; DPB2 ; DQA1 ; DQB1
  ; DRA
  (* DRB1; DRB3 ; DRB4     See hack below. *)
  ; E ; F ; G ; HFE ; H ; J ; K ; L ; MICA ; MICB
  (* ; P *)
  ; TAP1 ; TAP2 ; T ; V ; W ; Y
  ]

(*
LOCI="../foreign/IMGTHLA/alignments/A_gen.txt
../foreign/IMGTHLA/alignments/B_gen.txt
../foreign/IMGTHLA/alignments/C_gen.txt
../foreign/IMGTHLA/alignments/E_gen.txt
../foreign/IMGTHLA/alignments/F_gen.txt
../foreign/IMGTHLA/alignments/G_gen.txt
../foreign/IMGTHLA/alignments/HFE_gen.txt
../foreign/IMGTHLA/alignments/H_gen.txt
../foreign/IMGTHLA/alignments/J_gen.txt
../foreign/IMGTHLA/alignments/K_gen.txt
../foreign/IMGTHLA/alignments/L_gen.txt
../foreign/IMGTHLA/alignments/MICA_gen.txt
../foreign/IMGTHLA/alignments/P_gen.txt
../foreign/IMGTHLA/alignments/TAP1_gen.txt
../foreign/IMGTHLA/alignments/TAP2_gen.txt
../foreign/IMGTHLA/alignments/T_gen.txt
../foreign/IMGTHLA/alignments/V_gen.txt
../foreign/IMGTHLA/alignments/W_gen.txt"
*)

let locus_to_nuc_fname l =
  Nomenclature.show_locus l
  |> sprintf "%s_nuc"
  |> Common.to_alignment_file

let load_nuc l =
  let fname = locus_to_nuc_fname l in
  MSA.Parser.from_file fname

let nuc_mps =
  List.map loci_to_check ~f:load_nuc

let drb_nuc_mp =
  MSA.Parser.from_file (Common.to_alignment_file "DRB_nuc")

let p_gen_mp =
  MSA.Parser.from_file (Common.to_alignment_file "P_gen")

let funwrap = function
  | Ok o -> o
  | Error e -> failwithf "%s" e

let mp_to_seq mp =
  let open MSA in
  let open Parser in
  let no_sequences = List.filter ~f:(fun a -> not (is_sequence a)) in
  let als =
    List.map mp.alt_elems ~f:(fun { allele; seq; _} ->
      let seqs = allele_sequences ~reference:mp.ref_elems ~allele:seq |> funwrap in
      let rels = Segments.distances ~reference:mp.ref_elems ~allele:seq |> funwrap in
      ( allele, List.map2 seqs rels ~f:(fun s r -> (s, r))))
  in
  let rp =
    let allele = no_sequences mp.ref_elems in
    let s = allele_sequences ~reference:mp.ref_elems ~allele |> funwrap in
    let r = Segments.distances ~reference:mp.ref_elems ~allele |> funwrap in
    ( mp.reference, List.map2 s r ~f:(fun s r -> (s, r)))
  in
  List.sort ~cmp:compare (rp :: als)

let p_nuc_mp =
  let open MSA in
  let open Parser in
  let just_xons s =
    Boundaries.grouped s
    |> funwrap
    |> List.filter ~f:(function (Boundaries.{ label = Exon _; _ }, _) -> true | _ -> false)
    |> Boundaries.ungrouped
  in
  { p_gen_mp with ref_elems = just_xons p_gen_mp.ref_elems
                ; alt_elems = List.map p_gen_mp.alt_elems ~f:(fun a ->
                                { a with seq = just_xons a.seq })
  }

let mp_and_seqs =
  let open MSA.Parser in
  drb_nuc_mp :: p_nuc_mp :: nuc_mps
  |> List.sort ~cmp:(fun mp1 mp2 -> compare mp1.reference mp2.reference)
  |> List.map ~f:(fun mp -> mp, mp_to_seq mp)

let print_exon_sequences mp_and_seq =
  let open MSA in
  let open Parser in
  List.iter mp_and_seq ~f:(fun (mp, seqs, _distances) ->
    let ref_lst = assoc_exn mp.reference seqs in
    printf "%s\n" mp.reference;
    List.iter ref_lst ~f:(fun (bm, seqs) ->
        printf "\t%s:%d\n"
          (Gene_region.to_string bm.MSA.Boundaries.label)
          (String.length seqs)))

module IntSet = Set.Make (struct type t = int let compare = compare end)

let sizes seqs =
  let number_of_exons =
    List.hd_exn seqs |> snd |> List.length
  in
  let init =
    Array.init number_of_exons ~f:(fun _ -> IntSet.empty) |> Array.to_list
  in
  List.fold_left seqs ~init
    ~f:(fun sl (_, sls) ->
        let lengths = List.map sls ~f:(fun (_bm, s) -> String.length s) in
        List.map2 lengths sl ~f:IntSet.add)

module Naive = struct

  let mismatches_between ~end_ s1 o1 s2 o2 =
    let rec loop sum i =
      if i = end_ then
        sum
      else
        let c1 = String.get s1 ~index:(o1 + i) in
        let c2 = String.get s2 ~index:(o2 + i) in
        if c1 = c2 then
          loop sum (i + 1)
        else
          loop (sum + 1) (i + 1)
    in
    loop 0 0

  let mismatches_over ~smaller ~smaller_len ~larger ~larger_len =
    let last_large_start_offset = larger_len - smaller_len in
    let mb lo = mismatches_between ~end_:smaller_len smaller 0 larger lo in
    let rec loop mn sp =
      if sp > last_large_start_offset then
        mn
      else
        loop (min mn (mb sp)) (sp + 1)
    in
    loop (mb 0) 0

  let mismatches s1 s2 =
    let n1 = String.length s1 in
    let n2 = String.length s2 in
    if n1 < n2 then
      let m = mismatches_over ~smaller:s1 ~smaller_len:n1 ~larger:s2 ~larger_len:n2 in
      (n2 - n1) + m
    else if n1 = n2 then
      mismatches_between ~end_:(n1 - 1) s1 0 s2 0
    else (* n1 > n2 *)
      let m = mismatches_over ~smaller:s2 ~smaller_len:n2 ~larger:s1 ~larger_len:n1 in
      (n1 - n2) + m


  let search ?m ?n predicate extreme_value a1 a2 =
    let m = Option.value m ~default:(Array.length a1) in
    let n = Option.value n ~default:(Array.length a2) in
    let rec loop i j mm p =
      if i = m then
        mm, p
      else if j = n then
        loop (i + 1) 0 mm p
      else
        let al1, s1 = a1.(i) in
        let al2, s2 = a2.(j) in
        (*printf "comparing %s vs %s\n%!" al1 al2;*)
        let msms = mismatches s1 s2 in
        if predicate msms mm then
          loop i (j + 1) msms (al1, al2)
        else
          loop i (j + 1) mm p
    in
    loop 0 0 extreme_value ("","")

  let minimum = search (fun x y -> x < y) max_int
  let maximum = search (fun x y -> x > y) min_int

end (* Naive *)

let general_class_1_search_loci =
  Nomenclature.[ A ; B ; C ; E ; F ; G ; H ; J ; K ; L ; Y ]

let with_classII_search_loci =
  general_class_1_search_loci @
    Nomenclature.[ DMA ; DMB ; DOA ; DOB ; DPA1 ; DPA2 ; DPB1 ; DPB2 ; DQA1
                 ; DQB1 ; DRA ; DRB1 ]


(* Note that Mas.Parser starts the exon counting at 1. *)
let exons_to_check =
  [ 2; 3]

module Levenshtein = struct

  let min3 (a : int) (b : int) c =
    min a (min b c)

  let staged s =
    let l1 = String.length s in
    let n1 = l1 + 1 in
    let arr =
      Array.init 2 ~f:(fun _i ->
          Array.init n1 ~f:(fun _j -> 0))
    in
    let first_row = Array.init n1 ~f:(fun j -> j) in
    let clear () =
      Array.blit ~src:first_row ~src_pos:0 ~dst:arr.(0) ~dst_pos:0 ~len:n1;
      arr.(1).(0) <- 1
    in
    fun t ->
      clear ();
      let l2 = String.length t in
      let rec loop v i j =
        if i >= l2 then begin
          v
        end else if j >= l1 then begin
          arr.(i mod 2).(0) <- (i + 1);
          loop v (i + 1) 0
        end else
          let c_i = String.get_exn t ~index:(i) in
          let c_j = String.get_exn s ~index:(j) in
          let cp = if c_j = c_i then 0 else 1 in
          let ip1 = i + 1 in
          let jp1 = j + 1 in
          let u = min3
                    (arr.(i mod 2).(j) + cp)
                    (arr.(i mod 2).(jp1) + 1)
                    (arr.(ip1 mod 2).(j) + 1)
          in
          arr.(ip1 mod 2).(jp1) <- u;
          loop u i jp1
      in
      loop min_int 0 0

end (* Levenshtein *)

let all_distances a =
  let n = Array.length a in
  Array.mapi a ~f:(fun i (name, seq1) ->
    let sl1 = Levenshtein.staged seq1 in
    let dm =
      Array.map (Array.sub a ~pos:(i + 1) ~len:(n-i-1)) ~f:(fun (n2, seq2) ->
        let d = sl1 seq2 in
        printf "distance for %s vs %s: %d\n%!" name n2 d;
        d)
    in
    name, dm)

let spec_xon exon seqs =
  let open MSA in
  let e = Gene_region.Exon exon in
  List.filter_map seqs ~f:(fun (allele, exon_list) ->
    match Nomenclature.parse allele with
    | Ok (_, (_, None)) ->  (* Doesn't have a qualifier *)
        let xo =
          List.find_map exon_list ~f:(fun ((bm, s), r) ->
            if bm.Boundaries.label = e then begin
              match r.Segments.relationship with
              | Segments.Full _    -> Some s
              | Segments.Missing
              | Segments.Partial _ -> None
            end else
              None)
        in
        begin match xo with
        | None    -> printf "dropping %s because it is missing or partial exon: %d\n" allele exon; None
        | Some "" -> printf "dropping %s because the exon sequence is empty.\n" allele; None
        | Some s  -> Some (allele, s)
        end
    | _                 ->
      printf "dropping %s because it has an extension\n" allele;
      None)
  |> Array.of_list

let just_that_exon exon loci =
  let open MSA in
  List.map loci ~f:(fun locus ->
      let _, seqs = List.find_exn mp_and_seqs ~f:(fun (mp, _) -> mp.Parser.locus = locus) in
      spec_xon exon seqs)

let minimum_of arr1 arr2 =
  let n1 = Array.length arr1 in
  let n2 = Array.length arr2 in
  let d_min = ref max_int in
  let a1 = ref (fst arr1.(0)) in
  let a2 = ref (fst arr2.(0)) in
  for i = 0 to n1 - 1 do
    let sl = Levenshtein.staged (snd arr1.(i)) in
    for j = 0 to n2 - 1 do
      let d = sl (snd arr2.(j)) in
      if d < !d_min then begin
        d_min := d;
        a1 := fst arr1.(i);
        a2 := fst arr2.(j)
      end
    done
  done;
  !d_min, !a1, !a2

let max_dist allele ma =
  match Array.findi ma ~f:(fun (a,_) -> a = allele) with
  | None   -> invalid_argf "Couldn't find %s" allele
  | Some i ->
      let _ca, cd = ma.(i) in
      let ba = ref "" in
      let bd = ref min_int in
      for j = 0 to i - 1 do
        let al, da = ma.(j) in
        let d = da.(i-j-1) in
        if d > !bd then begin
          bd := d;
          ba := al
        end
      done;
      for j = 0 to Array.length cd - 1 do
        if cd.(j) > !bd then begin
          bd := cd.(j);
          ba := fst ma.(i + j + 1)
        end
      done;
      !ba, !bd

let grand_search seq_lst =
  let all_distances = List.map seq_lst ~f:(fun arr ->
    let alleles = Array.map arr ~f:fst in
    let ad = all_distances arr in
    alleles, ad, arr)
  in
  List.iteri all_distances ~f:(fun i (_alleles, ad, arr) ->
    let others = List.filteri all_distances ~f:(fun j _ -> i <> j) in
    List.iter others ~f:(fun (_, bd, brr) ->
      let min_d, a_a, b_a = minimum_of arr brr in
      let a_alt, max_dist_from_a = max_dist a_a ad in
      let b_alt, max_dist_from_b = max_dist b_a bd in
      printf "%s\t%s\t%d\t%s\t%d\t%s\t%d\n"
        a_a b_a min_d a_alt max_dist_from_a b_alt max_dist_from_b))

let grand_search_par seq_lst =
  let all_distances =
    Parmap.parmap
      (fun arr ->
        let alleles = Array.map arr ~f:fst in
        let ad = all_distances arr in
        alleles, ad, arr)
      (Parmap.L seq_lst)
  in
  let print_me =
    Parmap.parmapi
      (fun i (_alleles, ad, arr) ->
        let others = List.filteri all_distances ~f:(fun j _ -> i <> j) in
        List.map others ~f:(fun (_, bd, brr) ->
          let min_d, a_a, b_a = minimum_of arr brr in
          let a_alt, max_dist_from_a = max_dist a_a ad in
          let b_alt, max_dist_from_b = max_dist b_a bd in
          sprintf "%s\t%s\t%d\t%s\t%d\t%s\t%d\n"
            a_a b_a min_d a_alt max_dist_from_a b_alt max_dist_from_b))
      (Parmap.L all_distances)
  in
  List.iter print_me ~f:(List.iter ~f:(print_line stdout))


let () =
  if not (!Sys.interactive) then begin
    let exon =
      try if Array.length Sys.argv > 1 then int_of_string (Sys.argv.(1)) else 2
      with _ -> 2
    in
    printf "for exon %d\n" exon;
    printf "allele1\tallele2\tmin distance\tallele1_alt\tmax distance\tallele2_alt\tmax distance\n";
    (*grand_search_par (just_that_exon exon general_class_1_search_loci)*)
    grand_search_par (just_that_exon exon with_classII_search_loci)
  end
