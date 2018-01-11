(* Old and deprecated of better levenshtein search string. *)
open Util

let (//) = Filename.concat
let imgthla_dir =
  try Sys.getenv "IMGTHLA_DIR"
  with _ -> "../foreign/IMGTHLA"

let to_alignment_file f = imgthla_dir // "alignments" // (f ^ ".txt")

let anuc = MSA.Parser.from_file (to_alignment_file "A_nuc")
let bnuc = MSA.Parser.from_file (to_alignment_file "B_nuc")
let cnuc = MSA.Parser.from_file (to_alignment_file "C_nuc")

let mp_to_seq mp =
  let open MSA in
  let no_sequences = List.filter ~f:(fun a -> not (is_sequence a)) in
  let als =
    List.map mp.Parser.alt_elems ~f:(fun { Parser.allele; seq; _} ->
      ( allele
      , allele_sequences ~reference:mp.Parser.ref_elems ~allele:seq))
  in
  let rp =
    ( mp.Parser.reference
    , allele_sequences ~reference:mp.Parser.ref_elems
        ~allele:(no_sequences mp.Parser.ref_elems))
  in
  rp :: als

let alst = mp_to_seq anuc
let blst = mp_to_seq bnuc
let clst = mp_to_seq cnuc

let reference_length r l n =
  List.Assoc.get r l
  |> Option.value_exn ~msg:"no reference"
  |> fun l -> List.nth_exn l n
  |> fun (_, s) -> String.length s

let a2_rl = reference_length "A*01:01:01:01" alst 1 
let () = printf "reference length a2_rl: %d\n%!" a2_rl
let a3_rl = reference_length "A*01:01:01:01" alst 2 
let () = printf "reference length a3_rl: %d\n%!" a3_rl

let b2_rl = reference_length "B*07:02:01:01" blst 1
let () = printf "reference length b2_rl: %d\n%!" b2_rl
let b3_rl = reference_length "B*07:02:01:01" blst 2 
let () = printf "reference length b3_rl: %d\n%!" b3_rl

let c2_rl = reference_length "C*01:02:01:01" clst 1
let () = printf "reference length c2_rl: %d\n%!" c2_rl
let c3_rl = reference_length "C*01:02:01:01" clst 2
let () = printf "reference length c3_rl: %d\n%!" c3_rl


let spec_xon rl l n =
  List.filter_map l ~f:(fun (a, ll) ->
    List.nth ll n |> Option.bind ~f:(fun (_, s) ->
      if String.length s = rl then begin
        match Nomenclature.parse a with
        | Ok (_, (_, None)) -> Some (a, s)
        | _                 -> None
      end else
        None))

let a2 = spec_xon a2_rl alst 1 |> Array.of_list
let a3 = spec_xon a3_rl alst 2 |> Array.of_list

let b2 = spec_xon b2_rl blst 1 |> Array.of_list
let b3 = spec_xon b3_rl blst 2 |> Array.of_list

let c2 = spec_xon c2_rl clst 1 |> Array.of_list
let c3 = spec_xon c3_rl clst 2 |> Array.of_list

let mismatches s1 s2 =
  let n1 = String.length s1 in
  let n2 = String.length s2 in
  let m = min n1 n2 in
  let x = max n1 n2 - m in
  let rec loop index s =
    if index = m then s else
      if String.get s1 ~index = String.get s2 ~index then
        loop (index + 1) s
      else
        loop (index + 1) (s + 1)
  in
  loop 0 x

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

let disp (d, (a, b)) =
  printf "%s - %s : %d\n%!" a b d

let () =
  printf "maximum distances\n%!";
  maximum a2 a2 |> disp;
  maximum b2 b2 |> disp;
  maximum c2 c2 |> disp;
  maximum a3 a3 |> disp;
  maximum b3 b3 |> disp;
  maximum c3 c3 |> disp;

  printf "minimum distances\n%!";
  minimum a2 b2 |> disp;
  minimum a2 c2 |> disp;
  minimum b2 c2 |> disp;
  minimum a3 b3 |> disp;
  minimum a3 c3 |> disp;
  minimum b3 c3 |> disp

(* Output for  5800ab43fe4afd66825db1e574176370558a94bf 3.29.0.1 IMGT: 
reference length a2_rl: 270
reference length a3_rl: 276
reference length b2_rl: 270
reference length b3_rl: 276
reference length c2_rl: 270
reference length c3_rl: 276
maximum distances
A*02:81 - A*11:262 : 26
B*73:02 - B*15:95 : 36
C*04:01:10 - C*16:85 : 38
A*01:67:02 - A*02:339 : 27
B*45:03 - B*07:64 : 29
C*05:152 - C*03:296 : 21
minimum distances
A*25:01:10 - B*07:81 : 17
A*03:01:59 - C*14:27 : 20
B*44:151 - C*16:85 : 0
A*32:89 - B*27:90:04 : 9
A*32:30:02 - C*16:67 : 13
B*35:205 - C*15:15 : 2
  *)

let search_1 ?n c z s1 a2 =
  let n = Option.value n ~default:(Array.length a2) in
  let rec loop i mm p =
    if i = n then
      mm, p
    else
      let al2, s2 = a2.(i) in
      (*printf "comparing %s vs %s\n%!" al1 al2;*)
      let msms = mismatches s1 s2 in
      if c msms mm then
        loop (i + 1) msms al2
      else
        loop (i + 1) mm p
  in
  loop 0 z ""

let maximum1 = search_1 (fun x y -> x > y) min_int

let max_of allele xarr =
  let _, s =
    Array.findi xarr ~f:(fun (a, _) -> a = allele)
    |> Option.value_exn ~msg:"missing allele"
    |> fun i -> xarr.(i)
  in
  let () = printf "%s: %s\n%!" allele s in
  maximum1 s xarr
  
let disp1 s (d, a) =
  printf "for %s -- %s: %d\n%!" s a d

let () =
  max_of "A*25:01:10" a2 |> disp1 "A*25:01:10";
  max_of "B*07:81" b2 |> disp1 "B*07:81";
  max_of "A*03:01:59" a2 |> disp1 "A*03:01:59";
  max_of "C*14:27" c2 |> disp1 "C*14:27";
  max_of "B*44:151" b2 |> disp1 "B*44:151";
  max_of "C*16:85" c2 |> disp1 "C*16:85";
  max_of "A*32:89" a3 |> disp1 "A*32:89";
  max_of "B*27:90:04" b3 |> disp1 "B*27:90:04";
  max_of "A*32:30:02" a3 |> disp1 "A*32:30:02";
  max_of "C*16:67" c3 |> disp1 "C*16:67";
  max_of "B*35:205" b3 |> disp1 "B*35:205";
  max_of "C*15:15" c3 |> disp1 "C*15:15"
