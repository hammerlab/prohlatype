(** Decode HLA allele names. *)

open Util

type suffix =
  | N
  | L
  | S
  | C
  | A
  | Q

let suffix_opt_of_char = function
  | 'N' -> Some N
  | 'L' -> Some L
  | 'S' -> Some S
  | 'C' -> Some C
  | 'A' -> Some A
  | 'Q' -> Some Q
  | x   -> None

let suffix_to_string = function
  | N -> "N"
  | L -> "L"
  | S -> "S"
  | C -> "C"
  | A -> "A"
  | Q -> "Q"

let int_of_string ?msg s =
  let msg = Option.value msg ~default:(sprintf "int_of_string parsing: %s" s) in
  match Int.of_string s with
  | None   -> Error msg
  | Some i -> Ok i

let trim_suffix s =
  let index = String.length s - 1 in
  let c = String.get_exn s ~index in
  match suffix_opt_of_char c with
  | Some su -> Ok (String.take s ~index, Some su)
  | None    ->
    let msg = sprintf "Allele %s doesn't end on a known suffix or integer: %c" s c in
    match int_of_string ~msg (String.drop s ~index) with
    | Error m -> Error m
    | Ok _    -> Ok (s, None)

let parse_ints lst =
  let rec loop acc = function
    | []     -> Ok (List.rev acc)
    | h :: t ->
        begin
          match int_of_string h with
          | Ok x    -> loop (x :: acc) t
          | Error e -> Error e
        end
  in
  loop [] lst

type resolution =
  | One of int
  | Two of int * int
  | Three of int * int * int
  | Four of int * int * int * int

(* This logic compares alleles such that the sorted order is the same as in the
   alignment files. Without doing too much work the resulting order has a
   tremendous impact on performance: For example, ParPHMM forward passes are at
   least twice as fast because of the resulting partition_maps are more
   compressed! *)
let compare_by_resolution (r1,_) (r2,_) =
  let to_quad = function
    | One a           -> [ a; -1; -1; -1 ]
    | Two (a,b)       -> [ a;  b; -1; -1 ]
    | Three (a,b,c)   -> [ a;  b;  c; -1 ]
    | Four (a,b,c,d)  -> [ a;  b;  c;  d ]
  in
  Pervasives.compare (to_quad r1) (to_quad r2)

let parse_resolution s =
  trim_suffix s >>=
    fun (without_suffix, suffix_opt) ->
      String.split ~on:(`Character ':') without_suffix
      |> parse_ints  >>= function
          | []            -> Error (sprintf "Empty allele name: %s" without_suffix)
          | [ a ]         -> Ok (One a, suffix_opt)
          | [ a; b ]      -> Ok (Two (a, b), suffix_opt)
          | [ a; b; c]    -> Ok (Three (a, b, c), suffix_opt)
          | [ a; b; c; d] -> Ok (Four (a, b, c, d), suffix_opt)
          | lst           -> Error (sprintf "parsed more than 4 ints: %s" without_suffix)

(* This list was manually generated on 2017-10-12 and is unfortunately
   hard-coded here. Much of the current work has been focused on making the
   inference feasible for A,B,C (especially the imputation) that I would
   strongly recommend that any extension require a thoughtful experimentation
   with results on other genes. *)
type locus =
  | A
  | B
  | C
  | DMA
  | DMB
  | DOA
  | DOB
  | DPA1
  | DPA2
  | DPB1
  | DPB2
  | DQA1
  | DQB1
  | DRA
  | DRB1
  | DRB3
  | DRB4
  | E
  | F
  | G
  | HFE
  | H
  | J
  | K
  | L
  | MICA
  | MICB
  | P
  | TAP1
  | TAP2
  | T
  | V
  | W
  | Y
  [@@deriving show]

let parse_locus = function
  | "A"     -> Ok A
  | "B"     -> Ok B
  | "C"     -> Ok C
  | "DMA"   -> Ok DMA
  | "DMB"   -> Ok DMB
  | "DOA"   -> Ok DOA
  | "DOB"   -> Ok DOB
  | "DPA1"  -> Ok DPA1
  | "DPA2"  -> Ok DPA2
  | "DPB1"  -> Ok DPB1
  | "DPB2"  -> Ok DPB2
  | "DQA1"  -> Ok DQA1
  | "DQB1"  -> Ok DQB1
  | "DRA"   -> Ok DRA
  | "DRB1"  -> Ok DRB1
  | "DRB3"  -> Ok DRB3
  | "DRB4"  -> Ok DRB4
  | "E"     -> Ok E
  | "F"     -> Ok F
  | "G"     -> Ok G
  | "HFE"   -> Ok HFE
  | "H"     -> Ok H
  | "J"     -> Ok J
  | "K"     -> Ok K
  | "L"     -> Ok L
  | "MICA"  -> Ok MICA
  | "MICB"  -> Ok MICB
  | "P"     -> Ok P
  | "TAP1"  -> Ok TAP1
  | "TAP2"  -> Ok TAP2
  | "T"     -> Ok T
  | "V"     -> Ok V
  | "W"     -> Ok W
  | "Y"     -> Ok Y
  | l       -> error "Unrecognized locus: %s, please send a pull request!" l

type t = resolution * suffix option

let parse : string -> (locus * t, string) result =
  fun s ->
    match String.split s ~on:(`Character '*') with
    | [ l; n] -> parse_locus l >>= fun l ->
                    parse_resolution n >>= fun p -> Ok (l, p)
    | _ :: [] -> error "Did not find the '*' separator in %s" s
    | _       -> error "Found too many '*' separators in %s" s

let parse_to_resolution_exn s =
  parse s |> unwrap_ok |> snd

let resolution_to_string ?locus =
  let ls = Option.value_map ~default:"" ~f:(fun l -> show_locus l ^ "*") locus in
  function
  | One a             -> sprintf "%s%02d" ls a
  | Two (a, b)        -> sprintf "%s%02d:%02d" ls a b
  | Three (a, b, c)   -> sprintf "%s%02d:%02d:%02d" ls a b c
  | Four (a, b, c, d) -> sprintf "%s%02d:%02d:%02d:%02d" ls a b c d

let resolution_and_suffix_opt_to_string ?locus (r,so) =
  sprintf "%s%s" (resolution_to_string ?locus r)
    (Option.value_map so ~default:"" ~f:suffix_to_string)

let two_matches_full nr1 nr2 =
  match (fst nr1) with                     (* Match on just the name, ignore suffix *)
  | Two (r1, r2) ->
      begin
        match (fst nr2) with                               (* Also, ignore suffix. *)
        | One _               -> false
        | Two (s1, s2)
        | Three (s1, s2, _)
        | Four (s1, s2, _, _) -> r1 = s1 && r2 = s2
      end
  | fnr ->
      invalid_argf "Not Two resolution: %s " (resolution_to_string fnr)

module Trie (*: sig

  type t
  val empty : t
  val add : resolution -> t -> t
  val nearest : resolution -> t -> resolution

end *) = struct

  let insert_sorted ~cmp ~eq e lst =
    let rec loop acc = function
      | []            -> List.rev (e :: acc)
      | (h :: t) as l -> let r = cmp h in
                        if r < 0 then
                          loop (h :: acc) t
                        else if r = 0 then
                          (List.rev ((eq h) :: acc)) @ t
                        else
                          (List.rev (e :: acc)) @ l
    in
    loop [] lst

  let find_closest ~distance = function
    | []      -> invalid_argf "Empty!"
    | h :: t  ->
        let rec loop a = function
          | []     -> a
          | h :: t -> if distance a < distance h then a else loop h t
        in
        loop h t

  (*
  find_closest ~distance:(fun x -> abs (x - 3)) [1;2;3;4] = 3
  find_closest (fun x -> abs (x - 3)) [1;2;3;4;5] = 3
  find_closest (fun x -> abs (x - 3)) [1;2;4;5] = 4
  find_closest (fun x -> abs (x - 3)) [1;4;5] = 4
  find_closest (fun x -> abs (x - 3)) [1;2;5] = 2
  *)

  let empty = []

  type node =
    | Leaf of suffix option
    | Level of (int * node) list

  let cmp_fst_against a = fun (f,_) -> f - a
  let insert_by_fst a = insert_sorted ~cmp:(cmp_fst_against a)

  let add n so l =
    let insert_final_level s v so l =
      let eq _ = invalid_argf "Duplicate %s digit: %d" s v in
      (insert_by_fst v) ~eq (v, Leaf so) l
    in
    let insert_above_final s v w so l =
      let eq = function
        | (v, Leaf so)  -> invalid_argf "Nil already specified: %d" v
        | (v, Level tl) -> v, Level (insert_final_level s w so tl)
      in
      (insert_by_fst v) ~eq (eq (v, Level [])) l
    in
    let insert_two_above_final s v w x so l =
      let eq = function
        | (v, Leaf so)  -> invalid_argf "Nil already specified: %d" v
        | (v, Level tl) -> v, Level (insert_above_final s w x so tl)
      in
      (insert_by_fst v) ~eq (eq (v, Level [])) l
    in
    let insert_three_above_final s v w x y so l =
      let eq = function
        | (v, Leaf so)  -> invalid_argf "Nil already specified: %d" v
        | (v, Level tl) -> v, Level (insert_two_above_final s w x y so tl)
      in
      (insert_by_fst v) ~eq (eq (v, Level [])) l
    in
    match n with
    | One a             -> insert_final_level "2nd" a so l
    | Two (a, b)        -> insert_above_final "4th" a b so l
    | Three (a, b, c)   -> insert_two_above_final "6th" a b c so l
    | Four (a, b, c, d) -> insert_three_above_final "8th"  a b c d so l

  let nearest ?(distance=fun x y -> abs (x - y)) n l =
    let dst x (y, _) = distance x y in
    let find_or_first = function
      | None   -> List.hd_exn
      | Some v -> find_closest ~distance:(dst v)
    in
    let lookup a b c d =
      match find_or_first a l with
      | (a, Leaf so)          -> One a, so
      | (a, Level l2)         ->
        match find_or_first b l2 with
        | (b, Leaf so)        -> Two (a, b), so
        | (b, Level l3)       ->
          match find_or_first c l3 with
          | (c, Leaf so)      -> Three (a, b, c), so
          | (c, Level l4)     ->
            match find_or_first d l4 with
            | (d, Leaf so)    -> Four (a, b, c, d), so
            | (d, Level _)    -> invalid_argf "Found more than 4 levels for %d:%d:%d:%d"
                                    a b c d
    in
    match n with
    | One a             -> lookup (Some a)  None     None     None
    | Two (a, b)        -> lookup (Some a) (Some b)  None     None
    | Three (a, b, c)   -> lookup (Some a) (Some b) (Some c)  None
    | Four (a, b, c, d) -> lookup (Some a) (Some b) (Some c) (Some d)

end
