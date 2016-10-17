(** Encode HLA allele's *)

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

let parse_resolution s =
  trim_suffix s >>=
    fun (without_suffix, suffix_opt) ->
      String.split ~on:(`Character ':') without_suffix
      |> parse_ints  >>= function
          | []            -> Error (sprintf "Empty allele name: %s" without_suffix)
          | [ a ]         -> Ok (One a)
          | [ a; b ]      -> Ok (Two (a, b))
          | [ a; b; c]    -> Ok (Three (a, b, c))
          | [ a; b; c; d] -> Ok (Four (a, b, c, d))
          | lst           -> Error (sprintf "parsed more than 4 ints: %s" without_suffix)

let parse s =
  match String.split s ~on:(`Character '*') with
  | [ g; n] -> parse_resolution n >>= fun p -> Ok (g, p)
  | _ :: [] -> error "Did not find the '*' separator in %s" s
  | _       -> error "Found too many '*' separators in %s" s

let resolution_to_string ?gene =
  let gene_str = match gene with | None -> "" | Some g -> g ^ "*" in
  function
  | One a             -> sprintf "%s%02d" gene_str a
  | Two (a, b)        -> sprintf "%s%02d:%02d" gene_str a b
  | Three (a, b, c)   -> sprintf "%s%02d:%02d:%02d" gene_str a b c
  | Four (a, b, c, d) -> sprintf "%s%02d:%02d:%02d:%02d" gene_str a b c d
