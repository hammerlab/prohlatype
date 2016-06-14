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
  
let parse_numbers s =
  trim_suffix s >>=
    fun (without_suffix, suffix_opt) ->
      String.split ~on:(`Character ':') without_suffix
      |> parse_ints  >>= function
          | []            -> Error (sprintf "Empty allele name: %s" without_suffix)
          | [ a ]         -> Ok (a, 0, 0, 0)
          | [ a; b ]      -> Ok (a, b, 0, 0)
          | [ a; b; c]    -> Ok (a, b, c, 0)
          | [ a; b; c; d] -> Ok (a, b, c, d)
          | lst           -> Error (sprintf "parsed more than 4 ints: %s" without_suffix)
     
