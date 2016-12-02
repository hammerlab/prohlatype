
include Nonstd
module String = Sosa.Native_string

let invalid_argf ?(prefix="") fmt =
  ksprintf invalid_arg ("%s" ^^ fmt) prefix

let id x = x

let error_bind ~f oe =
  match oe with
  | Error e -> Error e
  | Ok o -> f o

let (>>=) oe f = error_bind ~f oe

let error_map ~f oe =
  match oe with
  | Error e -> Error e
  | Ok o    -> Ok (f o)

let (>>|) oe ~f = error_map ~f oe

let error fmt = ksprintf (fun s -> Error s) fmt

let unwrap_ok = function
  | Ok o    -> o
  | Error _ -> invalid_arg "Not Ok in unwrap_ok."

let unwrap_error = function
  | Ok _    -> invalid_argf "Not Error in unwrap_error."
  | Error e -> e

let short_seq s =
  let n = String.length s in
  if n > 10 then
    sprintf "%s...%s"
      (String.slice_exn ~finish:4 s) (String.slice_exn ~start:(n-3) s)
  else
    s

let index_string s index =
  let (b, a) = String.split_at s ~index in
  b ^ "." ^ a

let _pair_of_empty_strings = String.empty, String.empty
(** Compare two strings and display vertical bars for mismatches. *)
let manual_comp_display ?width ?(labels=_pair_of_empty_strings) s1 s2 =
  let msm = ref 0 in
  let mismatch_string =
    String.mapi s1 ~f:(fun index c1 ->
      match String.get s2 ~index with
      | None                 -> incr msm; 'X'
      | Some c2 when c1 = c2 -> ' '
      | Some _ (*c1 <> c2*)  -> incr msm; '|')
  in
  let ms = string_of_int !msm in
  let n  = String.length ms in
  let top, bottom = labels in
  let label_length = max (max (String.length top) (String.length bottom)) n in
  let top_row = sprintf "%-*s%s" label_length top s1 in
  let mid_row = sprintf "%*s%s" label_length ms mismatch_string in
  let bot_row = sprintf "%-*s%s" label_length bottom s2 in
  let index = Option.value width ~default:max_int in
  let rec loop acc t m b =
    let t, tl = String.split_at t ~index in
    let m, ml = String.split_at m ~index in
    let b, bl = String.split_at b ~index in
    let nacc = (sprintf "%s\n%s\n%s" t m b) :: acc in
    if String.is_empty tl && String.is_empty ml && String.is_empty bl then
      String.concat ~sep:"\n" (List.rev nacc)
    else
      loop nacc tl ml bl
  in
  loop [] top_row mid_row bot_row

let insert_chars ?(every=120) ?(token=';') ics s =
  String.to_character_list s
  |> List.fold_left ~init:(0,[]) ~f:(fun (i, acc) c ->
      if i > every && c = token then
        (0, ics @ (c :: acc))
      else
        (i + 1, c :: acc))
  |> snd
  |> List.rev
  |> String.of_character_list

let complement = function
  | 'A' -> 'T'
  | 'C' -> 'G'
  | 'G' -> 'C'
  | 'T' -> 'A'
  |  c  ->  c   (* could complain, but there are 'N' 's *)

let reverse_complement s =
  String.fold s ~init:[] ~f:(fun l c -> complement c :: l)
  |> String.of_character_list

let list_fold_ok lst ~f ~init =
  let rec loop acc = function
    | []      -> Ok acc
    | h :: t  -> f acc h >>= fun a -> loop a t
  in
  loop init lst

let opt_is_true = function
  | Some true  -> true
  | Some false -> false
  | None       -> false

let list_map_consecutives f lst =
  let rec loop acc = function
    | []
    | _ :: []     -> List.rev acc
    | a :: b :: t -> loop (f a b :: acc) (b :: t)
  in
  loop [] lst
