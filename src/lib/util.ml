include MoreLabels
module String = Sosa.Native_bytes

module NList = struct
  include Nonstd.List

  let map_snd lst ~f = map lst ~f:(fun (k, v) -> k, f v)

  let map2_snd l1 l2 ~f =
    map2 l1 l2 ~f:(fun (k1, v1) (k2, v2) ->
      assert (k1 = k2);
      (k1, f v1 v2))

end (* NList *)

module NArray = struct
  include Nonstd.Array

  let rev a =
    let n = Array.length a in
    Array.init n (fun i -> Array.unsafe_get a (n - i - 1))

  let findi a ~f =
    let n = Array.length a in
    let rec loop i =
      if i >= n then None
      else if f a.(i) then Some i
      else loop (i + 1)
    in
    loop 0

end (* NArray *)
include (Nonstd : module type of Nonstd with module List := NList
                                         and module Array := NArray)

module List = NList
module Array = NArray

let invalid_argf ?(prefix="") fmt =
  ksprintf invalid_arg ("%s" ^^ fmt) prefix

let failwithf fmt =
  ksprintf failwith fmt

let id x = x

let result_bind ~f oe =
  match oe with
  | Error e -> Error e
  | Ok o -> f o

let (>>=) oe f = result_bind ~f oe

let result_map ~f oe =
  match oe with
  | Error e -> Error e
  | Ok o    -> Ok (f o)

let (>>|) oe ~f = result_map ~f oe

let error_map ~f oe =
  match oe with
  | Error e -> Error (f e)
  | Ok o    -> Ok o

let error_bind ~f oe =
  match oe with
  | Error e -> f e
  | Ok o    -> Ok o

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
let manual_comp_display ?(msm_offset=0) ?width ?(labels=_pair_of_empty_strings) s1 s2 =
  let msm = ref msm_offset in
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

let mismatch_indices s1 s2 =
  String.fold2_exn s1 s2 ~init:(0, []) ~f:(fun (i, a) c1 c2 ->
    if c1 = c2 then (i + 1, a) else (i + 1, i :: a))
  |> snd
  |> List.rev

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

module StringSet = Set.Make (struct
  type t = string [@@deriving ord]
end)

module StringMap = Map.Make (struct
  type t = string [@@deriving ord]
end)

let string_set_of_list lst =
  List.fold_left lst ~init:StringSet.empty
    ~f:(fun s e -> StringSet.add e s)

let string_map_of_assoc asc =
  List.fold_left asc ~init:StringMap.empty
    ~f:(fun acc (key, data) -> StringMap.add ~key ~data acc)

let remove_and_assoc el list =
  let rec loop acc = function
    | []                      -> raise Not_found
    | (e, v) :: t when e = el -> v, (List.rev acc @ t)
    | h :: t                  -> loop (h :: acc) t
  in
  loop [] list

let group_by_assoc l =
  let insert assoc (k, v) =
    match List.Assoc.remove_and_get k assoc with
    | None              -> (k,[v]) :: assoc
    | Some (cv, rassoc) -> (k, v ::cv) :: rassoc
  in
  List.fold_left ~init:[] ~f:insert l

let log_likelihood ?(alph_size=4) ?(er=0.01) ~len mismatches =
  let lmp = log (er /. (float (alph_size - 1))) in
  let lcp = log (1. -. er) in
  let c = (float len) -. mismatches in
  c *. lcp +. mismatches *. lmp

let likelihood ?alph_size ?er ~len m =
  exp (log_likelihood ?alph_size ?er ~len m)

type too_short =
  | TooShort of { desired: int ; actual: int}
  [@@deriving show]

let manual_phred_llhd_lst s1 s2 probability_of_error =
  String.fold2_exn s1 s2 ~init:(0, []) ~f:(fun (i, a) c1 c2 ->
    if c1 = c2 then
      (i + 1, (`m (log1p (-. probability_of_error.(i)))) :: a)
    else
      (i + 1, (`X (log (probability_of_error.(i) /. 3.)) :: a)))
  |> snd
  |> List.rev

let manual_phred_llhd s1 s2 probability_of_error =
  manual_phred_llhd_lst s1 s2 probability_of_error
  |> List.fold_left ~init:0. ~f:(fun s -> function | `m p | `X p -> s +. p)

let time s f =
  let n = Sys.time () in
  let r = f () in
  printf "%s total running time in seconds: %f\n%!" s (Sys.time () -. n);
  r

type 'a single_or_paired =
  | Single of 'a
  | Paired of ('a * 'a)


