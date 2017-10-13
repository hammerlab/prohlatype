(*include MoreLabels 
module String = Sosa.Native_bytes *)
(*
module NList = struct
  include Nonstd.List

  let map_snd lst ~f = map lst ~f:(fun (k, v) -> k, f v)

  let map2_snd l1 l2 ~f =
    map2 l1 l2 ~f:(fun (k1, v1) (k2, v2) ->
      assert (k1 = k2);
      (k1, f v1 v2))

end (* NList *)
*)
(*
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
*)

(*
let invalid_argf ?(prefix="") fmt =
  Printf.ksprintf invalid_arg ("%s" ^^ fmt) prefix

let failwithf fmt =
  Printf.ksprintf failwith fmt
  *)
(*let id x = x *)

(*
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

let error fmt = Printf.ksprintf (fun s -> Error s) fmt

let unwrap_ok = function
  | Ok o    -> o
  | Error _ -> invalid_arg "Not Ok in unwrap_ok."

let unwrap_error = function
  | Ok _    -> invalid_argf "Not Error in unwrap_error."
  | Error e -> e
  *)

open Core_kernel.Std

let short_seq s =
  let n = String.length s in
  if n > 10 then
    Printf.sprintf "%s...%s" (String.prefix s 4) (String.suffix s 3)
  else
    s

let string_split_at s ~index =
  let sindex = String.length s - index in
  if sindex < 0 then
    String.prefix s index, ""
  else
    String.prefix s index, String.suffix s sindex

let index_string s index =
  let (b, a) = string_split_at s ~index in
  b ^ "." ^ a

let _pair_of_empty_strings = "", ""

(** Compare two strings and display vertical bars for mismatches. *)
let manual_comp_display ?(msm_offset=0) ?width ?(labels=_pair_of_empty_strings) s1 s2 =
  let msm = ref msm_offset in
  let mismatch_string =
    String.mapi s1 ~f:(fun index c1 ->
      try
        let c2 = String.get s2 index in
        if Char.equal c2 c1 then
          ' '
        else begin
          msm := !msm + 1;
          '|'
        end
      with (Invalid_argument _) ->
        msm := !msm + 1;
        'X')
      (*match String.get s2 index with
      | None                 -> incr msm; 'X'
      | Some c2 when c1 = c2 -> ' '
      | Some _ (*c1 <> c2*)  -> incr msm; '|' *)
  in
  let ms = Int.to_string !msm in
  let n  = String.length ms in
  let top, bottom = labels in
  let label_length = max (max (String.length top) (String.length bottom)) n in
  let top_row = Printf.sprintf "%-*s%s" label_length top s1 in
  let mid_row = Printf.sprintf "%*s%s" label_length ms mismatch_string in
  let bot_row = Printf.sprintf "%-*s%s" label_length bottom s2 in
  let index = Option.value width ~default:1000 in
  let rec loop acc t m b =
    let t, tl = string_split_at t ~index in
    let m, ml = string_split_at m ~index in
    let b, bl = string_split_at b ~index in
    let nacc = (Printf.sprintf "%s\n%s\n%s" t m b) :: acc in
    if String.is_empty tl && String.is_empty ml && String.is_empty bl then
      String.concat ~sep:"\n" (List.rev nacc)
    else
      loop nacc tl ml bl
  in
  loop [] top_row mid_row bot_row

let mismatch_indices s1 s2 =
  String.foldi s1 ~init:(0, []) ~f:(fun index (i, a) c1 ->
    let c2 = String.get s2 index in
    if Char.equal c1 c2 then (i + 1, a) else (i + 1, i :: a))
  |> snd
  |> List.rev

let insert_chars ?(every=120) ?(token=';') ics s =
  String.to_list s
  |> List.fold_left ~init:(0,[]) ~f:(fun (i, acc) c ->
      if i > every && Char.equal c token then
        (0, ics @ (c :: acc))
      else
        (i + 1, c :: acc))
  |> snd
  |> List.rev
  |> String.of_char_list

let complement = function
  | 'A' -> 'T'
  | 'C' -> 'G'
  | 'G' -> 'C'
  | 'T' -> 'A'
  |  c  ->  c   (* could complain, but there are 'N' 's *)

let reverse_complement s =
  String.fold s ~init:[] ~f:(fun l c -> complement c :: l)
  |> String.of_char_list

(*
let list_fold_ok lst ~f ~init =
  let rec loop acc = function
    | []      -> Ok acc
    | h :: t  -> f acc h >>= fun a -> loop a t
  in 
  loop init lst *)

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

(*
module StringSet = struct
  include Set.M (String)

module StringMap = Map.M (String)
*)

let string_set_of_list lst =
  List.fold_left lst ~init:(Set.empty String.comparator)
    ~f:(fun s e -> Set.add s e)

let string_map_of_assoc asc =
  List.fold_left asc ~init:(Map.empty String.comparator) 
    ~f:(fun acc (key, data) -> Map.add ~key ~data acc)

let remove_and_assoc el list =
  let rec loop acc = function
    | []                      -> raise Not_found
    | (e, v) :: t when e = el -> v, (List.rev acc @ t)
    | h :: t                  -> loop (h :: acc) t
  in
  loop [] list

(*
let assoc v l =
  Option.value_exn ~msg:"Not found" (List.Assoc.find v l)
  *)
let list_assoc_remove_and_get el list =
  let rec loop acc = function
    | []                      -> None
    | (e, v) :: t when e = el -> Some (v, (List.rev acc @ t))
    | h :: t                  -> loop (h :: acc) t
  in
  loop [] list

let group_by_assoc l =
  let insert assoc (k, v) =
    match list_assoc_remove_and_get k assoc with
    | None              -> (k,[v]) :: assoc
    | Some (cv, rassoc) -> (k, v ::cv) :: rassoc
  in
  List.fold_left ~init:[] ~f:insert l

let string_of_list ?(show_empty=false) ~sep ~f = function
  | [] -> if show_empty then "[]" else ""
  | l  -> String.concat ~sep (List.map ~f l)

let log_likelihood ?(alphabet_size=4) ?(er=0.01) ~len mismatches =
  let bs = Float.of_int (alphabet_size - 1) in
  let lmp = Float.(log (er / bs)) in
  let lcp = Float.(log (1. - er)) in
  let c = Float.(of_int len - mismatches) in
  Float.(c * lcp + mismatches * lmp)

let likelihood ?alphabet_size ?er ~len m =
  Float.exp (log_likelihood ?alphabet_size ?er ~len m)

type too_short =
  | TooShort of { desired: int ; actual: int}
  [@@deriving show]

let manual_phred_llhd_lst s1 s2 probability_of_error =
  String.foldi s1 ~init:(0, []) ~f:(fun index (i, a) c1 ->
    let c2 = String.get s2 index in
    if Char.equal c1 c2 then
      (i + 1, (`m Float.(log1p (-probability_of_error.(i)))) :: a)
    else
      (i + 1, (`X Float.(log (probability_of_error.(i) / 3.0)) :: a)))
  |> snd
  |> List.rev

let manual_phred_llhd s1 s2 probability_of_error =
  manual_phred_llhd_lst s1 s2 probability_of_error
  |> List.fold_left ~init:0.
      ~f:(fun s -> function | `m p | `X p -> Float.(s + p))

let time s f =
  let n = Sys.time () in
  try
    let r = f () in
    printf "%s, total running time in seconds: %f\n%!" s (Sys.time () -. n);
    r
  with e ->
    printf "%s, failed in seconds: %f\n%!" s (Sys.time () -. n);
    raise e

let gc_between s f =
  let open Gc.Stat in
  let before = Gc.stat () in
  let r = f () in
  let after = Gc.stat () in
  printf "%s Gc change: \n\
    \t { minor_words : %f;\n\
    \t   promoted_words : %f;\n\
    \t   major_words : %f;\n\
    \t   minor_collections : %d;\n\
    \t   major_collections : %d;\n\
    \t   heap_words : %d;\n\
    \t   heap_chunks : %d;\n\
    \t   live_words : %d;\n\
    \t   live_blocks : %d;\n\
    \t   free_words : %d;\n\
    \t   free_blocks : %d;\n\
    \t   largest_free : %d;\n\
    \t   fragments : %d;\n\
    \t   compactions : %d;\n\
    \t   top_heap_words : %d;\n\
    \t   stack_size : %d;\n\
    \t }\n"
      s
      (after.minor_words -. before.minor_words)
      (after.promoted_words -. before.promoted_words)
      (after.major_words -. before.major_words)
      (after.minor_collections - before.minor_collections)
      (after.major_collections - before.major_collections)
      (after.heap_words - before.heap_words)
      (after.heap_chunks - before.heap_chunks)
      (after.live_words - before.live_words)
      (after.live_blocks - before.live_blocks)
      (after.free_words - before.free_words)
      (after.free_blocks - before.free_blocks)
      (after.largest_free - before.largest_free)
      (after.fragments - before.fragments)
      (after.compactions - before.compactions)
      (after.top_heap_words - before.top_heap_words)
      (after.stack_size - before.stack_size);
  r

type 'a single_or_paired =
  | Single of 'a
  | Paired of ('a * 'a)

(* Maintain a sorted association list of the top n items. *)
let topn p k a i lst =
  let rec loop added n lst =
    if n >= k then
      []
    else
      match lst with
      | []         -> if added then [] else [a,i]
      | (u,j) :: t -> if p a u && not added then
                        (a,i) :: loop true (n + 1) lst
                      else
                        (u,j) :: loop added (n + 1)  t
  in
  loop false 0 lst

let insert_sorted p a i l =
  let rec loop lst = match lst with
    | []         -> [a, i]
    | (u,j) :: t -> if p a u then
                      (a, i) :: lst
                    else
                      (u, j) :: loop t
  in
  loop l
