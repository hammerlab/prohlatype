
open Util
module N = Nomenclature

let lower_resolution = min_int

(* TODO:
   The current "compression" only takes the first value encountered.
   It might be useful to parameterize further upon the aggregation logic. *)

module MC (M : sig
  include Map.S
  val to_key : N.resolution -> key
end) = struct

  let add mp (nr, value) =
    let key = M.to_key nr in
    match M.find key mp with
    | _current_value      -> mp
    | exception Not_found -> M.add key value mp

  let compress lst =
    List.fold_left lst ~init:M.empty ~f:add
    |> M.bindings

end

module M1 = Map.Make(struct type t = int let compare = compare end)
module MOne = MC(struct
  include M1
  let to_key = function
    | N.One a | N.Two (a,_) | N.Three (a, _, _) | N.Four (a, _, _, _) -> a
end)

let compress_one lst =
  MOne.compress lst
  |> List.map ~f:(fun (a, v) -> N.One a, v)

module M2 = Map.Make(struct type t = int * int let compare = compare end)
module MTwo = MC(struct
  include M2
  let to_key = function
    | N.One a             -> a, lower_resolution
    | N.Two (a, b)
    | N.Three (a, b, _)
    | N.Four (a, b, _, _) -> a, b
  end)

let compress_two lst =
  MTwo.compress lst
  |> List.map ~f:(fun ((a, b), v) ->
      if b = lower_resolution then
        N.One a, v
      else
        N.Two (a, b), v)

module M3 = Map.Make(struct type t = int * int * int let compare = compare end)
module MThree = MC(struct
  include M3
  let to_key = function
    | N.One a             -> a, lower_resolution, lower_resolution
    | N.Two (a, b)        -> a, b,                lower_resolution
    | N.Three (a, b, c)
    | N.Four (a, b, c, _) -> a, b,                c
  end)

let compress_three lst =
  MThree.compress lst
  |> List.map ~f:(fun ((a, b, c), v) ->
      if c <> lower_resolution then
        N.Three (a, b, c), v
      else if b <> lower_resolution then
        N.Two (a, b), v
      else
        N.One a, v)

let compress ?expected_gene ~level lst =
  (* TODO: Change the expected_gene check into a result. *)
  let gene, expected_check =
    match expected_gene with
    | None    ->
        let allele_str = fst (List.hd_exn lst) in
        begin match N.parse allele_str with
        | Error e       -> invalid_argf "parsing error for %s: %s" allele_str e
        | Ok (gene, _)  -> gene
        end
        , fun _ -> ()
    | Some e  ->
        e
        , begin fun gene ->
            if gene = e then () else
              invalid_argf "Found different gene %s than expected: %s" gene e
          end
  in
  (* TODO: How to handle the suffix in this case generally?
     The current usage of Compress is for outputing and comparison against
     other results, so this might not be necessary. *)
  let as_nomenclatures =
    List.map lst ~f:(fun (allele_str, value) ->
      match N.parse allele_str with
      | Error e              -> invalid_argf "parsing error for %s: %s" allele_str e
      | Ok (gene, (nm, _so)) -> expected_check gene; (nm, value))
  in
  let compressed =
    match level with
    | `One    -> compress_one as_nomenclatures
    | `Two    -> compress_two as_nomenclatures
    | `Three  -> compress_three as_nomenclatures
  in
  List.map compressed ~f:(fun (nm, v) -> N.resolution_to_string ~gene nm, v)
