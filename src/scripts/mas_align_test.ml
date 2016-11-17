
open Util
open Common
open MoreLabels

let test f =
  let r = Mas_parser.from_file f in
  let reference = r.Mas_parser.ref_elems in
  List.map r.Mas_parser.alt_elems
    ~f:(fun (a, allele) ->
        printf "testing %s\n" a;
        a, Mas_parser.apply ~reference ~allele)

let load_fasta f =
  List.map (Fasta.all f) ~f:(fun (hdr, s) ->
    let allele = List.nth_exn (String.split ~on:(`Character ' ') hdr) 1 in
    allele, s)

module Sm = Map.Make (struct type t = string let compare = compare end)

let merge mp fa = 
  let to_map =
    List.fold_left ~init:Sm.empty ~f:(fun m (key,data) -> Sm.add ~key ~data m)
  in
  let ms = to_map mp in
  let fs = to_map fa in
  let mgd =
    Sm.merge ms fs ~f:(fun allele align_seq_o fasta_o ->
      match (align_seq_o, fasta_o) with
      | None, None      -> assert false
      | None, Some f    -> Some (`JustFasta f)
      | Some a, None    -> Some (`JustAlign a)
      | Some a, Some f  -> Some (if a = f then `Same a else `Diff (a, f)))
    |> Sm.bindings
  in
  let just_fasta, rest = List.partition_map mgd 
    ~f:(function | (a, `JustFasta s) -> `Fst (a,s) | r -> `Snd r)
  in
  let return = ref 0 in
  if just_fasta <> [] then begin
    printf "Found only in Fasta:\n";
    List.iter just_fasta ~f:(fun (allele,s) -> printf "%s: %s\n" allele s);
    return := 1
  end;
  let just_align, rest = List.partition_map rest
    ~f:(function | (a, `JustAlign s) -> `Fst (a, s) | r -> `Snd r)
  in
  if just_align <> [] then begin
    printf "Found only in Align:\n";
    List.iter just_align ~f:(fun (allele,s) -> printf "%s: %s\n" allele s);
    return := 1
  end;
  let just_diff,rest = List.partition_map rest
    ~f:(function | (a, `Diff d) -> `Fst (a, d) | r -> `Snd r)
  in
  if just_diff <> [] then begin
    printf "Different\n";
    List.iter just_diff ~f:(fun (allele, (a, f)) ->
      printf "%s:\n %s\n" allele (manual_comp_display a f));
    return := 1
  end;
  !return

let () = 
  if !Sys.interactive then
    ()
  else
    let n = Array.length Sys.argv in
    let file =
      if n <= 1 then "A_gen"
      else if n = 2 then Sys.argv.(1)
      else failwith "At most 1 arg"
    in
    exit (merge (test (to_alignment_file file)) 
                (load_fasta (to_fasta_file file)))

