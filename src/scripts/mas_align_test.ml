(* Test that the alignment parsing and then applying gives the same sequences
   as the fasta's

   There are known (2016-11-18) differences:

    B_nuc:    B*44:02:01:02S
    C_nuc:    C*04:09N
    DRB_nuc:  DRB4*01:03:01:02N

*)
open Util
open Common

let apply_allele a r =
  let open Mas_parser in
  allele_sequence ~reference:r.ref_elems ~allele:(List.assoc a r.alt_elems) ()

let test f =
  let r = Mas_parser.from_file f in
  let reference = r.Mas_parser.ref_elems in
  (r.Mas_parser.reference, Mas_parser.reference_sequence r)
  :: List.map r.Mas_parser.alt_elems
      ~f:(fun (a, allele) ->
          printf "testing %s\n" a;
          a, Mas_parser.allele_sequence ~reference ~allele ())

let load_fasta f =
  List.map (Fasta.all f) ~f:(fun (hdr, s) ->
    let allele = List.nth_exn (String.split ~on:(`Character ' ') hdr) 1 in
    allele, s)

let merge ?(known=[||]) mp fa =
  let to_map =
    List.fold_left ~init:StringMap.empty
      ~f:(fun m (key,data) -> StringMap.add ~key ~data m)
  in
  let ms = to_map mp in
  let fs = to_map fa in
  let mgd =
    StringMap.merge ms fs ~f:(fun allele align_seq_o fasta_o ->
      match (align_seq_o, fasta_o) with
      | None, None      -> assert false
      | None, Some f    -> Some (`JustFasta f)
      | Some a, None    -> Some (`JustAlign a)
      | Some a, Some f  -> Some (if a = f then `Same a else `Diff (a, f)))
    |> StringMap.bindings
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
    let all_known =
      List.fold_left just_diff ~init:true ~f:(fun all_ok (allele, (a, f)) ->
        let k = Array.exists ~f:((=) allele) known in
        printf "known: %b %s:\n %s\n" k allele
          (manual_comp_display ~labels:("align: ","fasta: ") a f);
        k)
    in
    if all_known then return := 0 else return := 1
  end;
  !return

let () =
  if !Sys.interactive then
    ()
  else
    let n = Array.length Sys.argv in
    let file, known =
      if n <= 1 then "A_gen", [||]
      else if n = 2 then Sys.argv.(1), [||]
      else Sys.argv.(1), Array.sub Sys.argv ~pos:2 ~len:(n - 2)
    in
    let r =
      merge ~known
        (test (to_alignment_file file))
        (load_fasta (to_fasta_file file))
    in
    if r = 0 then
      printf "all tests passed!\n"
    else
      exit r
