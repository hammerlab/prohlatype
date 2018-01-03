

module SeqSet = Set.Make (
  struct
    type t = string MSA.alignment_element
    let compare = compare
  end)

let just_by_label s =
  List.map s ~f:(fun (bm, l) -> bm.MSA.Boundaries.label, l)

let is_full grpd =
  let open MSA in
  let _fbm, fbl = List.hd_exn grpd in
  match List.last grpd with
  | None -> false
  | Some (_, lbl) ->
      List.exists fbl ~f:is_start &&
      List.exists lbl ~f:is_end

let of_mp ?(skip_incomplete=true) ?(skip_gaps=true) mp =
  let open MSA in
  let ref_grps = Boundaries.grouped mp.Parser.ref_elems in
  let inits =
    List.map ref_grps ~f:(fun (b, l) -> b, (b, SeqSet.of_list l, SeqSet.empty))
    |> just_by_label
  in
  List.fold_left mp.Parser.alt_elems ~init:(inits, [])
    ~f:(fun (acc_lst, alst) {Parser.seq; allele } ->
        let seq_grps = just_by_label (Boundaries.grouped seq) in
        if skip_incomplete && not (is_full seq_grps) then begin
          (*printf "skipping not full %s\n" allele; *)
          acc_lst, alst
        end else
          merge_assoc acc_lst seq_grps ~f:(fun b ms ->
            match ms with
            | `First s                        ->
                Some s
            | `Second _                       ->
                invalid_argf "Reference doesn't have %s boundary!"
                  (boundary_label_to_string b)
            | `Both ((rb, ref_seq, acc_set), lst) ->
                let nacc_set =
                  List.fold_left lst ~init:acc_set ~f:(fun acc st ->
                    if SeqSet.mem st ref_seq then
                      acc
                    else if skip_gaps && is_gap st then
                      acc
                    else
                      SeqSet.add st acc)
                in
                Some (rb, ref_seq, nacc_set))
          , allele :: alst)

let describe l =
  let open MSA in
  List.map l ~f:(fun (_l, (b, _rs, acc)) ->
    let num_variations = SeqSet.cardinal acc in
    b
    , num_variations
    , (float num_variations) /.  (float b.Boundaries.seq_length))

let full ?skip_incomplete ?skip_gaps f =
  let s, alst =
    Common.to_alignment_file f
    |> MSA.Parser.from_file
    |> of_mp ?skip_incomplete ?skip_gaps
  in
  describe s, (List.length alst)

let pretty_full ?skip_incomplete ?skip_gaps f =
  let bl, num_alleles = full ?skip_incomplete ?skip_gaps f in
  printf "For %s, number of alleles %d\n" f num_alleles;
  let open MSA in
  List.iter bl ~f:(fun (bm, _tl, p) ->
    printf "%10s\t%.2f\n"
      (boundary_label_to_string bm.Boundaries.label)
      (100. *. p))



