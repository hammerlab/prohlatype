
let reads_from_fastq file =
  let ic = open_in file in
  let li = ref [] in
  try
    let rec loop i =
      let line = input_line ic in
      if i mod 4 = 1 then li := line :: !li;
      loop (i + 1)
    in
    loop 0
  with End_of_file ->
    close_in ic;
    !li


