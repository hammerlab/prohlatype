
let reads = reads_from_fastq "/Users/leonidrozenberg/Documents/projects/hlatyping/upenn/opti/merged/120013_TGACCA/120013_TGACCA_2fin.fastq" ;;

let greads =
  List.filter reads ~f:(fun r -> match String.index_of_character r 'N' with | Some _ -> false | _ -> true) 
  |> Array.of_list ;;
