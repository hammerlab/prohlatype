
let cargs = { To_graph.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
            ; To_graph.which = None
            } ;;
let aset, g3 = To_graph.(construct_from_file { cargs with which = Some (NumberOfAlts 3)}) ;;
let kt = Ref_graph.kmer_list ~k:5 g3 ;;
let kmt = kt.Ref_graph.full ;;

