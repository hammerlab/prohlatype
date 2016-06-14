
let cargs = { To_graph_ge.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
            ; To_graph_ge.which = None
            } ;;
let aset, g3 = To_graph_ge.(construct_from_file { cargs with which = Some (NumberOfAlts 3)}) ;;
let kt = Ref_graph.kmer_list_ge ~k:5 g3 ;;
let kmt = kt.Ref_graph.full ;;

