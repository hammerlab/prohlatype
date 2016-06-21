
open Prohlatype
let cargs = { Cache.alignment_file = "../foreign/IMGTHLA/alignments/A_nuc.txt"
            ; Cache.which = None
            } ;;
let k = 5 ;;
let g, aindex, idx = Cache.graph_and_two_index { Cache.k = k; Cache.g = cargs };;

