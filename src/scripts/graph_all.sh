#!/usr/bin/env bash

set -e    # stop if failure

FILES="../foreign/IMGTHLA/alignments/A_gen.txt
../foreign/IMGTHLA/alignments/A_nuc.txt
../foreign/IMGTHLA/alignments/B_gen.txt
../foreign/IMGTHLA/alignments/B_nuc.txt
../foreign/IMGTHLA/alignments/C_gen.txt
../foreign/IMGTHLA/alignments/C_nuc.txt
../foreign/IMGTHLA/alignments/DMA_gen.txt
../foreign/IMGTHLA/alignments/DMA_nuc.txt
../foreign/IMGTHLA/alignments/DMB_gen.txt
../foreign/IMGTHLA/alignments/DMB_nuc.txt
../foreign/IMGTHLA/alignments/DOA_gen.txt
../foreign/IMGTHLA/alignments/DOA_nuc.txt
../foreign/IMGTHLA/alignments/DOB_gen.txt
../foreign/IMGTHLA/alignments/DOB_nuc.txt
../foreign/IMGTHLA/alignments/DPA1_gen.txt
../foreign/IMGTHLA/alignments/DPA1_nuc.txt
../foreign/IMGTHLA/alignments/DPB1_gen.txt
../foreign/IMGTHLA/alignments/DPB1_nuc.txt
../foreign/IMGTHLA/alignments/DPB2_gen.txt
../foreign/IMGTHLA/alignments/DPB2_nuc.txt
../foreign/IMGTHLA/alignments/DQA1_gen.txt
../foreign/IMGTHLA/alignments/DQA1_nuc.txt
../foreign/IMGTHLA/alignments/DQB1_gen.txt
../foreign/IMGTHLA/alignments/DQB1_nuc.txt
../foreign/IMGTHLA/alignments/DRA_gen.txt
../foreign/IMGTHLA/alignments/DRA_nuc.txt
../foreign/IMGTHLA/alignments/DRB1_gen.txt
../foreign/IMGTHLA/alignments/DRB3_gen.txt
../foreign/IMGTHLA/alignments/DRB4_gen.txt
../foreign/IMGTHLA/alignments/DRB_nuc.txt
../foreign/IMGTHLA/alignments/E_gen.txt
../foreign/IMGTHLA/alignments/E_nuc.txt
../foreign/IMGTHLA/alignments/F_gen.txt
../foreign/IMGTHLA/alignments/F_nuc.txt
../foreign/IMGTHLA/alignments/G_gen.txt
../foreign/IMGTHLA/alignments/G_nuc.txt
../foreign/IMGTHLA/alignments/HFE_gen.txt
../foreign/IMGTHLA/alignments/HFE_nuc.txt
../foreign/IMGTHLA/alignments/H_gen.txt
../foreign/IMGTHLA/alignments/H_nuc.txt
../foreign/IMGTHLA/alignments/J_gen.txt
../foreign/IMGTHLA/alignments/J_nuc.txt
../foreign/IMGTHLA/alignments/K_gen.txt
../foreign/IMGTHLA/alignments/K_nuc.txt
../foreign/IMGTHLA/alignments/L_gen.txt
../foreign/IMGTHLA/alignments/L_nuc.txt
../foreign/IMGTHLA/alignments/MICA_gen.txt
../foreign/IMGTHLA/alignments/MICA_nuc.txt
../foreign/IMGTHLA/alignments/MICB_gen.txt
../foreign/IMGTHLA/alignments/MICB_nuc.txt
../foreign/IMGTHLA/alignments/P_gen.txt
../foreign/IMGTHLA/alignments/TAP1_gen.txt
../foreign/IMGTHLA/alignments/TAP1_nuc.txt
../foreign/IMGTHLA/alignments/TAP2_gen.txt
../foreign/IMGTHLA/alignments/TAP2_nuc.txt
../foreign/IMGTHLA/alignments/V_gen.txt
../foreign/IMGTHLA/alignments/V_nuc.txt
../foreign/IMGTHLA/alignments/Y_gen.txt
../foreign/IMGTHLA/alignments/Y_nuc.txt"

function graph() {
  FILE=$1
  S=`basename $FILE`
  O=graph_2016_10_17/${S%.txt}

  ./mhc2gpdf.native --no-cache --no-open -f $FILE -o $O
}

export -f graph

parallel --shuf --workdir `pwd` -j -1 --env graph graph ::: $FILES
