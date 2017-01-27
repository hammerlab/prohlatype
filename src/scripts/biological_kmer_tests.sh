#!/usr/bin/env bash

set -e

echo A_gen    && ./biological_kmers.native A_gen
echo A_nuc    && ./biological_kmers.native A_nuc      "A*01:11N"
echo B_gen    && ./biological_kmers.native B_gen
echo B_nuc    && ./biological_kmers.native B_nuc      "B*44:02:01:02S" "B*07:44N"
echo C_gen    && ./biological_kmers.native C_gen
echo C_nuc    && ./biological_kmers.native C_nuc      "C*04:09N"
echo DRB_nuc  && ./biological_kmers.native DRB_nuc    "DRB4*01:03:01:02N"
echo DRB1_gen && ./biological_kmers.native DRB1_gen
echo DMA_gen  && ./biological_kmers.native DMA_gen
echo DMA_nuc  && ./biological_kmers.native DMA_nuc
echo DMB_gen  && ./biological_kmers.native DMB_gen
echo DMB_nuc  && ./biological_kmers.native DMB_nuc
echo DOA_gen  && ./biological_kmers.native DOA_gen    "DOA*01:01:04:02" "DOA*01:01:02:01"
echo DOA_nuc  && ./biological_kmers.native DOA_nuc
echo DOB_gen  && ./biological_kmers.native DOB_gen    "DOB*01:04:01:02"
echo DOB_nuc  && ./biological_kmers.native DOB_nuc
echo DPA1_gen && ./biological_kmers.native DPA1_gen
echo DPA1_nuc && ./biological_kmers.native DPA1_nuc
echo DPB1_gen && ./biological_kmers.native DPB1_gen
echo DPB1_nuc && ./biological_kmers.native DPB1_nuc
echo DPB2_gen && ./biological_kmers.native DPB2_gen   "DPB2*01:01:01"
echo DPB2_nuc && ./biological_kmers.native DPB2_nuc
echo DQA1_gen && ./biological_kmers.native DQA1_gen
echo DQA1_nuc && ./biological_kmers.native DQA1_nuc
echo DQB1_gen && ./biological_kmers.native DQB1_gen
echo DQB1_nuc && ./biological_kmers.native DQB1_nuc   "DQB1*06:79:01"
echo DRA_gen  && ./biological_kmers.native DRA_gen
echo DRA_nuc  && ./biological_kmers.native DRA_nuc
echo DRB1_gen && ./biological_kmers.native DRB1_gen
echo DRB3_gen && ./biological_kmers.native DRB3_gen
echo DRB4_gen && ./biological_kmers.native DRB4_gen
echo DRB_nuc  && ./biological_kmers.native DRB_nuc    "DRB4*01:03:01:02N"
echo E_gen    && ./biological_kmers.native E_gen
echo E_nuc    && ./biological_kmers.native E_nuc
echo F_gen    && ./biological_kmers.native F_gen
echo F_nuc    && ./biological_kmers.native F_nuc
echo G_gen    && ./biological_kmers.native G_gen
echo G_nuc    && ./biological_kmers.native G_nuc
echo HFE_gen  && ./biological_kmers.native HFE_gen
echo HFE_nuc  && ./biological_kmers.native HFE_nuc

#Bug?  AAGCACAGAC
#echo H_gen    && ./biological_kmers.native H_gen

echo H_nuc    && ./biological_kmers.native H_nuc

# Bug? ACCACCCTCT
#echo J_gen    && ./biological_kmers.native J_gen
# Bug? ACCACCCTCT
#echo J_nuc    && ./biological_kmers.native J_nuc

echo K_gen    && ./biological_kmers.native K_gen
echo K_nuc    && ./biological_kmers.native K_nuc
echo L_gen    && ./biological_kmers.native L_gen
echo L_nuc    && ./biological_kmers.native L_nuc
echo MICA_gen && ./biological_kmers.native MICA_gen
echo MICA_nuc && ./biological_kmers.native MICA_nuc
echo MICB_gen && ./biological_kmers.native MICB_gen
echo MICB_nuc && ./biological_kmers.native MICB_nuc
echo P_gen    && ./biological_kmers.native P_gen
echo TAP1_gen && ./biological_kmers.native TAP1_gen
echo TAP1_nuc && ./biological_kmers.native TAP1_nuc
echo TAP2_gen && ./biological_kmers.native TAP2_gen
echo TAP2_nuc && ./biological_kmers.native TAP2_nuc
echo V_gen    && ./biological_kmers.native V_gen
echo V_nuc    && ./biological_kmers.native V_nuc

# No fasta's for comparison
# echo Y_gen &&  ./biological_kmers.native Y_gen
# echo Y_nuc &&  ./biological_kmers.native Y_nuc
