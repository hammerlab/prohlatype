#!/usr/bin/env bash

set -e

echo A_gen    && ./biological_kmers.exe A_gen
# They really messed up this fasta file
echo A_nuc    && ./biological_kmers.exe A_nuc      "A*01:11N" "fA*01:11N"
echo B_gen    && ./biological_kmers.exe B_gen
echo B_nuc    && ./biological_kmers.exe B_nuc      "B*44:02:01:02S" "B*07:44N" "fB*07:44N"
echo C_gen    && ./biological_kmers.exe C_gen
echo C_nuc    && ./biological_kmers.exe C_nuc      "C*04:09N"
#echo DRB_nuc  && ./biological_kmers.exe DRB_nuc    "DRB4*01:03:01:02N"
echo DRB1_gen && ./biological_kmers.exe DRB1_gen   "DRB1*15:02:01:01" "DRB1*12:01:01:01"
echo DMA_gen  && ./biological_kmers.exe DMA_gen
echo DMA_nuc  && ./biological_kmers.exe DMA_nuc
echo DMB_gen  && ./biological_kmers.exe DMB_gen
echo DMB_nuc  && ./biological_kmers.exe DMB_nuc
echo DOA_gen  && ./biological_kmers.exe DOA_gen    "DOA*01:01:04:02" "DOA*01:01:02:01"
echo DOA_nuc  && ./biological_kmers.exe DOA_nuc
echo DOB_gen  && ./biological_kmers.exe DOB_gen    "DOB*01:04:01:02"
echo DOB_nuc  && ./biological_kmers.exe DOB_nuc
echo DPA1_gen && ./biological_kmers.exe DPA1_gen
echo DPA1_nuc && ./biological_kmers.exe DPA1_nuc
echo DPB1_gen && ./biological_kmers.exe DPB1_gen   "DPB1*13:01:02"
echo DPB1_nuc && ./biological_kmers.exe DPB1_nuc
echo DPB2_gen && ./biological_kmers.exe DPB2_gen   "DPB2*01:01:01"
echo DPB2_nuc && ./biological_kmers.exe DPB2_nuc
echo DQA1_gen && ./biological_kmers.exe DQA1_gen
echo DQA1_nuc && ./biological_kmers.exe DQA1_nuc
echo DQB1_gen && ./biological_kmers.exe DQB1_gen
echo DQB1_nuc && ./biological_kmers.exe DQB1_nuc   "DQB1*06:79:01"
echo DRA_gen  && ./biological_kmers.exe DRA_gen
echo DRA_nuc  && ./biological_kmers.exe DRA_nuc
echo DRB1_gen && ./biological_kmers.exe DRB1_gen   "DRB1*12:01:01:01" "DRB1*15:02:01:01"
echo DRB3_gen && ./biological_kmers.exe DRB3_gen
echo DRB4_gen && ./biological_kmers.exe DRB4_gen
#echo DRB_nuc  && ./biological_kmers.exe DRB_nuc    "DRB4*01:03:01:02N"
echo E_gen    && ./biological_kmers.exe E_gen
echo E_nuc    && ./biological_kmers.exe E_nuc
echo F_gen    && ./biological_kmers.exe F_gen
echo F_nuc    && ./biological_kmers.exe F_nuc
echo G_gen    && ./biological_kmers.exe G_gen      "G*01:01:01:08"
echo G_nuc    && ./biological_kmers.exe G_nuc
echo HFE_gen  && ./biological_kmers.exe HFE_gen
echo HFE_nuc  && ./biological_kmers.exe HFE_nuc

#Bug?  AAGCACAGAC
#echo H_gen    && ./biological_kmers.exe H_gen

echo H_nuc    && ./biological_kmers.exe H_nuc

# Bug? ACCACCCTCT
#echo J_gen    && ./biological_kmers.exe J_gen
# Bug? ACCACCCTCT
#echo J_nuc    && ./biological_kmers.exe J_nuc

echo K_gen    && ./biological_kmers.exe K_gen
echo K_nuc    && ./biological_kmers.exe K_nuc
echo L_gen    && ./biological_kmers.exe L_gen
echo L_nuc    && ./biological_kmers.exe L_nuc
echo MICA_gen && ./biological_kmers.exe MICA_gen
echo MICA_nuc && ./biological_kmers.exe MICA_nuc
echo MICB_gen && ./biological_kmers.exe MICB_gen
echo MICB_nuc && ./biological_kmers.exe MICB_nuc
echo P_gen    && ./biological_kmers.exe P_gen
echo TAP1_gen && ./biological_kmers.exe TAP1_gen
echo TAP1_nuc && ./biological_kmers.exe TAP1_nuc
echo TAP2_gen && ./biological_kmers.exe TAP2_gen
echo TAP2_nuc && ./biological_kmers.exe TAP2_nuc
echo V_gen    && ./biological_kmers.exe V_gen
echo V_nuc    && ./biological_kmers.exe V_nuc

# No fasta's for comparison
# echo Y_gen &&  ./biological_kmers.exe Y_gen
# echo Y_nuc &&  ./biological_kmers.exe Y_nuc
