#!/usr/bin/env bash

set -e

echo A_gen && ./mas_align_test.exe A_gen
echo A_nuc && ./mas_align_test.exe A_nuc
echo B_gen && ./mas_align_test.exe B_gen
echo B_nuc && ./mas_align_test.exe B_nuc
echo C_gen && ./mas_align_test.exe C_gen
echo C_nuc && ./mas_align_test.exe C_nuc "C*04:09N"
