#!/usr/bin/env bash

set -e

echo A_gen && ./mas_align.exe A_gen
echo A_nuc && ./mas_align.exe A_nuc
echo B_gen && ./mas_align.exe B_gen
echo B_nuc && ./mas_align.exe B_nuc
echo C_gen && ./mas_align.exe C_gen
echo C_nuc && ./mas_align.exe C_nuc "C*04:09N"
