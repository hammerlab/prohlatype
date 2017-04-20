#!/usr/bin/env bash

set -e

echo A_gen && ./mas_align_test.native A_gen
echo A_nuc && ./mas_align_test.native A_nuc
echo B_gen && ./mas_align_test.native B_gen
echo B_nuc && ./mas_align_test.native B_nuc
echo C_gen && ./mas_align_test.native C_gen
echo C_nuc && ./mas_align_test.native C_nuc "C*04:09N"
