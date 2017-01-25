#!/usr/bin/env bash

set -e

echo A_gen && ./mas_align_test.native A_gen
echo A_nuc && ./mas_align_test.native A_nuc "A*01:11N"
echo B_gen && ./mas_align_test.native B_gen
echo B_nuc && ./mas_align_test.native B_nuc "B*44:02:01:02S" "B*07:44N"
echo C_gen && ./mas_align_test.native C_gen
echo C_nuc && ./mas_align_test.native C_nuc "C*04:09N"
