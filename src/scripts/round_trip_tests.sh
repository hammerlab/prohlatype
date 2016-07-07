#!/usr/bin/env bash

set -e
set -x

./round_trip.native 
./round_trip.native A_gen "A*01:01:01:01"
./round_trip.native B_gen "B*07:02:01"
./round_trip.native B_nuc "B*07:02:01"
./round_trip.native C_gen "C*01:02:01"
./round_trip.native C_nuc "C*01:02:01"
./round_trip.native DRB_nuc "DRB1*01:01:01"
./round_trip.native DRB1_gen "DRB1*01:02:01"
