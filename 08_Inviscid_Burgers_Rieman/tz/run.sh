#!/bin/bash
set -xe

clang++ burgers_riemann.cpp -o a -DNDEBUG -O3 -std=c++17
time ./a
python3 plotting.py
rm a
rm solution_riemann.csv
