#!/bin/bash
set -xe

clang++ weno_dirichlet.cpp -o a -DNDEBUG -O3 -std=c++17
time ./a
clang++ weno_periodic.cpp -o a -DNDEBUG -O3 -std=c++17
time ./a
python3 plotting.py
rm a solution_d.csv solution_p.csv
