#!/bin/bash
set -x

clang++ ftcs.cpp -o ftcs -DNDEBUG -O3 -std=c++17
time ./ftcs
python3 plotting2.py
rm ftcs
