#!/bin/bash
set -xe

clang++ burger_ftcs.cpp -o a -DNDEBUG -O3 -std=c++17
time ./a
python3 plotting2.py
rm a
