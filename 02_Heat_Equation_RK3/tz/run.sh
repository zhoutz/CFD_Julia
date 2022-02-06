#!/bin/bash
set -x

clang++ rk3.cpp -o rk3 -DNDEBUG -O3 -std=c++17
time ./rk3
python3 plotting2.py
rm rk3
