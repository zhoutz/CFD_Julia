#!/bin/bash
set -x

clang++ cn.cpp -o cn -DNDEBUG -O3 -std=c++17
time ./cn
python3 plotting2.py
rm cn
