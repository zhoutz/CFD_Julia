#!/bin/bash
set -x

clang++ icp.cpp -o icp -DNDEBUG -O3 -std=c++17
time ./icp
python3 plotting2.py
rm icp
