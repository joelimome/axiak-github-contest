#!/bin/bash

DIR=`dirname "$0"`
cd "$DIR"
cd src/py
./svd.py < ../../dat/mapped_test.txt > ../../dat/output_winners.dat
cd ../../dat
../src/pl/genoutput.pl < output_winners.dat > results-nomap.dat
cd ../src/cpp/map
./map -r < ../../../dat/results-nomap.dat > ../../../results.txt
