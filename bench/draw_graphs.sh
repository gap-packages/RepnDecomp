#!/bin/bash

./plot.py mymethod_symmetric symmetric_mymethod_{kronecker,naive,orbit_sum}.csv log
./plot.py serre_symmetric symmetric_serre_{naive,kronecker}.csv log
./plot.py random random_*.csv nolog
