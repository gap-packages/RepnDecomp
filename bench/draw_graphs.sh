#!/bin/bash

./plot.py mymethod_symmetric symmetric_mymethod_{kronecker,naive,orbit_sum}.csv log
./plot.py serre_symmetric symmetric_serre_{naive,kronecker}.csv log
./plot.py mymethod_cyclic cyclic_mymethod_{naive,orbit_sum,kronecker}.csv log
./plot.py serre_cyclic cyclic_serre_{naive,kronecker}.csv log
./plot.py mymethod_tensor tensor_mymethod_{naive,orbit_sum,kronecker}.csv log
./plot.py serre_tensor tensor_serre_{naive,kronecker}.csv log
./plot.py mymethod_cyclic_naive cyclic_mymethod_{,parallel_}naive.csv log
./plot.py mymethod_cyclic_kronecker cyclic_mymethod_{,parallel_}kronecker.csv log
./plot.py serre_cyclic_naive cyclic_serre_{,parallel_}naive.csv log
./plot.py serre_cyclic_kronecker cyclic_serre_{,parallel_}kronecker.csv log
./plot.py vs_random random_{serre_{naive,kronecker},mymethod_parallel_kronecker}.csv log
./plot.py vs_cyclic cyclic_{serre_{naive,kronecker},mymethod_parallel_kronecker}.csv log
./plot.py vs_symmetric symmetric_{serre_{naive,kronecker},mymethod_parallel_kronecker}.csv log
