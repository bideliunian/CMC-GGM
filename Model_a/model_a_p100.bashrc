#!/bin/bash

for args in `seq 1 50`;
do
qsub model_a_p100_gauss.PBS -v "args=$args"
qsub model_a_p100_copula.PBS -v "args=$args"
qsub model_a_p100_cmcopula.PBS -v "args=$args"
done
