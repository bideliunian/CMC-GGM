#!/bin/bash

for args in `seq 1 50`;
do
qsub model_c_p50_cmcopula.PBS -v "args=$args"
qsub model_c_p100_cmcopula.PBS -v "args=$args"
done
