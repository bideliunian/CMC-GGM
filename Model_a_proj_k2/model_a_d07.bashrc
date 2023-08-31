#!/bin/bash

for args in `seq 1 50`;
do
qsub model_a_d07_exp.PBS -v "args=$args"
qsub model_a_d07_power.PBS -v "args=$args"
done
