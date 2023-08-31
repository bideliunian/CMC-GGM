#!/bin/bash

for args in `seq 1 50`;
do
qsub model_b_p50.PBS -v "args=$args"
done
