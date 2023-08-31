#!/bin/bash

for args in `seq 1 50`;
do
qsub model_a_p100.PBS -v "args=$args"
done
