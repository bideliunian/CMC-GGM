#!/bin/bash

for args in `seq 1 50`;
do
qsub model_a_all.PBS -v "args=$args"
done
