#!/bin/bash

for args in `seq 1 100`;
do
qsub model_c.PBS -v "args=$args"
done
