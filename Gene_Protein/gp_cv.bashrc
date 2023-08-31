#!/bin/bash

for args in `seq 1 10`;
do
qsub gp_cv.PBS -v "args=$args"
done
