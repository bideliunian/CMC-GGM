#!/bin/bash

for args in `seq 1 10`;
do
qsub D79_cv.PBS -v "args=$args"
done
