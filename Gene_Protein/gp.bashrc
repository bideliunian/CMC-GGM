#!/bin/bash

for args in `seq 1 50`;
do
qsub gp.PBS -v "args=$args"
done
