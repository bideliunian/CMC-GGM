#!/bin/bash

for args in `seq 1 20`;
do
qsub eeg_tune.PBS -v "args=$args"
done
