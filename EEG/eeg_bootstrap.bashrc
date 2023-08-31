#!/bin/bash

for args in `seq 1 50`;
do
qsub eeg_bootstrap.PBS -v "args=$args"
done
