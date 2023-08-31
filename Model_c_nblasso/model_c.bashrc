#!/bin/bash


qsub model_c_p50_gauss.PBS 
qsub model_c_p50_copula.PBS 
qsub model_c_p50_cmcopula.PBS 

qsub model_c_p100_gauss.PBS 
qsub model_c_p100_copula.PBS 
qsub model_c_p100_cmcopula.PBS 