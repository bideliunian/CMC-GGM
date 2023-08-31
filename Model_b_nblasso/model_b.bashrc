#!/bin/bash


qsub model_b_p50_gauss.PBS 
qsub model_b_p50_copula.PBS 
qsub model_b_p50_cmcopula.PBS 

qsub model_b_p100_gauss.PBS 
qsub model_b_p100_copula.PBS 
qsub model_b_p100_cmcopula.PBS 
