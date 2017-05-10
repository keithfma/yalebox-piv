#!/bin/bash -l
#
#  Submission script for experiment XXX PIV
#
#$ -N XXX_piv
#$ -P glaciermod
#$ -j y
#$ -pe omp 16
#$ -l h_rt=72:00:00

# define parameters
standalone=XXX/piv_series_standalone
param_file=XXX/XXX_piv_param.mat

# run
$standalone $param_file
