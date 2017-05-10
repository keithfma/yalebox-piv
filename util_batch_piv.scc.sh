#!/bin/bash -l
#
#  Submission script for experiment XXX PIV
#
#$ -N XXX_piv
#$ -P glaciermod
#$ -j y
#$ -pe omp 16
#$ -l h_rt=72:00:00

module load mcr/9.0.1_2016a 

standalone=XXX/piv_series_standalone
param_file=XXX/XXX_piv_param.mat
version=$(git -C $(dirname $standalone) rev-parse HEAD)

mcr $standalone $param_file $version
