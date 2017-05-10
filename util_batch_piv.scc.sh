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

version=$(git -C $(dirname $standalone) rev-parse HEAD)
standalone_dir=/projectnb/glaciermod/yalebox-piv/standalone/$version/
standalone_file=piv_series_standalone
param_dir=/projectnb/glaciermod/faultless_wedge/experiment/proc/
param_file=XXX_piv_param.mat

mcr $standalone_dir/$standalone_file $param_dir/$param_file $version
