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

yalebox_dir=/projectnb/glaciermod/yalebox-piv/
version=$(git -C $yalebox_dir rev-parse HEAD)
standalone_dir=$yalebox_dir/standalone/$version/
standalone_file=$standalone_dir/piv_series_standalone
param_dir=/projectnb/glaciermod/faultless_wedge/experiment/proc/
param_file=$param_dir/XXX_piv_param.mat

mcr /$standalone_file $param_file $version
